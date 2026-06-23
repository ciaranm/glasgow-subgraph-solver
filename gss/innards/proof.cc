#include <gss/innards/proof.hh>

#include <algorithm>
#include <fstream>
#include <iterator>
#include <map>
#include <memory>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <utility>

using namespace gss;
using namespace gss::innards;

using std::copy;
using std::endl;
using std::find;
using std::function;
using std::istreambuf_iterator;
using std::make_unique;
using std::map;
using std::max;
using std::min;
using std::move;
using std::ofstream;
using std::optional;
using std::ostream;
using std::ostreambuf_iterator;
using std::pair;
using std::set;
using std::string;
using std::stringstream;
using std::to_string;
using std::tuple;
using std::unique_ptr;
using std::unordered_map;
using std::vector;

ProofError::ProofError(const string & m) noexcept :
    _message("Proof error: " + m)
{
}

auto ProofError::what() const noexcept -> const char *
{
    return _message.c_str();
}

struct Proof::Imp
{
    string opb_filename, log_filename;
    stringstream model_stream, model_prelude_stream;
    unique_ptr<ostream> proof_stream;
    bool super_extra_verbose = false;

    map<pair<long, long>, string> variable_mappings;
    map<long, string> binary_variable_mappings;
    map<tuple<long, long, long>, string> connected_variable_mappings;
    map<tuple<long, long, long, long>, string> connected_variable_mappings_aux;
    map<long, string> at_least_one_value_constraints, at_most_one_value_constraints, injectivity_constraints;
    map<pair<long, long>, string> locally_injective_constraints;
    bool locally_injective = false;
    // The adjacency-line proof state: labels, the numeric proof-line ids of supplemental
    // adjacency lines (so transient ones re-derived for a degree/NDS check can be deleted,
    // see weaken_supplemental_adjacency), and permitted target sets. Grouped so the
    // derivations consuming it can migrate to the middle layer (adjacency.labels/.ids/.permitted).
    AdjacencyProofLines adjacency;
    map<pair<long, long>, long> eliminations;
    map<pair<long, long>, string> non_edge_constraints;
    long objective_line = 0;
    stringstream objective_sum;

    unordered_map<string, string> cached_proof_lines;

    long nb_constraints = 0;
    long proof_line = 0;
    int largest_level_set = 0;
    int active_level = 0;

    bool clique_encoding = false;
    bool doing_mcs_by_clique = false;

    bool doing_hom_colour_proof = false;
    NamedVertex hom_colour_proof_p, hom_colour_proof_t;
    vector<NamedVertex> p_clique;
    map<int, NamedVertex> t_clique_neighbourhood;
    map<pair<pair<NamedVertex, NamedVertex>, pair<NamedVertex, NamedVertex>>, long> clique_for_hom_non_edge_constraints;

    vector<pair<int, int>> zero_in_proof_objectives;
};

Proof::Proof(const ProofOptions & options) :
    _imp(make_unique<Imp>())
{
    _imp->opb_filename = options.opb_file;
    _imp->log_filename = options.log_file;
    _imp->super_extra_verbose = options.super_extra_verbose;
}

Proof::Proof(Proof &&) = default;

Proof::~Proof() = default;

auto Proof::operator=(Proof &&) -> Proof & = default;

auto Proof::create_cp_variable(int pattern_vertex, int target_size,
    const function<auto(int)->string> & pattern_name,
    const function<auto(int)->string> & target_name) -> void
{
    for (int i = 0; i < target_size; ++i)
        _imp->variable_mappings.emplace(pair{pattern_vertex, i}, pattern_name(pattern_vertex) + "_" + target_name(i));

    _imp->model_stream << "* vertex " << pattern_vertex << " domain\n";
    auto al1_label = "@al1" + pattern_name(pattern_vertex);
    stringstream al1_constraint;
    al1_constraint << al1_label << " ";
    for (int i = 0; i < target_size; ++i)
        al1_constraint << "1 x" << _imp->variable_mappings[{pattern_vertex, i}] << " ";
    al1_constraint << ">= 1";
    _imp->model_stream << al1_constraint.str() << " ;\n";
    ++_imp->nb_constraints;
    _imp->at_least_one_value_constraints.emplace(pattern_vertex, al1_label);

    auto am1_label = "@am1" + pattern_name(pattern_vertex);
    stringstream am1_constraint;
    am1_constraint << am1_label << " ";
    for (int i = 0; i < target_size; ++i)
        am1_constraint << "-1 x" << _imp->variable_mappings[{pattern_vertex, i}] << " ";
    am1_constraint << ">= -1";
    _imp->model_stream << am1_constraint.str() << " ;\n";
    ++_imp->nb_constraints;
    _imp->at_most_one_value_constraints.emplace(pattern_vertex, am1_label);
}

auto Proof::create_injectivity_constraints(int pattern_size, int target_size,
    const function<auto(int)->string> & target_name) -> void
{
    for (int v = 0; v < target_size; ++v) {
        _imp->model_stream << "* injectivity on value " << v << '\n';
        auto inj_label = "@inj" + target_name(v);
        _imp->model_stream << inj_label << " ";

        for (int p = 0; p < pattern_size; ++p) {
            auto x = _imp->variable_mappings.find(pair{p, v});
            if (x != _imp->variable_mappings.end())
                _imp->model_stream << "-1 x" << x->second << " ";
        }
        _imp->model_stream << ">= -1 ;\n";
        ++_imp->nb_constraints;
        _imp->injectivity_constraints.emplace(v, inj_label);
    }
}

auto Proof::create_locally_injective_constraints(int pattern_size, int target_size,
    const function<auto(int, int)->bool> & adjacent,
    const function<auto(int)->string> & pattern_name,
    const function<auto(int)->string> & target_name) -> void
{
    _imp->locally_injective = true;

    for (int v = 0; v < pattern_size; ++v) {
        // the neighbourhood of v: if v has a self-loop then v is its own neighbour, so
        // local injectivity also forces phi(v) to differ from its neighbours' images.
        vector<int> neighbours;
        for (int u = 0; u < pattern_size; ++u)
            if (adjacent(v, u))
                neighbours.push_back(u);

        // at most one neighbour maps to a given target only bites for |N(v)| >= 2
        if (neighbours.size() < 2)
            continue;

        for (int t = 0; t < target_size; ++t) {
            _imp->model_stream << "* local injectivity on neighbourhood of " << v << " for value " << t << '\n';
            auto label = "@linj" + pattern_name(v) + "_" + target_name(t);
            _imp->model_stream << label;
            for (auto & u : neighbours)
                _imp->model_stream << " -1 x" << _imp->variable_mappings[pair{u, t}];
            _imp->model_stream << " >= -1 ;\n";
            ++_imp->nb_constraints;
            _imp->locally_injective_constraints.emplace(pair{v, t}, label);
        }
    }
}

auto Proof::create_forbidden_assignment_constraint(int p, int t) -> void
{
    auto & var_name = _imp->variable_mappings[pair{p, t}];
    _imp->model_stream << "* forbidden assignment\n";
    _imp->model_stream << "@forb" << var_name << " 1 ~x" << var_name << " >= 1 ;\n";
    ++_imp->nb_constraints;
    _imp->eliminations.emplace(pair{p, t}, _imp->nb_constraints);
}

auto Proof::start_adjacency_constraints_for(int p, int t) -> void
{
    _imp->model_stream << "* adjacency " << p << " maps to " << t << '\n';
}

auto Proof::create_adjacency_constraint(const NamedVertex & p, const NamedVertex & q, const NamedVertex & t,
    const vector<int> & uu, bool) -> void
{
    auto adj_label = "@adj" + p.second + "_" + t.second + "_" + q.second;

    _imp->model_stream << adj_label << " ";
    _imp->model_stream << "1 ~x" << _imp->variable_mappings[pair{p.first, t.first}];
    for (auto & u : uu)
        _imp->model_stream << " 1 x" << _imp->variable_mappings[pair{q.first, u}];
    _imp->model_stream << " >= 1 ;\n";
    ++_imp->nb_constraints;
    _imp->adjacency.labels.emplace(tuple{0, p.first, q.first, t.first}, adj_label);
    _imp->adjacency.permitted.emplace(tuple{0, p.first, q.first, t.first}, vector<long>(uu.begin(), uu.end()));
}

auto Proof::emit_preserved_assignment_variables() -> void
{
    // List exactly the assignment variables as the projected (preserved) set,
    // so VeriPB counts solutions in terms of the high-level mapping. The line
    // goes into the prelude, which finalise_model writes immediately after the
    // header and before any constraint.
    _imp->model_prelude_stream << "preserved:";
    for (auto & [_, name] : _imp->variable_mappings)
        _imp->model_prelude_stream << " x" << name;
    _imp->model_prelude_stream << " ;\n";
}

auto Proof::finalise_model() -> void
{
    unique_ptr<ostream> f = make_unique<ofstream>(_imp->opb_filename);

    *f << "* #variable= " << (_imp->variable_mappings.size() + _imp->binary_variable_mappings.size() + _imp->connected_variable_mappings.size() + _imp->connected_variable_mappings_aux.size())
       << " #constraint= " << _imp->nb_constraints << ";\n";
    copy(istreambuf_iterator<char>{_imp->model_prelude_stream}, istreambuf_iterator<char>{}, ostreambuf_iterator<char>{*f});
    _imp->model_prelude_stream.clear();
    copy(istreambuf_iterator<char>{_imp->model_stream}, istreambuf_iterator<char>{}, ostreambuf_iterator<char>{*f});
    _imp->model_stream.clear();

    if (! *f)
        throw ProofError{"Error writing opb file to '" + _imp->opb_filename + "'"};

    _imp->proof_stream = make_unique<ofstream>(_imp->log_filename);

    *_imp->proof_stream << "pseudo-Boolean proof version 3.0\n";

    *_imp->proof_stream << "f " << _imp->nb_constraints << " ;\n";
    _imp->proof_line += _imp->nb_constraints;

    if (! *_imp->proof_stream)
        throw ProofError{"Error writing proof file to '" + _imp->log_filename + "'"};
}

auto Proof::finish_unsat_proof() -> void
{
    *_imp->proof_stream << "% asserting that we've proved unsat\n";
    *_imp->proof_stream << "@unsatconc rup >= 1 ;\n";
    ++_imp->proof_line;
    *_imp->proof_stream << "output NONE;\n"
                        << "conclusion UNSAT : -1;\n"
                        << "end pseudo-Boolean proof;\n";
}

auto Proof::finish_sat_proof() -> void
{
    *_imp->proof_stream << "output NONE;\n"
        << "conclusion SAT;\n"
        << "end pseudo-Boolean proof;\n";
}

auto Proof::finish_enumeration_proof(const loooong & number_of_solutions, bool complete) -> void
{
    if (complete) {
        // Every solution has been logged with solx and excluded, so the
        // remaining problem is unsatisfiable: assert the contradiction and
        // conclude a complete enumeration of exactly this many solutions.
        *_imp->proof_stream << "% asserting that we've enumerated every solution\n";
        *_imp->proof_stream << "rup >= 1 ;\n";
        ++_imp->proof_line;
        *_imp->proof_stream << "output NONE;\n"
                            << "conclusion ENUMERATION_COMPLETE " << number_of_solutions << " : -1;\n"
                            << "end pseudo-Boolean proof;\n";
    }
    else {
        // We stopped early (timeout or solution limit): we have witnessed this
        // many solutions but make no claim that there are no others.
        *_imp->proof_stream << "output NONE;\n"
                            << "conclusion ENUMERATION_PARTIAL " << number_of_solutions << ";\n"
                            << "end pseudo-Boolean proof;\n";
    }
}

auto Proof::finish_unknown_proof() -> void
{
    *_imp->proof_stream << "output NONE;\n"
        << "conclusion NONE;\n"
        << "end pseudo-Boolean proof;\n";
}

auto Proof::finish_optimisation_proof(int size) -> void
{
    *_imp->proof_stream << "rup" << _imp->objective_sum.str() << " >= " << size << ";\n";
    *_imp->proof_stream << "output NONE;\n"
        << "conclusion BOUNDS " << size << " " << size << ";\n"
        << "end pseudo-Boolean proof;\n";
}

auto Proof::failure_due_to_pattern_bigger_than_target() -> void
{
    *_imp->proof_stream << "% failure due to the pattern being bigger than the target\n";

    // we get a hall violator by adding up all of the things
    *_imp->proof_stream << "@ptbig pol";
    bool first = true;

    for (auto & [_, label] : _imp->at_least_one_value_constraints) {
        if (first) {
            *_imp->proof_stream << " " << label;
            first = false;
        }
        else
            *_imp->proof_stream << " " << label << " +";
    }

    for (auto & [_, label] : _imp->injectivity_constraints)
        *_imp->proof_stream << " " << label << " +";
    *_imp->proof_stream << " ;\n";
    ++_imp->proof_line;
}

auto Proof::need_elimination(int p, int t) -> void
{
    if (! _imp->eliminations.contains(pair{p, t})) {
        auto & var_name = _imp->variable_mappings[pair{p, t}];
        *_imp->proof_stream << "setlvl 0;\n";
        *_imp->proof_stream << "@elimnds" << var_name << " rup 1 ~x" << var_name << " >= 1 ;\n";
        _imp->eliminations[pair{p, t}] = ++_imp->proof_line;
        *_imp->proof_stream << "setlvl " << _imp->active_level << ";\n";
    }
}

auto Proof::incompatible_by_degrees(
    int g,
    const NamedVertex & p,
    const vector<int> & n_p,
    const NamedVertex & t,
    const vector<int> & n_t) -> void
{
    *_imp->proof_stream << "% cannot map " << p.second << " to " << t.second << " due to degrees in graph pairs " << g << '\n';

    auto & var_deg = _imp->variable_mappings[pair{p.first, t.first}];
    bool first_time_deg = ! _imp->eliminations.contains(pair{p.first, t.first});
    if (first_time_deg)
        *_imp->proof_stream << "@elimdegpol" << var_deg << " ";
    else
        *_imp->proof_stream << "@reelimdegpol" << g << "_" << var_deg << " ";
    *_imp->proof_stream << "pol";
    bool first = true;
    for (auto & n : n_p) {
        // due to loops or labels, it might not be possible to map n to t.first
        if (_imp->adjacency.labels.count(tuple{g, p.first, n, t.first})) {
            if (first) {
                first = false;
                *_imp->proof_stream << " " << _imp->adjacency.labels.at(tuple{g, p.first, n, t.first});
            }
            else
                *_imp->proof_stream << " " << _imp->adjacency.labels.at(tuple{g, p.first, n, t.first}) << " +";
        }
    }

    // if I map p to t, I have to map the neighbours of p to distinct neighbours of t.
    // Under full injectivity that distinctness is the global injectivity on each value;
    // under local injectivity it is the neighbourhood-injectivity of p (phi|N(p) is
    // injective), which is exactly what the degree pigeonhole needs.
    for (auto & n : n_t)
        *_imp->proof_stream << " " << (_imp->locally_injective ? _imp->locally_injective_constraints[pair{p.first, n}] : _imp->injectivity_constraints[n]) << " +";

    *_imp->proof_stream << " s ;\n";
    ++_imp->proof_line;

    if (first_time_deg)
        *_imp->proof_stream << "@elimdeg" << var_deg << " ";
    else
        *_imp->proof_stream << "@reelimdeg" << g << "_" << var_deg << " ";
    *_imp->proof_stream << "ia 1 ~x" << var_deg << " >= 1 : " << _imp->proof_line << " ;\n";
    ++_imp->proof_line;
    _imp->eliminations.emplace(pair{p.first, t.first}, _imp->proof_line);

    *_imp->proof_stream << "del id " << _imp->proof_line - 1 << " ;\n";
}

auto Proof::incompatible_by_nds(
    int g,
    const NamedVertex & p,
    const NamedVertex & t,
    const vector<int> & p_subsequence,
    const vector<int> & t_subsequence,
    const vector<int> & t_remaining) -> void
{
    *_imp->proof_stream << "% cannot map " << p.second << " to " << t.second << " due to nds in graph pairs " << g << '\n';

    for (auto & n : p_subsequence)
        for (auto & u : t_remaining)
            need_elimination(n, u);
    for (auto & n : p_subsequence)
        need_elimination(n, t_subsequence.back());

    auto & var_nds = _imp->variable_mappings[pair{p.first, t.first}];
    bool first_time_nds = ! _imp->eliminations.contains(pair{p.first, t.first});
    if (first_time_nds)
        *_imp->proof_stream << "@elimndspol" << var_nds << " ";
    else
        *_imp->proof_stream << "@reelimndspol" << g << "_" << var_nds << " ";
    // summing up horizontally
    *_imp->proof_stream << "pol";
    bool first = true;
    for (auto & n : p_subsequence) {
        // due to loops or labels, it might not be possible to map n to t.first
        if (_imp->adjacency.labels.count(tuple{g, p.first, n, t.first})) {
            if (first) {
                first = false;
                *_imp->proof_stream << " " << _imp->adjacency.labels.at(tuple{g, p.first, n, t.first});
            }
            else
                *_imp->proof_stream << " " << _imp->adjacency.labels.at(tuple{g, p.first, n, t.first}) << " +";
        }
    }

    // injectivity in the square: each column of the square holds at most one of p's
    // neighbours. Under full injectivity that is the global injectivity on the value;
    // under local injectivity it is the neighbourhood-injectivity of p (at most one
    // neighbour of p maps to t), exactly as in the degree pigeonhole above.
    for (auto & t : t_subsequence) {
        if (t != t_subsequence.back())
            *_imp->proof_stream << " " << (_imp->locally_injective ? _imp->locally_injective_constraints.at(pair{p.first, t}) : _imp->injectivity_constraints.at(t)) << " +";
    }

    // block to the right of the failing square
    for (auto & n : p_subsequence) {
        for (auto & u : t_remaining) {
            /* n -> t is already eliminated by degree or loop */
            *_imp->proof_stream << " " << _imp->eliminations[pair{n, u}] << " +";
        }
    }

    // final column
    for (auto & n : p_subsequence) {
        /* n -> t is already eliminated by degree or loop */
        *_imp->proof_stream << " " << _imp->eliminations[pair{n, t_subsequence.back()}] << " +";
    }

    *_imp->proof_stream << " s ;\n";
    ++_imp->proof_line;

    if (first_time_nds)
        *_imp->proof_stream << "@elimndsconc" << var_nds << " ";
    else
        *_imp->proof_stream << "@reelimndsconc" << g << "_" << var_nds << " ";
    *_imp->proof_stream << "ia 1 ~x" << var_nds << " >= 1 : " << _imp->proof_line << " ;\n";
    ++_imp->proof_line;

    *_imp->proof_stream << "del id " << _imp->proof_line - 1 << " ;\n";
}

auto Proof::incompatible_by_loops(
    const NamedVertex & p,
    const NamedVertex & t) -> void
{
    // may be requested both up front (so the unit is available to later derivations and
    // search propagations) and again during domain initialisation: only emit it once.
    if (_imp->eliminations.contains(pair{p.first, t.first}))
        return;
    auto & var_loop = _imp->variable_mappings[pair{p.first, t.first}];
    *_imp->proof_stream << "% cannot map " << p.second << " to " << t.second << " due to loop\n";
    *_imp->proof_stream << "@loop" << var_loop << " rup 1 ~x" << var_loop << " >= 1 ;\n";
    _imp->eliminations.emplace(pair{p.first, t.first}, ++_imp->proof_line);
}

auto Proof::initial_domain_is_empty(int p, const string & where) -> void
{
    *_imp->proof_stream << "% failure due to domain " << p << " being empty at " << where << '\n';
}

auto Proof::emit_hall_set_or_violator(const vector<NamedVertex> & lhs, const vector<NamedVertex> & rhs) -> void
{
    *_imp->proof_stream << "% hall set or violator {";
    for (auto & l : lhs)
        *_imp->proof_stream << " " << l.second;
    *_imp->proof_stream << " } / {";
    for (auto & r : rhs)
        *_imp->proof_stream << " " << r.second;
    *_imp->proof_stream << " }\n";

    *_imp->proof_stream << "@hall" << (_imp->proof_line + 1) << " pol";
    bool first = true;
    for (auto & l : lhs) {
        if (first) {
            first = false;
            *_imp->proof_stream << " " << _imp->at_least_one_value_constraints.at(l.first);
        }
        else
            *_imp->proof_stream << " " << _imp->at_least_one_value_constraints.at(l.first) << " +";
    }
    for (auto & r : rhs)
        *_imp->proof_stream << " " << _imp->injectivity_constraints.at(r.first) << " +";
    *_imp->proof_stream << " ;\n";
    ++_imp->proof_line;
}

auto Proof::root_propagation_failed() -> void
{
    *_imp->proof_stream << "% root node propagation failed\n";
}

auto Proof::guessing(int depth, const NamedVertex & branch_v, const NamedVertex & val) -> void
{
    *_imp->proof_stream << "% [" << depth << "] guessing " << branch_v.second << "=" << val.second << '\n';
}

auto Proof::propagation_failure(const vector<pair<int, int>> & decisions, const NamedVertex & branch_v, const NamedVertex & val) -> void
{
    *_imp->proof_stream << "% [" << decisions.size() << "] propagation failure on " << branch_v.second << "=" << val.second << '\n';
    *_imp->proof_stream << "@prop" << (_imp->proof_line + 1) << " rup";
    for (auto & [var, val] : decisions)
        *_imp->proof_stream << " 1 ~x" << _imp->variable_mappings[pair{var, val}];
    *_imp->proof_stream << " >= 1 ;\n";
    ++_imp->proof_line;
}

auto Proof::incorrect_guess(const vector<pair<int, int>> & decisions, bool failure) -> void
{
    if (failure)
        *_imp->proof_stream << "% [" << decisions.size() << "] incorrect guess\n";
    else
        *_imp->proof_stream << "% [" << decisions.size() << "] backtracking\n";

    *_imp->proof_stream << "@guess" << (_imp->proof_line + 1) << " rup";
    for (auto & [var, val] : decisions)
        *_imp->proof_stream << " 1 ~x" << _imp->variable_mappings[pair{var, val}];
    *_imp->proof_stream << " >= 1 ;\n";
    ++_imp->proof_line;
}

auto Proof::out_of_guesses(const vector<pair<int, int>> &) -> void
{
}

auto Proof::unit_propagating(const NamedVertex & var, const NamedVertex & val) -> void
{
    *_imp->proof_stream << "% unit propagating " << var.second << "=" << val.second << '\n';
}

auto Proof::start_level(int l) -> void
{
    *_imp->proof_stream << "setlvl " << l << ";\n";
    _imp->largest_level_set = max(_imp->largest_level_set, l);
    _imp->active_level = l;
}

auto Proof::back_up_to_level(int l) -> void
{
    *_imp->proof_stream << "setlvl " << l << ";\n";
    _imp->largest_level_set = max(_imp->largest_level_set, l);
    _imp->active_level = l;
}

auto Proof::forget_level(int l) -> void
{
    if (_imp->largest_level_set >= l)
        *_imp->proof_stream << "wiplvl " << l << ";\n";
}

auto Proof::back_up_to_top() -> void
{
    *_imp->proof_stream << "setlvl " << 0 << ";\n";
    _imp->active_level = 0;
}

auto Proof::post_restart_nogood(const vector<pair<int, int>> & decisions) -> void
{
    *_imp->proof_stream << "% [" << decisions.size() << "] restart nogood\n";
    *_imp->proof_stream << "@nogood" << (_imp->proof_line + 1) << " rup";
    for (auto & [var, val] : decisions)
        *_imp->proof_stream << " 1 ~x" << _imp->variable_mappings[pair{var, val}];
    *_imp->proof_stream << " >= 1 ;\n";
    ++_imp->proof_line;
}

auto Proof::post_solution(const vector<pair<NamedVertex, NamedVertex>> & decisions) -> void
{
    *_imp->proof_stream << "% found solution";
    for (auto & [var, val] : decisions)
        *_imp->proof_stream << " " << var.second << "=" << val.second;
    *_imp->proof_stream << ";\n";

    // Emit the solution-excluding (solx) rule at the top proof level. The
    // blocking constraint it introduces must persist for the rest of the proof
    // (so the solution count stays sound); if we logged it at the current deep
    // search level, the wiplvl that cleans up this subtree on backtrack would
    // delete it, weakening the guarantee.
    if (0 != _imp->active_level)
        *_imp->proof_stream << "setlvl 0;\n";

    *_imp->proof_stream << "solx";
    for (auto & [var, val] : decisions)
        *_imp->proof_stream << " x" << _imp->variable_mappings[pair{var.first, val.first}];
    *_imp->proof_stream << ";\n";
    ++_imp->proof_line;

    if (0 != _imp->active_level)
        *_imp->proof_stream << "setlvl " << _imp->active_level << ";\n";
}

auto Proof::post_solution(const vector<int> & solution) -> void
{
    *_imp->proof_stream << "solx";
    for (auto & v : solution)
        *_imp->proof_stream << " x" << _imp->binary_variable_mappings[v];
    *_imp->proof_stream << ";\n";
    ++_imp->proof_line;
}

auto Proof::new_incumbent(const vector<pair<int, bool>> & solution) -> void
{
    *_imp->proof_stream << "soli";
    for (auto & [v, t] : solution)
        *_imp->proof_stream << " " << (t ? "" : "~") << "x" << _imp->binary_variable_mappings[v];
    for (auto & [v, w] : _imp->zero_in_proof_objectives)
        *_imp->proof_stream << " ~"
                            << "x" << _imp->variable_mappings[pair{v, w}];
    *_imp->proof_stream << ";\n";
    _imp->objective_line = ++_imp->proof_line;
}

auto Proof::new_incumbent(const vector<tuple<NamedVertex, NamedVertex, bool>> & decisions) -> void
{
    *_imp->proof_stream << "soli";
    for (auto & [var, val, t] : decisions)
        *_imp->proof_stream << " " << (t ? "" : "~") << "x" << _imp->variable_mappings[pair{var.first, val.first}];
    *_imp->proof_stream << ";\n";
    _imp->objective_line = ++_imp->proof_line;
}

auto Proof::adjacency_proof_lines() -> AdjacencyProofLines &
{
    return _imp->adjacency;
}

auto Proof::emit_proof_line(const string & line) -> long
{
    *_imp->proof_stream << line << "\n";
    return ++_imp->proof_line;
}

auto Proof::emit_proof_directive(const string & line) -> void
{
    *_imp->proof_stream << line << "\n";
}

auto Proof::current_proof_line() const -> long
{
    return _imp->proof_line;
}

auto Proof::variable_name(int p, int t) const -> const string &
{
    return _imp->variable_mappings.at(pair<long, long>{p, t});
}

auto Proof::has_variable_mapping(int p, int t) const -> bool
{
    return _imp->variable_mappings.contains(pair<long, long>{p, t});
}

auto Proof::is_locally_injective() const -> bool
{
    return _imp->locally_injective;
}

auto Proof::emit_model_constraint(const string & line) -> void
{
    _imp->model_stream << line << "\n";
    ++_imp->nb_constraints;
}

auto Proof::emit_model_comment(const string & line) -> void
{
    _imp->model_stream << line << "\n";
}

auto Proof::injectivity_label(int t) const -> const string &
{
    return _imp->injectivity_constraints.at(t);
}

auto Proof::locally_injective_label(int p, int t) const -> const string &
{
    return _imp->locally_injective_constraints.at(pair<long, long>{p, t});
}

auto Proof::at_most_one_value_label(int p) const -> const string &
{
    return _imp->at_most_one_value_constraints.at(p);
}

auto Proof::cached_proof_line(const string & key) const -> optional<string>
{
    auto it = _imp->cached_proof_lines.find(key);
    if (it == _imp->cached_proof_lines.end())
        return {};
    return it->second;
}

auto Proof::cache_proof_line(const string & key, const string & label) -> void
{
    _imp->cached_proof_lines.emplace(key, label);
}

auto Proof::create_binary_variable(int vertex,
    const function<auto(int)->string> & name) -> void
{
    _imp->binary_variable_mappings.emplace(vertex, name(vertex));
}

auto Proof::create_objective(int n, optional<int> d) -> void
{
    if (d) {
        _imp->model_stream << "* objective\n";
        for (int v = 0; v < n; ++v)
            _imp->model_stream << "1 x" << _imp->binary_variable_mappings[v] << " ";
        _imp->model_stream << ">= " << *d << ";\n";
        _imp->objective_line = ++_imp->nb_constraints;
    }
    else {
        _imp->model_prelude_stream << "min:";
        for (int v = 0; v < n; ++v)
            _imp->model_prelude_stream << " 1 ~x" << _imp->binary_variable_mappings[v];
        _imp->model_prelude_stream << " ;\n";

        for (int v = 0; v < n; ++v)
            _imp->objective_sum << " 1 ~x" << _imp->binary_variable_mappings[v];
    }
}

auto Proof::create_non_edge_constraint(const NamedVertex & p, const NamedVertex & q) -> void
{
    auto name = string{"@noedge"} + (p.first < q.first ? p.second : q.second) + "_" + (p.first < q.first ? q.second : p.second);

    _imp->model_stream << name << " -1 x" << _imp->binary_variable_mappings[p.first] << " -1 x" << _imp->binary_variable_mappings[q.first] << " >= -1 ;\n";

    ++_imp->nb_constraints;
    _imp->non_edge_constraints.emplace(pair{p.first, q.first}, name);
    _imp->non_edge_constraints.emplace(pair{q.first, p.first}, name);
}

auto Proof::create_null_decision_bound(int p, int t, optional<int> d) -> void
{
    if (d) {
        _imp->model_stream << "* objective\n";
        for (int v = 0; v < p; ++v)
            _imp->model_stream << " 1 x" << _imp->variable_mappings[pair{v, t}] << " ";
        _imp->model_stream << ">= " << *d << " ;\n";
        _imp->objective_line = ++_imp->nb_constraints;
    }
    else {
        _imp->model_prelude_stream << "min:";
        for (int v = 0; v < p; ++v)
            _imp->model_prelude_stream << " 1 x" << _imp->variable_mappings[pair{v, t}] << " ";
        _imp->model_prelude_stream << " ;\n";

        for (int v = 0; v < p; ++v)
            _imp->objective_sum << " 1 x" << _imp->variable_mappings[pair{v, t}] << " ";
    }
}

auto Proof::backtrack_from_binary_variables(const vector<int> & v) -> void
{
    if (! _imp->doing_hom_colour_proof) {
        *_imp->proof_stream << "@binback" << (_imp->proof_line + 1) << " rup";
        for (auto & w : v)
            *_imp->proof_stream << " 1 ~x" << _imp->binary_variable_mappings[w];
        *_imp->proof_stream << " >= 1 ;\n";
        ++_imp->proof_line;
    }
    else {
        *_imp->proof_stream << "% backtrack shenanigans, depth " << v.size() << '\n';
        function<auto(unsigned, const vector<pair<int, int>> &)->void> f;
        f = [&](unsigned d, const vector<pair<int, int>> & trail) -> void {
            if (d == v.size()) {
                *_imp->proof_stream << "@binback" << (_imp->proof_line + 1) << " rup 1 ~x" << _imp->variable_mappings[pair{_imp->hom_colour_proof_p.first, _imp->hom_colour_proof_t.first}];
                for (auto & t : trail)
                    *_imp->proof_stream << " 1 ~x" << _imp->variable_mappings[t];
                *_imp->proof_stream << " >= 1 ;\n";
                ++_imp->proof_line;
            }
            else {
                for (auto & p : _imp->p_clique) {
                    vector<pair<int, int>> new_trail{trail};
                    new_trail.emplace_back(pair{p.first, _imp->t_clique_neighbourhood.find(v[d])->second.first});
                    f(d + 1, new_trail);
                }
            }
        };
        f(0, {});
    }
}

auto Proof::colour_bound(const vector<vector<int>> & ccs) -> void
{
    *_imp->proof_stream << "% bound, ccs";
    for (auto & cc : ccs) {
        *_imp->proof_stream << " [";
        for (auto & c : cc)
            *_imp->proof_stream << " " << c;
        *_imp->proof_stream << " ]";
    }
    *_imp->proof_stream << '\n';

    auto how_many_summands = 0u;
    auto do_one_cc = [&](const auto & cc, const auto & non_edge_constraint) {
        if (cc.size() > 2) {
            *_imp->proof_stream << " " << non_edge_constraint(cc[0], cc[1]);

            for (unsigned i = 2; i < cc.size(); ++i) {
                *_imp->proof_stream << " " << i << " *";
                for (unsigned j = 0; j < i; ++j)
                    *_imp->proof_stream << " " << non_edge_constraint(cc[i], cc[j]) << " +";
                *_imp->proof_stream << " " << (i + 1) << " d";
            }

            ++how_many_summands;
        }
        else if (cc.size() == 2) {
            *_imp->proof_stream << " " << non_edge_constraint(cc[0], cc[1]);
            ++how_many_summands;
        }
    };

    *_imp->proof_stream << "@colpol" << (_imp->proof_line + 1) << " pol ";

    for (auto & cc : ccs) {
        if (cc.size() == 1)
            continue;

        if (_imp->doing_hom_colour_proof) {
            vector<pair<NamedVertex, NamedVertex>> bigger_cc;
            for (auto & c : cc)
                for (auto & v : _imp->p_clique)
                    bigger_cc.push_back(pair{v, _imp->t_clique_neighbourhood.find(c)->second});

            do_one_cc(bigger_cc, [&](const pair<NamedVertex, NamedVertex> & a, const pair<NamedVertex, NamedVertex> & b) -> long {
                return _imp->clique_for_hom_non_edge_constraints[pair{a, b}];
            });
        }
        else
            do_one_cc(cc, [&](int a, int b) -> string { return _imp->non_edge_constraints[pair{a, b}]; });
    }

    *_imp->proof_stream << " " << _imp->objective_line;
    for (unsigned n = 0; n < how_many_summands; ++n)
        *_imp->proof_stream << " +";
    *_imp->proof_stream << ";\n";
    ++_imp->proof_line;
}

auto Proof::prepare_hom_clique_proof(const NamedVertex & p, const NamedVertex & t, unsigned size) -> void
{
    *_imp->proof_stream << "% clique of size " << size << " around neighbourhood of " << p.second << " but not " << t.second << '\n';
    *_imp->proof_stream << "setlvl 1;\n";
    _imp->doing_hom_colour_proof = true;
    _imp->hom_colour_proof_p = p;
    _imp->hom_colour_proof_t = t;
}

auto Proof::start_hom_clique_proof(const NamedVertex & p, vector<NamedVertex> && p_clique, const NamedVertex & t, map<int, NamedVertex> && t_clique_neighbourhood) -> void
{
    _imp->p_clique = move(p_clique);
    _imp->t_clique_neighbourhood = move(t_clique_neighbourhood);

    *_imp->proof_stream << "% hom clique objective\n";
    vector<long> to_sum;
    for (auto & q : _imp->p_clique) {
        *_imp->proof_stream << "@hombd" << q.second << "_" << t.second << "_" << (_imp->proof_line + 1) << " rup 1 ~x" << _imp->variable_mappings[pair{p.first, t.first}];
        for (auto & u : _imp->t_clique_neighbourhood)
            *_imp->proof_stream << " 1 x" << _imp->variable_mappings[pair{q.first, u.second.first}];
        *_imp->proof_stream << " >= 1 ;\n";
        to_sum.push_back(++_imp->proof_line);
    }

    *_imp->proof_stream << "@hompol" << (_imp->proof_line + 1) << " pol";
    bool first = true;
    for (auto & t : to_sum) {
        *_imp->proof_stream << " " << t;
        if (! first)
            *_imp->proof_stream << " +";
        first = false;
    }
    *_imp->proof_stream << ";\n";
    _imp->objective_line = ++_imp->proof_line;

    *_imp->proof_stream << "% hom clique non edges for injectivity\n";

    for (auto & p : _imp->p_clique)
        for (auto & q : _imp->p_clique)
            if (p != q) {
                for (auto & [_, t] : _imp->t_clique_neighbourhood) {
                    *_imp->proof_stream << "@hominj" << p.second << "_" << q.second << "_" << t.second << "_" << (_imp->proof_line + 1) << " rup 1 ~x" << _imp->variable_mappings[pair{p.first, t.first}] << " 1 ~x" << _imp->variable_mappings[pair{q.first, t.first}] << " >= 1 ;\n";
                    ++_imp->proof_line;
                    _imp->clique_for_hom_non_edge_constraints.emplace(pair{pair{p, t}, pair{q, t}}, _imp->proof_line);
                    _imp->clique_for_hom_non_edge_constraints.emplace(pair{pair{q, t}, pair{p, t}}, _imp->proof_line);
                }
            }

    *_imp->proof_stream << "% hom clique non edges for variables\n";

    for (auto & p : _imp->p_clique)
        for (auto & [_, t] : _imp->t_clique_neighbourhood) {
            for (auto & [_, u] : _imp->t_clique_neighbourhood) {
                if (t != u) {
                    *_imp->proof_stream << "@homdom" << p.second << "_" << t.second << "_" << u.second << "_" << (_imp->proof_line + 1) << " rup 1 ~x" << _imp->variable_mappings[pair{p.first, t.first}] << " 1 ~x" << _imp->variable_mappings[pair{p.first, u.first}] << " >= 1 ;\n";
                    ++_imp->proof_line;
                    _imp->clique_for_hom_non_edge_constraints.emplace(pair{pair{p, t}, pair{p, u}}, _imp->proof_line);
                    _imp->clique_for_hom_non_edge_constraints.emplace(pair{pair{p, u}, pair{p, t}}, _imp->proof_line);
                }
            }
        }
}

auto Proof::finish_hom_clique_proof(const NamedVertex & p, const NamedVertex & t, unsigned size) -> void
{
    *_imp->proof_stream << "% end clique of size " << size << " around neighbourhood of " << p.second << " but not " << t.second << '\n';
    *_imp->proof_stream << "setlvl 0;\n";
    *_imp->proof_stream << "@homfin" << p.second << "_" << t.second << " rup 1 ~x" << _imp->variable_mappings[pair{p.first, t.first}] << " >= 1 ;\n";
    *_imp->proof_stream << "wiplvl 1;\n";
    ++_imp->proof_line;
    _imp->doing_hom_colour_proof = false;
    _imp->clique_for_hom_non_edge_constraints.clear();
}

auto Proof::add_hom_clique_non_edge(
    const NamedVertex & pp,
    const NamedVertex & tt,
    const vector<NamedVertex> & p_clique,
    const NamedVertex & t,
    const NamedVertex & u) -> void
{
    *_imp->proof_stream << "% hom clique non edges for " << t.second << " " << u.second << '\n';
    for (auto & p : p_clique) {
        for (auto & q : p_clique) {
            if (p != q) {
                *_imp->proof_stream << "@homcross" << (_imp->proof_line + 1) << " rup 1 ~x" << _imp->variable_mappings[pair{pp.first, tt.first}]
                                    << " 1 ~x" << _imp->variable_mappings[pair{p.first, t.first}]
                                    << " 1 ~x" << _imp->variable_mappings[pair{q.first, u.first}] << " >= 1 ;\n";
                ++_imp->proof_line;
                _imp->clique_for_hom_non_edge_constraints.emplace(pair{pair{p, t}, pair{q, u}}, _imp->proof_line);
                _imp->clique_for_hom_non_edge_constraints.emplace(pair{pair{q, u}, pair{p, t}}, _imp->proof_line);
            }
        }
    }
}

auto Proof::mcs_bound(
    const vector<pair<set<int>, set<int>>> & partitions) -> void
{
    *_imp->proof_stream << "% failed bound\n";

    vector<string> to_sum;
    for (auto & [l, r] : partitions) {
        if (r.size() >= l.size())
            continue;

        *_imp->proof_stream << "@mcspart" << (_imp->proof_line + 1) << " pol";
        bool first = true;
        for (auto & v : l) {
            *_imp->proof_stream << " " << _imp->at_least_one_value_constraints[v];
            if (first)
                first = false;
            else
                *_imp->proof_stream << " +";
        }
        for (auto & v : r)
            *_imp->proof_stream << " " << _imp->injectivity_constraints[v] << " +";

        *_imp->proof_stream << ";\n";
        to_sum.push_back(to_string(++_imp->proof_line));
    }

    if (! to_sum.empty()) {
        *_imp->proof_stream << "@mcsfin" << (_imp->proof_line + 1) << " pol " << _imp->objective_line;
        for (auto & t : to_sum)
            *_imp->proof_stream << " " << t << " +";
        *_imp->proof_stream << ";\n";
        ++_imp->proof_line;
    }
}

auto Proof::create_connected_constraints(int p, int t, const function<auto(int, int)->bool> & adj) -> void
{
    _imp->model_stream << "* selected vertices must be connected, walk 1\n";
    int mapped_to_null = t;

    for (int v = 0; v < p; ++v)
        for (int w = 0; w < v; ++w) {
            string n = "conn1_" + to_string(v) + "_" + to_string(w);
            _imp->connected_variable_mappings.emplace(tuple{1, v, w}, n);
            if (! adj(v, w)) {
                // v not adjacent to w, so the walk does not exist
                _imp->model_stream << "1 ~x" << n << " >= 1 ;\n";
                ++_imp->nb_constraints;
            }
            else {
                // v = null -> the walk does not exist
                _imp->model_stream << "1 ~x" << n << " 1 ~x" << _imp->variable_mappings[pair{v, mapped_to_null}] << " >= 1 ;\n";
                // w = null -> the walk does not exist
                _imp->model_stream << "1 ~x" << n << " 1 ~x" << _imp->variable_mappings[pair{w, mapped_to_null}] << " >= 1 ;\n";
                // either v = null, or w = null, or the walk exists
                _imp->model_stream << "1 x" << n << " 1 x" << _imp->variable_mappings[pair{v, mapped_to_null}]
                                   << " 1 x" << _imp->variable_mappings[pair{w, mapped_to_null}] << " >= 1 ;\n";
                _imp->nb_constraints += 3;
            }
        }

    int last_k = 0;
    for (int k = 2 ; ; k *= 2) {
        last_k = k;
        _imp->model_stream << "* selected vertices must be connected, walk " << k << '\n';
        for (int v = 0; v < p; ++v)
            for (int w = 0; w < v; ++w) {
                string n = "conn" + to_string(k) + "_" + to_string(v) + "_" + to_string(w);
                _imp->connected_variable_mappings.emplace(tuple{k, v, w}, n);

                vector<string> ors;
                for (int u = 0; u < p; ++u) {
                    if (v != w && v != u && u != w) {
                        string m = "conn" + to_string(k) + "_" + to_string(v) + "_" + to_string(w) + "_via_" + to_string(u);
                        ors.push_back(m);
                        _imp->connected_variable_mappings_aux.emplace(tuple{k, v, w, u}, m);
                        // either the first half walk exists, or the via term is false
                        _imp->model_stream << "1 x" << _imp->connected_variable_mappings[tuple{k / 2, max(u, v), min(u, v)}]
                                           << " 1 ~x" << m << " >= 1 ;\n";
                        // either the second half walk exists, or the via term is false
                        _imp->model_stream << "1 x" << _imp->connected_variable_mappings[tuple{k / 2, max(u, w), min(u, w)}]
                                           << " 1 ~x" << m << " >= 1 ;\n";
                        // one of the half walks is false, or the via term must be true
                        _imp->model_stream << "1 x" << m
                                           << " 1 ~x" << _imp->connected_variable_mappings[tuple{k / 2, max(v, u), min(v, u)}]
                                           << " 1 ~x" << _imp->connected_variable_mappings[tuple{k / 2, max(u, w), min(u, w)}] << " >= 1 ;\n";
                        _imp->nb_constraints += 3;
                    }
                }

                // one of the vias must be true, or a shorter walk exists, or the entry is false
                _imp->model_stream << "1 ~x" << n;
                for (auto & o : ors)
                    _imp->model_stream << " 1 x" << o;
                _imp->model_stream << " 1 x" << _imp->connected_variable_mappings[tuple{k / 2, v, w}];
                _imp->model_stream << " >= 1 ;\n";
                ++_imp->nb_constraints;

                // if the entry is false, then all of the vias must be false and the shorter walk must be false
                for (auto & o : ors) {
                    _imp->model_stream << "1 x" << n << " 1 ~x" << o << " >= 1 ;\n";
                    ++_imp->nb_constraints;
                }
                _imp->model_stream << "1 x" << n << " 1 ~x" << _imp->connected_variable_mappings[tuple{k / 2, v, w}] << " >= 1 ;\n";
                ++_imp->nb_constraints;
            }

        if (k >= min(p, t))
            break;
    }

    _imp->model_stream << "* if two vertices are used, they must be connected\n";
    for (int v = 0; v < p; ++v)
        for (int w = 0; w < v; ++w) {
            _imp->model_stream << "1 x" << _imp->variable_mappings[pair{v, mapped_to_null}]
                               << " 1 x" << _imp->variable_mappings[pair{w, mapped_to_null}];
            auto var = _imp->connected_variable_mappings.find(tuple{last_k, v, w});
            if (var != _imp->connected_variable_mappings.end())
                _imp->model_stream << " 1 x" << _imp->connected_variable_mappings[tuple{last_k, v, w}];
            _imp->model_stream << " >= 1 ;\n";
            ++_imp->nb_constraints;
        }
}

auto Proof::not_connected_in_underlying_graph(const vector<int> & x, int y) -> void
{
    *_imp->proof_stream << "@notconn" << y << "_" << (_imp->proof_line + 1) << " rup 1 ~x" << _imp->binary_variable_mappings[y];
    for (auto & v : x)
        *_imp->proof_stream << " 1 ~x" << _imp->binary_variable_mappings[v];
    *_imp->proof_stream << " >= 1 ;\n";
    ++_imp->proof_line;
}

auto Proof::has_clique_model() const -> bool
{
    return _imp->clique_encoding;
}

auto Proof::create_clique_encoding(
    const vector<pair<int, int>> & enc,
    const vector<pair<int, int>> & zero_in_proof_objectives) -> void
{
    _imp->clique_encoding = true;
    for (unsigned i = 0; i < enc.size(); ++i)
        _imp->binary_variable_mappings.emplace(i, _imp->variable_mappings[enc[i]]);

    _imp->zero_in_proof_objectives = zero_in_proof_objectives;
    _imp->doing_mcs_by_clique = true;
}

auto Proof::create_clique_nonedge(int v, int w) -> void
{
    *_imp->proof_stream << "@cliqedge" << min(v, w) << "_" << max(v, w) << " rup 1 ~x" << _imp->binary_variable_mappings[v]
                        << " 1 ~x" << _imp->binary_variable_mappings[w] << " >= 1 ;\n";
    ++_imp->proof_line;
    _imp->non_edge_constraints.emplace(pair{v, w}, to_string(_imp->proof_line));
    _imp->non_edge_constraints.emplace(pair{w, v}, to_string(_imp->proof_line));
}

auto Proof::super_extra_verbose() const -> bool
{
    return _imp->super_extra_verbose;
}

auto Proof::show_domains(const string & s, const vector<pair<NamedVertex, vector<NamedVertex>>> & domains) -> void
{
    *_imp->proof_stream << "% " << s << ", domains follow\n";
    for (auto & [p, ts] : domains) {
        *_imp->proof_stream << "%    " << p.second << " size " << ts.size() << " = {";
        for (auto & t : ts)
            *_imp->proof_stream << " " << t.second;
        *_imp->proof_stream << " }\n";
    }
}

auto Proof::propagated(const NamedVertex & p, const NamedVertex & t, int g, int n_values, const NamedVertex & q) -> void
{
    *_imp->proof_stream << "% adjacency propagation from " << p.second << " -> " << t.second << " in graph pairs " << g << " deleted " << n_values << " values from " << q.second << '\n';
}
