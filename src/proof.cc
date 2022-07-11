/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "proof.hh"

#include <algorithm>
#include <iterator>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <tuple>
#include <utility>

#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/stream.hpp>

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
using std::vector;

using boost::iostreams::bzip2_compressor;
using boost::iostreams::file_sink;
using boost::iostreams::filtering_ostream;

namespace
{
    auto make_compressed_ostream(const string & fn) -> unique_ptr<ostream>
    {
        auto out = make_unique<filtering_ostream>();
        out->push(bzip2_compressor());
        out->push(file_sink(fn));
        return out;
    }
}

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
    bool friendly_names;
    bool bz2 = false;
    bool super_extra_verbose = false;

    map<pair<long, long>, string> variable_mappings;
    map<long, string> binary_variable_mappings;
    map<tuple<long, long, long>, string> connected_variable_mappings;
    map<tuple<long, long, long, long>, string> connected_variable_mappings_aux;
    map<long, long> at_least_one_value_constraints, at_most_one_value_constraints, injectivity_constraints;
    map<tuple<long, long, long, long>, long> adjacency_lines;
    map<pair<long, long>, long> eliminations;
    map<pair<long, long>, long> non_edge_constraints;
    long objective_line = 0;

    long nb_constraints = 0;
    long proof_line = 0;
    int largest_level_set = 0;

    bool clique_encoding = false;

    bool doing_hom_colour_proof = false;
    NamedVertex hom_colour_proof_p, hom_colour_proof_t;
    vector<NamedVertex> p_clique;
    map<int, NamedVertex> t_clique_neighbourhood;
    map<pair<pair<NamedVertex, NamedVertex>, pair<NamedVertex, NamedVertex> >, long> clique_for_hom_non_edge_constraints;

    vector<pair<int, int> > zero_in_proof_objectives;
};

Proof::Proof(const string & opb_file, const string & log_file, bool f, bool b, bool s) :
    _imp(new Imp)
{
    _imp->opb_filename = opb_file;
    _imp->log_filename = log_file;
    _imp->friendly_names = f;
    _imp->bz2 = b;
    _imp->super_extra_verbose = s;
}

Proof::Proof(Proof &&) = default;

Proof::~Proof() = default;

auto Proof::operator= (Proof &&) -> Proof & = default;

auto Proof::create_cp_variable(int pattern_vertex, int target_size,
        const function<auto (int) -> string> & pattern_name,
        const function<auto (int) -> string> & target_name) -> void
{
    for (int i = 0 ; i < target_size ; ++i)
        if (_imp->friendly_names)
            _imp->variable_mappings.emplace(pair{ pattern_vertex, i }, pattern_name(pattern_vertex) + "_" + target_name(i));
        else
            _imp->variable_mappings.emplace(pair{ pattern_vertex, i }, to_string(_imp->variable_mappings.size() + 1));

    _imp->model_stream << "* vertex " << pattern_vertex << " domain" << endl;
    for (int i = 0 ; i < target_size ; ++i)
        _imp->model_stream << "1 x" << _imp->variable_mappings[{ pattern_vertex, i }] << " ";
    _imp->model_stream << ">= 1 ;" << endl;
    _imp->at_least_one_value_constraints.emplace(pattern_vertex, ++_imp->nb_constraints);

    for (int i = 0 ; i < target_size ; ++i)
        _imp->model_stream << "-1 x" << _imp->variable_mappings[{ pattern_vertex, i }] << " ";
    _imp->model_stream << ">= -1 ;" << endl;
    _imp->at_most_one_value_constraints.emplace(pattern_vertex, ++_imp->nb_constraints);
}

auto Proof::create_injectivity_constraints(int pattern_size, int target_size) -> void
{
    for (int v = 0 ; v < target_size ; ++v) {
        _imp->model_stream << "* injectivity on value " << v << endl;

        for (int p = 0 ; p < pattern_size ; ++p) {
            auto x = _imp->variable_mappings.find(pair{ p, v });
            if (x != _imp->variable_mappings.end())
                _imp->model_stream << "-1 x" << x->second << " ";
        }
        _imp->model_stream << ">= -1 ;" << endl;
        _imp->injectivity_constraints.emplace(v, ++_imp->nb_constraints);
    }
}

auto Proof::create_forbidden_assignment_constraint(int p, int t) -> void
{
    _imp->model_stream << "* forbidden assignment" << endl;
    _imp->model_stream << "1 ~x" << _imp->variable_mappings[pair{ p, t }] << " >= 1 ;" << endl;
    ++_imp->nb_constraints;
    _imp->eliminations.emplace(pair{ p, t }, _imp->nb_constraints);
}

auto Proof::start_adjacency_constraints_for(int p, int t) -> void
{
    _imp->model_stream << "* adjacency " << p << " maps to " << t << endl;
}

auto Proof::create_adjacency_constraint(int p, int q, int t, const vector<int> & uu, bool) -> void
{
    _imp->model_stream << "1 ~x" << _imp->variable_mappings[pair{ p, t }];
    for (auto & u : uu)
        _imp->model_stream << " 1 x" << _imp->variable_mappings[pair{ q, u }];
    _imp->model_stream << " >= 1 ;" << endl;
    _imp->adjacency_lines.emplace(tuple{ 0, p, q, t }, ++_imp->nb_constraints);
}

auto Proof::finalise_model() -> void
{
    unique_ptr<ostream> f = (_imp->bz2 ? make_compressed_ostream(_imp->opb_filename + ".bz2") : make_unique<ofstream>(_imp->opb_filename));

    *f << "* #variable= " << (_imp->variable_mappings.size() + _imp->binary_variable_mappings.size()
            + _imp->connected_variable_mappings.size() + _imp->connected_variable_mappings_aux.size())
        << " #constraint= " << _imp->nb_constraints << endl;
    copy(istreambuf_iterator<char>{ _imp->model_prelude_stream }, istreambuf_iterator<char>{}, ostreambuf_iterator<char>{ *f });
    _imp->model_prelude_stream.clear();
    copy(istreambuf_iterator<char>{ _imp->model_stream }, istreambuf_iterator<char>{}, ostreambuf_iterator<char>{ *f });
    _imp->model_stream.clear();

    if (! *f)
        throw ProofError{ "Error writing opb file to '" + _imp->opb_filename + "'" };

    _imp->proof_stream = (_imp->bz2 ? make_compressed_ostream(_imp->log_filename + ".bz2") : make_unique<ofstream>(_imp->log_filename));

    *_imp->proof_stream << "pseudo-Boolean proof version 1.0" << endl;

    *_imp->proof_stream << "f " << _imp->nb_constraints << " 0" << endl;
    _imp->proof_line += _imp->nb_constraints;

    if (! *_imp->proof_stream)
        throw ProofError{ "Error writing proof file to '" + _imp->log_filename + "'" };
}

auto Proof::finish_unsat_proof() -> void
{
    *_imp->proof_stream << "* asserting that we've proved unsat" << endl;
    *_imp->proof_stream << "u >= 1 ;" << endl;
    ++_imp->proof_line;
    *_imp->proof_stream << "c " << _imp->proof_line << " 0" << endl;
}

auto Proof::failure_due_to_pattern_bigger_than_target() -> void
{
    *_imp->proof_stream << "* failure due to the pattern being bigger than the target" << endl;

    // we get a hall violator by adding up all of the things
    *_imp->proof_stream << "p";
    bool first = true;

    for (auto & [ _, line ] : _imp->at_least_one_value_constraints) {
        if (first) {
            *_imp->proof_stream << " " << line;
            first = false;
        }
        else
            *_imp->proof_stream << " " << line << " +";
    }

    for (auto & [ _, line ] : _imp->injectivity_constraints)
        *_imp->proof_stream << " " << line << " +";
    *_imp->proof_stream << " 0" << endl;
    ++_imp->proof_line;
}

auto Proof::incompatible_by_degrees(
        int g,
        const NamedVertex & p,
        const vector<int> & n_p,
        const NamedVertex & t,
        const vector<int> & n_t) -> void
{
    *_imp->proof_stream << "* cannot map " << p.second << " to " << t.second << " due to degrees in graph pairs " << g << endl;

    *_imp->proof_stream << "p";
    bool first = true;
    for (auto & n : n_p) {
        // due to loops or labels, it might not be possible to map n to t.first
        if (_imp->adjacency_lines.count(tuple{ g, p.first, n, t.first })) {
            if (first) {
                first = false;
                *_imp->proof_stream << " " << _imp->adjacency_lines[tuple{ g, p.first, n, t.first }];
            }
            else
                *_imp->proof_stream << " " << _imp->adjacency_lines[tuple{ g, p.first, n, t.first }] << " +";
        }
    }

    // if I map p to t, I have to map the neighbours of p to neighbours of t
    for (auto & n : n_t)
        *_imp->proof_stream << " " << _imp->injectivity_constraints[n] << " +";

    *_imp->proof_stream << " 0" << endl;
    ++_imp->proof_line;

    *_imp->proof_stream << "j " << _imp->proof_line << " 1 ~x" << _imp->variable_mappings[pair{ p.first, t.first }] << " >= 1 ;" << endl;
    ++_imp->proof_line;
    _imp->eliminations.emplace(pair{ p.first, t.first }, _imp->proof_line);

    *_imp->proof_stream << "d " << _imp->proof_line - 1 << " 0" << endl;
}

auto Proof::incompatible_by_nds(
        int g,
        const NamedVertex & p,
        const NamedVertex & t,
        const vector<int> & p_subsequence,
        const vector<int> & t_subsequence,
        const vector<int> & t_remaining) -> void
{
    *_imp->proof_stream << "* cannot map " << p.second << " to " << t.second << " due to nds in graph pairs " << g << endl;

    // summing up horizontally
    *_imp->proof_stream << "p";
    bool first = true;
    for (auto & n : p_subsequence) {
        // due to loops or labels, it might not be possible to map n to t.first
        if (_imp->adjacency_lines.count(tuple{ g, p.first, n, t.first })) {
            if (first) {
                first = false;
                *_imp->proof_stream << " " << _imp->adjacency_lines[tuple{ g, p.first, n, t.first }];
            }
            else
                *_imp->proof_stream << " " << _imp->adjacency_lines[tuple{ g, p.first, n, t.first }] << " +";
        }
    }

    // injectivity in the square
    for (auto & t : t_subsequence) {
        if (t != t_subsequence.back())
            *_imp->proof_stream << " " << _imp->injectivity_constraints.find(t)->second << " +";
    }

    // block to the right of the failing square
    for (auto & n : p_subsequence) {
        for (auto & u : t_remaining) {
            /* n -> t is already eliminated by degree or loop */
            *_imp->proof_stream << " " << _imp->eliminations[pair{ n, u }] << " +";
        }
    }

    // final column
    for (auto & n : p_subsequence) {
        /* n -> t is already eliminated by degree or loop */
        *_imp->proof_stream << " " << _imp->eliminations[pair{ n, t_subsequence.back() }] << " +";
    }

    *_imp->proof_stream << " 0" << endl;
    ++_imp->proof_line;

    *_imp->proof_stream << "j " << _imp->proof_line << " 1 ~x" << _imp->variable_mappings[pair{ p.first, t.first }] << " >= 1 ;" << endl;
    ++_imp->proof_line;

    *_imp->proof_stream << "d " << _imp->proof_line - 1 << " 0" << endl;
}

auto Proof::initial_domain_is_empty(int p) -> void
{
    *_imp->proof_stream << "* failure due to domain " << p << " being empty" << endl;
}

auto Proof::emit_hall_set_or_violator(const vector<NamedVertex> & lhs, const vector<NamedVertex> & rhs) -> void
{
    *_imp->proof_stream << "* hall set or violator {";
    for (auto & l : lhs)
        *_imp->proof_stream << " " << l.second;
    *_imp->proof_stream << " } / {";
    for (auto & r : rhs)
        *_imp->proof_stream << " " << r.second;
    *_imp->proof_stream << " }" << endl;
    *_imp->proof_stream << "p";
    bool first = true;
    for (auto & l : lhs) {
        if (first) {
            first = false;
            *_imp->proof_stream << " " << _imp->at_least_one_value_constraints[l.first];
        }
        else
            *_imp->proof_stream << " " << _imp->at_least_one_value_constraints[l.first] << " +";
    }
    for (auto & r : rhs)
        *_imp->proof_stream << " " << _imp->injectivity_constraints[r.first] << " +";
    *_imp->proof_stream << " 0" << endl;
    ++_imp->proof_line;
}

auto Proof::root_propagation_failed() -> void
{
    *_imp->proof_stream << "* root node propagation failed" << endl;
}

auto Proof::guessing(int depth, const NamedVertex & branch_v, const NamedVertex & val) -> void
{
    *_imp->proof_stream << "* [" << depth << "] guessing " << branch_v.second << "=" << val.second << endl;
}

auto Proof::propagation_failure(const vector<pair<int, int> > & decisions, const NamedVertex & branch_v, const NamedVertex & val) -> void
{
    *_imp->proof_stream << "* [" << decisions.size() << "] propagation failure on " << branch_v.second << "=" << val.second << endl;
    *_imp->proof_stream << "u ";
    for (auto & [ var, val ] : decisions)
        *_imp->proof_stream << " 1 ~x" << _imp->variable_mappings[pair{ var, val }];
    *_imp->proof_stream << " >= 1 ;" << endl;
    ++_imp->proof_line;
}

auto Proof::incorrect_guess(const vector<pair<int, int> > & decisions, bool failure) -> void
{
    if (failure)
        *_imp->proof_stream << "* [" << decisions.size() << "] incorrect guess" << endl;
    else
        *_imp->proof_stream << "* [" << decisions.size() << "] backtracking" << endl;

    *_imp->proof_stream << "u";
    for (auto & [ var, val ] : decisions)
        *_imp->proof_stream << " 1 ~x" << _imp->variable_mappings[pair{ var, val }];
    *_imp->proof_stream << " >= 1 ;" << endl;
    ++_imp->proof_line;
}

auto Proof::out_of_guesses(const vector<pair<int, int> > &) -> void
{
}

auto Proof::unit_propagating(const NamedVertex & var, const NamedVertex & val) -> void
{
    *_imp->proof_stream << "* unit propagating " << var.second << "=" << val.second << endl;
}

auto Proof::start_level(int l) -> void
{
    *_imp->proof_stream << "# " << l << endl;
    _imp->largest_level_set = max(_imp->largest_level_set, l);
}

auto Proof::back_up_to_level(int l) -> void
{
    *_imp->proof_stream << "# " << l << endl;
    _imp->largest_level_set = max(_imp->largest_level_set, l);
}

auto Proof::forget_level(int l) -> void
{
    if (_imp->largest_level_set >= l)
        *_imp->proof_stream << "w " << l << endl;
}

auto Proof::back_up_to_top() -> void
{
    *_imp->proof_stream << "# " << 0 << endl;
}

auto Proof::post_restart_nogood(const vector<pair<int, int> > & decisions) -> void
{
    *_imp->proof_stream << "* [" << decisions.size() << "] restart nogood" << endl;
    *_imp->proof_stream << "u";
    for (auto & [ var, val ] : decisions)
        *_imp->proof_stream << " 1 ~x" << _imp->variable_mappings[pair{ var, val }];
    *_imp->proof_stream << " >= 1 ;" << endl;
    ++_imp->proof_line;
}

auto Proof::post_solution(const vector<pair<NamedVertex, NamedVertex> > & decisions) -> void
{
    *_imp->proof_stream << "* found solution";
    for (auto & [ var, val ] : decisions)
        *_imp->proof_stream << " " << var.second << "=" << val.second;
    *_imp->proof_stream << endl;

    *_imp->proof_stream << "v";
    for (auto & [ var, val ] : decisions)
        *_imp->proof_stream << " x" << _imp->variable_mappings[pair{ var.first, val.first }];
    *_imp->proof_stream << endl;
    ++_imp->proof_line;
}

auto Proof::post_solution(const vector<int> & solution) -> void
{
    *_imp->proof_stream << "v";
    for (auto & v : solution)
        *_imp->proof_stream << " x" << _imp->binary_variable_mappings[v];
    *_imp->proof_stream << endl;
    ++_imp->proof_line;
}

auto Proof::new_incumbent(const vector<pair<int, bool> > & solution) -> void
{
    *_imp->proof_stream << "o";
    for (auto & [ v, t ] : solution)
        *_imp->proof_stream << " " << (t ? "" : "~") << "x" << _imp->binary_variable_mappings[v];
    for (auto & [ v, w ] : _imp->zero_in_proof_objectives)
        *_imp->proof_stream << " ~" << "x" << _imp->variable_mappings[pair{ v, w }];
    *_imp->proof_stream << endl;
    _imp->objective_line = ++_imp->proof_line;
}

auto Proof::new_incumbent(const vector<tuple<NamedVertex, NamedVertex, bool> > & decisions) -> void
{
    *_imp->proof_stream << "o";
    for (auto & [ var, val, t ] : decisions)
        *_imp->proof_stream << " " << (t ? "" : "~") << "x" << _imp->variable_mappings[pair{ var.first, val.first }];
    *_imp->proof_stream << endl;
    _imp->objective_line = ++_imp->proof_line;
}

auto Proof::create_exact_path_graphs(
        int g,
        const NamedVertex & p,
        const NamedVertex & q,
        const vector<NamedVertex> & between_p_and_q,
        const NamedVertex & t,
        const vector<NamedVertex> & n_t,
        const vector<pair<NamedVertex, vector<NamedVertex> > > & two_away_from_t,
        const vector<NamedVertex> & d_n_t
        ) -> void
{
    *_imp->proof_stream << "* adjacency " << p.second << " maps to " << t.second <<
        " in G^[" << g << "x2] so " << q.second << " maps to one of..." << endl;

    *_imp->proof_stream << "# 1" << endl;

    *_imp->proof_stream << "p";

    // if p maps to t then things in between_p_and_q have to go to one of these...
    bool first = true;
    for (auto & b : between_p_and_q) {
        *_imp->proof_stream << " " << _imp->adjacency_lines[tuple{ 0, p.first, b.first, t.first }];
        if (! first)
            *_imp->proof_stream << " +";
        first = false;
    }

    // now go two hops out: cancel between_p_and_q things with where can q go
    for (auto & b : between_p_and_q) {
        for (auto & w : n_t) {
            // due to loops or labels, it might not be possible to map to w
            auto i = _imp->adjacency_lines.find(tuple{ 0, b.first, q.first, w.first });
            if (i != _imp->adjacency_lines.end())
                *_imp->proof_stream << " " << i->second << " +";
        }
    }

    *_imp->proof_stream << " 0" << endl;
    ++_imp->proof_line;

    // first tidy-up step: if p maps to t then q maps to something a two-walk away from t
    *_imp->proof_stream << "j " << _imp->proof_line << " 1 ~x" << _imp->variable_mappings[pair{ p.first, t.first }];
    for (auto & u : two_away_from_t)
        *_imp->proof_stream << " 1 x" << _imp->variable_mappings[pair{ q.first, u.first.first }];
    *_imp->proof_stream << " >= 1 ;" << endl;
    ++_imp->proof_line;

    // if p maps to t then q does not map to t
    *_imp->proof_stream << "p " << _imp->proof_line << " " << _imp->injectivity_constraints[t.first] << " + 0" << endl;
    ++_imp->proof_line;

    // and cancel out stray extras from injectivity
    *_imp->proof_stream << "j " << _imp->proof_line << " 1 ~x" << _imp->variable_mappings[pair{ p.first, t.first }];
    for (auto & u : two_away_from_t)
        if (u.first != t)
            *_imp->proof_stream << " 1 x" << _imp->variable_mappings[pair{ q.first, u.first.first }];
    *_imp->proof_stream << " >= 1 ;" << endl;
    ++_imp->proof_line;

    vector<long> things_to_add_up;
    things_to_add_up.push_back(_imp->proof_line);

    // cancel out anything that is two away from t, but by insufficiently many paths
    for (auto & u : two_away_from_t) {
        if ((u.first == t) || (d_n_t.end() != find(d_n_t.begin(), d_n_t.end(), u.first)))
            continue;

        *_imp->proof_stream << "p";
        bool first = true;
        for (auto & b : between_p_and_q) {
            *_imp->proof_stream << " " << _imp->adjacency_lines[tuple{ 0, p.first, b.first, t.first }];
            if (! first)
                *_imp->proof_stream << " +";
            first = false;
            *_imp->proof_stream << " " << _imp->adjacency_lines[tuple{ 0, q.first, b.first, u.first.first }] << " +";
            *_imp->proof_stream << " " << _imp->at_most_one_value_constraints[b.first] << " +";
        }

        for (auto & z : u.second)
            *_imp->proof_stream << " " << _imp->injectivity_constraints[z.first] << " +";

        *_imp->proof_stream << " 0" << endl;
        ++_imp->proof_line;

        // want: ~x_p_t + ~x_q_u >= 1
        *_imp->proof_stream << "j " << _imp->proof_line << " 1 ~x" << _imp->variable_mappings[pair{ p.first, t.first }]
            << " 1 ~x" << _imp->variable_mappings[pair{ q.first, u.first.first }] << " >= 1 ;" << endl;
        things_to_add_up.push_back(++_imp->proof_line);
    }

    // do the getting rid of
    if (things_to_add_up.size() > 1) {
        bool first = true;
        *_imp->proof_stream << "p";
        for (auto & t : things_to_add_up) {
            *_imp->proof_stream << " " << t;
            if (! first)
                *_imp->proof_stream << " +";
            first = false;
        }
        *_imp->proof_stream << " 0" << endl;
        ++_imp->proof_line;
    }

    *_imp->proof_stream << "# 0" << endl;

    // and finally, tidy up to get what we wanted
    *_imp->proof_stream << "j " << _imp->proof_line << " 1 ~x" << _imp->variable_mappings[pair{ p.first, t.first }];
    for (auto & u : d_n_t)
        if (u != t)
            *_imp->proof_stream << " 1 x" << _imp->variable_mappings[pair{ q.first, u.first }];
    *_imp->proof_stream << " >= 1 ;" << endl;
    ++_imp->proof_line;

    _imp->adjacency_lines.emplace(tuple{ g, p.first, q.first, t.first }, _imp->proof_line);

    *_imp->proof_stream << "w 1" << endl;
}

auto Proof::hack_in_shape_graph(
        int g,
        const NamedVertex & p,
        const NamedVertex & q,
        const NamedVertex & t,
        const std::vector<NamedVertex> & n_t
        ) -> void
{
    *_imp->proof_stream << "* adjacency " << p.second << " maps to " << t.second <<
        " in shape graph " << g << " so " << q.second << " maps to one of..." << endl;
    *_imp->proof_stream << "a 1 ~x" << _imp->variable_mappings[pair{ p.first, t.first }];
    for (auto & u : n_t)
        *_imp->proof_stream << " 1 x" << _imp->variable_mappings[pair{ q.first, u.first }];
    *_imp->proof_stream << " >= 1 ;" << endl;
    ++_imp->proof_line;

    _imp->adjacency_lines.emplace(tuple{ g, p.first, q.first, t.first }, _imp->proof_line);
}

auto Proof::create_distance3_graphs_but_actually_distance_1(
        int g,
        const NamedVertex & p,
        const NamedVertex & q,
        const NamedVertex & t,
        const vector<NamedVertex> & d3_from_t) -> void
{
    *_imp->proof_stream << "* adjacency " << p.second << " maps to " << t.second <<
        " in G^3 so by adjacency, " << q.second << " maps to one of..." << endl;

    *_imp->proof_stream << "j " << _imp->adjacency_lines[tuple{ 0, p.first, q.first, t.first }]
        << " 1 ~x" << _imp->variable_mappings[pair{ p.first, t.first }];
    for (auto & u : d3_from_t)
        *_imp->proof_stream << " 1 x" << _imp->variable_mappings[pair{ q.first, u.first }];
    *_imp->proof_stream << " >= 1 ;" << endl;
    ++_imp->proof_line;

    _imp->adjacency_lines.emplace(tuple{ g, p.first, q.first, t.first }, _imp->proof_line);
}

auto Proof::create_distance3_graphs_but_actually_distance_2(
        int g,
        const NamedVertex & p,
        const NamedVertex & q,
        const NamedVertex & path_from_p_to_q,
        const NamedVertex & t,
        const vector<NamedVertex> & d1_from_t,
        const vector<NamedVertex> & d2_from_t,
        const vector<NamedVertex> & d3_from_t
        ) -> void
{
    *_imp->proof_stream << "* adjacency " << p.second << " maps to " << t.second <<
        " in G^3 so using vertex " << path_from_p_to_q.second << ", " << q.second << " maps to one of..." << endl;

    *_imp->proof_stream << "# 1" << endl;

    *_imp->proof_stream << "p";

    // if p maps to t then the first thing on the path from p to q has to go to one of...
    *_imp->proof_stream << " " << _imp->adjacency_lines[tuple{ 0, p.first, path_from_p_to_q.first, t.first }];
    // so the second thing on the path from p to q has to go to one of...
    for (auto & u : d1_from_t)
        *_imp->proof_stream << " " << _imp->adjacency_lines[tuple{ 0, path_from_p_to_q.first, q.first, u.first }] << " +";

    *_imp->proof_stream << " 0" << endl;
    ++_imp->proof_line;

    // tidy up
    *_imp->proof_stream << "j " << _imp->proof_line << " 1 ~x" << _imp->variable_mappings[pair{ p.first, t.first }];
    for (auto & u : d2_from_t)
        *_imp->proof_stream << " 1 x" << _imp->variable_mappings[pair{ q.first, u.first }];
    *_imp->proof_stream << " >= 1 ;" << endl;
    ++_imp->proof_line;

    *_imp->proof_stream << "# 0" << endl;

    *_imp->proof_stream << "j " << _imp->proof_line << " 1 ~x" << _imp->variable_mappings[pair{ p.first, t.first }];
    for (auto & u : d3_from_t)
        *_imp->proof_stream << " 1 x" << _imp->variable_mappings[pair{ q.first, u.first }];
    *_imp->proof_stream << " >= 1 ;" << endl;
    ++_imp->proof_line;

    _imp->adjacency_lines.emplace(tuple{ g, p.first, q.first, t.first }, _imp->proof_line);
}

auto Proof::create_distance3_graphs(
        int g,
        const NamedVertex & p,
        const NamedVertex & q,
        const NamedVertex & path_from_p_to_q_1,
        const NamedVertex & path_from_p_to_q_2,
        const NamedVertex & t,
        const vector<NamedVertex> & d1_from_t,
        const vector<NamedVertex> & d2_from_t,
        const vector<NamedVertex> & d3_from_t
        ) -> void
{
    *_imp->proof_stream << "* adjacency " << p.second << " maps to " << t.second <<
        " in G^3 so using path " << path_from_p_to_q_1.second << " -- " << path_from_p_to_q_2.second << ", " << q.second << " maps to one of..." << endl;

    *_imp->proof_stream << "# 1" << endl;

    *_imp->proof_stream << "p";

    // if p maps to t then the first thing on the path from p to q has to go to one of...
    *_imp->proof_stream << " " << _imp->adjacency_lines[tuple{ 0, p.first, path_from_p_to_q_1.first, t.first }];
    // so the second thing on the path from p to q has to go to one of...
    for (auto & u : d1_from_t)
        *_imp->proof_stream << " " << _imp->adjacency_lines[tuple{ 0, path_from_p_to_q_1.first, path_from_p_to_q_2.first, u.first }] << " +";

    *_imp->proof_stream << " 0" << endl;
    ++_imp->proof_line;

    // tidy up
    *_imp->proof_stream << "j " << _imp->proof_line << " 1 ~x" << _imp->variable_mappings[pair{ p.first, t.first }];
    for (auto & u : d2_from_t)
        *_imp->proof_stream << " 1 x" << _imp->variable_mappings[pair{ path_from_p_to_q_2.first, u.first }];
    *_imp->proof_stream << " >= 1 ;" << endl;
    ++_imp->proof_line;

    *_imp->proof_stream << "p " << _imp->proof_line;
    for (auto & u : d2_from_t)
        *_imp->proof_stream << " " << _imp->adjacency_lines[tuple{ 0, path_from_p_to_q_2.first, q.first, u.first }] << " +";
    *_imp->proof_stream << " 0" << endl;
    ++_imp->proof_line;

    *_imp->proof_stream << "# 0" << endl;

    *_imp->proof_stream << "j " << _imp->proof_line << " 1 ~x" << _imp->variable_mappings[pair{ p.first, t.first }];
    for (auto & u : d3_from_t)
        *_imp->proof_stream << " 1 x" << _imp->variable_mappings[pair{ q.first, u.first }];
    *_imp->proof_stream << " >= 1 ;" << endl;
    ++_imp->proof_line;

    _imp->adjacency_lines.emplace(tuple{ g, p.first, q.first, t.first }, _imp->proof_line);
}

auto Proof::create_binary_variable(int vertex,
                const function<auto (int) -> string> & name) -> void
{
    if (_imp->friendly_names)
        _imp->binary_variable_mappings.emplace(vertex, name(vertex));
    else
        _imp->binary_variable_mappings.emplace(vertex, to_string(_imp->binary_variable_mappings.size() + 1));
}

auto Proof::create_objective(int n, optional<int> d) -> void
{
    if (d) {
        _imp->model_stream << "* objective" << endl;
        for (int v = 0 ; v < n ; ++ v)
            _imp->model_stream << "1 x" << _imp->binary_variable_mappings[v] << " ";
        _imp->model_stream << ">= " << *d << ";" << endl;
        _imp->objective_line = ++_imp->nb_constraints;
    }
    else {
        _imp->model_prelude_stream << "min:";
        for (int v = 0 ; v < n ; ++ v)
            _imp->model_prelude_stream << " -1 x" << _imp->binary_variable_mappings[v];
        _imp->model_prelude_stream << " ;" << endl;
    }
}

auto Proof::create_non_edge_constraint(int p, int q) -> void
{
    _imp->model_stream << "-1 x" << _imp->binary_variable_mappings[p] << " -1 x" << _imp->binary_variable_mappings[q] << " >= -1 ;" << endl;

    ++_imp->nb_constraints;
    _imp->non_edge_constraints.emplace(pair{ p, q }, _imp->nb_constraints);
    _imp->non_edge_constraints.emplace(pair{ q, p }, _imp->nb_constraints);
}

auto Proof::create_non_null_decision_bound(int p, int t, optional<int> d) -> void
{
    if (d) {
        _imp->model_stream << "* objective" << endl;
        for (int v = 0 ; v < p ; ++v)
            for (int w = 0 ; w < t ; ++w)
                _imp->model_stream << "1 x" << _imp->variable_mappings[pair{ v, w }] << " ";
        _imp->model_stream << ">= " << *d << " ;" << endl;
        _imp->objective_line = ++_imp->nb_constraints;
    }
    else {
        _imp->model_prelude_stream << "min:";
        for (int v = 0 ; v < p ; ++v)
            for (int w = 0 ; w < t ; ++w)
                _imp->model_prelude_stream << " -1 x" << _imp->variable_mappings[pair{ v, w }] << " ";
        _imp->model_prelude_stream << " ;" << endl;
    }
}

auto Proof::backtrack_from_binary_variables(const vector<int> & v) -> void
{
    if (! _imp->doing_hom_colour_proof) {
        *_imp->proof_stream << "u";
        for (auto & w : v)
            *_imp->proof_stream << " 1 ~x" << _imp->binary_variable_mappings[w];
        *_imp->proof_stream << " >= 1 ;" << endl;
        ++_imp->proof_line;
    }
    else {
        *_imp->proof_stream << "* backtrack shenanigans, depth " << v.size() << endl;
        function<auto (unsigned, const vector<pair<int, int> > &) -> void> f;
        f = [&] (unsigned d, const vector<pair<int, int> > & trail) -> void {
            if (d == v.size()) {
                *_imp->proof_stream << "u 1 ~x" << _imp->variable_mappings[pair{ _imp->hom_colour_proof_p.first, _imp->hom_colour_proof_t.first }];
                for (auto & t : trail)
                    *_imp->proof_stream << " 1 ~x" << _imp->variable_mappings[t];
                *_imp->proof_stream << " >= 1 ;" << endl;
                ++_imp->proof_line;
            }
            else {
                for (auto & p : _imp->p_clique) {
                    vector<pair<int, int> > new_trail{ trail };
                    new_trail.emplace_back(pair{ p.first, _imp->t_clique_neighbourhood.find(v[d])->second.first });
                    f(d + 1, new_trail);
                }
            }
        };
        f(0, {});
    }
}

auto Proof::colour_bound(const vector<vector<int> > & ccs) -> void
{
    *_imp->proof_stream << "* bound, ccs";
    for (auto & cc : ccs) {
        *_imp->proof_stream << " [";
        for (auto & c : cc)
            *_imp->proof_stream << " " << c;
        *_imp->proof_stream << " ]";
    }
    *_imp->proof_stream << endl;

    vector<long> to_sum;
    auto do_one_cc = [&] (const auto & cc, const auto & non_edge_constraint) {
        if (cc.size() > 2) {
            *_imp->proof_stream << "p " << non_edge_constraint(cc[0], cc[1]);

            for (unsigned i = 2 ; i < cc.size() ; ++i) {
                *_imp->proof_stream << " " << i << " *";
                for (unsigned j = 0 ; j < i ; ++j)
                    *_imp->proof_stream << " " << non_edge_constraint(cc[i], cc[j]) << " +";
                *_imp->proof_stream << " " << (i + 1) << " d";
            }

            *_imp->proof_stream << endl;
            to_sum.push_back(++_imp->proof_line);
        }
        else if (cc.size() == 2) {
            to_sum.push_back(non_edge_constraint(cc[0], cc[1]));
        }
    };

    for (auto & cc : ccs) {
        if (_imp->doing_hom_colour_proof) {
            vector<pair<NamedVertex, NamedVertex> > bigger_cc;
            for (auto & c : cc)
                for (auto & v : _imp->p_clique)
                    bigger_cc.push_back(pair{ v, _imp->t_clique_neighbourhood.find(c)->second });

            *_imp->proof_stream << "* colour class [";
            for (auto & c : bigger_cc)
                *_imp->proof_stream << " " << c.first.second << "/" << c.second.second;
            *_imp->proof_stream << " ]" << endl;

            do_one_cc(bigger_cc, [&] (const pair<NamedVertex, NamedVertex> & a, const pair<NamedVertex, NamedVertex> & b) -> long {
                    return _imp->clique_for_hom_non_edge_constraints[pair{ a, b }];
                    });
        }
        else
            do_one_cc(cc, [&] (int a, int b) -> long { return _imp->non_edge_constraints[pair{ a, b }]; });

        *_imp->proof_stream << "p " << _imp->objective_line;
        for (auto & t : to_sum)
            *_imp->proof_stream << " " << t << " +";
        *_imp->proof_stream << endl;
        ++_imp->proof_line;
    }
}

auto Proof::prepare_hom_clique_proof(const NamedVertex & p, const NamedVertex & t, unsigned size) -> void
{
    *_imp->proof_stream << "* clique of size " << size << " around neighbourhood of " << p.second << " but not " << t.second << endl;
    *_imp->proof_stream << "# 1" << endl;
    _imp->doing_hom_colour_proof = true;
    _imp->hom_colour_proof_p = p;
    _imp->hom_colour_proof_t = t;
}

auto Proof::start_hom_clique_proof(const NamedVertex & p, vector<NamedVertex> && p_clique, const NamedVertex & t, map<int, NamedVertex> && t_clique_neighbourhood) -> void
{
    _imp->p_clique = move(p_clique);
    _imp->t_clique_neighbourhood = move(t_clique_neighbourhood);

    *_imp->proof_stream << "* hom clique objective" << endl;
    vector<long> to_sum;
    for (auto & q : _imp->p_clique) {
        *_imp->proof_stream << "u 1 ~x" << _imp->variable_mappings[pair{ p.first, t.first }];
        for (auto & u : _imp->t_clique_neighbourhood)
            *_imp->proof_stream << " 1 x" << _imp->variable_mappings[pair{ q.first, u.second.first }];
        *_imp->proof_stream << " >= 1 ;" << endl;
        to_sum.push_back(++_imp->proof_line);
    }

    *_imp->proof_stream << "p";
    bool first = true;
    for (auto & t : to_sum) {
        *_imp->proof_stream << " " << t;
        if (! first)
            *_imp->proof_stream << " +";
        first = false;
    }
    *_imp->proof_stream << endl;
    _imp->objective_line = ++_imp->proof_line;

    *_imp->proof_stream << "* hom clique non edges for injectivity" << endl;

    for (auto & p : _imp->p_clique)
        for (auto & q : _imp->p_clique)
            if (p != q) {
                for (auto & [ _, t ] : _imp->t_clique_neighbourhood) {
                    *_imp->proof_stream << "u 1 ~x" << _imp->variable_mappings[pair{ p.first, t.first }] << " 1 ~x" << _imp->variable_mappings[pair{ q.first, t.first }] << " >= 1 ;" << endl;
                    ++_imp->proof_line;
                    _imp->clique_for_hom_non_edge_constraints.emplace(pair{ pair{ p, t }, pair{ q, t } }, _imp->proof_line);
                    _imp->clique_for_hom_non_edge_constraints.emplace(pair{ pair{ q, t }, pair{ p, t } }, _imp->proof_line);
                }
            }

    *_imp->proof_stream << "* hom clique non edges for variables" << endl;

    for (auto & p : _imp->p_clique)
        for (auto & [ _, t ] : _imp->t_clique_neighbourhood) {
            for (auto & [ _, u ] : _imp->t_clique_neighbourhood) {
                if (t != u) {
                    *_imp->proof_stream << "u 1 ~x" << _imp->variable_mappings[pair{ p.first, t.first }] << " 1 ~x" << _imp->variable_mappings[pair{ p.first, u.first }] << " >= 1 ;" << endl;
                    ++_imp->proof_line;
                    _imp->clique_for_hom_non_edge_constraints.emplace(pair{ pair{ p, t }, pair{ p, u } }, _imp->proof_line);
                    _imp->clique_for_hom_non_edge_constraints.emplace(pair{ pair{ p, u }, pair{ p, t } }, _imp->proof_line);
                }
            }
        }
}

auto Proof::finish_hom_clique_proof(const NamedVertex & p, const NamedVertex & t, unsigned size) -> void
{
    *_imp->proof_stream << "* end clique of size " << size << " around neighbourhood of " << p.second << " but not " << t.second << endl;
    *_imp->proof_stream << "# 0" << endl;
    *_imp->proof_stream << "u 1 ~x" << _imp->variable_mappings[pair{ p.first, t.first }] << " >= 1 ;" << endl;
    *_imp->proof_stream << "w 1" << endl;
    ++_imp->proof_line;
    _imp->doing_hom_colour_proof = false;
    _imp->clique_for_hom_non_edge_constraints.clear();
}

auto Proof::add_hom_clique_non_edge(
        const NamedVertex & pp,
        const NamedVertex & tt,
        const std::vector<NamedVertex> & p_clique,
        const NamedVertex & t,
        const NamedVertex & u) -> void
{
    *_imp->proof_stream << "* hom clique non edges for " << t.second << " " << u.second << endl;
    for (auto & p : p_clique) {
        for (auto & q : p_clique) {
            if (p != q) {
                *_imp->proof_stream << "u 1 ~x" << _imp->variable_mappings[pair{ pp.first, tt.first }]
                    << " 1 ~x" << _imp->variable_mappings[pair{ p.first, t.first }]
                    << " 1 ~x" << _imp->variable_mappings[pair{ q.first, u.first }] << " >= 1 ;" << endl;
                ++_imp->proof_line;
                _imp->clique_for_hom_non_edge_constraints.emplace(pair{ pair{ p, t }, pair{ q, u } }, _imp->proof_line);
                _imp->clique_for_hom_non_edge_constraints.emplace(pair{ pair{ q, u }, pair{ p, t } }, _imp->proof_line);
            }
        }
    }
}

auto Proof::mcs_bound(
        const vector<pair<set<int>, set<int> > > & partitions) -> void
{
    *_imp->proof_stream << "* failed bound" << endl;

    vector<string> to_sum;
    for (auto & [ l, r ] : partitions) {
        if (r.size() >= l.size())
            continue;

        *_imp->proof_stream << "p";
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

        *_imp->proof_stream << endl;
        to_sum.push_back(to_string(++_imp->proof_line));
    }

    if (! to_sum.empty()) {
        *_imp->proof_stream << "p " << _imp->objective_line;
        for (auto & t : to_sum)
            *_imp->proof_stream << " " << t << " +";
        *_imp->proof_stream << endl;
        ++_imp->proof_line;
    }
}

auto Proof::rewrite_mcs_objective(int pattern_size) -> void
{
    *_imp->proof_stream << "* get the objective function to talk about nulls, not non-nulls" << endl;
    *_imp->proof_stream << "p " << _imp->objective_line;
    for (int v = 0 ; v < pattern_size ; ++v)
        *_imp->proof_stream << " " << _imp->at_most_one_value_constraints[v] << " +";
    *_imp->proof_stream << endl;
    _imp->objective_line = ++_imp->proof_line;
}

auto Proof::create_connected_constraints(int p, int t, const function<auto (int, int) -> bool> & adj) -> void
{
    _imp->model_stream << "* selected vertices must be connected, walk 1" << endl;
    int mapped_to_null = t;
    int cnum = _imp->variable_mappings.size();

    for (int v = 0 ; v < p ; ++v)
        for (int w = 0 ; w < v ; ++w) {
            string n = _imp->friendly_names ? ("conn1_" + to_string(v) + "_" + to_string(w)) : to_string(++cnum);
            _imp->connected_variable_mappings.emplace(tuple{ 1, v, w }, n);
            if (! adj(v, w)) {
                // v not adjacent to w, so the walk does not exist
                _imp->model_stream << "1 ~x" << n << " >= 1 ;" << endl;
                ++_imp->nb_constraints;
            }
            else {
                // v = null -> the walk does not exist
                _imp->model_stream << "1 ~x" << n << " 1 ~x" << _imp->variable_mappings[pair{ v, mapped_to_null }] << " >= 1 ;" << endl;
                // w = null -> the walk does not exist
                _imp->model_stream << "1 ~x" << n << " 1 ~x" << _imp->variable_mappings[pair{ w, mapped_to_null }] << " >= 1 ;" << endl;
                // either v = null, or w = null, or the walk exists
                _imp->model_stream << "1 x" << n << " 1 x" << _imp->variable_mappings[pair{ v, mapped_to_null }]
                    << " 1 x" << _imp->variable_mappings[pair{ w, mapped_to_null }] << " >= 1 ;" << endl;
                _imp->nb_constraints += 3;
            }
        }

    int last_k = 0;
    for (int k = 2 ; k < 2 * t ; k *= 2) {
        last_k = k;
        _imp->model_stream << "* selected vertices must be connected, walk " << k << endl;
        for (int v = 0 ; v < p ; ++v)
            for (int w = 0 ; w < v ; ++w) {
                string n = _imp->friendly_names ? ("conn" + to_string(k) + "_" + to_string(v) + "_" + to_string(w)) : to_string(++cnum);
                _imp->connected_variable_mappings.emplace(tuple{ k, v, w }, n);

                vector<string> ors;
                for (int u = 0 ; u < p ; ++u) {
                    if (v != w && v != u && u != w) {
                        string m = _imp->friendly_names ?
                            ("conn" + to_string(k) + "_" + to_string(v) + "_" + to_string(w) + "_via_" + to_string(u)) :
                            to_string(++cnum);
                        ors.push_back(m);
                        _imp->connected_variable_mappings_aux.emplace(tuple{ k, v, w, u }, m);
                        // either the first half walk exists, or the via term is false
                        _imp->model_stream << "1 x" << _imp->connected_variable_mappings[tuple{ k / 2, max(u, v), min(u, v) }]
                            << " 1 ~x" << m << " >= 1 ;" << endl;
                        // either the second half walk exists, or the via term is false
                        _imp->model_stream << "1 x" << _imp->connected_variable_mappings[tuple{ k / 2, max(u, w), min(u, w) }]
                            << " 1 ~x" << m << " >= 1 ;" << endl;
                        // one of the half walks is false, or the via term must be true
                        _imp->model_stream << "1 x" << m
                            << " 1 ~x" << _imp->connected_variable_mappings[tuple{ k / 2, max(v, u), min(v, u) }]
                            << " 1 ~x" << _imp->connected_variable_mappings[tuple{ k / 2, max(u, w), min(u, w) }] << " >= 1 ;" << endl;
                        _imp->nb_constraints += 3;
                    }
                }

                // one of the vias must be true, or a shorter walk exists, or the entry is false
                _imp->model_stream << "1 ~x" << n;
                for (auto & o : ors)
                    _imp->model_stream << " 1 x" << o;
                _imp->model_stream << " 1 x" << _imp->connected_variable_mappings[tuple{ k / 2, v, w }];
                _imp->model_stream << " >= 1 ;" << endl;
                ++_imp->nb_constraints;

                // if the entry is false, then all of the vias must be false and the shorter walk must be false
                for (auto & o : ors) {
                    _imp->model_stream << "1 x" << n << " 1 ~x" << o << " >= 1 ;" << endl;
                    ++_imp->nb_constraints;
                }
                _imp->model_stream << "1 x" << n << " 1 ~x" << _imp->connected_variable_mappings[tuple{ k / 2, v, w }] << " >= 1 ;" << endl;
                ++_imp->nb_constraints;
            }
    }

    _imp->model_stream << "* if two vertices are used, they must be connected" << endl;
    for (int v = 0 ; v < p ; ++v)
        for (int w = 0 ; w < v ; ++w) {
            _imp->model_stream << "1 x" << _imp->variable_mappings[pair{ v, mapped_to_null }]
                << " 1 x" << _imp->variable_mappings[pair{ w, mapped_to_null }];
            auto var = _imp->connected_variable_mappings.find(tuple{ last_k, v, w });
            if (var != _imp->connected_variable_mappings.end())
                _imp->model_stream << " 1 x" << _imp->connected_variable_mappings[tuple{ last_k, v, w }];
           _imp->model_stream << " >= 1 ;" << endl;
            ++_imp->nb_constraints;
        }
}

auto Proof::not_connected_in_underlying_graph(const std::vector<int> & x, int y) -> void
{
    *_imp->proof_stream << "u 1 ~x" << _imp->binary_variable_mappings[y];
    for (auto & v : x)
        *_imp->proof_stream << " 1 ~x" << _imp->binary_variable_mappings[v];
    *_imp->proof_stream << " >= 1 ;" << endl;
    ++_imp->proof_line;
}

auto Proof::has_clique_model() const -> bool
{
    return _imp->clique_encoding;
}

auto Proof::create_clique_encoding(
        const vector<pair<int, int> > & enc,
        const vector<pair<int, int> > & zero_in_proof_objectives) -> void
{
    _imp->clique_encoding = true;
    for (unsigned i = 0 ; i < enc.size() ; ++i)
        _imp->binary_variable_mappings.emplace(i, _imp->variable_mappings[enc[i]]);

    _imp->zero_in_proof_objectives = zero_in_proof_objectives;
}

auto Proof::create_clique_nonedge(int v, int w) -> void
{
    *_imp->proof_stream << "u 1 ~x" << _imp->binary_variable_mappings[v]
        << " 1 ~x" << _imp->binary_variable_mappings[w] << " >= 1 ;" << endl;
    ++_imp->proof_line;
    _imp->non_edge_constraints.emplace(pair{ v, w }, _imp->proof_line);
    _imp->non_edge_constraints.emplace(pair{ w, v }, _imp->proof_line);
}

auto Proof::super_extra_verbose() const -> bool
{
    return _imp->super_extra_verbose;
}

auto Proof::show_domains(const string & s, const std::vector<std::pair<NamedVertex, std::vector<NamedVertex> > > & domains) -> void
{
    *_imp->proof_stream << "* " << s << ", domains follow" << endl;
    for (auto & [ p, ts ] : domains) {
        *_imp->proof_stream << "*    " << p.second << " size " << ts.size() << " = {";
        for (auto & t : ts)
            *_imp->proof_stream << " " << t.second;
        *_imp->proof_stream << " }" << endl;
    }
}

auto Proof::propagated(const NamedVertex & p, const NamedVertex & t, int g, int n_values, const NamedVertex & q) -> void
{
    *_imp->proof_stream << "* adjacency propagation from " << p.second << " -> " << t.second << " in graph pairs " << g << " deleted " << n_values << " values from " << q.second << endl;
}

