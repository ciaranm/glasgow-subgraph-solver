#include <gss/clique.hh>
#include <gss/innards/homomorphism_proofs.hh>

#include <algorithm>
#include <chrono>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

using namespace gss;
using namespace gss::innards;

using std::make_optional;
using std::make_unique;
using std::map;
using std::move;
using std::nullopt;
using std::optional;
using std::pair;
using std::set;
using std::shared_ptr;
using std::string;
using std::vector;

using std::chrono::steady_clock;

HomomorphismProofs::HomomorphismProofs(const shared_ptr<Proof> & proof, const InputGraph & pattern, const InputGraph & target) :
    _proof(proof)
{
    for (int v = 0; v < pattern.size(); ++v)
        _pattern_names.push_back(string{pattern.vertex_name(v)});
    for (int v = 0; v < target.size(); ++v)
        _target_names.push_back(string{target.vertex_name(v)});
}

auto HomomorphismProofs::pattern_vertex(int v) const -> NamedVertex
{
    if (v < 0 || unsigned(v) >= _pattern_names.size())
        throw ProofError{"Oops, there's a bug: v out of range in pattern"};
    return pair{v, _pattern_names[v]};
}

auto HomomorphismProofs::target_vertex(int v) const -> NamedVertex
{
    if (v < 0 || unsigned(v) >= _target_names.size())
        throw ProofError{"Oops, there's a bug: v out of range in target"};
    return pair{v, _target_names[v]};
}

auto HomomorphismProofs::register_supplemental(const std::tuple<long, long, long, long> & key, int p, int t,
    std::function<void()> emit) -> void
{
    if (_pending_supplementals.contains(key) || _adjacency.labels.contains(key))
        return;
    _pending_supplementals.emplace(key, move(emit));
    _pending_by_antecedent[pair<int, int>{p, t}].push_back(key);
}

auto HomomorphismProofs::materialise_one(const std::tuple<long, long, long, long> & key) -> bool
{
    auto it = _pending_supplementals.find(key);
    if (it == _pending_supplementals.end())
        return false;
    auto emit = move(it->second);
    _pending_supplementals.erase(it);
    emit();
    return true;
}

auto HomomorphismProofs::materialise_adjacency_for(int p, int t) -> void
{
    auto it = _pending_by_antecedent.find(pair<int, int>{p, t});
    if (it == _pending_by_antecedent.end())
        return;
    auto keys = move(it->second);
    _pending_by_antecedent.erase(it);
    // each emit leaves the proof at level 0 (where the @label persists); restore the
    // search's active level once for the whole batch rather than once per supplemental.
    bool emitted_any = false;
    for (auto & key : keys)
        emitted_any = materialise_one(key) || emitted_any;
    if (emitted_any && _proof->active_level() != 0)
        _proof->emit_proof_directive("setlvl " + std::to_string(_proof->active_level()) + ";");
}

auto HomomorphismProofs::has_pending_supplementals() const -> bool
{
    return ! _pending_supplementals.empty();
}

auto HomomorphismProofs::create_locally_injective_constraints(const InputGraph & pattern, const InputGraph & target) -> void
{
    _locally_injective = true;

    for (int v = 0; v < pattern.size(); ++v) {
        // the neighbourhood of v: if v has a self-loop then v is its own neighbour, so
        // local injectivity also forces phi(v) to differ from its neighbours' images.
        std::vector<int> neighbours;
        for (int u = 0; u < pattern.size(); ++u)
            if (pattern.adjacent(v, u))
                neighbours.push_back(u);

        // at most one neighbour maps to a given target only bites for |N(v)| >= 2
        if (neighbours.size() < 2)
            continue;

        for (int t = 0; t < target.size(); ++t) {
            _proof->emit_model_comment("* local injectivity on neighbourhood of " + std::to_string(v) + " for value " + std::to_string(t));
            auto label = "@linj" + _pattern_names[v] + "_" + _target_names[t];
            std::string line = label;
            for (auto & u : neighbours)
                line += " -1 x" + _proof->variable_name(u, t);
            line += " >= -1 ;";
            _proof->emit_model_constraint(line);
            _locally_injective_constraints.emplace(std::pair{v, t}, label);
        }
    }
}

auto HomomorphismProofs::is_locally_injective() const -> bool
{
    return _locally_injective;
}

auto HomomorphismProofs::locally_injective_label(int p, int t) const -> const std::string &
{
    return _locally_injective_constraints.at(std::pair{p, t});
}

auto HomomorphismProofs::failure_due_to_pattern_bigger_than_target() -> void
{
    _proof->emit_proof_directive("% failure due to the pattern being bigger than the target");

    // we get a hall violator by adding up all of the things
    std::string pol = "pol";
    bool first = true;
    for (auto & [_, label] : _proof->at_least_one_value_labels()) {
        if (first) {
            pol += " " + label;
            first = false;
        }
        else
            pol += " " + label + " +";
    }

    for (auto & [_, label] : _proof->injectivity_labels())
        pol += " " + label + " +";
    pol += " ;";
    _proof->emit_proof_line(pol);
}

auto HomomorphismProofs::guessing(int depth, int p, int t) -> void
{
    _proof->emit_proof_directive("% [" + std::to_string(depth) + "] guessing " + _pattern_names[p] + "=" + _target_names[t]);
}

auto HomomorphismProofs::unit_propagating(int p, int t) -> void
{
    _proof->emit_proof_directive("% unit propagating " + _pattern_names[p] + "=" + _target_names[t]);
}

auto HomomorphismProofs::propagation_failure(const std::vector<std::pair<int, int>> & decisions, int p, int t) -> void
{
    _proof->emit_proof_directive("% [" + std::to_string(decisions.size()) + "] propagation failure on " +
        _pattern_names[p] + "=" + _target_names[t]);
    std::string line = "rup ";
    for (auto & [var, val] : decisions)
        line += " 1 ~x" + _proof->variable_name(var, val);
    line += " >= 1 ;";
    _proof->emit_proof_line(line);
}

auto HomomorphismProofs::propagated(int p, int t, int g, int n_values, int q) -> void
{
    _proof->emit_proof_directive("% adjacency propagation from " + _pattern_names[p] + " -> " + _target_names[t] +
        " in graph pairs " + std::to_string(g) + " deleted " + std::to_string(n_values) + " values from " + _pattern_names[q]);
}

auto HomomorphismProofs::show_domains(const std::string & where, const std::vector<std::pair<int, std::vector<int>>> & domains) -> void
{
    _proof->emit_proof_directive("% " + where + ", domains follow");
    for (auto & [p, ts] : domains) {
        std::string line = "%    " + _pattern_names[p] + " size " + std::to_string(ts.size()) + " = {";
        for (auto & t : ts)
            line += " " + _target_names[t];
        line += " }";
        _proof->emit_proof_directive(line);
    }
}

auto HomomorphismProofs::initial_domain_is_empty(int p, const std::string & where) -> void
{
    _proof->emit_proof_directive("% failure due to domain " + std::to_string(p) + " being empty at " + where);
}

auto HomomorphismProofs::emit_hall_set_or_violator(const std::vector<int> & lhs, const std::vector<int> & rhs) -> void
{
    std::string comment = "% hall set or violator {";
    for (auto & l : lhs)
        comment += " " + _pattern_names[l];
    comment += " } / {";
    for (auto & r : rhs)
        comment += " " + _target_names[r];
    comment += " }";
    _proof->emit_proof_directive(comment);

    std::string pol = "pol";
    bool first = true;
    for (auto & l : lhs) {
        if (first) {
            first = false;
            pol += " " + _proof->at_least_one_value_label(l);
        }
        else
            pol += " " + _proof->at_least_one_value_label(l) + " +";
    }
    for (auto & r : rhs)
        pol += " " + _proof->injectivity_label(r) + " +";
    pol += " ;";
    _proof->emit_proof_line(pol);
}

auto HomomorphismProofs::need_elimination(int p, int t) -> void
{
    if (! _eliminations.contains(std::pair{p, t})) {
        _proof->emit_proof_directive("setlvl 0;");
        _eliminations[std::pair{p, t}] = _proof->emit_proof_line("rup 1 ~x" + _proof->variable_name(p, t) + " >= 1 ;");
        _proof->emit_proof_directive("setlvl " + std::to_string(_proof->active_level()) + ";");
    }
}

auto HomomorphismProofs::incompatible_by_loops(int p, int t) -> void
{
    // may be requested both up front (so the unit is available to later derivations and
    // search propagations) and again during domain initialisation: only emit it once.
    if (_eliminations.contains(std::pair{p, t}))
        return;
    _proof->emit_proof_directive("% cannot map " + _pattern_names[p] + " to " + _target_names[t] + " due to loop");
    _eliminations.emplace(std::pair{p, t}, _proof->emit_proof_line("rup 1 ~x" + _proof->variable_name(p, t) + " >= 1 ;"));
}

auto HomomorphismProofs::incompatible_by_degrees(int g, int p, const std::vector<int> & n_p, int t, const std::vector<int> & n_t) -> void
{
    auto & adjacency = _adjacency;
    _proof->emit_proof_directive("% cannot map " + _pattern_names[p] + " to " + _target_names[t] +
        " due to degrees in graph pairs " + std::to_string(g));

    std::string pol = "pol";
    bool first = true;
    for (auto & n : n_p) {
        // due to loops or labels, it might not be possible to map n to t
        if (adjacency.labels.count(std::tuple<long, long, long, long>{g, p, n, t})) {
            if (first) {
                first = false;
                pol += " " + adjacency.labels.at(std::tuple<long, long, long, long>{g, p, n, t});
            }
            else
                pol += " " + adjacency.labels.at(std::tuple<long, long, long, long>{g, p, n, t}) + " +";
        }
    }

    // if I map p to t, I have to map the neighbours of p to distinct neighbours of t.
    // Under full injectivity that distinctness is the global injectivity on each value;
    // under local injectivity it is the neighbourhood-injectivity of p (phi|N(p) is
    // injective), which is exactly what the degree pigeonhole needs.
    for (auto & n : n_t)
        pol += " " + (is_locally_injective() ? locally_injective_label(p, n) : _proof->injectivity_label(n)) + " +";

    pol += " s ;";
    auto sum_line = _proof->emit_proof_line(pol);

    _proof->emit_proof_line("ia 1 ~x" + _proof->variable_name(p, t) + " >= 1 : " + std::to_string(sum_line) + " ;");
    auto ia_line = _proof->current_proof_line();
    _eliminations.emplace(std::pair{p, t}, ia_line);

    _proof->emit_proof_directive("del id " + std::to_string(ia_line - 1) + " ;");
}

auto HomomorphismProofs::incompatible_by_nds(int g, int p, int t, const std::vector<int> & p_subsequence,
    const std::vector<int> & t_subsequence, const std::vector<int> & t_remaining) -> void
{
    auto & adjacency = _adjacency;
    _proof->emit_proof_directive("% cannot map " + _pattern_names[p] + " to " + _target_names[t] +
        " due to nds in graph pairs " + std::to_string(g));

    for (auto & n : p_subsequence)
        for (auto & u : t_remaining)
            need_elimination(n, u);
    for (auto & n : p_subsequence)
        need_elimination(n, t_subsequence.back());

    // summing up horizontally
    std::string pol = "pol";
    bool first = true;
    for (auto & n : p_subsequence) {
        // due to loops or labels, it might not be possible to map n to t
        if (adjacency.labels.count(std::tuple<long, long, long, long>{g, p, n, t})) {
            if (first) {
                first = false;
                pol += " " + adjacency.labels.at(std::tuple<long, long, long, long>{g, p, n, t});
            }
            else
                pol += " " + adjacency.labels.at(std::tuple<long, long, long, long>{g, p, n, t}) + " +";
        }
    }

    // injectivity in the square: each column of the square holds at most one of p's
    // neighbours. Under full injectivity that is the global injectivity on the value;
    // under local injectivity it is the neighbourhood-injectivity of p (at most one
    // neighbour of p maps to t), exactly as in the degree pigeonhole above.
    for (auto & tsub : t_subsequence) {
        if (tsub != t_subsequence.back())
            pol += " " + (is_locally_injective() ? locally_injective_label(p, tsub) : _proof->injectivity_label(tsub)) + " +";
    }

    // block to the right of the failing square
    for (auto & n : p_subsequence) {
        for (auto & u : t_remaining) {
            /* n -> u is already eliminated by degree or loop */
            pol += " " + std::to_string(_eliminations[std::pair{n, u}]) + " +";
        }
    }

    // final column
    for (auto & n : p_subsequence) {
        /* n -> t_subsequence.back() is already eliminated by degree or loop */
        pol += " " + std::to_string(_eliminations[std::pair{n, t_subsequence.back()}]) + " +";
    }

    pol += " s ;";
    auto sum_line = _proof->emit_proof_line(pol);

    _proof->emit_proof_line("ia 1 ~x" + _proof->variable_name(p, t) + " >= 1 : " + std::to_string(sum_line) + " ;");
    auto ia_line = _proof->current_proof_line();

    _proof->emit_proof_directive("del id " + std::to_string(ia_line - 1) + " ;");
}

auto HomomorphismProofs::emit_adjacency_constraint(int p, int q, int t, const std::vector<int> & permitted) -> void
{
    std::string adj_label = "@adj" + _pattern_names[p] + "_" + _target_names[t] + "_" + _pattern_names[q];
    std::string line = adj_label + " 1 ~x" + _proof->variable_name(p, t);
    for (auto & u : permitted)
        line += " 1 x" + _proof->variable_name(q, u);
    line += " >= 1 ;";
    _proof->emit_model_constraint(line);

    auto & adjacency = _adjacency;
    adjacency.labels.emplace(std::tuple<long, long, long, long>{0, p, q, t}, adj_label);
    adjacency.permitted.emplace(std::tuple<long, long, long, long>{0, p, q, t},
        std::vector<long>(permitted.begin(), permitted.end()));
}

auto HomomorphismProofs::emit_exact_path_graph(int g, int p, int q, const std::vector<int> & between_p_and_q,
    int t, const std::vector<int> & n_t, const std::vector<std::pair<int, std::vector<int>>> & two_away_from_t,
    const std::vector<int> & d_n_t) -> void
{
    auto & adjacency = _adjacency;

    // tidy up to get what we wanted. do this first so we can check for duplicates
    std::string tidied_up = "1 ~x" + _proof->variable_name(p, t);
    for (auto & u : d_n_t)
        if (u != t)
            tidied_up += " 1 x" + _proof->variable_name(q, u);
    tidied_up += " >= 1 :";

    if (auto cached = _proof->cached_proof_line(tidied_up)) {
        adjacency.labels.emplace(std::tuple<long, long, long, long>{g, p, q, t}, *cached);
        return;
    }

    _proof->emit_proof_directive("% adjacency " + _pattern_names[p] + " maps to " + _target_names[t] +
        " in G^[" + std::to_string(g) + "x2] so " + _pattern_names[q] + " maps to one of...");

    // Scratch one level above the search: wiplvl wipes every level >= its argument (that is
    // how forget_level discards a search subtree), so a hardcoded level 1 would wipe the whole
    // search trail when materialising mid-search. The @label is emitted at level 0 below so it
    // persists; the caller restores the search level.
    int scratch = _proof->active_level() + 1;
    _proof->emit_proof_directive("setlvl " + std::to_string(scratch) + ";");

    // if p maps to t then things in between_p_and_q have to go to one of these, and then go
    // two hops out cancelling between_p_and_q things with where q can go
    std::string pol = "pol";
    bool first = true;
    for (auto & b : between_p_and_q) {
        pol += " " + adjacency.labels.at(std::tuple<long, long, long, long>{0, p, b, t});
        if (! first)
            pol += " s +";
        first = false;
    }
    for (auto & b : between_p_and_q)
        for (auto & w : n_t)
            // due to loops or labels, it might not be possible to map to w
            if (adjacency.labels.contains(std::tuple<long, long, long, long>{0, b, q, w}))
                pol += " " + adjacency.labels.at(std::tuple<long, long, long, long>{0, b, q, w}) + " +";
    pol += " s ;";
    _proof->emit_proof_line(pol);

    // first tidy-up step: if p maps to t then q maps to something a two-walk away from t. The
    // adjacency constraints summed above are the loop-cancelled forms, so plain implication
    // addition closes it.
    {
        std::string line = "ia 1 ~x" + _proof->variable_name(p, t);
        for (auto & u : two_away_from_t)
            line += " 1 x" + _proof->variable_name(q, u.first);
        line += " >= 1 : " + std::to_string(_proof->current_proof_line()) + " ;";
        _proof->emit_proof_line(line);
    }

    // if p maps to t then q does not map to t. Under full injectivity that is the global
    // injectivity on t; under local injectivity p and q share a common neighbour b (that is
    // what between_p_and_q holds), so the neighbourhood-injectivity of b forbids them both
    // mapping to t. Either way the constraint cancels the "q maps to t" term.
    {
        const std::string & inj = is_locally_injective()
            ? locally_injective_label(between_p_and_q.front(), t)
            : _proof->injectivity_label(t);
        _proof->emit_proof_line("pol " + std::to_string(_proof->current_proof_line()) + " " + inj + " + s ;");
    }

    // and cancel out stray extras from injectivity
    {
        std::string line = "ia 1 ~x" + _proof->variable_name(p, t);
        for (auto & u : two_away_from_t)
            if (u.first != t)
                line += " 1 x" + _proof->variable_name(q, u.first);
        line += " >= 1 : " + std::to_string(_proof->current_proof_line()) + " ;";
        _proof->emit_proof_line(line);
    }

    std::vector<long> things_to_add_up;
    things_to_add_up.push_back(_proof->current_proof_line());

    // cancel out anything that is two away from t, but by insufficiently many paths
    for (auto & u : two_away_from_t) {
        if ((u.first == t) || (std::find(d_n_t.begin(), d_n_t.end(), u.first) != d_n_t.end()))
            continue;

        std::string pol2 = "pol";
        bool first2 = true;
        for (auto & b : between_p_and_q) {
            pol2 += " " + adjacency.labels.at(std::tuple<long, long, long, long>{0, p, b, t});
            if (! first2)
                pol2 += " +";
            first2 = false;
            pol2 += " " + adjacency.labels.at(std::tuple<long, long, long, long>{0, q, b, u.first}) + " +";
            pol2 += " " + _proof->at_most_one_value_label(b) + " +";
        }
        // the between-vertices must map to distinct common neighbours of t and u (the z's):
        // global injectivity on each z, or under local injectivity the neighbourhood-
        // injectivity of p (the between-vertices are all neighbours of p) -- the same pigeonhole.
        for (auto & z : u.second)
            pol2 += " " + (is_locally_injective() ? locally_injective_label(p, z) : _proof->injectivity_label(z)) + " +";
        pol2 += " s ;";
        _proof->emit_proof_line(pol2);

        // want: ~x_p_t + ~x_q_u >= 1
        std::string line = "ia 1 ~x" + _proof->variable_name(p, t) + " 1 ~x" + _proof->variable_name(q, u.first) +
            " >= 1 : " + std::to_string(_proof->current_proof_line()) + " ;";
        things_to_add_up.push_back(_proof->emit_proof_line(line));
    }

    // do the getting rid of
    if (things_to_add_up.size() > 1) {
        std::string pol3 = "pol";
        bool first3 = true;
        for (auto & line_id : things_to_add_up) {
            pol3 += " " + std::to_string(line_id);
            if (! first3)
                pol3 += " +";
            first3 = false;
        }
        pol3 += " s ;";
        _proof->emit_proof_line(pol3);
    }

    _proof->emit_proof_directive("setlvl 0;");
    std::string adj_label = "@g" + std::to_string(g) + "adj" + _pattern_names[p] + "_" + _target_names[t] + "_" + _pattern_names[q];
    _proof->emit_proof_line(adj_label + " ia " + tidied_up + " " + std::to_string(_proof->current_proof_line()) + " ;");
    adjacency.labels.emplace(std::tuple<long, long, long, long>{g, p, q, t}, adj_label);
    _proof->cache_proof_line(tidied_up, adj_label);
    _proof->wipe_level(scratch);
}

auto HomomorphismProofs::prove_exact_path_graphs(const ProcessedGraphsData & graphs, unsigned max_graphs,
    const std::vector<std::pair<int, unsigned>> & exact_path_index_and_slot, unsigned exact_path_1_slot,
    bool elide_subsumed) -> std::set<std::pair<int, int>>
{
    const unsigned pattern_size = _pattern_names.size();
    const unsigned target_size = _target_names.size();

    std::set<std::pair<int, int>> covered;

    for (unsigned p = 0; p < pattern_size; ++p) {
        for (unsigned q = 0; q < pattern_size; ++q) {
            if (p == q)
                continue;

            // the exact-path indices (with their slot) whose pattern edge (p,q) holds.
            vector<pair<int, unsigned>> emit_for;
            for (auto & [gg, ss] : exact_path_index_and_slot)
                if (graphs.pattern_graph_rows[p * max_graphs + ss].test(q))
                    emit_for.emplace_back(gg, ss);
            if (emit_for.empty())
                continue;

            // Eliding: the highest index has the smallest target set and so subsumes every
            // lower index's constraint for this head, so emit only it and record (p,q) as
            // covered (distance3, which is wider still, then skips it). Otherwise emit them
            // all (the index list is ascending, so the last entry is the highest index).
            if (elide_subsumed) {
                covered.emplace(p, q);
                emit_for = {emit_for.back()};
                // remember which slot's constraint we kept, so a degree/NDS check on a lower
                // (elided) exact-path graph or on distance3 can weaken from it on demand.
                _kept_supplemental_slot[pair{int(p), int(q)}] = emit_for.front().second;
            }

            for (auto & [g, slot] : emit_for) {
                auto n_p_q = graphs.pattern_graph_rows[p * max_graphs + 0];
                n_p_q &= graphs.pattern_graph_rows[q * max_graphs + 0];
                vector<int> between_p_and_q;
                for (auto v = n_p_q.find_first(); v != decltype(n_p_q)::npos; v = n_p_q.find_first()) {
                    n_p_q.reset(v);
                    between_p_and_q.push_back(int(v));
                    if (between_p_and_q.size() >= unsigned(g))
                        break;
                }

                for (unsigned t = 0; t < target_size; ++t) {
                    vector<int> n_t, d_n_t;
                    vector<pair<int, vector<int>>> two_away_from_t;
                    auto n_t_row = graphs.target_graph_rows[t * max_graphs + 0];
                    for (auto w = n_t_row.find_first(); w != decltype(n_t_row)::npos; w = n_t_row.find_first()) {
                        n_t_row.reset(w);
                        n_t.push_back(int(w));
                    }

                    auto nd_t = graphs.target_graph_rows[t * max_graphs + slot];
                    for (auto w = nd_t.find_first(); w != decltype(nd_t)::npos; w = nd_t.find_first()) {
                        nd_t.reset(w);
                        d_n_t.push_back(int(w));
                    }

                    auto n2_t = graphs.target_graph_rows[t * max_graphs + exact_path_1_slot];
                    for (auto w = n2_t.find_first(); w != decltype(n2_t)::npos; w = n2_t.find_first()) {
                        n2_t.reset(w);
                        auto n_t_w = graphs.target_graph_rows[w * max_graphs + 0];
                        n_t_w &= graphs.target_graph_rows[t * max_graphs + 0];
                        vector<int> n_t_w_idx;
                        for (auto x = n_t_w.find_first(); x != decltype(n_t_w)::npos; x = n_t_w.find_first()) {
                            n_t_w.reset(x);
                            n_t_w_idx.push_back(int(x));
                        }
                        two_away_from_t.emplace_back(int(w), n_t_w_idx);
                    }

                    int slot_i = int(slot), p_i = int(p), q_i = int(q), t_i = int(t);
                    register_supplemental(std::tuple<long, long, long, long>{slot_i, p_i, q_i, t_i}, p_i, t_i,
                        [this, slot_i, p_i, q_i, t_i, between_p_and_q, n_t, two_away_from_t, d_n_t]() {
                            emit_exact_path_graph(slot_i, p_i, q_i, between_p_and_q, t_i, n_t, two_away_from_t, d_n_t);
                        });
                }
            }
        }
    }

    return covered;
}

auto HomomorphismProofs::emit_distance3_graph_distance_1(int g, int p, int q, int t,
    const std::vector<int> & d3_from_t) -> void
{
    auto & adjacency = _adjacency;
    _proof->emit_proof_directive("% adjacency " + _pattern_names[p] + " maps to " + _target_names[t] +
        " in G^3 so by adjacency, " + _pattern_names[q] + " maps to one of...");

    // single-line derivation: emit the @label at the top level so it persists across search
    // backtracking (the caller restores the active level for the batch).
    _proof->emit_proof_directive("setlvl 0;");
    std::string adj_label = "@d3adj" + _pattern_names[p] + "_" + _target_names[t] + "_" + _pattern_names[q];
    std::string line = adj_label + " ia 1 ~x" + _proof->variable_name(p, t);
    for (auto & u : d3_from_t)
        line += " 1 x" + _proof->variable_name(q, u);
    line += " >= 1 : " + adjacency.labels.at(std::tuple<long, long, long, long>{0, p, q, t}) + " ;";
    _proof->emit_proof_line(line);

    adjacency.labels.emplace(std::tuple<long, long, long, long>{g, p, q, t}, adj_label);
}

auto HomomorphismProofs::emit_distance3_graph_distance_2(int g, int p, int q, int path1, int t,
    const std::vector<int> & d1_from_t, const std::vector<int> & d2_from_t,
    const std::vector<int> & d3_from_t) -> void
{
    auto & adjacency = _adjacency;
    _proof->emit_proof_directive("% adjacency " + _pattern_names[p] + " maps to " + _target_names[t] +
        " in G^3 so using vertex " + _pattern_names[path1] + ", " + _pattern_names[q] + " maps to one of...");

    // scratch one level above the search (see emit_exact_path_graph).
    int scratch = _proof->active_level() + 1;
    _proof->emit_proof_directive("setlvl " + std::to_string(scratch) + ";");

    // if p maps to t then the first thing on the path from p to q has to go to one of, so the
    // second thing on the path from p to q has to go to one of...
    std::string pol = "pol " + adjacency.labels.at(std::tuple<long, long, long, long>{0, p, path1, t});
    for (auto & u : d1_from_t)
        pol += " " + adjacency.labels.at(std::tuple<long, long, long, long>{0, path1, q, u}) + " +";
    pol += " ;";
    _proof->emit_proof_line(pol);

    // tidy up
    std::string ia = "ia 1 ~x" + _proof->variable_name(p, t);
    for (auto & u : d2_from_t)
        ia += " 1 x" + _proof->variable_name(q, u);
    ia += " >= 1 : " + std::to_string(_proof->current_proof_line()) + " ;";
    _proof->emit_proof_line(ia);

    _proof->emit_proof_directive("setlvl 0;");

    std::string adj_label = "@d3adj" + _pattern_names[p] + "_" + _target_names[t] + "_" + _pattern_names[q];
    std::string line = adj_label + " ia 1 ~x" + _proof->variable_name(p, t);
    for (auto & u : d3_from_t)
        line += " 1 x" + _proof->variable_name(q, u);
    line += " >= 1 : " + std::to_string(_proof->current_proof_line()) + " ;";
    _proof->emit_proof_line(line);

    adjacency.labels.emplace(std::tuple<long, long, long, long>{g, p, q, t}, adj_label);
    // self-clean the scratch (lazy ordering means we can't rely on a later supplemental's
    // wiplvl); the @label is at level 0, and the caller restores the search's active level.
    _proof->wipe_level(scratch);
}

auto HomomorphismProofs::emit_distance3_graph(int g, int p, int q, int path1, int path2, int t,
    const std::vector<int> & d1_from_t, const std::vector<int> & d2_from_t,
    const std::vector<int> & d3_from_t) -> void
{
    auto & adjacency = _adjacency;
    _proof->emit_proof_directive("% adjacency " + _pattern_names[p] + " maps to " + _target_names[t] +
        " in G^3 so using path " + _pattern_names[path1] + " -- " + _pattern_names[path2] + ", " +
        _pattern_names[q] + " maps to one of...");

    // scratch one level above the search (see emit_exact_path_graph).
    int scratch = _proof->active_level() + 1;
    _proof->emit_proof_directive("setlvl " + std::to_string(scratch) + ";");

    // if p maps to t then the first thing on the path from p to q has to go to one of, so the
    // second thing on the path from p to q has to go to one of...
    std::string pol = "pol " + adjacency.labels.at(std::tuple<long, long, long, long>{0, p, path1, t});
    for (auto & u : d1_from_t)
        pol += " " + adjacency.labels.at(std::tuple<long, long, long, long>{0, path1, path2, u}) + " +";
    pol += " ;";
    _proof->emit_proof_line(pol);

    // tidy up
    std::string ia = "ia 1 ~x" + _proof->variable_name(p, t);
    for (auto & u : d2_from_t)
        ia += " 1 x" + _proof->variable_name(path2, u);
    ia += " >= 1 : " + std::to_string(_proof->current_proof_line()) + " ;";
    _proof->emit_proof_line(ia);

    std::string pol2 = "pol " + std::to_string(_proof->current_proof_line());
    for (auto & u : d2_from_t)
        pol2 += " " + adjacency.labels.at(std::tuple<long, long, long, long>{0, path2, q, u}) + " s +";
    pol2 += " ;";
    _proof->emit_proof_line(pol2);

    _proof->emit_proof_directive("setlvl 0;");

    std::string adj_label = "@d3adj" + _pattern_names[p] + "_" + _target_names[t] + "_" + _pattern_names[q];
    std::string line = adj_label + " ia 1 ~x" + _proof->variable_name(p, t);
    for (auto & u : d3_from_t)
        line += " 1 x" + _proof->variable_name(q, u);
    line += " >= 1 : " + std::to_string(_proof->current_proof_line()) + " ;";
    _proof->emit_proof_line(line);

    adjacency.labels.emplace(std::tuple<long, long, long, long>{g, p, q, t}, adj_label);
    // self-clean the scratch (lazy ordering means we can't rely on a later supplemental's
    // wiplvl); the @label is at level 0, and the caller restores the search's active level.
    _proof->wipe_level(scratch);
}

auto HomomorphismProofs::prove_distance3_graphs(const ProcessedGraphsData & graphs, unsigned max_graphs, unsigned slot,
    const std::set<std::pair<int, int>> & covered_by_exact_path) -> void
{
    const unsigned pattern_size = _pattern_names.size();
    const unsigned target_size = _target_names.size();

    for (unsigned p = 0; p < pattern_size; ++p) {
        for (unsigned q = 0; q < pattern_size; ++q) {
            // only do this if they're actually adjacent
            if (p == q || ! graphs.pattern_graph_rows[p * max_graphs + slot].test(q))
                continue;

            // a distance-3 constraint's target set (within distance 3 of t) is a superset of
            // any exact-path set (within two paths of t) for the same head, so if exact-path
            // already covered (p,q) the distance-3 form is subsumed -- skip it.
            if (covered_by_exact_path.contains(pair{int(p), int(q)}))
                continue;

            // this (p,q) is not exact-path-covered, so the distance-3 graph holds its kept
            // (strongest) constraint -- record the slot for on-demand weakening if needed.
            _kept_supplemental_slot[pair{int(p), int(q)}] = slot;

            bool actually_adjacent = false;
            optional<int> path_from_p_to_q_1 = nullopt, path_from_p_to_q_2 = nullopt;

            auto n_p = graphs.pattern_graph_rows[p * max_graphs + 0];

            // are they actually distance 1 apart?
            if (n_p.test(q))
                actually_adjacent = true;
            else {
                auto n_q = graphs.pattern_graph_rows[q * max_graphs + 0];

                auto n_p_q = n_p;
                n_p_q &= n_q;
                n_p_q.reset(p);
                n_p_q.reset(q);

                if (n_p_q.any()) {
                    // they're actually distance 2 apart
                    path_from_p_to_q_1 = int(n_p_q.find_first());
                }
                else {
                    // find a path of length 3
                    n_p.reset(p);
                    n_p.reset(q);
                    for (auto v = n_p.find_first(); v != decltype(n_p)::npos && ! path_from_p_to_q_1; v = n_p.find_first()) {
                        n_p.reset(v);
                        auto n_v = graphs.pattern_graph_rows[v * max_graphs + 0];
                        n_v.reset(v);
                        n_v.reset(p);
                        n_v.reset(q);
                        for (auto w = n_v.find_first(); w != decltype(n_v)::npos && ! path_from_p_to_q_1; w = n_v.find_first()) {
                            n_v.reset(w);
                            if (graphs.pattern_graph_rows[w * max_graphs + 0].test(q)) {
                                path_from_p_to_q_1 = int(v);
                                path_from_p_to_q_2 = int(w);
                            }
                        }
                    }
                }

                if (! path_from_p_to_q_1)
                    throw ProofError{"Oops, there's a bug: missing path from " + _pattern_names[p] + " to " + _pattern_names[q]};
            }

            for (unsigned t = 0; t < target_size; ++t) {
                vector<int> d1_from_t, d2_from_t, d3_from_t;
                set<int> d2_from_t_set, d3_from_t_set;
                auto n_t = graphs.target_graph_rows[t * max_graphs + 0];
                n_t.set(t);
                for (auto v = n_t.find_first(); v != decltype(n_t)::npos; v = n_t.find_first()) {
                    n_t.reset(v);
                    d1_from_t.push_back(int(v));
                    auto n_v = graphs.target_graph_rows[v * max_graphs + 0];
                    n_v.set(v);
                    for (auto w = n_v.find_first(); w != decltype(n_v)::npos; w = n_v.find_first()) {
                        n_v.reset(w);
                        d2_from_t_set.insert(int(w));
                        auto n_w = graphs.target_graph_rows[w * max_graphs + 0];
                        n_w.set(w);
                        for (auto x = n_w.find_first(); x != decltype(n_w)::npos; x = n_w.find_first()) {
                            n_w.reset(x);
                            d3_from_t_set.insert(int(x));
                        }
                    }
                }

                d2_from_t.assign(d2_from_t_set.begin(), d2_from_t_set.end());
                d3_from_t.assign(d3_from_t_set.begin(), d3_from_t_set.end());

                int slot_i = int(slot), p_i = int(p), q_i = int(q), t_i = int(t);
                bool adj = actually_adjacent;
                int path1 = path_from_p_to_q_1.value_or(-1);
                optional<int> path2 = path_from_p_to_q_2;
                register_supplemental(std::tuple<long, long, long, long>{slot_i, p_i, q_i, t_i}, p_i, t_i,
                    [this, slot_i, p_i, q_i, t_i, adj, path1, path2, d1_from_t, d2_from_t, d3_from_t]() {
                        if (adj)
                            emit_distance3_graph_distance_1(slot_i, p_i, q_i, t_i, d3_from_t);
                        else if (path2)
                            emit_distance3_graph(slot_i, p_i, q_i, path1, *path2, t_i, d1_from_t, d2_from_t, d3_from_t);
                        else
                            emit_distance3_graph_distance_2(slot_i, p_i, q_i, path1, t_i, d1_from_t, d2_from_t, d3_from_t);
                    });
            }
        }
    }
}

auto HomomorphismProofs::emit_shape_graph(int g, int p, int q, int t, const std::vector<int> & n_t) -> void
{
    _proof->emit_proof_directive("% adjacency " + _pattern_names[p] + " maps to " + _target_names[t] +
        " in shape graph " + std::to_string(g) + " so " + _pattern_names[q] + " maps to one of...");
    // single-line assertion: emit the @label at the top level so it persists across search
    // backtracking (the caller restores the active level for the batch).
    _proof->emit_proof_directive("setlvl 0;");
    std::string adj_label = "@g" + std::to_string(g) + "adj" + _pattern_names[p] + "_" + _target_names[t] + "_" + _pattern_names[q];
    std::string line = adj_label + " a 1 ~x" + _proof->variable_name(p, t);
    for (auto & u : n_t)
        line += " 1 x" + _proof->variable_name(q, u);
    line += " >= 1 ;";
    _proof->emit_proof_line(line);

    _adjacency.labels.emplace(std::tuple<long, long, long, long>{g, p, q, t}, adj_label);
}

auto HomomorphismProofs::prove_extra_shape(const ProcessedGraphsData & graphs, unsigned max_graphs, unsigned slot) -> void
{
    const unsigned pattern_size = _pattern_names.size();
    const unsigned target_size = _target_names.size();

    for (unsigned p = 0; p < pattern_size; ++p) {
        for (unsigned q = 0; q < pattern_size; ++q) {
            // only do this if they're actually adjacent
            if (! graphs.pattern_graph_rows[p * max_graphs + slot].test(q))
                continue;

            for (unsigned t = 0; t < target_size; ++t) {
                vector<int> n_t;
                auto n_t_row = graphs.target_graph_rows[t * max_graphs + slot];
                for (auto v = n_t_row.find_first(); v != decltype(n_t_row)::npos; v = n_t_row.find_first()) {
                    n_t_row.reset(v);
                    n_t.push_back(int(v));
                }
                int slot_i = int(slot), p_i = int(p), q_i = int(q), t_i = int(t);
                register_supplemental(std::tuple<long, long, long, long>{slot_i, p_i, q_i, t_i}, p_i, t_i,
                    [this, slot_i, p_i, q_i, t_i, n_t]() {
                        emit_shape_graph(slot_i, p_i, q_i, t_i, n_t);
                    });
            }
        }
    }
}

auto HomomorphismProofs::emit_model(const InputGraph & pattern, const InputGraph & target, const HomomorphismParams & params) -> void
{
    // set up our model file, with a set of OPB variables for each CP variable
    for (int n = 0; n < pattern.size(); ++n) {
        _proof->create_cp_variable(
            n, target.size(),
            [&](int v) { return pattern.vertex_name(v); },
            [&](int v) { return target.vertex_name(v); });
    }

    // generate constraints for injectivity
    if (params.injectivity == Injectivity::Injective)
        _proof->create_injectivity_constraints(pattern.size(), target.size(),
            [&](int v) { return target.vertex_name(v); });
    else if (params.injectivity == Injectivity::LocallyInjective)
        // local injectivity: for each pattern vertex and each target, at most one of
        // that vertex's neighbours may map there (so phi restricted to a neighbourhood
        // is injective). The neighbourhood analogue of the injectivity constraints.
        create_locally_injective_constraints(pattern, target);

    // generate edge constraints, and also handle loops here
    for (int p = 0; p < pattern.size(); ++p) {
        for (int t = 0; t < target.size(); ++t) {
            // it's simpler to always have the adjacency constraints, even
            // if the assignment is forbidden
            _proof->emit_model_comment("* adjacency " + std::to_string(p) + " maps to " + std::to_string(t));

            // if p can be mapped to t, then each neighbour of p...
            for (int q = 0; q < pattern.size(); ++q)
                if (pattern.adjacent(p, q)) {
                    // ... must be mapped to a neighbour of t. A target self-loop
                    // (u == t) is kept in the sum, so a loop-preserving mapping
                    // satisfies the constraint (this matches the verified CakePB
                    // encoding; see issue #49).
                    vector<int> permitted;
                    for (int u = 0; u < target.size(); ++u)
                        if (target.adjacent(t, u))
                            permitted.push_back(u);
                    emit_adjacency_constraint(p, q, t, permitted);
                }

            // same for non-adjacency for induced
            if (params.induced) {
                for (int q = 0; q < pattern.size(); ++q)
                    if (q != p && ! pattern.adjacent(p, q)) {
                        // ... must be mapped to a non-neighbour of t. t itself counts as
                        // a non-neighbour exactly when it has no self-loop, so the
                        // permitted set is just the non-neighbours of t (the same test
                        // the q == p case below uses). Under full injectivity q cannot
                        // share t with p anyway, so whether t is in the set is moot; but
                        // under local injectivity p and q may both map to a loopless t,
                        // and that is a legitimate induced non-edge (t is not adjacent to
                        // itself), so t must stay in the set or the model wrongly rejects
                        // it.
                        vector<int> permitted;
                        for (int u = 0; u < target.size(); ++u)
                            if (! target.adjacent(t, u))
                                permitted.push_back(u);
                        emit_adjacency_constraint(p, q, t, permitted);
                    }

                // the q == p case of non-edge preservation: a non-loopy pattern
                // vertex cannot map to a loopy target, since induced isomorphism
                // requires loop(p) == loop(t). The loop above skips q == p, and the
                // edge loop only constrains a loopy p, so without this the model
                // admits invalid induced mappings (issue #56). p -> t then forces p
                // onto a non-neighbour of t; with t a neighbour of itself and
                // at-most-one-value, that is a contradiction.
                if (! pattern.adjacent(p, p) && target.adjacent(t, t)) {
                    vector<int> permitted;
                    for (int u = 0; u < target.size(); ++u)
                        if (! target.adjacent(t, u))
                            permitted.push_back(u);
                    emit_adjacency_constraint(p, p, t, permitted);
                }
            }
        }
    }

    // declare the projected set (the assignment variables) so the proof's
    // solution count is in terms of the high-level mapping
    _proof->emit_preserved_assignment_variables();

    // output the model file
    _proof->finalise_model();
}

auto HomomorphismProofs::emit_reified_model(const InputGraph & pattern, const InputGraph & target, const HomomorphismParams & params) -> void
{
    // Names by index, since vertex names may contain anything.
    auto pattern_name = [](int v) { return "p" + std::to_string(v); };
    auto target_name = [](int v) { return "t" + std::to_string(v); };

    // As reification and the model do: a pattern without labels of a kind ignores the
    // target's, and if either graph is directed, an undirected edge is two arcs.
    bool use_vertex_labels = pattern.has_vertex_labels();
    bool use_edge_labels = pattern.has_edge_labels();
    bool directed = pattern.directed() || target.directed();

    // The target's edges, by ordered endpoint pair, then label, with their costs.
    map<pair<int, int>, map<string, long long>> target_edges;
    target.for_each_edge_and_cost([&](int f, int t, std::string_view l, optional<long long> c) {
        target_edges[{f, t}][string{l}] = c.value_or(0);
    });

    // The cheapest target edge from x to y that a pattern edge with this label may use.
    auto edge_cost = [&](int x, int y, const string & label) -> optional<long long> {
        auto it = target_edges.find({x, y});
        if (it == target_edges.end())
            return nullopt;
        if (use_edge_labels) {
            auto l = it->second.find(label);
            if (l == it->second.end())
                return nullopt;
            return l->second;
        }
        optional<long long> best;
        for (auto & [_, c] : it->second)
            if ((! best) || c < *best)
                best = c;
        return best;
    };

    vector<vector<int>> values(pattern.size());
    for (int p = 0; p < pattern.size(); ++p) {
        for (int t = 0; t < target.size(); ++t)
            if ((! use_vertex_labels) || pattern.vertex_label(p) == target.vertex_label(t))
                values[p].push_back(t);
        _proof->create_cp_variable(p, values[p], pattern_name, target_name);
    }

    _proof->create_injectivity_constraints(pattern.size(), target.size(), target_name);

    vector<pair<string, long long>> objective;
    map<pair<int, int>, long long> unary_cost;
    if (target.has_vertex_costs())
        for (int p = 0; p < pattern.size(); ++p)
            for (auto t : values[p])
                unary_cost[{p, t}] += target.vertex_cost(t);

    // The pattern's edges: a loop is a requirement on one vertex's image, and every other
    // edge belongs to the pair of its endpoints.
    struct PatternEdge
    {
        int from, to;
        string label;
    };
    map<pair<int, int>, vector<PatternEdge>> pairs;
    vector<PatternEdge> loops;
    pattern.for_each_edge([&](int f, int t, std::string_view l) {
        if ((! directed) && t < f)
            return;
        if (f == t)
            loops.push_back(PatternEdge{f, t, string{l}});
        else
            pairs[{std::min(f, t), std::max(f, t)}].push_back(PatternEdge{f, t, string{l}});
    });

    for (auto & loop : loops)
        for (auto t : values[loop.from]) {
            auto c = edge_cost(t, t, loop.label);
            if (! c) {
                _proof->emit_model_comment("* no loop " + loop.label + " on " + target_name(t));
                _proof->emit_model_constraint("1 ~x" + _proof->variable_name(loop.from, t) + " >= 1 ;");
            }
            else
                unary_cost[{loop.from, t}] += *c;
        }

    for (auto & [key, c] : unary_cost)
        if (c != 0)
            objective.emplace_back("x" + _proof->variable_name(key.first, key.second), c);

    long extra_variables = 0;
    for (auto & [ab, edges] : pairs) {
        auto [a, b] = ab;
        _proof->emit_model_comment("* pair " + pattern_name(a) + " " + pattern_name(b));

        // z(a, b, x, y) for every x, y whose images carry all of the pair's edges.
        map<int, vector<string>> by_x, by_y;
        for (auto x : values[a])
            for (auto y : values[b]) {
                if (x == y)
                    continue;
                long long cost = 0;
                bool ok = true;
                for (auto & e : edges) {
                    auto c = (e.from == a) ? edge_cost(x, y, e.label) : edge_cost(y, x, e.label);
                    if (! c) {
                        ok = false;
                        break;
                    }
                    cost += *c;
                }
                if (! ok)
                    continue;

                auto z = "z" + pattern_name(a) + "_" + pattern_name(b) + "_" + target_name(x) + "_" + target_name(y);
                ++extra_variables;
                by_x[x].push_back(z);
                by_y[y].push_back(z);
                if (cost != 0)
                    objective.emplace_back(z, cost);
            }

        // sum_y z(a, b, x, y) = x(a, x) and sum_x z(a, b, x, y) = x(b, y), each as two
        // inequalities with labels, for derivations to cite.
        auto link = [&](int p, int t, const vector<string> & zs, const string & side) {
            auto label = "@lnk" + pattern_name(a) + "_" + pattern_name(b) + "_" + side + target_name(t);
            string ge = label + "ge", le = label + "le";
            for (auto & z : zs) {
                ge += " 1 " + z;
                le += " -1 " + z;
            }
            ge += " -1 x" + _proof->variable_name(p, t) + " >= 0 ;";
            le += " 1 x" + _proof->variable_name(p, t) + " >= 0 ;";
            _proof->emit_model_constraint(ge);
            _proof->emit_model_constraint(le);
        };
        for (auto x : values[a])
            link(a, x, by_x[x], "a");
        for (auto y : values[b])
            link(b, y, by_y[y], "b");
    }

    _proof->declare_extra_variables(extra_variables);
    if (params.minimise_cost)
        _proof->create_weighted_objective(objective);

    _proof->emit_preserved_assignment_variables();
    _proof->finalise_model();
}

auto HomomorphismProofs::cost_bound(const std::vector<std::pair<int, int>> & decisions, const CostBoundCertificate & certificate, bool failed) -> void
{
    _proof->emit_proof_directive("% cost bound " + std::string(failed ? "failed" : "removed values") + " at depth " + std::to_string(decisions.size()));

    switch (certificate.kind) {
    case CostBoundCertificate::Kind::None:
        // Anything done for lack of support follows by propagation over the linking
        // equalities.
        break;

    case CostBoundCertificate::Kind::HallViolator: {
        vector<int> rows(certificate.hall_rows.begin(), certificate.hall_rows.end());
        vector<int> columns(certificate.hall_columns.begin(), certificate.hall_columns.end());
        emit_hall_set_or_violator(rows, columns);
        break;
    }

    case CostBoundCertificate::Kind::Bound: {
        // The objective-improving constraint, plus each certificate constraint times its
        // multiplier, using the other half of an equality for a negative one.
        std::string pol = "pol " + std::to_string(_proof->objective_line());
        auto term = [&](const std::string & label, long long multiplier) {
            pol += " " + label;
            if (multiplier != 1)
                pol += " " + std::to_string(multiplier) + " *";
            pol += " +";
        };
        for (auto & [p, m] : certificate.exactly_one)
            if (m > 0)
                term(_proof->at_least_one_value_label(p), m);
            else
                term(_proof->at_most_one_value_label(p), -m);
        for (auto & [t, m] : certificate.at_most_one)
            term(_proof->injectivity_label(t), m);
        for (auto & l : certificate.links) {
            auto label = "@lnkp" + std::to_string(l.a) + "_p" + std::to_string(l.b) + "_" + (l.on_a ? "a" : "b") + "t" + std::to_string(l.value);
            if (l.multiplier > 0)
                term(label + "ge", l.multiplier);
            else
                term(label + "le", -l.multiplier);
        }
        pol += " ;";
        _proof->emit_proof_line(pol);
        break;
    }
    }

    // Nothing more is needed. At this node the derived constraint conflicts, or
    // propagates exactly the values the bound removed, so the search's own nogoods follow
    // from it by RUP, as do later steps relying on the removals. (Checking each removal
    // with its own RUP here found no failures, and cost only size.)
}

auto HomomorphismProofs::derive_loop_fixed_adjacencies() -> void
{
    // Derive the loop-cancelled form of each loopy adjacency constraint, so the degree,
    // supplemental-graph and distance-3 derivations can sum them into pols without a stray
    // "maps to the loopy target" term (issue #56).
    //
    // This is deferred out of emit_model (the OPB stays complete and is emitted up front,
    // but these are PBP *derivations*) until the search step runs it, just before the first
    // @adj-citing pol. A cheap concluding step (pattern-too-big, target-loop, clique) then
    // pays nothing for it: its refutation cites injectivity, never adjacency. For any
    // instance that does reach search nothing is emitted to the proof in between, so the
    // proof is byte-identical -- this only removes the derivations on an early conclusion
    // (e.g. the induced pattern-bigger-than-target case, where they were all dead).
    //
    // Mechanism: before issue #49 the adjacency constraint left the target's self-loop term
    // out, so it could be summed into a pol cleanly. We now keep that term (so a loop->loop
    // mapping satisfies the model), but it then appears as a stray "q maps to the loopy
    // target" term in every pol. For each adjacency constraint over a loopy target t, rewrite
    // its @adj label to the loop-cancelled version -- ~x_p_t together with the neighbours of t
    // other than t -- which follows from the constraint plus injectivity on t (mapping p to t
    // forbids q from also mapping to t). VeriPB lets a proof line reassign an existing label,
    // so every later @adj reference in a pol picks up the loop-cancelled form; the original
    // loop-bearing constraint stays in the database (by number) so solutions still satisfy the
    // model. The loop-cancelled form relies on global injectivity on t, which local injectivity
    // does not give -- but under local injectivity the pol-summing filters that would need it
    // are disabled anyway (issue #58), so skip the relabelling entirely.
    if (is_locally_injective())
        return;

    auto & adjacency = _adjacency;
    for (auto & [key, label] : adjacency.labels) {
        auto & [g, p, q, t] = key;
        // a pattern self-loop edge (p == q) has its loop term pinned by at-most-one rather
        // than injectivity, and is not summed into the supplemental/degree pols; skip it.
        if (p == q)
            continue;
        auto pit = adjacency.permitted.find(key);
        if (pit == adjacency.permitted.end())
            continue;
        if (std::find(pit->second.begin(), pit->second.end(), t) == pit->second.end())
            continue; // t is not a neighbour of itself: no loop term to cancel
        std::string line = label + " rup 1 ~x" + _proof->variable_name(p, t);
        for (auto & u : pit->second)
            if (u != t)
                line += " 1 x" + _proof->variable_name(q, u);
        line += " >= 1 ;";
        _proof->emit_proof_line(line);
        // (The original loop-bearing constraint is now redundant, but it is left in place:
        // deleting it needs `del id <number>`, and that number does not correspond to the
        // same constraint in CakePB's independently-numbered OPB, so the deletion breaks the
        // verified-pipeline elaboration. The extra constraint is cheap.)
    }
}

auto HomomorphismProofs::prove_no_clique(const ProcessedGraphsData & graphs, unsigned max_graphs, unsigned pattern_size,
    unsigned target_size, const HomomorphismParams & params, unsigned g, int p, int tt) -> void
{
    vector<NamedVertex> p_clique;
    map<int, NamedVertex> t_clique_neighbourhood;
    unsigned decide_size;

    {
        vector<int> include(pattern_size, -1), invinclude(pattern_size, 0);
        int count = 0;
        for (int w = 0; w < int(pattern_size); ++w)
            if (w != p && graphs.pattern_graph_rows[w * max_graphs + g].test(p)) {
                include[w] = count;
                invinclude[count] = w;
                ++count;
            }

        InputGraph gv(count, false, false);
        for (unsigned f = 0; f < pattern_size; ++f)
            if (include[f] != -1)
                for (unsigned t = 0; t < pattern_size; ++t)
                    if (f != t && include[t] != -1 && graphs.pattern_graph_rows[f * max_graphs + g].test(t))
                        gv.add_edge(include[f], include[t]);

        CliqueParams clique_params;
        clique_params.timeout = params.timeout;
        clique_params.start_time = steady_clock::now();
        clique_params.restarts_schedule = make_unique<NoRestartsSchedule>();
        auto result = solve_clique_problem(gv, clique_params);
        for (auto & v : result.clique)
            p_clique.push_back(pattern_vertex(invinclude[v]));
        decide_size = result.clique.size();
    }

    {
        vector<int> include(target_size, -1), invinclude(target_size, 0);
        int count = 0;
        for (int w = 0; w < int(target_size); ++w)
            if (w != tt && graphs.target_graph_rows[w * max_graphs + g].test(tt)) {
                t_clique_neighbourhood.emplace(count, target_vertex(w));
                include[w] = count;
                invinclude[count] = w;
                ++count;
            }

        _proof->prepare_hom_clique_proof(pattern_vertex(p), target_vertex(tt), decide_size);

        InputGraph gv(count, false, false);
        for (unsigned f = 0; f < target_size; ++f)
            if (include[f] != -1)
                for (unsigned t = 0; t < target_size; ++t) {
                    if (f != t && include[t] != -1) {
                        if (graphs.target_graph_rows[f * max_graphs + g].test(t))
                            gv.add_edge(include[f], include[t]);
                        else if (f < t)
                            _proof->add_hom_clique_non_edge(
                                pattern_vertex(p), target_vertex(tt),
                                p_clique, target_vertex(f), target_vertex(t));
                    }
                }

        _proof->start_hom_clique_proof(pattern_vertex(p), move(p_clique), target_vertex(tt), move(t_clique_neighbourhood));

        CliqueParams clique_params;
        clique_params.timeout = params.timeout;
        clique_params.start_time = steady_clock::now();
        clique_params.decide = make_optional(decide_size);
        clique_params.restarts_schedule = make_unique<NoRestartsSchedule>();
        clique_params.extend_proof = _proof;
        clique_params.proof_is_for_hom = true;

        auto result = solve_clique_problem(gv, clique_params);
        if (result.complete && ! result.clique.empty())
            throw ProofError{"Oops, found a clique that shound't exist"};
        _proof->finish_hom_clique_proof(pattern_vertex(p), target_vertex(tt), decide_size);
    }
}

auto HomomorphismProofs::ensure_supplemental_adjacency(const ProcessedGraphsData & graphs, unsigned max_graphs,
    int g, int p, int q, int t) -> void
{
    // present already (it is the kept constraint, the original graph, or elision is off): nothing to do.
    auto & adjacency = _adjacency;
    // with lazy emission the kept constraint for this head may still be pending; if (g,p,q,t)
    // is itself the kept one, materialising it makes it present and we are done.
    materialise_one(std::tuple<long, long, long, long>{g, p, q, t});
    if (adjacency.labels.contains(std::tuple<long, long, long, long>{g, p, q, t}))
        return;

    // no kept supplemental constraint for this head means we never emitted (or could emit) a
    // supplemental adjacency line here -- e.g. the original graph (g 0), or a distance-2 /
    // k4 graph, which carry no adjacency lines. The degree/NDS pigeonhole already tolerates a
    // missing term in those cases, so leave it missing (matches the no-elision behaviour).
    auto kept = _kept_supplemental_slot.find(pair{p, q});
    if (kept == _kept_supplemental_slot.end())
        return;

    // otherwise it was elided in favour of a stronger same-head one (in graph from_g, a
    // narrower target set). The wider graph-g constraint follows from the narrower by a single
    // implication step, so derive it that way, citing the kept constraint by its label.
    int from_g = int(kept->second);
    // the kept (stronger) constraint we weaken from may itself still be pending; emit it first.
    materialise_one(std::tuple<long, long, long, long>{from_g, p, q, t});
    auto from_label = adjacency.labels.at(std::tuple<long, long, long, long>{from_g, p, q, t});
    std::string adj_label = "@g" + std::to_string(g) + "adj" + _pattern_names[p] + "_" + _target_names[t] + "_" + _pattern_names[q];
    std::string line = adj_label + " ia 1 ~x" + _proof->variable_name(p, t);
    auto row = graphs.target_graph_rows[t * max_graphs + g];
    for (auto u = row.find_first(); u != decltype(row)::npos; u = row.find_first()) {
        row.reset(u);
        if (int(u) != t)
            line += " 1 x" + _proof->variable_name(q, int(u));
    }
    line += " >= 1 : " + from_label + " ;";
    auto id = _proof->emit_proof_line(line);
    adjacency.labels.emplace(std::tuple<long, long, long, long>{g, p, q, t}, adj_label);
    adjacency.ids.emplace(std::tuple<long, long, long, long>{g, p, q, t}, id);
    _pending_transient_adjacencies.emplace_back(g, p, q, t);
}

auto HomomorphismProofs::forget_transient_supplemental_adjacencies() -> void
{
    auto & adjacency = _adjacency;
    for (auto & [g, p, q, t] : _pending_transient_adjacencies) {
        std::tuple<long, long, long, long> key{g, p, q, t};
        _proof->emit_proof_directive("del id " + std::to_string(adjacency.ids.at(key)) + " ;");
        adjacency.labels.erase(key);
        adjacency.ids.erase(key);
    }
    _pending_transient_adjacencies.clear();
}
