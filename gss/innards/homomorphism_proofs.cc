#include <gss/clique.hh>
#include <gss/innards/homomorphism_proofs.hh>

#include <algorithm>
#include <chrono>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
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

auto HomomorphismProofs::emit_exact_path_graph(int g, int p, int q, const std::vector<int> & between_p_and_q,
    int t, const std::vector<int> & n_t, const std::vector<std::pair<int, std::vector<int>>> & two_away_from_t,
    const std::vector<int> & d_n_t) -> void
{
    auto & adjacency = _proof->adjacency_proof_lines();

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

    _proof->emit_proof_directive("setlvl 1;");

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
        const std::string & inj = _proof->is_locally_injective()
            ? _proof->locally_injective_label(between_p_and_q.front(), t)
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
            pol2 += " " + (_proof->is_locally_injective() ? _proof->locally_injective_label(p, z) : _proof->injectivity_label(z)) + " +";
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
    _proof->emit_proof_directive("wiplvl 1;");
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

                    emit_exact_path_graph(int(slot), int(p), int(q), between_p_and_q,
                        int(t), n_t, two_away_from_t, d_n_t);
                }
            }
        }
    }

    return covered;
}

auto HomomorphismProofs::emit_distance3_graph_distance_1(int g, int p, int q, int t,
    const std::vector<int> & d3_from_t) -> void
{
    auto & adjacency = _proof->adjacency_proof_lines();
    _proof->emit_proof_directive("% adjacency " + _pattern_names[p] + " maps to " + _target_names[t] +
        " in G^3 so by adjacency, " + _pattern_names[q] + " maps to one of...");

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
    auto & adjacency = _proof->adjacency_proof_lines();
    _proof->emit_proof_directive("% adjacency " + _pattern_names[p] + " maps to " + _target_names[t] +
        " in G^3 so using vertex " + _pattern_names[path1] + ", " + _pattern_names[q] + " maps to one of...");

    _proof->emit_proof_directive("setlvl 1;");

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
}

auto HomomorphismProofs::emit_distance3_graph(int g, int p, int q, int path1, int path2, int t,
    const std::vector<int> & d1_from_t, const std::vector<int> & d2_from_t,
    const std::vector<int> & d3_from_t) -> void
{
    auto & adjacency = _proof->adjacency_proof_lines();
    _proof->emit_proof_directive("% adjacency " + _pattern_names[p] + " maps to " + _target_names[t] +
        " in G^3 so using path " + _pattern_names[path1] + " -- " + _pattern_names[path2] + ", " +
        _pattern_names[q] + " maps to one of...");

    _proof->emit_proof_directive("setlvl 1;");

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

                if (actually_adjacent)
                    emit_distance3_graph_distance_1(int(slot), int(p), int(q), int(t), d3_from_t);
                else if (path_from_p_to_q_2)
                    emit_distance3_graph(int(slot), int(p), int(q), *path_from_p_to_q_1,
                        *path_from_p_to_q_2, int(t), d1_from_t, d2_from_t, d3_from_t);
                else
                    emit_distance3_graph_distance_2(int(slot), int(p), int(q), *path_from_p_to_q_1,
                        int(t), d1_from_t, d2_from_t, d3_from_t);
            }
        }
    }
}

auto HomomorphismProofs::prove_extra_shape(const ProcessedGraphsData & graphs, unsigned max_graphs, unsigned slot) -> void
{
    const unsigned pattern_size = _pattern_names.size();
    const unsigned target_size = _target_names.size();

    for (unsigned p = 0; p < pattern_size; ++p) {
        for (unsigned q = 0; q < pattern_size; ++q) {
            auto named_p = pattern_vertex(p);
            auto named_q = pattern_vertex(q);

            // only do this if they're actually adjacent
            if (! graphs.pattern_graph_rows[p * max_graphs + slot].test(q))
                continue;

            for (unsigned t = 0; t < target_size; ++t) {
                auto named_t = target_vertex(t);
                vector<NamedVertex> named_n_t;
                auto n_t = graphs.target_graph_rows[t * max_graphs + slot];
                for (auto v = n_t.find_first(); v != decltype(n_t)::npos; v = n_t.find_first()) {
                    n_t.reset(v);
                    named_n_t.push_back(target_vertex(v));
                }
                _proof->hack_in_shape_graph(slot, named_p, named_q, named_t, named_n_t);
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
        _proof->create_locally_injective_constraints(pattern.size(), target.size(),
            [&](int a, int b) { return pattern.adjacent(a, b); },
            [&](int v) { return pattern.vertex_name(v); },
            [&](int v) { return target.vertex_name(v); });

    // generate edge constraints, and also handle loops here
    for (int p = 0; p < pattern.size(); ++p) {
        for (int t = 0; t < target.size(); ++t) {
            // it's simpler to always have the adjacency constraints, even
            // if the assignment is forbidden
            _proof->start_adjacency_constraints_for(p, t);

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
                    _proof->create_adjacency_constraint(
                        NamedVertex{p, pattern.vertex_name(p)},
                        NamedVertex{q, pattern.vertex_name(q)},
                        NamedVertex{t, target.vertex_name(t)},
                        permitted, false);
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
                        _proof->create_adjacency_constraint(
                            NamedVertex{p, pattern.vertex_name(p)},
                            NamedVertex{q, pattern.vertex_name(q)},
                            NamedVertex{t, target.vertex_name(t)},
                            permitted, true);
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
                    _proof->create_adjacency_constraint(
                        NamedVertex{p, pattern.vertex_name(p)},
                        NamedVertex{p, pattern.vertex_name(p)},
                        NamedVertex{t, target.vertex_name(t)},
                        permitted, true);
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
    if (_proof->is_locally_injective())
        return;

    auto & adjacency = _proof->adjacency_proof_lines();
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
    auto & adjacency = _proof->adjacency_proof_lines();
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
    auto & adjacency = _proof->adjacency_proof_lines();
    for (auto & [g, p, q, t] : _pending_transient_adjacencies) {
        std::tuple<long, long, long, long> key{g, p, q, t};
        _proof->emit_proof_directive("del id " + std::to_string(adjacency.ids.at(key)) + " ;");
        adjacency.labels.erase(key);
        adjacency.ids.erase(key);
    }
    _pending_transient_adjacencies.clear();
}
