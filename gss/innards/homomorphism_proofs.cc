#include <gss/innards/homomorphism_proofs.hh>

#include <memory>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

using namespace gss;
using namespace gss::innards;

using std::nullopt;
using std::optional;
using std::pair;
using std::set;
using std::shared_ptr;
using std::string;
using std::vector;

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

auto HomomorphismProofs::prove_exact_path_graphs(const ProcessedGraphsData & graphs, unsigned max_graphs, int number_of_exact_path_graphs) -> void
{
    const unsigned pattern_size = _pattern_names.size();
    const unsigned target_size = _target_names.size();

    for (int g = 1; g <= number_of_exact_path_graphs; ++g) {
        for (unsigned p = 0; p < pattern_size; ++p) {
            for (unsigned q = 0; q < pattern_size; ++q) {
                if (p == q || ! graphs.pattern_graph_rows[p * max_graphs + g].test(q))
                    continue;

                auto named_p = pattern_vertex(p);
                auto named_q = pattern_vertex(q);

                auto n_p_q = graphs.pattern_graph_rows[p * max_graphs + 0];
                n_p_q &= graphs.pattern_graph_rows[q * max_graphs + 0];
                vector<NamedVertex> between_p_and_q;
                for (auto v = n_p_q.find_first(); v != decltype(n_p_q)::npos; v = n_p_q.find_first()) {
                    n_p_q.reset(v);
                    between_p_and_q.push_back(pattern_vertex(v));
                    if (between_p_and_q.size() >= unsigned(g))
                        break;
                }

                for (unsigned t = 0; t < target_size; ++t) {
                    auto named_t = target_vertex(t);

                    vector<NamedVertex> named_n_t, named_d_n_t;
                    vector<pair<NamedVertex, vector<NamedVertex>>> named_two_away_from_t;
                    auto n_t = graphs.target_graph_rows[t * max_graphs + 0];
                    for (auto w = n_t.find_first(); w != decltype(n_t)::npos; w = n_t.find_first()) {
                        n_t.reset(w);
                        named_n_t.push_back(target_vertex(w));
                    }

                    auto nd_t = graphs.target_graph_rows[t * max_graphs + g];
                    for (auto w = nd_t.find_first(); w != decltype(nd_t)::npos; w = nd_t.find_first()) {
                        nd_t.reset(w);
                        named_d_n_t.push_back(target_vertex(w));
                    }

                    auto n2_t = graphs.target_graph_rows[t * max_graphs + 1];
                    for (auto w = n2_t.find_first(); w != decltype(n2_t)::npos; w = n2_t.find_first()) {
                        n2_t.reset(w);
                        auto n_t_w = graphs.target_graph_rows[w * max_graphs + 0];
                        n_t_w &= graphs.target_graph_rows[t * max_graphs + 0];
                        vector<NamedVertex> named_n_t_w;
                        for (auto x = n_t_w.find_first(); x != decltype(n_t_w)::npos; x = n_t_w.find_first()) {
                            n_t_w.reset(x);
                            named_n_t_w.push_back(target_vertex(x));
                        }
                        named_two_away_from_t.emplace_back(target_vertex(w), named_n_t_w);
                    }

                    _proof->create_exact_path_graphs(g, named_p, named_q, between_p_and_q,
                        named_t, named_n_t, named_two_away_from_t, named_d_n_t);
                }
            }
        }
    }
}

auto HomomorphismProofs::prove_distance3_graphs(const ProcessedGraphsData & graphs, unsigned max_graphs, unsigned slot) -> void
{
    const unsigned pattern_size = _pattern_names.size();
    const unsigned target_size = _target_names.size();

    for (unsigned p = 0; p < pattern_size; ++p) {
        for (unsigned q = 0; q < pattern_size; ++q) {
            auto named_p = pattern_vertex(p);
            auto named_q = pattern_vertex(q);

            // only do this if they're actually adjacent
            if (p == q || ! graphs.pattern_graph_rows[p * max_graphs + slot].test(q))
                continue;

            bool actually_adjacent = false;
            optional<NamedVertex> path_from_p_to_q_1 = nullopt, path_from_p_to_q_2 = nullopt;

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
                    path_from_p_to_q_1 = pattern_vertex(n_p_q.find_first());
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
                                path_from_p_to_q_1 = pattern_vertex(v);
                                path_from_p_to_q_2 = pattern_vertex(w);
                            }
                        }
                    }
                }

                if (! path_from_p_to_q_1)
                    throw ProofError{"Oops, there's a bug: missing path from " + named_p.second + " to " + named_q.second};
            }

            for (unsigned t = 0; t < target_size; ++t) {
                auto named_t = target_vertex(t);

                vector<NamedVertex> d1_from_t, d2_from_t, d3_from_t;
                set<NamedVertex> d2_from_t_set, d3_from_t_set;
                auto n_t = graphs.target_graph_rows[t * max_graphs + 0];
                n_t.set(t);
                for (auto v = n_t.find_first(); v != decltype(n_t)::npos; v = n_t.find_first()) {
                    n_t.reset(v);
                    d1_from_t.push_back(target_vertex(v));
                    auto n_v = graphs.target_graph_rows[v * max_graphs + 0];
                    n_v.set(v);
                    for (auto w = n_v.find_first(); w != decltype(n_v)::npos; w = n_v.find_first()) {
                        n_v.reset(w);
                        d2_from_t_set.insert(target_vertex(w));
                        auto n_w = graphs.target_graph_rows[w * max_graphs + 0];
                        n_w.set(w);
                        for (auto x = n_w.find_first(); x != decltype(n_w)::npos; x = n_w.find_first()) {
                            n_w.reset(x);
                            d3_from_t_set.insert(target_vertex(x));
                        }
                    }
                }

                d2_from_t.assign(d2_from_t_set.begin(), d2_from_t_set.end());
                d3_from_t.assign(d3_from_t_set.begin(), d3_from_t_set.end());

                if (actually_adjacent)
                    _proof->create_distance3_graphs_but_actually_distance_1(slot, named_p, named_q, named_t, d3_from_t);
                else if (path_from_p_to_q_2)
                    _proof->create_distance3_graphs(slot, named_p, named_q, *path_from_p_to_q_1,
                        *path_from_p_to_q_2, named_t, d1_from_t, d2_from_t, d3_from_t);
                else
                    _proof->create_distance3_graphs_but_actually_distance_2(slot, named_p, named_q, *path_from_p_to_q_1,
                        named_t, d1_from_t, d2_from_t, d3_from_t);
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

    // derive the loop-cancelled form of each loopy adjacency constraint, so the
    // degree, supplemental-graph and distance-3 derivations can sum them into pols
    // without a stray "maps to the loopy target" term (issue #56).
    _proof->loop_fix_adjacencies();
}
