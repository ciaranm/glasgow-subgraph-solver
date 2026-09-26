#include <gss/configuration.hh>
#include <gss/formats/input_graph.hh>
#include <gss/homomorphism.hh>
#include <gss/innards/verify.hh>

#include <catch2/catch_test_macros.hpp>

#include <chrono>
#include <map>
#include <memory>
#include <optional>
#include <random>
#include <string>
#include <vector>

using namespace gss;
using namespace gss::innards;

using std::make_shared;
using std::make_unique;
using std::map;
using std::mt19937;
using std::optional;
using std::string;
using std::to_string;
using std::uniform_real_distribution;
using std::vector;

using std::chrono::operator""s;

namespace
{
    // The axes that change what reification and the cost bound do. Directedness is
    // per graph, since a directed target under an undirected pattern is exactly the case
    // where reification has to make both directed; loops are where reified edges get a
    // label of their own; a multigraph is what reification is for; and the three kinds
    // of cost exercise the unary-only path, the edge-vertex path, and both together.
    struct Family
    {
        bool pattern_directed, target_directed, loops, vertex_labels, edge_labels, multigraph;
        bool vertex_costs, edge_costs;

        [[nodiscard]] auto name() const -> string
        {
            string n = string(pattern_directed ? "directed" : "undirected") + " pattern, " + (target_directed ? "directed" : "undirected") + " target";
            n += loops ? ", loopy" : ", loopless";
            if (vertex_labels)
                n += ", vertex-labelled";
            if (edge_labels)
                n += ", edge-labelled";
            if (multigraph)
                n += ", multigraph";
            if (vertex_costs)
                n += ", vertex costs";
            if (edge_costs)
                n += ", edge costs";
            return n;
        }
    };

    auto all_families() -> vector<Family>
    {
        vector<Family> families;
        for (bool pattern_directed : {false, true})
            for (bool target_directed : {false, true})
                for (bool loops : {false, true})
                    for (bool vertex_labels : {false, true})
                        for (bool edge_labels : {false, true})
                            for (bool multigraph : {false, true})
                                for (int costs : {1, 2, 3}) {
                                    bool vertex_costs = costs & 1, edge_costs = costs & 2;
                                    // A multigraph without edge labels can have no parallel edges.
                                    if (multigraph && ! edge_labels)
                                        continue;
                                    // Mixed directedness means what reification says only when
                                    // there is a reification; without one the model takes its
                                    // directedness from the pattern, which is a documented gap
                                    // of its own and not what this tests.
                                    if (pattern_directed != target_directed && ! (multigraph || edge_costs))
                                        continue;
                                    families.push_back(Family{pattern_directed, target_directed, loops, vertex_labels, edge_labels, multigraph, vertex_costs, edge_costs});
                                }
        return families;
    }

    auto random_label(mt19937 & rng) -> string
    {
        return (rng() % 2) ? "x" : "y";
    }

    // Costs include negative ones, which nothing about the bound should mind.
    auto random_cost(mt19937 & rng) -> long long
    {
        return long(rng() % 13) - 3;
    }

    auto random_graph(int n, double edge_probability, double loop_probability, const Family & family,
        bool directed, bool costed, mt19937 & rng) -> InputGraph
    {
        InputGraph g{n, InputGraphProperties{.has_vertex_labels = family.vertex_labels, .has_edge_labels = family.edge_labels, .directed = directed, .multigraph = family.multigraph, .has_vertex_costs = costed && family.vertex_costs, .has_edge_costs = costed && family.edge_costs}};
        uniform_real_distribution<double> dist{0.0, 1.0};

        auto add = [&](int a, int b, const string & label) {
            if (g.has_edge_costs()) {
                auto c = random_cost(rng);
                directed ? g.add_directed_edge(a, b, label, c) : g.add_edge(a, b, label, c);
            }
            else
                directed ? g.add_directed_edge(a, b, label) : g.add_edge(a, b, label);
        };

        // A multigraph gets each label independently, so a pair may have both.
        auto labels_for_one_edge = [&]() -> vector<string> {
            if (! family.edge_labels)
                return {""};
            if (! family.multigraph)
                return {random_label(rng)};
            vector<string> result;
            for (auto l : {"x", "y"})
                if (dist(rng) < 0.6)
                    result.push_back(l);
            if (result.empty())
                result.push_back(random_label(rng));
            return result;
        };

        for (int v = 0; v < n; ++v) {
            if (family.vertex_labels)
                g.set_vertex_label(v, random_label(rng));
            if (g.has_vertex_costs())
                g.set_vertex_cost(v, random_cost(rng));

            if (family.loops && dist(rng) < loop_probability)
                for (auto & l : labels_for_one_edge())
                    add(v, v, l);

            for (int w = (directed ? 0 : v + 1); w < n; ++w)
                if (v != w && dist(rng) < edge_probability)
                    for (auto & l : labels_for_one_edge())
                        add(v, w, l);
        }

        return g;
    }

    auto describe(const InputGraph & g) -> string
    {
        string s;
        g.for_each_edge_and_cost([&](int f, int t, std::string_view label, optional<long long> c) {
            s += " " + to_string(f) + (g.directed() ? ">" : "-") + to_string(t);
            if (g.has_edge_labels())
                s += ":" + string{label};
            if (c)
                s += "$" + to_string(*c);
        });
        for (int v = 0; v < g.size(); ++v) {
            if (g.has_vertex_labels())
                s += " " + to_string(v) + "=" + string{g.vertex_label(v)};
            if (g.has_vertex_costs())
                s += " " + to_string(v) + "$" + to_string(g.vertex_cost(v));
        }
        return s;
    }

    // The oracle: over every injective mapping the verifier accepts, the cheapest, by a
    // cost computed straight from the graphs. Neither knows anything about reification
    // or the bound.
    auto brute_force_cheapest(const InputGraph & pattern, const InputGraph & target) -> optional<long long>
    {
        optional<long long> best;
        int np = pattern.size(), nt = target.size();
        vector<int> assignment(np, 0);
        for (;;) {
            map<int, int> mapping;
            for (int i = 0; i < np; ++i)
                mapping.emplace(i, assignment[i]);

            try {
                verify_homomorphism(pattern, target, true, false, false, mapping);
                auto cost = cost_of_mapping(pattern, target, mapping);
                if ((! best) || cost < *best)
                    best = cost;
            }
            catch (const BuggySolution &) {
            }

            int i = 0;
            for (; i < np; ++i) {
                if (++assignment[i] < nt)
                    break;
                assignment[i] = 0;
            }
            if (i == np)
                break;
        }
        return best;
    }

    auto make_params() -> HomomorphismParams
    {
        HomomorphismParams params;
        params.timeout = make_shared<Timeout>(0s);
        params.restarts_schedule = make_unique<NoRestartsSchedule>();
        params.injectivity = Injectivity::Injective;
        return params;
    }
}

TEST_CASE("random weighted instances: the cheapest mapping matches the brute-force oracle")
{
    mt19937 rng{0xc057};

    for (auto & family : all_families()) {
        for (int iter = 0; iter < 25; ++iter) {
            int np = 1 + int(rng() % 4); // 1..4
            int nt = 1 + int(rng() % 6); // 1..6
            double pattern_density = (rng() % 100) / 100.0;
            double target_density = 0.3 + (rng() % 70) / 100.0;

            auto pattern = random_graph(np, pattern_density, 0.4, family, family.pattern_directed, false, rng);
            auto target = random_graph(nt, target_density, 0.5, family, family.target_directed, true, rng);

            auto expected = brute_force_cheapest(pattern, target);

            INFO("family " << family.name() << ", iteration " << iter
                           << "\n  pattern:" << describe(pattern)
                           << "\n  target:" << describe(target)
                           << "\n  oracle=" << (expected ? to_string(*expected) : "none"));

            auto params = make_params();
            params.minimise_cost = true;
            auto result = solve_homomorphism_problem(pattern, target, params);

            // As in random_homomorphism_test, the pattern-bigger-than-target refutation
            // reports itself as incomplete, so only a found mapping has to be complete.
            CHECK((result.complete || ! expected));
            REQUIRE(result.cost.has_value() == expected.has_value());
            if (expected) {
                CHECK(*result.cost == *expected);
                CHECK_NOTHROW(verify_homomorphism(pattern, target, true, false, false, result.mapping));
                CHECK(cost_of_mapping(pattern, target, result.mapping) == *result.cost);
            }
            else
                CHECK(result.mapping.empty());

            // Deciding, not minimising, on the same graphs: a multigraph is reified either
            // way, and must find a mapping exactly when one exists.
            if (family.multigraph) {
                auto decide_params = make_params();
                auto decided = solve_homomorphism_problem(pattern, target, decide_params);
                CHECK((decided.complete || ! expected));
                CHECK(decided.mapping.empty() == ! expected.has_value());
                if (! decided.mapping.empty())
                    CHECK_NOTHROW(verify_homomorphism(pattern, target, true, false, false, decided.mapping));
            }
        }
    }
}

TEST_CASE("weighted: without pattern edge labels, the cheapest parallel edge is used")
{
    InputGraph pattern{2, false, false};
    pattern.add_edge(0, 1);

    InputGraph target{2, InputGraphProperties{.has_edge_labels = true, .multigraph = true, .has_edge_costs = true}};
    target.add_edge(0, 1, "expensive", 10);
    target.add_edge(0, 1, "cheap", 3);

    auto params = make_params();
    params.minimise_cost = true;
    auto result = solve_homomorphism_problem(pattern, target, params);
    REQUIRE(result.cost.has_value());
    CHECK(*result.cost == 3);
    CHECK(result.complete);
}

TEST_CASE("weighted: a loop only lands on a loop, and an edge only on an edge")
{
    // The pattern is one vertex with a loop. The target's only loop-free vertex has a
    // cheap edge to the other, and the other has an expensive loop: an undirected edge-
    // vertex for a loop has one neighbour, so without its own label a pattern loop could
    // land on any edge touching its image.
    InputGraph pattern{1, InputGraphProperties{.has_edge_labels = true, .multigraph = true}};
    pattern.add_edge(0, 0, "r");

    InputGraph target{2, InputGraphProperties{.has_edge_labels = true, .multigraph = true, .has_edge_costs = true}};
    target.add_edge(0, 1, "r", 1);
    target.add_edge(1, 1, "r", 50);

    auto params = make_params();
    params.minimise_cost = true;
    auto result = solve_homomorphism_problem(pattern, target, params);
    REQUIRE(result.cost.has_value());
    CHECK(*result.cost == 50);
    CHECK(result.mapping == VertexToVertexMapping{{0, 1}});
}

TEST_CASE("weighted: what cannot yet be combined with minimising cost is refused")
{
    auto costed_target = [] {
        InputGraph g{2, InputGraphProperties{.has_vertex_costs = true}};
        g.set_vertex_cost(0, 1);
        g.set_vertex_cost(1, 2);
        g.add_edge(0, 1);
        return g;
    };
    InputGraph pattern{2, false, false};
    pattern.add_edge(0, 1);
    auto target = costed_target();

    auto refused = [&](auto && adjust) {
        auto params = make_params();
        params.minimise_cost = true;
        adjust(params);
        CHECK_THROWS_AS(solve_homomorphism_problem(pattern, target, params), UnsupportedConfiguration);
    };

    refused([](HomomorphismParams & p) { p.injectivity = Injectivity::NonInjective; });
    refused([](HomomorphismParams & p) { p.injectivity = Injectivity::LocallyInjective; });
    refused([](HomomorphismParams & p) { p.induced = true; });
    refused([](HomomorphismParams & p) { p.count_solutions = true; });
    refused([](HomomorphismParams & p) { p.n_threads = 2; });
    refused([](HomomorphismParams & p) { p.staged = true; });
    refused([](HomomorphismParams & p) { p.restarts_schedule = make_unique<LubyRestartsSchedule>(LubyRestartsSchedule::default_multiplier); });
    refused([](HomomorphismParams & p) { p.pattern_less_constraints.emplace_back("0", "1"); });

    // Minimising needs something to minimise.
    InputGraph uncosted{2, false, false};
    uncosted.add_edge(0, 1);
    auto params = make_params();
    params.minimise_cost = true;
    CHECK_THROWS_AS(solve_homomorphism_problem(pattern, uncosted, params), UnsupportedConfiguration);

    // A cost on the pattern has no meaning here, minimising or not.
    InputGraph costed_pattern{2, InputGraphProperties{.has_vertex_costs = true}};
    costed_pattern.set_vertex_cost(0, 1);
    costed_pattern.set_vertex_cost(1, 1);
    costed_pattern.add_edge(0, 1);
    auto plain = make_params();
    CHECK_THROWS_AS(solve_homomorphism_problem(costed_pattern, target, plain), UnsupportedConfiguration);
}

TEST_CASE("weighted: what cannot yet be combined with a multigraph is refused")
{
    InputGraph pattern{2, InputGraphProperties{.has_edge_labels = true, .multigraph = true}};
    pattern.add_edge(0, 1, "a");
    pattern.add_edge(0, 1, "b");
    InputGraph target{2, InputGraphProperties{.has_edge_labels = true, .multigraph = true}};
    target.add_edge(0, 1, "a");
    target.add_edge(0, 1, "b");

    auto refused = [&](auto && adjust) {
        auto params = make_params();
        adjust(params);
        CHECK_THROWS_AS(solve_homomorphism_problem(pattern, target, params), UnsupportedConfiguration);
    };

    refused([](HomomorphismParams & p) { p.injectivity = Injectivity::NonInjective; });
    refused([](HomomorphismParams & p) { p.induced = true; });
    refused([](HomomorphismParams & p) { p.count_solutions = true; });

    // And the plain injective decision problem works.
    auto params = make_params();
    auto result = solve_homomorphism_problem(pattern, target, params);
    CHECK(result.mapping.size() == 2);
}
