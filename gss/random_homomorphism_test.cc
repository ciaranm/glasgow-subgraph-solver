#include <gss/formats/input_graph.hh>
#include <gss/homomorphism.hh>
#include <gss/innards/verify.hh>

#include <catch2/catch_test_macros.hpp>

#include <chrono>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <vector>

using namespace gss;
using namespace gss::innards;

using std::make_shared;
using std::make_unique;
using std::map;
using std::mt19937;
using std::set;
using std::string;
using std::to_string;
using std::uniform_real_distribution;
using std::vector;

using std::chrono::operator""s;

namespace
{
    // The shape of instance to generate. These four axes are what gate different code paths
    // in the model and the searcher: directedness and edge labels each select a different
    // instantiation of the forward-checking template, loops turn off whole filters through
    // the traits layer, and vertex labels are filtered for at domain-initialisation time.
    //
    // A family applies the same shape to the pattern and to the target. The model takes its
    // directedness and its label-ness from the *pattern* -- HomomorphismModel's constructor
    // branches on pattern.directed() and pattern.has_edge_labels() -- so a directed target
    // under an undirected pattern is not the problem either file appears to describe. That
    // is a gap to document rather than to sweep; see dev_docs/option-compatibility.md.
    struct Family
    {
        bool loops, directed, vertex_labels, edge_labels;

        [[nodiscard]] auto name() const -> string
        {
            string n = loops ? "loopy" : "loopless";
            n += directed ? " directed" : " undirected";
            if (vertex_labels)
                n += " vertex-labelled";
            if (edge_labels)
                n += " edge-labelled";
            if (! vertex_labels && ! edge_labels)
                n += " unlabelled";
            return n;
        }
    };

    auto all_families() -> vector<Family>
    {
        vector<Family> families;
        for (bool loops : {false, true})
            for (bool directed : {false, true})
                for (bool vertex_labels : {false, true})
                    for (bool edge_labels : {false, true})
                        families.push_back(Family{loops, directed, vertex_labels, edge_labels});
        return families;
    }

    // Two labels is the smallest alphabet that can both match and mismatch, and keeping it
    // that small is what makes a labelled instance interesting: with a large alphabet almost
    // every instance is unsatisfiable for the same boring reason, and the searcher's label
    // checks never have to discriminate.
    auto random_label(mt19937 & rng) -> string
    {
        return (rng() % 2) ? "x" : "y";
    }

    // A random graph of the given family. Loops and edges are generated separately because a
    // loop is the interesting case: it is an edge that carries a label like any other, but
    // one the solver keeps out of its adjacency rows.
    auto random_graph(int n, double edge_probability, double loop_probability, const Family & family, mt19937 & rng) -> InputGraph
    {
        InputGraph g{n, family.vertex_labels, family.edge_labels, family.directed};
        uniform_real_distribution<double> dist{0.0, 1.0};

        auto add = [&](int a, int b) {
            // add_edge keeps an undirected graph undirected however its edges are labelled,
            // which is the distinction #86 turned on; add_directed_edge requires the graph
            // to have been declared directed, which is why the family drives the constructor.
            if (family.directed)
                g.add_directed_edge(a, b, family.edge_labels ? random_label(rng) : "");
            else if (family.edge_labels)
                g.add_edge(a, b, random_label(rng));
            else
                g.add_edge(a, b);
        };

        for (int v = 0; v < n; ++v) {
            if (family.vertex_labels)
                g.set_vertex_label(v, random_label(rng));

            if (family.loops && dist(rng) < loop_probability)
                add(v, v);

            // undirected: each unordered pair once. directed: each ordered pair
            // independently, so that one-way arcs actually turn up rather than every edge
            // being reciprocated.
            for (int w = (family.directed ? 0 : v + 1); w < n; ++w)
                if (v != w && dist(rng) < edge_probability)
                    add(v, w);
        }

        return g;
    }

    auto describe(const InputGraph & g) -> string
    {
        string s;
        g.for_each_edge([&](int f, int t, std::string_view label) {
            s += " " + to_string(f) + (g.directed() ? ">" : "-") + to_string(t);
            if (g.has_edge_labels())
                s += ":" + string{label};
        });
        if (g.has_vertex_labels())
            for (int v = 0; v < g.size(); ++v)
                s += " " + to_string(v) + "=" + string{g.vertex_label(v)};
        return s;
    }

    // The brute-force oracle: every complete mapping pattern -> target that the solver's own
    // verifier accepts. verify_homomorphism does no search, so this is an independent
    // specification of "what counts as a solution" for these options.
    auto brute_force_solutions(const InputGraph & pattern, const InputGraph & target,
        bool injective, bool locally_injective, bool induced) -> set<map<int, int>>
    {
        set<map<int, int>> solutions;
        int np = pattern.size(), nt = target.size();
        if (np == 0 || nt == 0)
            return solutions;

        // iterate over all nt^np functions as a mixed-radix counter
        vector<int> assignment(np, 0);
        for (;;) {
            map<int, int> mapping;
            for (int i = 0; i < np; ++i)
                mapping.emplace(i, assignment[i]);

            try {
                verify_homomorphism(pattern, target, injective, locally_injective, induced, mapping);
                solutions.insert(move(mapping));
            }
            catch (const BuggySolution &) {
                // not a valid mapping for these options
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

        return solutions;
    }

    auto make_params(Injectivity injectivity, bool induced) -> HomomorphismParams
    {
        HomomorphismParams params;
        params.timeout = make_shared<Timeout>(0s);
        params.restarts_schedule = make_unique<NoRestartsSchedule>();
        params.injectivity = injectivity;
        params.induced = induced;
        return params;
    }
}

// The oracle is exponential (nt^np), so the instances stay tiny -- but there are only a
// handful of problem-defining axes, so rather than sampling them this runs every one of them
// on every instance: three injectivity modes x induced or not, over all sixteen combinations
// of loops x directed x vertex labels x edge labels. Filtering options are *not* varied here;
// that is what the metamorphic sweep in option_sweep_test.cc is for, which needs no oracle
// and so can afford instances large enough for the filters to fire.
TEST_CASE("random instances: solver enumeration matches the brute-force oracle")
{
    // Fixed seed so failures reproduce; the family name and iteration index identify the
    // instance.
    mt19937 rng{0x5eed};

    for (auto & family : all_families()) {
        for (int iter = 0; iter < 60; ++iter) {
            int np = 1 + int(rng() % 4); // 1..4
            int nt = 1 + int(rng() % 5); // 1..5
            double pattern_density = (rng() % 100) / 100.0;
            double target_density = (rng() % 100) / 100.0;
            double loop_probability = (rng() % 3 == 0) ? 0.25 : 0.5;

            auto pattern = random_graph(np, pattern_density, loop_probability, family, rng);
            auto target = random_graph(nt, target_density, loop_probability, family, rng);

            for (auto injectivity : {Injectivity::Injective, Injectivity::LocallyInjective, Injectivity::NonInjective}) {
                for (bool induced : {false, true}) {
                    bool injective = (injectivity == Injectivity::Injective);
                    bool locally_injective = (injectivity == Injectivity::LocallyInjective);

                    auto expected = brute_force_solutions(pattern, target, injective, locally_injective, induced);

                    INFO("family " << family.name() << ", iteration " << iter
                                   << ", |pattern|=" << np << ", |target|=" << nt
                                   << ", induced=" << induced << ", injectivity=" << int(injectivity)
                                   << "\n  pattern:" << describe(pattern)
                                   << "\n  target:" << describe(target)
                                   << "\n  oracle=" << expected.size());

                    auto params = make_params(injectivity, induced);
                    params.count_solutions = true;

                    set<map<int, int>> got;
                    params.enumerate_callback = [&](const VertexToVertexMapping & mapping) {
                        got.insert(mapping);
                        return true;
                    };

                    auto result = solve_homomorphism_problem(pattern, target, params);
                    CHECK(got == expected);
                    CHECK(result.solution_count == int(expected.size()));
                    // a non-empty enumeration must have run to completion; the only
                    // incomplete case is the trivial injective pattern-bigger-than-target
                    // shortcut (0 solutions)
                    CHECK((result.complete || expected.empty()));

                    // Staged counting must match the oracle exactly, including under the
                    // no-restarts schedule counting uses: a tiny first-round budget forces
                    // the Stage-1 -> Stage-2 transition mid-enumeration, and the
                    // restart-resumption nogoods must stop Stage 2 re-counting the solutions
                    // Stage 1 already found. (This is the regression guard for the
                    // might_have_watches fix: with the watch machinery disabled under
                    // staging, the transition posts no nogoods and Stage 2 re-explores the
                    // whole tree, inflating the count.) The distinct-set check catches a
                    // missed solution; the solution_count check catches a double-counted one.
                    set<map<int, int>> staged_got;
                    params.staged = true;
                    params.staged_first_round_backtracks = 1;
                    params.enumerate_callback = [&](const VertexToVertexMapping & mapping) {
                        staged_got.insert(mapping);
                        return true;
                    };
                    auto staged_count = solve_homomorphism_problem(pattern, target, params);
                    CHECK(staged_got == expected);
                    CHECK(staged_count.solution_count == int(expected.size()));
                    CHECK((staged_count.complete || expected.empty()));

                    // Decision mode takes a different route through the pipeline: the
                    // pattern-bigger-than-target refutation, the target-loop shortcut and the
                    // clique reduction all conclude without searching, and each is only sound
                    // under conditions of its own (#93, #94). So it is checked against the
                    // oracle too, rather than only against counting -- satisfiability, and
                    // that whatever mapping comes back really is one of the solutions.
                    params.count_solutions = false;
                    params.enumerate_callback = {};
                    for (bool staged : {false, true}) {
                        params.staged = staged;
                        auto decision = solve_homomorphism_problem(pattern, target, params);
                        INFO("decision, staged=" << staged << ", satisfiable=" << (! decision.mapping.empty()));
                        CHECK((! decision.mapping.empty()) == (! expected.empty()));
                        CHECK((decision.mapping.empty() || expected.contains(decision.mapping)));
                    }
                }
            }
        }
    }
}
