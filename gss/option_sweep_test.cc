#include <gss/configuration.hh>
#include <gss/formats/input_graph.hh>
#include <gss/homomorphism.hh>
#include <gss/innards/verify.hh>

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <list>
#include <map>
#include <memory>
#include <optional>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace gss;
using namespace gss::innards;

using std::function;
using std::ifstream;
using std::make_shared;
using std::make_unique;
using std::map;
using std::mt19937;
using std::ofstream;
using std::pair;
using std::string;
using std::to_string;
using std::tuple;
using std::uniform_real_distribution;
using std::vector;

using std::chrono::operator""s;

// The metamorphic half of the option coverage proposed in issue #84. Its reference is not a
// brute-force oracle -- there isn't one, at these sizes -- but the *same instance solved with
// the filtering switched off*. Every option here is supposed to prune the search without
// changing which mappings exist, so:
//
//   - satisfiability is preserved by every configuration, always;
//   - the solution count is preserved too, but is only observable when counting;
//   - the particular mapping handed back is *not* preserved -- clique detection returns the
//     clique solver's mapping and the value-ordering heuristics change which solution
//     surfaces first -- so a returned mapping is checked against the verifier instead of
//     against the baseline's.
//
// The structure is a pairwise covering array over the option columns, multiplied by an
// enumerated set of instance families. The multiplication is the point: pairwise over
// options alone would not have caught #58, which needs local injectivity *and* loops *and*
// supplemental graphs together. Making graph structure a multiplier rather than another
// column gives strength 2 across options and full strength across structure, hence strength
// 3 for every (option, option, graph property) triple -- the shape of every soundness
// condition in homomorphism_traits.cc.
//
// What is deliberately not swept, and why, is in dev_docs/option-compatibility.md.

namespace
{
    // --- option columns -------------------------------------------------------------
    //
    // Value 0 of every column is the solver's own default, so the all-zeros row is the
    // default configuration. The first three columns define the *problem* rather than the
    // filtering: the baseline for a row copies those and turns everything else off.
    constexpr int semantic_columns = 3;

    struct Column
    {
        const char * name;
        int arity;
        function<auto(HomomorphismParams &, int)->void> apply;
    };

    auto columns() -> const vector<Column> &
    {
        static const vector<Column> cols{
            {"injectivity", 3, [](HomomorphismParams & p, int v) {
                 p.injectivity = (v == 0 ? Injectivity::Injective
                         : v == 1        ? Injectivity::LocallyInjective
                                         : Injectivity::NonInjective);
             }},
            {"induced", 2, [](HomomorphismParams & p, int v) { p.induced = (v == 1); }},
            {"count", 2, [](HomomorphismParams & p, int v) { p.count_solutions = (v == 1); }},
            {"clique-detection", 2, [](HomomorphismParams & p, int v) { p.clique_detection = (v == 0); }},
            {"supplementals", 2, [](HomomorphismParams & p, int v) { p.no_supplementals = (v == 1); }},
            {"exact-paths", 3, [](HomomorphismParams & p, int v) { p.number_of_exact_path_graphs = (v == 0 ? 4 : v == 1 ? 1
                                                                                                                        : 0); }},
            {"distance3", 2, [](HomomorphismParams & p, int v) { p.distance3 = (v == 1); }},
            {"k4", 2, [](HomomorphismParams & p, int v) { p.k4 = (v == 1); }},
            {"nds", 2, [](HomomorphismParams & p, int v) { p.no_nds = (v == 1); }},
            {"cliques", 2, [](HomomorphismParams & p, int v) { p.clique_size_constraints = (v == 1); }},
            {"cliques-on-supplementals", 2, [](HomomorphismParams & p, int v) { p.clique_size_constraints_on_supplementals = (v == 1); }},
            {"staged", 2, [](HomomorphismParams & p, int v) {
                 p.staged = (v == 1);
                 // a budget small enough to force the Stage-1 -> Stage-2 transition on
                 // instances this size, which is the whole point of exercising it here
                 p.staged_first_round_backtracks = 1;
             }},
            {"restarts", 3, [](HomomorphismParams & p, int v) {
                 // No timed restarts: they are wall-clock dependent, so a cell's outcome
                 // would not be reproducible.
                 if (v == 0)
                     p.restarts_schedule = make_unique<NoRestartsSchedule>();
                 else if (v == 1)
                     p.restarts_schedule = make_unique<LubyRestartsSchedule>(LubyRestartsSchedule::default_multiplier);
                 else
                     p.restarts_schedule = make_unique<GeometricRestartsSchedule>(
                         GeometricRestartsSchedule::default_initial_value, GeometricRestartsSchedule::default_multiplier);
             }},
            {"nogoods", 2, [](HomomorphismParams & p, int v) {
                 if (v == 1)
                     p.nogood_size_limit = 0;
             }},
            {"value-ordering", 5, [](HomomorphismParams & p, int v) {
                 p.value_ordering_heuristic = (v == 0 ? ValueOrdering::Biased
                         : v == 1                     ? ValueOrdering::None
                         : v == 2                     ? ValueOrdering::Degree
                         : v == 3                     ? ValueOrdering::AntiDegree
                                                      : ValueOrdering::Random);
             }},
        };
        return cols;
    }

    using Row = vector<int>;

    auto encode(const Row & row) -> string
    {
        string s;
        for (auto v : row)
            s += char('0' + v);
        return s;
    }

    auto describe(const Row & row) -> string
    {
        string s;
        for (unsigned c = 0; c < columns().size(); ++c)
            s += string{columns()[c].name} + "=" + to_string(row[c]) + " ";
        return s;
    }

    auto make_params(const Row & row) -> HomomorphismParams
    {
        // No timeout, so every solve here runs to completion and the two answers being
        // compared are both exact. That is why nothing below compares `complete`: the one
        // configuration that reports an incomplete search is the trivial
        // pattern-bigger-than-target refutation, which is exact anyway.
        HomomorphismParams params;
        params.timeout = make_shared<Timeout>(0s);
        params.restarts_schedule = make_unique<NoRestartsSchedule>();
        for (unsigned c = 0; c < columns().size(); ++c)
            columns()[c].apply(params, row[c]);
        return params;
    }

    // The reference configuration: the same problem, with nothing filtering it. Everything a
    // row can switch on is off here, so a disagreement is always the row's doing.
    auto make_baseline_params(const Row & row) -> HomomorphismParams
    {
        HomomorphismParams params;
        params.timeout = make_shared<Timeout>(0s);
        params.restarts_schedule = make_unique<NoRestartsSchedule>();
        for (int c = 0; c < semantic_columns; ++c)
            columns()[c].apply(params, row[c]);
        params.clique_detection = false;
        params.no_supplementals = true;
        params.no_nds = true;
        return params;
    }

    // --- the covering array ---------------------------------------------------------

    // A pairwise covering array, built greedily: take the first pair no row covers yet, fix
    // those two columns, then fill the rest choosing whichever value covers most of what is
    // still missing. Deterministic, so the rows -- and hence the golden file -- only move
    // when a column is added or changed. Typically within half again of optimal, which for
    // these arities is around two dozen rows.
    auto pairwise_covering_array(const vector<Column> & cols) -> vector<Row>
    {
        unsigned n = cols.size();
        // covered[i][j][vi * arity_j + vj] for i < j
        vector<vector<vector<bool>>> covered(n, vector<vector<bool>>(n));
        for (unsigned i = 0; i < n; ++i)
            for (unsigned j = i + 1; j < n; ++j)
                covered[i][j].resize(cols[i].arity * cols[j].arity, false);

        auto mark = [&](const Row & row) {
            for (unsigned i = 0; i < n; ++i)
                for (unsigned j = i + 1; j < n; ++j)
                    covered[i][j][row[i] * cols[j].arity + row[j]] = true;
        };

        auto first_uncovered = [&]() -> std::optional<tuple<unsigned, unsigned, int, int>> {
            for (unsigned i = 0; i < n; ++i)
                for (unsigned j = i + 1; j < n; ++j)
                    for (int vi = 0; vi < cols[i].arity; ++vi)
                        for (int vj = 0; vj < cols[j].arity; ++vj)
                            if (! covered[i][j][vi * cols[j].arity + vj])
                                return tuple{i, j, vi, vj};
            return std::nullopt;
        };

        vector<Row> rows;
        while (auto seed = first_uncovered()) {
            auto [si, sj, svi, svj] = *seed;
            Row row(n, -1);
            row[si] = svi;
            row[sj] = svj;

            for (unsigned c = 0; c < n; ++c) {
                if (row[c] != -1)
                    continue;
                int best_value = 0, best_gain = -1;
                for (int v = 0; v < cols[c].arity; ++v) {
                    int gain = 0;
                    for (unsigned o = 0; o < n; ++o) {
                        if (o == c || row[o] == -1)
                            continue;
                        auto [i, j] = std::minmax(c, o);
                        int vi = (i == c ? v : row[o]), vj = (j == c ? v : row[o]);
                        if (! covered[i][j][vi * cols[j].arity + vj])
                            ++gain;
                    }
                    if (gain > best_gain) {
                        best_gain = gain;
                        best_value = v;
                    }
                }
                row[c] = best_value;
            }

            mark(row);
            rows.push_back(row);
        }

        return rows;
    }

    // Pairwise is a floor, not a ceiling. Both #58 and #91 live in a triple -- an injectivity
    // mode, a graph property and one filter -- and the graph property is a family, so what
    // has to be added here is each filter on its own against each injectivity mode. Without
    // these, whether a given (injectivity, filter) pair is ever tried in isolation is an
    // accident of how the greedy array came out.
    auto dangerous_triple_rows(const vector<Column> & cols) -> vector<Row>
    {
        vector<Row> rows;
        vector<pair<string, int>> one_filter_at_a_time{
            {"clique-detection", 0}, // on, which is the default
            {"supplementals", 0}, // on
            {"nds", 0}, // on
            {"distance3", 1},
            {"k4", 1},
            {"cliques", 1},
            {"cliques-on-supplementals", 1},
        };

        for (int injectivity = 0; injectivity < 3; ++injectivity)
            for (int counting = 0; counting < 2; ++counting)
                for (auto & [name, value] : one_filter_at_a_time) {
                    // everything off, then exactly one filter on
                    Row row(cols.size(), 0);
                    row[0] = injectivity;
                    row[2] = counting;
                    for (unsigned c = semantic_columns; c < cols.size(); ++c) {
                        // the "off" value of each filter column: 1 for the ones whose default
                        // is on, 0 for the rest
                        string column_name = cols[c].name;
                        row[c] = (column_name == "clique-detection" || column_name == "supplementals" || column_name == "nds") ? 1 : 0;
                    }
                    for (unsigned c = semantic_columns; c < cols.size(); ++c)
                        if (name == cols[c].name)
                            row[c] = value;
                    // clique-size constraints on supplementals need the supplementals, and
                    // need the base flag too, or the cell is vacuous by construction
                    if (name == "cliques-on-supplementals") {
                        for (unsigned c = semantic_columns; c < cols.size(); ++c) {
                            string column_name = cols[c].name;
                            if (column_name == "cliques")
                                row[c] = 1;
                            else if (column_name == "supplementals")
                                row[c] = 0;
                        }
                    }
                    rows.push_back(row);
                }

        return rows;
    }

    auto sweep_rows() -> vector<Row>
    {
        auto rows = pairwise_covering_array(columns());
        for (auto & row : dangerous_triple_rows(columns()))
            if (rows.end() == find(rows.begin(), rows.end(), row))
                rows.push_back(row);
        return rows;
    }

    // --- instance families ----------------------------------------------------------

    using Instance = pair<InputGraph, InputGraph>;

    struct Family
    {
        string name;
        // Instances alternate between two regimes rather than being drawn from one
        // distribution, because the two things the sweep wants are in tension. An even
        // index is filter-friendly: a dense pattern into a sparse target, mostly
        // unsatisfiable, where everything that can prune gets to prune. An odd index is
        // solution-friendly: the other way round, so that the solution counts being
        // compared are large enough to disagree in interesting ways. Alternating puts one
        // of each in every cell.
        function<auto(mt19937 &, int index)->Instance> generate;
    };

    auto random_label(mt19937 & rng) -> string
    {
        return (rng() % 2) ? "x" : "y";
    }

    auto random_graph(int n, double edge_probability, double loop_probability,
        bool directed, bool vertex_labels, bool edge_labels, mt19937 & rng) -> InputGraph
    {
        InputGraph g{n, vertex_labels, edge_labels, directed};
        uniform_real_distribution<double> dist{0.0, 1.0};

        auto add = [&](int a, int b) {
            if (directed)
                g.add_directed_edge(a, b, edge_labels ? random_label(rng) : "");
            else if (edge_labels)
                g.add_edge(a, b, random_label(rng));
            else
                g.add_edge(a, b);
        };

        for (int v = 0; v < n; ++v) {
            if (vertex_labels)
                g.set_vertex_label(v, random_label(rng));
            if (loop_probability > 0.0 && dist(rng) < loop_probability)
                add(v, v);
            for (int w = (directed ? 0 : v + 1); w < n; ++w)
                if (v != w && dist(rng) < edge_probability)
                    add(v, w);
        }

        return g;
    }

    // A failing cell is no use without the instance that failed, and these are generated,
    // not checked in: this is what a reproduction gets built from.
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

    auto clique(int n) -> InputGraph
    {
        InputGraph g{n, false, false};
        for (int v = 0; v < n; ++v)
            for (int w = v + 1; w < n; ++w)
                g.add_edge(v, w);
        return g;
    }

    // Sizes are kept small enough that a non-injective count stays in the hundreds rather
    // than the millions: the sweep runs thousands of solves, and an instance that takes a
    // second would make it unrunnable. They are still far larger than the oracle test's,
    // which is the point of not needing an oracle -- the filters actually fire here.
    auto all_families() -> vector<Family>
    {
        vector<Family> families;

        for (bool loops : {false, true})
            for (bool directed : {false, true})
                for (bool vertex_labels : {false, true})
                    for (bool edge_labels : {false, true}) {
                        string name = (loops ? "loopy" : "loopless");
                        name += directed ? "-directed" : "-undirected";
                        name += vertex_labels ? "-vlabels" : "";
                        name += edge_labels ? "-elabels" : "";
                        if (! vertex_labels && ! edge_labels)
                            name += "-unlabelled";
                        families.push_back(Family{name, [=](mt19937 & rng, int index) -> Instance {
                                                      bool filter_friendly = (0 == index % 2);
                                                      int np = filter_friendly ? 5 + int(rng() % 2) : 4 + int(rng() % 2);
                                                      int nt = filter_friendly ? 9 + int(rng() % 3) : 7 + int(rng() % 2);
                                                      double pattern_density = (filter_friendly ? 0.45 : 0.4) + (rng() % 15) / 100.0;
                                                      double target_density = (filter_friendly ? 0.2 : 0.45) + (rng() % 15) / 100.0;
                                                      // The target is loopier than the pattern on purpose: a pattern
                                                      // loop needs somewhere to go, and an instance nothing can satisfy
                                                      // for a structural reason never reaches a filter at all.
                                                      return Instance{
                                                          random_graph(np, pattern_density, loops ? 0.15 : 0.0, directed, vertex_labels, edge_labels, rng),
                                                          random_graph(nt, target_density, loops ? 0.4 : 0.0, directed, vertex_labels, edge_labels, rng)};
                                                  }});
                    }

        // Three shapes that random generation essentially never produces, each of which
        // reaches a piece of the pipeline that nothing else does.
        families.push_back(Family{"clique-pattern", [](mt19937 & rng, int) -> Instance {
                                      // fires CliqueShortcutStep, and gives the clique-size filters something to bite
                                      // on: a target this sparse has plenty of vertices in no clique that big
                                      int k = 4 + int(rng() % 2);
                                      return Instance{clique(k), random_graph(10, 0.4 + (rng() % 15) / 100.0, 0.0, false, false, false, rng)};
                                  }});

        families.push_back(Family{"pattern-bigger-than-target", [](mt19937 & rng, int) -> Instance {
                                      // fires PatternBiggerThanTargetStep under injectivity, and is a plain hard
                                      // instance without it
                                      return Instance{
                                          random_graph(8, 0.4, 0.0, false, false, false, rng),
                                          random_graph(5, 0.5, 0.0, false, false, false, rng)};
                                  }});

        families.push_back(Family{"equal-size", [](mt19937 & rng, int) -> Instance {
                                      // equal sizes make degree_and_nds_are_exact() true under induced injectivity,
                                      // which is a different filter from the usual inequality
                                      int n = 5 + int(rng() % 2);
                                      double density = 0.35 + (rng() % 20) / 100.0;
                                      return Instance{
                                          random_graph(n, density, 0.0, false, false, false, rng),
                                          random_graph(n, density, 0.0, false, false, false, rng)};
                                  }});

        return families;
    }

    // --- reading the answer back ----------------------------------------------------

    struct Answer
    {
        bool satisfiable = false;
        loooong count{0};
        VertexToVertexMapping mapping;
        bool any_optional_filter_fired = false;
    };

    // Did any filter that an option can turn off remove anything? Read back out of the
    // extra stats, which is where the instrumentation surfaces it: everything except the
    // original graph's adjacency, the vertex labels and the loops, none of which are
    // optional filtering.
    auto parse_activations(const std::list<string> & extra_stats) -> bool
    {
        unsigned long long total = 0;
        for (auto & line : extra_stats) {
            if (! (line.starts_with("filter_activations_initial =") || line.starts_with("filter_activations_search =")))
                continue;
            std::istringstream tokens{line.substr(line.find('=') + 1)};
            string token;
            while (tokens >> token) {
                auto colon = token.find(':');
                auto key = token.substr(0, colon);
                if (key == "vertex_labels" || key == "loops" || key == "adjacency" || key == "edge_labels")
                    continue;
                // degree is a list over graph pairs, supplemental a list over the
                // supplemental ones; either way every number in it is optional filtering,
                // except the original graph's degree column, which is not optional at all
                // (the injectivity mode decides whether it runs)
                std::istringstream values{token.substr(colon + 1)};
                string value;
                while (getline(values, value, ','))
                    if (value != "-")
                        total += std::stoull(value);
            }
        }
        return total != 0;
    }

    auto solve(const InputGraph & pattern, const InputGraph & target, const HomomorphismParams & params) -> Answer
    {
        auto result = solve_homomorphism_problem(pattern, target, params);
        return Answer{
            params.count_solutions ? (result.solution_count > loooong{0}) : (! result.mapping.empty()),
            result.solution_count,
            result.mapping,
            parse_activations(result.extra_stats)};
    }

    // --- the golden table -----------------------------------------------------------

    enum class Outcome
    {
        Ok,
        Unsupported,
        Diverges
    };

    auto outcome_name(Outcome outcome) -> string
    {
        switch (outcome) {
        case Outcome::Ok: return "ok";
        case Outcome::Unsupported: return "unsupported";
        case Outcome::Diverges: return "diverges";
        }
        return "?";
    }

    struct Cell
    {
        Outcome outcome = Outcome::Ok;
        bool active = false;
    };

    auto golden_header() -> vector<string>
    {
        vector<string> header{
            "# Outcome of every (configuration, instance family) cell of the option sweep.",
            "# Generated by option_sweep_test; after reviewing a diff, regenerate with",
            "#     GSS_SWEEP_REGENERATE=1 ./build/option_sweep_test",
            "# A silent divergence fails the build whatever this file says; what this file is",
            "# for is making a newly unsupported or newly vacuous cell show up as a reviewable",
            "# diff rather than as nothing at all. See dev_docs/option-compatibility.md.",
            "#",
            "# config is one digit per option column, value 0 being the solver's default:",
            "#   "};
        string names;
        for (auto & column : columns())
            names += string{column.name} + " ";
        header.back() += names;
        header.push_back("# outcome: ok | unsupported (UnsupportedConfiguration thrown) | diverges");
        header.push_back("# activation: active if any optional filter removed anything in any instance of the cell");
        header.push_back("config\tfamily\toutcome\tactivation");
        return header;
    }
}

TEST_CASE("option sweep: filtering options do not change the answer")
{
    auto rows = sweep_rows();
    auto families = all_families();

    // Four instances per cell -- two of each regime -- keeps the whole sweep to a fraction of
    // a second. More is useful when hunting, and can only ever find more, so
    // GSS_SWEEP_INSTANCES raises it; in that mode the activation column of the golden table is
    // read as a floor rather than as an exact value, since extra instances can only wake
    // filters up.
    constexpr int default_instances_per_cell = 4;
    int instances_per_cell = default_instances_per_cell;
    if (auto env = std::getenv("GSS_SWEEP_INSTANCES"))
        instances_per_cell = std::max(1, atoi(env));
    const bool deep = (instances_per_cell != default_instances_per_cell);

    // Instances belong to the family, not to the configuration, so that every row sees the
    // same ones and a divergence can be compared across rows. Seeded per family so that
    // adding a family does not renumber the others.
    vector<vector<Instance>> instances(families.size());
    for (unsigned f = 0; f < families.size(); ++f) {
        mt19937 rng{0xc0ffee + f};
        for (int i = 0; i < instances_per_cell; ++i)
            instances[f].push_back(families[f].generate(rng, i));
    }

    // The baseline depends only on the problem, not on the filtering, so it is computed once
    // per (family, instance, semantics) and shared by every row that asks the same question.
    map<tuple<unsigned, int, int, int, int>, Answer> baselines;
    auto baseline_for = [&](unsigned f, int i, const Row & row) -> const Answer & {
        auto key = tuple{f, i, row[0], row[1], row[2]};
        auto found = baselines.find(key);
        if (found == baselines.end())
            found = baselines.emplace(key, solve(instances[f][i].first, instances[f][i].second, make_baseline_params(row))).first;
        return found->second;
    };

    map<string, Cell> cells;

    for (auto & row : rows) {
        auto config = encode(row);
        bool counting = (row[2] == 1);
        bool injective = (row[0] == 0), locally_injective = (row[0] == 1), induced = (row[1] == 1);

        for (unsigned f = 0; f < families.size(); ++f) {
            Cell cell;

            for (int i = 0; i < instances_per_cell; ++i) {
                const auto & [pattern, target] = instances[f][i];
                const auto & baseline = baseline_for(f, i, row);

                auto params = make_params(row);
                // The filter-friendly instances are solved with the instrumentation on, which
                // is where the cell's activation reading comes from -- they are the ones where
                // a filter having nothing to do means something. The rest go through the
                // ordinary path, because recording selects a different instantiation of the
                // propagation template and the sweep should mostly be testing the real one;
                // the first instance is solved both ways to pin down that the two agree.
                params.record_filter_activations = (0 == i % 2);

                Answer answer;
                try {
                    answer = solve(pattern, target, params);
                }
                catch (const UnsupportedConfiguration &) {
                    cell.outcome = Outcome::Unsupported;
                    break;
                }

                cell.active = cell.active || answer.any_optional_filter_fired;

                std::ostringstream where;
                where << "config " << config << " (" << describe(row) << ")"
                      << "\n  family " << families[f].name << ", instance " << i
                      << "\n  pattern:" << describe(pattern)
                      << "\n  target: " << describe(target)
                      << "\n  baseline: satisfiable=" << baseline.satisfiable << " count=" << baseline.count
                      << "\n  this row: satisfiable=" << answer.satisfiable << " count=" << answer.count;
                INFO(where.str());

                bool ok = (answer.satisfiable == baseline.satisfiable);
                CHECK(answer.satisfiable == baseline.satisfiable);
                if (counting) {
                    ok = ok && (answer.count == baseline.count);
                    CHECK(answer.count == baseline.count);
                }

                // The mapping itself is not preserved, so it is checked against the verifier
                // rather than against the baseline's mapping.
                if (! answer.mapping.empty()) {
                    try {
                        verify_homomorphism(pattern, target, injective, locally_injective, induced, answer.mapping);
                    }
                    catch (const BuggySolution & e) {
                        ok = false;
                        UNSCOPED_INFO("returned mapping does not verify: " << e.what());
                        CHECK(false);
                    }
                }

                if (i == 0) {
                    // recording must not change what comes back
                    params.record_filter_activations = false;
                    auto unrecorded = solve(pattern, target, params);
                    ok = ok && unrecorded.satisfiable == answer.satisfiable && unrecorded.count == answer.count;
                    CHECK(unrecorded.satisfiable == answer.satisfiable);
                    CHECK(unrecorded.count == answer.count);
                }

                if (! ok)
                    cell.outcome = Outcome::Diverges;
            }

            // one cell per (configuration, family), so a duplicated row would silently
            // discard a result rather than showing up as one
            REQUIRE(cells.emplace(config + "\t" + families[f].name, cell).second);
        }
    }

    // --- compare against, or regenerate, the golden table ---------------------------

    vector<string> lines = golden_header();
    for (auto & [key, cell] : cells)
        lines.push_back(key + "\t" + outcome_name(cell.outcome) + "\t" + (cell.active ? "active" : "vacuous"));

    if (std::getenv("GSS_SWEEP_REGENERATE")) {
        ofstream out{GSS_OPTION_SWEEP_GOLDEN};
        for (auto & line : lines)
            out << line << "\n";
        REQUIRE(out.good());
        WARN("regenerated " << GSS_OPTION_SWEEP_GOLDEN << " (" << cells.size() << " cells)");
        return;
    }

    ifstream in{GSS_OPTION_SWEEP_GOLDEN};
    INFO("golden table " << GSS_OPTION_SWEEP_GOLDEN);
    REQUIRE(in);
    vector<string> expected;
    for (string line; getline(in, line);)
        if (! line.empty())
            expected.push_back(line);

    // Compared as a table rather than line by line, so that a genuinely new cell reads as an
    // addition rather than shifting everything after it.
    map<string, string> expected_cells, got_cells;
    for (auto & line : expected)
        if (! line.starts_with("#") && ! line.starts_with("config\t"))
            expected_cells.emplace(line.substr(0, line.find('\t', line.find('\t') + 1)), line);
    for (auto & line : lines)
        if (! line.starts_with("#") && ! line.starts_with("config\t"))
            got_cells.emplace(line.substr(0, line.find('\t', line.find('\t') + 1)), line);

    for (auto & [key, line] : got_cells) {
        INFO("cell " << key << "\n  in the golden table: " << (expected_cells.count(key) ? expected_cells[key] : string{"(absent)"}) << "\n  this run:            " << line << "\n  if this is intended, review it and regenerate with GSS_SWEEP_REGENERATE=1");
        CHECK(expected_cells.count(key));
        if (! expected_cells.count(key))
            continue;
        if (! deep)
            CHECK(expected_cells[key] == line);
        else {
            // outcome must match exactly; activation only has to be at least what the table
            // records, because the extra instances can wake a filter up but never silence one
            auto field = [](const string & l, int n) {
                size_t from = 0;
                for (int i = 0; i < n; ++i)
                    from = l.find('\t', from) + 1;
                auto to = l.find('\t', from);
                return l.substr(from, to == string::npos ? string::npos : to - from);
            };
            CHECK(field(expected_cells[key], 2) == field(line, 2));
            if (field(expected_cells[key], 3) == "active")
                CHECK(field(line, 3) == "active");
        }
    }

    for (auto & [key, line] : expected_cells) {
        INFO("cell " << key << " is in the golden table but was not produced by this run: " << line);
        CHECK(got_cells.count(key));
    }
}
