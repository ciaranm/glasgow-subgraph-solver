#include <gss/clique.hh>
#include <gss/innards/clique_size_constraints.hh>
#include <gss/innards/homomorphism_proofs.hh>

#include <chrono>
#include <memory>
#include <optional>
#include <string>
#include <vector>

using namespace gss;
using namespace gss::innards;

using std::make_unique;
using std::max;
using std::nullopt;
using std::optional;
using std::shared_ptr;
using std::to_string;
using std::vector;

using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::steady_clock;

namespace
{
    auto find_clique(
        const shared_ptr<Timeout> & timeout,
        unsigned size,
        const vector<SVOBitset> & rows,
        unsigned g,
        unsigned max_graphs,
        unsigned v,
        optional<int> largest_if_target,
        vector<int> & best_knowns,
        std::list<std::string> & build_times,
        std::list<std::string> & solve_times,
        std::list<std::string> & find_nodes,
        std::list<std::string> & prove_nodes) -> int
    {
        if (largest_if_target && (best_knowns[v] >= *largest_if_target))
            return best_knowns[v];

        auto start_time = steady_clock::now();

        CliqueParams params;
        params.timeout = timeout;
        params.start_time = steady_clock::now();
        params.decide = nullopt;
        params.restarts_schedule = make_unique<NoRestartsSchedule>();
        if (largest_if_target)
            params.stop_after_finding = *largest_if_target - 1;

        vector<int> include(size, -1), invinclude(size, 0);
        int count = 0;
        for (unsigned w = 0; w < size; ++w)
            if (w != v && rows[w * max_graphs + g].test(v)) {
                include[w] = count;
                invinclude[count] = w;
                ++count;
            }

        InputGraph gv(count, false, false);
        for (unsigned f = 0; f < size; ++f)
            if (include[f] != -1)
                for (unsigned t = 0; t < size; ++t)
                    if (f != t && include[t] != -1 && rows[f * max_graphs + g].test(t))
                        gv.add_edge(include[f], include[t]);

        build_times.push_back(to_string(duration_cast<milliseconds>(steady_clock::now() - start_time).count()));

        start_time = steady_clock::now();

        auto result = solve_clique_problem(gv, params);

        solve_times.push_back(to_string(duration_cast<milliseconds>(steady_clock::now() - start_time).count()));
        find_nodes.push_back(to_string(result.find_nodes));
        prove_nodes.push_back(to_string(result.prove_nodes));

        best_knowns[v] = max<int>(best_knowns[v], result.clique.size() + 1);
        for (auto & w : result.clique)
            best_knowns[invinclude[w]] = max<int>(best_knowns[invinclude[w]], result.clique.size() + 1);

        return result.clique.size() + 1;
    }

    auto build_pattern_clique_sizes(CliqueSizeData & data, const ProcessedGraphsData & graphs,
        unsigned max_graphs, unsigned pattern_size, const shared_ptr<Timeout> & timeout) -> void
    {
        for (unsigned g = 0; g < data.max_graphs_for_clique_size_constraints; ++g) {
            for (unsigned v = 0; v < pattern_size; ++v) {
                auto c = find_clique(timeout, pattern_size, graphs.pattern_graph_rows, g, max_graphs, v, nullopt,
                    data.pattern_cliques_best_knowns[g], data.pattern_cliques_build_times, data.pattern_cliques_solve_times,
                    data.pattern_cliques_solve_find_nodes, data.pattern_cliques_solve_prove_nodes);
                data.pattern_cliques_sizes[g][v] = c;
                data.largest_pattern_clique[g] = max(data.largest_pattern_clique[g], c);
            }
            data.has_pattern_cliques_sizes = true;
        }
    }

    auto build_target_clique_size(CliqueSizeData & data, const ProcessedGraphsData & graphs,
        unsigned max_graphs, unsigned target_size, const shared_ptr<Timeout> & timeout, int v) -> void
    {
        if (0 == data.target_cliques_sizes[0][v])
            for (unsigned g = 0; g < data.max_graphs_for_clique_size_constraints; ++g) {
                data.target_cliques_sizes[g][v] = find_clique(timeout, target_size, graphs.target_graph_rows, g, max_graphs, v,
                    data.largest_pattern_clique[g], data.target_cliques_best_knowns[0], data.target_cliques_build_times,
                    data.target_cliques_solve_times, data.target_cliques_solve_find_nodes, data.target_cliques_solve_prove_nodes);
            }
    }
}

auto gss::innards::init_clique_size_data(CliqueSizeData & data, const HomomorphismParams & params,
    unsigned max_graphs, unsigned pattern_size, unsigned target_size) -> void
{
    if (! params.clique_size_constraints)
        return;

    data.max_graphs_for_clique_size_constraints = (params.clique_size_constraints_on_supplementals ? max_graphs : 1);
    for (unsigned g = 0; g < data.max_graphs_for_clique_size_constraints; ++g) {
        data.pattern_cliques_sizes.push_back(vector<int>(pattern_size, 0));
        data.target_cliques_sizes.push_back(vector<int>(target_size, 0));
        data.pattern_cliques_best_knowns.push_back(vector<int>(pattern_size, 0));
        data.target_cliques_best_knowns.push_back(vector<int>(target_size, 0));
    }
    data.largest_pattern_clique.resize(data.max_graphs_for_clique_size_constraints);
}

auto gss::innards::check_clique_compatibility(CliqueSizeData & data, const ProcessedGraphsData & graphs,
    unsigned max_graphs, unsigned pattern_size, unsigned target_size, const HomomorphismParams & params,
    HomomorphismProofs * proofs, int p, int t) -> bool
{
    if (! params.clique_size_constraints)
        return true;

    if (! data.has_pattern_cliques_sizes)
        build_pattern_clique_sizes(data, graphs, max_graphs, pattern_size, params.timeout);

    build_target_clique_size(data, graphs, max_graphs, target_size, params.timeout, t);

    for (unsigned g = 0; g < data.max_graphs_for_clique_size_constraints; ++g) {
        if (data.pattern_cliques_sizes[g][p] > data.target_cliques_sizes[g][t]) {
            if (proofs)
                proofs->prove_no_clique(graphs, max_graphs, pattern_size, target_size, params, g, p, t);
            return false;
        }
    }

    return true;
}
