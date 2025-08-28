#include "json_utils.hh"
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>
#include <string>

using namespace std::chrono;

using json = nlohmann::json;

namespace gss::utils
{

    json commandline_to_json(int argc, char ** argv)
    {
        std::ostringstream cmd;
        for (int i = 0; i < argc; ++i) {
            if (i > 0) cmd << " ";
            cmd << argv[i];
        }
        return cmd.str();
    }

    template <typename ResultT>
    json extra_stats_to_json(const ResultT & result)
    {
        json j = json::object();

        std::set<std::string> numeric_keys;

        if constexpr (std::is_same_v<ResultT, HomomorphismResult>) {
            numeric_keys = {"restarts", "shape_graphs", "search_time", "nogoods_size", "nogoods_lengths"};
        }
        else if constexpr (std::is_same_v<ResultT, CommonSubgraphResult>) {
            numeric_keys = {"search_time", "subgraph_nodes"};
        }
        else if constexpr (std::is_same_v<ResultT, CliqueResult>) {
            numeric_keys = {"nodes", "find_nodes", "prove_nodes"};
        }

        for (const auto & stat : result.extra_stats) {
            auto pos = stat.find('=');
            if (pos == std::string::npos) continue;

            std::string key = stat.substr(0, pos);
            std::string value = stat.substr(pos + 1);

            key.erase(key.find_last_not_of(" \t") + 1);
            value.erase(0, value.find_first_not_of(" \t"));

            if (numeric_keys.contains(key))
                j[key] = std::stoll(value);
            else
                j[key] = value;
        }
        return j;
    }

    json mapping_to_json(
        const HomomorphismResult & result,
        const InputGraph & pattern,
        const InputGraph & target)
    {
        json mapping_json = json::array();
        for (auto v : result.mapping) {
            mapping_json.push_back({{"pattern_vertex", pattern.vertex_name(v.first)},
                {"target_vertex", target.vertex_name(v.second)}});
        }
        return mapping_json;
    }

    template <typename ResultT, typename ParamsT>
    void make_solver_json(
        int argc,
        char ** argv,
        const string & pattern_file,
        const string & target_file,
        const InputGraph & pattern,
        const InputGraph & target,
        const ResultT & result,
        const ParamsT & params,
        const std::chrono::milliseconds & overall_time,
        const string & status,
        string & filename)
    {
        json j;

        j["commandline"] = commandline_to_json(argc, argv);

        auto started_at = system_clock::to_time_t(system_clock::now());
        std::ostringstream ts;
        ts << std::put_time(std::localtime(&started_at), "%F %T");
        j["started_at"] = ts.str();
        j["status"] = status;

        j["nodes"] = result.nodes;
        j["runtime"] = overall_time.count();

        if constexpr (std::is_same_v<ResultT, HomomorphismResult>) {
            j["pattern_file"] = pattern_file;
            j["pattern_properties"] = describe(pattern);
            j["pattern_vertices"] = pattern.size();
            j["pattern_directed_edges"] = pattern.number_of_directed_edges();

            j["target_file"] = target_file;
            j["target_properties"] = describe(target);
            j["target_vertices"] = target.size();
            j["target_directed_edges"] = target.number_of_directed_edges();
            j["propagations"] = result.propagations;
            if (params.count_solutions) j["solution_count"] = result.solution_count.to_string();
            if (params.n_threads) j["threads"] = params.n_threads;
            if (! result.extra_stats.empty()) j.update(extra_stats_to_json(result));
            if (! result.mapping.empty()) j["mapping"] = mapping_to_json(result, pattern, target);
            if (! result.all_mappings.empty()) j["all_mappings"] = result.all_mappings;
        }
        else if constexpr (std::is_same_v<ResultT, CommonSubgraphResult>) {
            j["first_file"] = pattern_file;
            j["first_properties"] = describe(pattern);
            j["first_vertices"] = pattern.size();
            j["first_directed_edges"] = pattern.number_of_directed_edges();

            j["second_file"] = target_file;
            j["second_properties"] = describe(target);
            j["second_vertices"] = target.size();
            j["second_directed_edges"] = target.number_of_directed_edges();

            j["nodes"] = result.nodes;
            if (params.count_solutions) j["solution_count"] = result.solution_count.to_string();
            if (! result.extra_stats.empty()) j.update(extra_stats_to_json(result));
            if (! result.mapping.empty()) j["mappings"] = result.mapping;
            if (! result.all_mappings.empty()) j["all_mappings"] = result.all_mappings;
        }

        if (! filename.ends_with(".json"))
            filename += ".json";
        std::ofstream out(filename);
        out << j.dump(4) << std::endl;
    }

    template void make_solver_json<HomomorphismResult, HomomorphismParams>(
        int, char **, const string &, const string &, const InputGraph &, const InputGraph &,
        const HomomorphismResult &, const HomomorphismParams &, const std::chrono::milliseconds &, const string &, string &);

    template void make_solver_json<CommonSubgraphResult, CommonSubgraphParams>(
        int, char **, const string &, const string &, const InputGraph &, const InputGraph &,
        const CommonSubgraphResult &, const CommonSubgraphParams &, const std::chrono::milliseconds &, const string &, string &);

    void make_solver_json(
        int argc,
        char ** argv,
        const string & pattern_file,
        const InputGraph & pattern,
        const CliqueResult & result,
        const CliqueParams & params,
        const std::chrono::milliseconds & overall_time,
        const string & status,
        string & filename)
    {
        json j;

        j["commandline"] = commandline_to_json(argc, argv);

        auto started_at = system_clock::to_time_t(system_clock::now());
        std::ostringstream ts;
        ts << std::put_time(std::localtime(&started_at), "%F %T");
        j["started_at"] = ts.str();
        j["status"] = status;

        j["pattern_file"] = pattern_file;
        j["pattern_properties"] = describe(pattern);
        j["pattern_vertices"] = pattern.size();
        j["pattern_directed_edges"] = pattern.number_of_directed_edges();

        j["clique"] = result.clique;
        j["omega"] = result.clique.size();
        j["nodes"] = result.nodes;
        j["prove_nodes"] = result.prove_nodes;
        j["find_nodes"] = result.find_nodes;

        j["runtime"] = overall_time.count();

        if (! result.extra_stats.empty()) {
            j.update(extra_stats_to_json(result));
        }

        if (! filename.ends_with(".json"))
            filename += ".json";
        std::ofstream out(filename);
        out << j.dump(4) << std::endl;
    }

    std::string describe(const InputGraph & g)
    {
        std::string shape_group;
        if (g.directed())
            shape_group = " directed";
        if (g.loopy())
            shape_group = " loopy";
        if (g.has_vertex_labels())
            shape_group = " vertex_labels";
        if (g.has_edge_labels())
            shape_group = " edge_labels";
        return shape_group;
    };
}
