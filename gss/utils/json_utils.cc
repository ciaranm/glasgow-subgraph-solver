#include "json_utils.hh"
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

    json extra_stats_to_json(const HomomorphismResult & result)
    {
        json j;
        const std::set<std::string> numeric_keys = {
            "restarts", "shape_graphs", "search_time", "nogoods_size", "nogoods_lengths"};

        for (const auto & stat : result.extra_stats) {
            auto pos = stat.find('=');
            if (pos == std::string::npos) continue;

            std::string key = stat.substr(0, pos);
            std::string value = stat.substr(pos + 1);

            key.erase(key.find_last_not_of(" \t") + 1);
            value.erase(0, value.find_first_not_of(" \t"));

            if (numeric_keys.count(key)) {
                j[key] = value.empty() ? nullptr : nlohmann::json(std::stoll(value));
            }
            else {
                j[key] = value;
            }
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

    json make_solver_json(
        int argc,
        char ** argv,
        const string & pattern_file,
        const string & target_file,
        const InputGraph & pattern,
        const InputGraph & target,
        const HomomorphismResult & result,
        const HomomorphismParams & params,
        const std::chrono::milliseconds & overall_time,
        const string & status)
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

        j["target_file"] = target_file;
        j["taget_properties"] = describe(target);
        j["target_vertices"] = target.size();
        j["target_directed_edges"] = target.number_of_directed_edges();

        j["nodes"] = result.nodes;
        j["propagations"] = result.propagations;

        if (params.count_solutions) {
            j["solution_count"] = params.count_solutions;
        }
        if (! result.extra_stats.empty()) {
            j.update(extra_stats_to_json(result));
        }
        if (! result.mapping.empty()) {
            j.update(mapping_to_json(result, pattern, target));
        }

        j["runtime"] = overall_time.count();

        return j;
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
