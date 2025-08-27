#ifndef GLASGOW_SUBGRAPH_SOLVER_JSON_UTILS_HH
#define GLASGOW_SUBGRAPH_SOLVER_JSON_UTILS_HH

#include <gss/homomorphism.hh>

#include <nlohmann/json.hpp>
#include <string>

using json = nlohmann::json;
using string = std::string;

namespace gss::utils
{

    /// Build commandline as a single string
    json commandline_to_json(int argc, char ** argv);

    /// Parse solver extra_stats into structured JSON
    json extra_stats_to_json(const HomomorphismResult & result);

    // Build json from result.mapping
    json mapping_to_json(
        const HomomorphismResult & result,
        const InputGraph & pattern,
        const InputGraph & target);

    /// Returns json for solver runs:
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
        const string & status);

    std::string describe(const InputGraph & g);
}

#endif // GLASGOW_SUBGRAPH_SOLVER_JSON_UTILS_HH
