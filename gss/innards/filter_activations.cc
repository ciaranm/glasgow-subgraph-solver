#include <gss/innards/filter_activations.hh>

#include <numeric>
#include <string>

using namespace gss;
using namespace gss::innards;

using std::accumulate;
using std::list;
using std::string;
using std::to_string;

namespace
{
    auto join(const std::vector<unsigned long long> & values, unsigned from) -> string
    {
        string result;
        for (unsigned g = from; g < values.size(); ++g) {
            if (g != from)
                result += ",";
            result += to_string(values[g]);
        }
        return result.empty() ? "-" : result;
    }
}

FilterActivations::FilterActivations(unsigned max_graphs) :
    degree(max_graphs, 0),
    search(max_graphs, 0)
{
}

auto FilterActivations::any_optional_filter_fired() const -> bool
{
    // Degree on the original graph counts: --no-nds leaves it on, but the locally
    // injective and non-injective modes switch it off entirely, which is exactly the kind
    // of silently-disabled filtering the sweep is looking for.
    return nds != 0 || cliques != 0 ||
        0 != accumulate(degree.begin(), degree.end(), 0ull) ||
        (search.size() > 1 && 0 != accumulate(search.begin() + 1, search.end(), 0ull));
}

auto FilterActivations::add_initial_stats(list<string> & stats) const -> void
{
    stats.emplace_back("filter_activations_initial = vertex_labels:" + to_string(vertex_labels) +
        " loops:" + to_string(loops) +
        " degree:" + join(degree, 0) +
        " nds:" + to_string(nds) +
        " cliques:" + to_string(cliques));
}

auto FilterActivations::add_search_stats(list<string> & stats) const -> void
{
    stats.emplace_back("filter_activations_search = adjacency:" + (search.empty() ? string{"-"} : to_string(search[0])) +
        " edge_labels:" + to_string(search_edge_labels) +
        " supplemental:" + join(search, 1));
}
