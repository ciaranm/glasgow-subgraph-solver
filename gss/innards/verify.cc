#include <gss/innards/verify.hh>
#include <map>

using namespace gss;
using namespace gss::innards;

using std::map;
using std::string;

BuggySolution::BuggySolution(const string & message) noexcept :
    _what(message)
{
}

auto BuggySolution::what() const noexcept -> const char *
{
    return _what.c_str();
}

auto gss::innards::verify_homomorphism(
    const InputGraph & pattern,
    const InputGraph & target,
    bool injective,
    bool locally_injective,
    bool induced,
    const map<int, int> & mapping) -> void
{
    // nothing to verify, if unsat
    if (mapping.empty())
        return;

    // totality and range
    for (int i = 0; i < pattern.size(); ++i) {
        if (! mapping.count(i))
            throw BuggySolution{"No mapping for vertex " + pattern.vertex_name(i)};
        else if (mapping.find(i)->second < 0 || mapping.find(i)->second >= target.size())
            throw BuggySolution{"Mapping " + pattern.vertex_name(i) + " -> " +
                target.vertex_name(mapping.find(i)->second) + " out of range"};
    }

    // no extra stuff
    for (auto & i : mapping)
        if (i.first < 0 || i.first >= pattern.size())
            throw BuggySolution{"Vertex " + pattern.vertex_name(i.first) + " out of range"};

    // labels
    if (pattern.has_vertex_labels())
        for (int i = 0; i < pattern.size(); ++i)
            if (pattern.vertex_label(i) != target.vertex_label(mapping.find(i)->second))
                throw BuggySolution{"Mismatched vertex label for assignment " + pattern.vertex_name(i) + " -> " +
                    target.vertex_name(mapping.find(i)->second)};

    // injectivity
    if (injective) {
        map<int, int> seen;
        for (auto & [i, j] : mapping) {
            if (! seen.emplace(j, i).second)
                throw BuggySolution{"Non-injective mapping: " + pattern.vertex_name(i) + " -> " +
                    target.vertex_name(mapping.find(i)->second) + " and " +
                    pattern.vertex_name(seen.find(j)->second) + " -> " + target.vertex_name(j)};
        }
    }

    // local injectivity
    if (locally_injective) {
        for (int v = 0; v < pattern.size(); ++v) {
            map<int, int> seen;
            for (auto & [i, j] : mapping) {
                if (pattern.adjacent(v, i) && ! seen.emplace(j, i).second)
                    throw BuggySolution{"Non locally-injective mapping: on neighbourhood of " + pattern.vertex_name(v) + ", " + pattern.vertex_name(i) + " -> " +
                        target.vertex_name(mapping.find(i)->second) + " and " +
                        pattern.vertex_name(seen.find(j)->second) + " -> " + target.vertex_name(j)};
            }
        }
    }

    // loops
    for (int i = 0; i < pattern.size(); ++i) {
        if (pattern.adjacent(i, i) && ! target.adjacent(mapping.find(i)->second, mapping.find(i)->second))
            throw BuggySolution{"Vertex " + pattern.vertex_name(i) + " has a loop but mapped vertex " +
                target.vertex_name(mapping.find(i)->second) + " does not"};
        else if (induced && target.adjacent(mapping.find(i)->second, mapping.find(i)->second) &&
            ! pattern.adjacent(i, i))
            throw BuggySolution{"Vertex " + pattern.vertex_name(i) + " has no loop but mapped vertex " +
                target.vertex_name(mapping.find(i)->second) + " does"};
    }

    // adjacency, non-adjacency, and edge labels
    for (auto & [i, t] : mapping) {
        for (auto & [j, u] : mapping) {
            if (pattern.adjacent(i, j) && ! target.adjacent(t, u))
                throw BuggySolution{"Edge " + pattern.vertex_name(i) + " -- " + pattern.vertex_name(j) +
                    " mapped to non-edge " + target.vertex_name(t) + " -/- " + target.vertex_name(u)};
            else if (induced && ! pattern.adjacent(i, j) && target.adjacent(t, u))
                throw BuggySolution{"Non-edge " + pattern.vertex_name(i) + " -/- " + pattern.vertex_name(j) +
                    " mapped to edge " + target.vertex_name(t) + " -- " + target.vertex_name(u)};
        }
    }
}
