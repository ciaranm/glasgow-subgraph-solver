/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "verify.hh"
#include <map>

using std::map;
using std::string;

BuggySolution::BuggySolution(const string & message) throw () :
    _what(message)
{
}

auto BuggySolution::what() const throw () -> const char *
{
    return _what.c_str();
}

auto verify(const std::pair<InputGraph, InputGraph> & graphs, const Params & params, const Result & result) -> void
{
    // nothing to verify, if unsat
    if (result.mapping.empty())
        return;

    // totality and range
    for (int i = 0 ; i < graphs.first.size() ; ++i) {
        if (! result.mapping.count(i))
            throw BuggySolution{ "No mapping for vertex " + graphs.first.vertex_name(i) };
        else if (result.mapping.find(i)->second < 0 || result.mapping.find(i)->second >= graphs.second.size())
            throw BuggySolution{ "Mapping " + graphs.first.vertex_name(i) + " -> " +
                graphs.second.vertex_name(result.mapping.find(i)->second) + " out of range" };
    }

    // no extra stuff
    for (auto & i : result.mapping)
        if (i.first < 0 || i.first >= graphs.first.size())
            throw BuggySolution{ "Vertex " + graphs.first.vertex_name(i.first)  + " out of range" };

    // labels
    if (graphs.first.has_vertex_labels())
        for (int i = 0 ; i < graphs.first.size() ; ++i)
            if (graphs.first.vertex_label(i) != graphs.second.vertex_label(result.mapping.find(i)->second))
                throw BuggySolution{ "Mismatched vertex label for assignment " + graphs.first.vertex_name(i) + " -> " +
                    graphs.second.vertex_name(result.mapping.find(i)->second) };

    // injectivity
    if (! params.noninjective) {
        map<int, int> seen;
        for (auto & [ i, j ] : result.mapping) {
            if (! seen.emplace(j, i).second)
                throw BuggySolution{ "Non-injective mapping: " + graphs.first.vertex_name(i) + " -> " +
                    graphs.second.vertex_name(result.mapping.find(i)->second) + " and " +
                    graphs.second.vertex_name(seen.find(j)->second) + " -> " + graphs.first.vertex_name(j) };
        }
    }

    // loops
    for (int i = 0 ; i < graphs.first.size() ; ++i) {
        if (graphs.first.adjacent(i, i) && ! graphs.second.adjacent(result.mapping.find(i)->second, result.mapping.find(i)->second))
            throw BuggySolution{ "Vertex " + graphs.first.vertex_name(i) + " has a loop but mapped vertex " +
                graphs.second.vertex_name(result.mapping.find(i)->second) + " does not" };
        else if (params.induced && graphs.second.adjacent(result.mapping.find(i)->second, result.mapping.find(i)->second &&
                    ! graphs.first.adjacent(i, i)))
            throw BuggySolution{ "Vertex " + graphs.first.vertex_name(i) + " has no loop but mapped vertex " +
                graphs.second.vertex_name(result.mapping.find(i)->second) + " does" };
    }

    // adjacency, non-adjacency, and edge labels
    for (auto & [ i, t ] : result.mapping) {
        for (auto & [ j, u ] : result.mapping) {
            if (graphs.first.adjacent(i, j) && ! graphs.second.adjacent(t, u))
                throw BuggySolution{ "Edge " + graphs.first.vertex_name(i) + " -- " + graphs.first.vertex_name(j) +
                    " mapped to non-edge " + graphs.second.vertex_name(t) + " -/- " + graphs.second.vertex_name(u) };
            else if (params.induced && ! graphs.first.adjacent(i, j) && graphs.second.adjacent(t, u))
                throw BuggySolution{ "Non-edge " + graphs.first.vertex_name(i) + " -/- " + graphs.first.vertex_name(j) +
                    " mapped to edge " + graphs.second.vertex_name(t) + " -- " + graphs.second.vertex_name(u) };
        }
    }
}

