/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "formats/csv.hh"
#include "formats/input_graph.hh"

#include <fstream>
#include <unordered_map>
#include <vector>

using std::ifstream;
using std::pair;
using std::string;
using std::unordered_map;
using std::vector;

auto read_csv(ifstream && infile, const string & filename) -> InputGraph
{
    InputGraph result{ 0, false, false };

    if (! infile)
        throw GraphFileError{ filename, "error opening file" };

    unordered_map<string, int> vertices;
    string line;

    vector<pair<int, int> > edges;

    while (getline(infile, line)) {
        auto pos = line.find(',');
        if (string::npos == pos)
            throw GraphFileError{ filename, "expected a comma but didn't get one" };
        string left = line.substr(0, pos), right = line.substr(pos + 1);
        int left_idx = vertices.emplace(left, vertices.size()).first->second;
        int right_idx = vertices.emplace(right, vertices.size()).first->second;
        edges.emplace_back(left_idx, right_idx);
    }

    result.resize(vertices.size());

    for (auto & e : edges)
        result.add_edge(e.first, e.second);

    for (auto & [v, l] : vertices)
        result.set_vertex_label(l, v);

    return result;
}

