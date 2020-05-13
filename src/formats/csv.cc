/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "formats/csv.hh"
#include "formats/input_graph.hh"

#include <fstream>
#include <optional>
#include <unordered_map>
#include <vector>

using std::ifstream;
using std::nullopt;
using std::optional;
using std::pair;
using std::string;
using std::unordered_map;
using std::vector;

namespace
{
    auto read_csv(ifstream && infile, const string & filename, const optional<unordered_map<string, string> > & rename_map) -> InputGraph
    {
        InputGraph result{ 0, false, false };

        if (! infile)
            throw GraphFileError{ filename, "error opening file", false };

        unordered_map<string, int> vertices;
        string line;

        vector<pair<int, int> > edges;

        while (getline(infile, line)) {
            auto pos = line.find(',');
            if (string::npos == pos)
                throw GraphFileError{ filename, "expected a comma but didn't get one", true };
            string left = line.substr(0, pos), right = line.substr(pos + 1);
            if (right.empty() && ! left.empty() && ! vertices.count(left))
                vertices.emplace(left, vertices.size());
            else {
                int left_idx = vertices.emplace(left, vertices.size()).first->second;
                int right_idx = vertices.emplace(right, vertices.size()).first->second;
                edges.emplace_back(left_idx, right_idx);
            }
        }

        result.resize(vertices.size());

        for (auto & e : edges)
            result.add_edge(e.first, e.second);

        auto rename = [&] (const string & s) -> string {
            if (rename_map) {
                auto r = rename_map->find(s);
                if (r == rename_map->end())
                    throw GraphFileError{ filename, "did not find a name for vertex '" + s + "'", true };
                return r->second;
            }
            else
                return s;
        };

        for (auto & [v, l] : vertices)
            result.set_vertex_name(l, rename(v));

        return result;
    }
}

auto read_csv(ifstream && infile, const string & filename) -> InputGraph
{
    return read_csv(move(infile), filename, nullopt);
}

auto read_csv_name(std::ifstream && infile, const std::string & filename, const std::string & name_map_filename) -> InputGraph
{
    ifstream name_map_file{ name_map_filename };
    if (! name_map_file)
        throw GraphFileError{ name_map_filename, "could not open rename map file", false };

    optional<unordered_map<string, string> > rename_map{ unordered_map<string, string>{ } };

    string line;
    while (getline(name_map_file, line)) {
        auto pos = line.find(',');
        if (string::npos == pos)
            throw GraphFileError{ filename, "expected a comma but didn't get one", true };
        string left = line.substr(0, pos), right = line.substr(pos + 1);
        rename_map->emplace(left, right);
    }

    return read_csv(move(infile), filename, rename_map);
}


