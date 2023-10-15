#include <gss/formats/csv.hh>
#include <gss/formats/input_graph.hh>

#include <fstream>
#include <istream>
#include <optional>
#include <tuple>
#include <unordered_map>
#include <vector>

using std::ifstream;
using std::istream;
using std::nullopt;
using std::optional;
using std::string;
using std::tuple;
using std::unordered_map;
using std::vector;

namespace
{
    auto read_csv(istream && infile, const string & filename, const optional<unordered_map<string, string>> & rename_map) -> InputGraph
    {
        if (! infile)
            throw GraphFileError{filename, "error opening file", false};

        unordered_map<string, int> vertices;
        unordered_map<string, string> vertex_labels;
        vector<tuple<int, int, string>> edges;
        bool seen_vertex_label = false, seen_edge_label = false, seen_directed_edge = false;

        string line;

        while (getline(infile, line)) {
            auto pos = line.find_first_of(",>");
            if (string::npos == pos)
                throw GraphFileError{filename, "expected a comma but didn't get one", true};

            string left = line.substr(0, pos), right = line.substr(pos + 1), label;
            char delim = line.at(pos);

            auto pos2 = right.find(',');
            if (string::npos != pos2) {
                label = right.substr(pos2 + 1);
                right = right.substr(0, pos2);
            }

            if (right.empty() && ! left.empty()) {
                if (! label.empty()) {
                    seen_vertex_label = true;
                    vertex_labels.emplace(left, label);
                }

                vertices.emplace(left, vertices.size());
            }
            else {
                int left_idx = vertices.emplace(left, vertices.size()).first->second;
                int right_idx = vertices.emplace(right, vertices.size()).first->second;

                if (! label.empty())
                    seen_edge_label = true;

                if (delim == '>') {
                    seen_directed_edge = true;
                    edges.emplace_back(left_idx, right_idx, label);
                }
                else {
                    edges.emplace_back(left_idx, right_idx, label);
                    edges.emplace_back(right_idx, left_idx, label);
                }
            }
        }

        InputGraph result{int(vertices.size()), seen_vertex_label, seen_edge_label};

        for (auto & [f, t, l] : edges)
            if (seen_directed_edge)
                result.add_directed_edge(f, t, l);
            else if (seen_edge_label) {
                result.add_directed_edge(f, t, l);
                result.add_directed_edge(t, f, l);
            }
            else
                result.add_edge(f, t);

        auto rename = [&](const string & s) -> string {
            if (rename_map) {
                auto r = rename_map->find(s);
                if (r == rename_map->end())
                    throw GraphFileError{filename, "did not find a name for vertex '" + s + "'", true};
                return r->second;
            }
            else
                return s;
        };

        for (auto & [v, l] : vertices)
            result.set_vertex_name(l, rename(v));

        if (seen_vertex_label)
            for (auto & [v, l] : vertices)
                result.set_vertex_label(l, vertex_labels[v]);

        return result;
    }
}

auto read_csv(istream && infile, const string & filename) -> InputGraph
{
    return read_csv(move(infile), filename, nullopt);
}

auto read_csv_name(std::istream && infile, const std::string & filename, const std::string & name_map_filename) -> InputGraph
{
    ifstream name_map_file{name_map_filename};
    if (! name_map_file)
        throw GraphFileError{name_map_filename, "could not open rename map file", false};

    optional<unordered_map<string, string>> rename_map{unordered_map<string, string>{}};

    string line;
    while (getline(name_map_file, line)) {
        auto pos = line.find(',');
        if (string::npos == pos)
            throw GraphFileError{filename, "expected a comma but didn't get one", true};
        string left = line.substr(0, pos), right = line.substr(pos + 1);
        rename_map->emplace(left, right);
    }

    return read_csv(move(infile), filename, rename_map);
}
