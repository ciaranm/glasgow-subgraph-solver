#include <gss/formats/csv.hh>
#include <gss/formats/input_graph.hh>

#include <fstream>
#include <istream>
#include <optional>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

using std::ifstream;
using std::istream;
using std::nullopt;
using std::optional;
using std::pair;
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
        optional<pair<string, string>> first_unlabelled_edge;

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
                else if (! first_unlabelled_edge)
                    first_unlabelled_edge = pair{left, right};

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

        // Labels are all or nothing, per element type: an element with no declared
        // label is not the same thing as one labelled with the empty string, and
        // since label matching is exact, quietly treating it as the latter would
        // silently constrain it to unlabelled target elements.

        // For vertices, report the lowest-numbered offender, so the message doesn't
        // depend on hash order. A vertex mentioned only by an edge line declares no
        // label either, so this has to look at every vertex rather than at the
        // declaration lines.
        if (seen_vertex_label) {
            optional<pair<int, string>> unlabelled;
            for (auto & [v, idx] : vertices)
                if (! vertex_labels.contains(v))
                    if ((! unlabelled) || idx < unlabelled->first)
                        unlabelled = pair{idx, v};

            if (unlabelled)
                throw GraphFileError{filename, "vertex '" + unlabelled->second + "' has no label, but other vertices do: vertex labels must be given for every vertex, or for none",
                    true};
        }

        // Every edge comes from a line of its own, so for edges this is the first
        // offending line, named as it was written.
        if (seen_edge_label && first_unlabelled_edge)
            throw GraphFileError{filename, "the edge between '" + first_unlabelled_edge->first + "' and '" + first_unlabelled_edge->second + "' has no label, but other edges do: edge labels must be given for every edge, or for none",
                true};

        InputGraph result{int(vertices.size()), seen_vertex_label, seen_edge_label, seen_directed_edge};

        // Note that the undirected case has both (f, t) and (t, f) in edges already,
        // so add_edge() is called once for each direction: harmless, and it keeps
        // labelled and unlabelled undirected edges on the same path. Using
        // add_directed_edge() here instead would mark the graph directed.
        for (auto & [f, t, l] : edges)
            if (seen_directed_edge)
                result.add_directed_edge(f, t, l);
            else
                result.add_edge(f, t, l);

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
                result.set_vertex_label(l, vertex_labels.at(v));

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
