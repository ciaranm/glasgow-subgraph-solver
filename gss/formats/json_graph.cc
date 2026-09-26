#include <gss/formats/json_graph.hh>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <iterator>
#include <limits>
#include <map>
#include <optional>
#include <ostream>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using std::find;
using std::get;
using std::istream;
using std::istreambuf_iterator;
using std::map;
using std::nullopt;
using std::optional;
using std::ostream;
using std::pair;
using std::sort;
using std::string;
using std::to_string;
using std::tuple;
using std::vector;

using json = nlohmann::json;

namespace
{
    // Bumped only when the format changes in a way an older reader would get wrong.
    // A reader refuses anything above what it knows rather than guessing at it.
    const long long known_version = 1;

    // How an edge names its endpoints. A file picks one and sticks to it: this is
    // the thing JSON buys us over every text format we read, since the vertex named
    // "1" and the vertex at index 1 are different tokens and so stop being a trap.
    enum class Addressing
    {
        by_index,
        by_name
    };

    auto addressing_name(Addressing a) -> string
    {
        return a == Addressing::by_index ? "index" : "name";
    }

    auto quoted_string(const string & s) -> string
    {
        return json(s).dump();
    }

    // Keys beginning "x-" are reserved for third-party annotation and ignored
    // wherever they appear. Anything else unrecognised is an error, so that a
    // misspelled key fails loudly instead of being silently dropped.
    auto check_keys(const json & object, const string & where, const vector<string> & allowed, const string & filename) -> void
    {
        for (const auto & item : object.items()) {
            const auto & key = item.key();
            if (key.starts_with("x-"))
                continue;

            if (find(allowed.begin(), allowed.end(), key) == allowed.end()) {
                string list;
                for (const auto & a : allowed)
                    list += (list.empty() ? "" : ", ") + quoted_string(a);
                throw GraphFileError{filename, "unrecognised key " + quoted_string(key) + " in " + where + " (this format allows " + list + " there, plus anything beginning \"x-\")",
                    true};
            }
        }
    }

    auto require(const json & object, const string & key, const string & where, const string & filename) -> const json &
    {
        auto it = object.find(key);
        if (it == object.end())
            throw GraphFileError{filename, "missing required key " + quoted_string(key) + " in " + where, true};
        return *it;
    }

    auto as_bool(const json & value, const string & what, const string & filename) -> bool
    {
        if (! value.is_boolean())
            throw GraphFileError{filename, what + " must be true or false", true};
        return value.get<bool>();
    }

    auto as_string(const json & value, const string & what, const string & filename) -> string
    {
        if (! value.is_string())
            throw GraphFileError{filename, what + " must be a string", true};
        return value.get<string>();
    }

    // Costs are 64-bit signed integers. JSON has no integer width, so anything that
    // would not fit is refused rather than wrapped.
    auto as_cost(const json & value, const string & what, const string & filename) -> long long
    {
        if (value.is_number_unsigned() && value.get<unsigned long long>() > (unsigned long long)std::numeric_limits<long long>::max())
            throw GraphFileError{filename, what + " does not fit in a 64-bit signed integer", true};
        if (! value.is_number_integer())
            throw GraphFileError{filename, what + " must be an integer", true};
        return value.get<long long>();
    }

    struct VertexData
    {
        optional<string> name;
        optional<string> label;
        optional<long long> cost;
    };

    auto read_vertices(const json & doc, const string & filename) -> vector<VertexData>
    {
        const auto & vertices = require(doc, "vertices", "the top-level object", filename);

        // A count means that many anonymous vertices, addressed by index.
        if (vertices.is_number_integer()) {
            auto count = vertices.get<long long>();
            if (count < 0)
                throw GraphFileError{filename, "\"vertices\" is negative (" + to_string(count) + ")", true};
            return vector<VertexData>(count);
        }

        if (! vertices.is_array())
            throw GraphFileError{filename, "\"vertices\" must be a count or an array of vertices", true};

        vector<VertexData> result;
        for (size_t v = 0; v < vertices.size(); ++v) {
            const auto & entry = vertices.at(v);
            auto where = "vertex " + to_string(v);

            // A bare string is sugar for {"name": ...}, exactly.
            if (entry.is_string()) {
                result.push_back(VertexData{entry.get<string>(), nullopt});
                continue;
            }

            if (! entry.is_object())
                throw GraphFileError{filename, where + " must be a name or an object", true};

            check_keys(entry, where, {"name", "label", "cost"}, filename);

            VertexData data;
            if (auto it = entry.find("name"); it != entry.end())
                data.name = as_string(*it, where + "'s \"name\"", filename);
            if (auto it = entry.find("label"); it != entry.end())
                data.label = as_string(*it, where + "'s \"label\"", filename);
            if (auto it = entry.find("cost"); it != entry.end())
                data.cost = as_cost(*it, where + "'s \"cost\"", filename);
            result.push_back(std::move(data));
        }

        return result;
    }

    // Labels are all or nothing per element type: absent and "" are different
    // things, "" being a real label, so a file that labels only some of its
    // vertices or only some of its edges cannot be read as meaning either one.
    auto check_all_or_nothing(size_t labelled, size_t total, const string & what, const string & filename) -> bool
    {
        if (labelled != 0 && labelled != total)
            throw GraphFileError{filename, to_string(labelled) + " of " + to_string(total) + " " + what + " carry a \"label\": labels are all or nothing (an absent label is not the same as \"\", which is itself a label)",
                true};
        return labelled != 0;
    }

    // Costs likewise: an absent cost is not a cost of 0, and reading it as one would
    // make an element free rather than making the file an error.
    auto check_costs_all_or_nothing(size_t costed, size_t total, const string & what, const string & filename) -> bool
    {
        if (costed != 0 && costed != total)
            throw GraphFileError{filename, to_string(costed) + " of " + to_string(total) + " " + what + " carry a \"cost\": costs are all or nothing (an absent cost is not the same as 0)",
                true};
        return costed != 0;
    }

    struct EdgeData
    {
        int from, to;
        optional<string> label;
        optional<long long> cost;
    };

    auto resolve_endpoint(const json & endpoint, const string & what, size_t vertex_count,
        const map<string, int> & by_name, optional<Addressing> & addressing, const string & filename) -> int
    {
        auto fix_addressing = [&](Addressing found) {
            if (! addressing)
                addressing = found;
            else if (*addressing != found)
                throw GraphFileError{filename, what + " addresses its vertex by " + addressing_name(found) + ", but an earlier endpoint used " + addressing_name(*addressing) + ": a file must address vertices one way throughout",
                    true};
        };

        if (endpoint.is_number_integer()) {
            fix_addressing(Addressing::by_index);
            auto index = endpoint.get<long long>();
            if (index < 0 || size_t(index) >= vertex_count)
                throw GraphFileError{filename, what + " is index " + to_string(index) + ", but the graph has " + to_string(vertex_count) + " vertices",
                    true};
            return int(index);
        }

        if (endpoint.is_string()) {
            fix_addressing(Addressing::by_name);
            auto name = endpoint.get<string>();
            auto it = by_name.find(name);
            if (it == by_name.end())
                throw GraphFileError{filename, what + " names vertex " + quoted_string(name) + ", which is not in \"vertices\"", true};
            return it->second;
        }

        throw GraphFileError{filename, what + " must be an integer index or a string name", true};
    }

    auto read_edges(const json & doc, size_t vertex_count, const map<string, int> & by_name,
        bool directed, bool multigraph, const string & filename) -> vector<EdgeData>
    {
        const auto & edges = require(doc, "edges", "the top-level object", filename);
        if (! edges.is_array())
            throw GraphFileError{filename, "\"edges\" must be an array", true};

        vector<EdgeData> result;
        optional<Addressing> addressing;

        // Remembers where each edge was first given, so a duplicate can point at it. A
        // multigraph's edges are unique by (from, to, label), and a simple graph's by
        // (from, to), which is the same thing with every label taken to be equal.
        map<tuple<int, int, string>, size_t> already_seen;

        for (size_t e = 0; e < edges.size(); ++e) {
            const auto & entry = edges.at(e);
            auto where = "edge " + to_string(e);

            const json * from = nullptr;
            const json * to = nullptr;
            optional<string> label;
            optional<long long> cost;

            if (entry.is_array()) {
                // The positional form is frozen at two or three elements for good:
                // every key added in future lives only in the object form, so this
                // spelling can never come to mean something else.
                if (entry.size() != 2 && entry.size() != 3)
                    throw GraphFileError{filename, where + " has " + to_string(entry.size()) + " elements; the array form of an edge is [from, to] or [from, to, label] and is fixed at that, so anything further has to use the object form",
                        true};
                from = &entry.at(0);
                to = &entry.at(1);
                if (entry.size() == 3)
                    label = as_string(entry.at(2), where + "'s label", filename);
            }
            else if (entry.is_object()) {
                check_keys(entry, where, {"from", "to", "label", "cost", "multiplicity"}, filename);
                from = &require(entry, "from", where, filename);
                to = &require(entry, "to", where, filename);
                if (auto it = entry.find("label"); it != entry.end())
                    label = as_string(*it, where + "'s \"label\"", filename);
                if (auto it = entry.find("cost"); it != entry.end())
                    cost = as_cost(*it, where + "'s \"cost\"", filename);

                // Multiplicity is how the format expresses several edges with the same
                // label between one pair, so it is a key this version knows about and
                // refuses rather than one it has never heard of. Parallel edges with
                // different labels need only "multigraph": true.
                if (auto it = entry.find("multiplicity"); it != entry.end()) {
                    if (! it->is_number_integer() || it->get<long long>() < 1)
                        throw GraphFileError{filename, where + "'s \"multiplicity\" must be an integer of at least 1", true};
                    if (it->get<long long>() != 1)
                        throw GraphFileError{filename, where + " has \"multiplicity\" " + to_string(it->get<long long>()) + ", and several edges with the same label between one pair are not supported by this build",
                            true};
                }
            }
            else
                throw GraphFileError{filename, where + " must be an array or an object", true};

            EdgeData data;
            data.from = resolve_endpoint(*from, where + "'s \"from\"", vertex_count, by_name, addressing, filename);
            data.to = resolve_endpoint(*to, where + "'s \"to\"", vertex_count, by_name, addressing, filename);
            data.label = std::move(label);
            data.cost = cost;

            // In an undirected graph [u, v] and [v, u] are the same edge, so writing
            // both is a duplicate rather than being quietly idempotent as it is in
            // CSV. A self-loop [v, v] is one edge either way, never doubled.
            auto key = directed || data.from <= data.to
                ? tuple{data.from, data.to, multigraph ? data.label.value_or("") : string{}}
                : tuple{data.to, data.from, multigraph ? data.label.value_or("") : string{}};
            auto [it, inserted] = already_seen.emplace(key, e);
            if (! inserted)
                throw GraphFileError{filename, where + " repeats the edge already given as edge " + to_string(it->second) + (multigraph ? " (in a multigraph an edge is its endpoints and its label)" : "") + (directed ? "" : " (in an undirected graph [u, v] and [v, u] are the same edge)"),
                    true};

            result.push_back(std::move(data));
        }

        return result;
    }
}

auto read_json_graph(istream && infile, const string & filename) -> InputGraph
{
    if (! infile)
        throw GraphFileError{filename, "error opening file", false};

    // Parsed from a string rather than the stream so that trailing rubbish after
    // the closing brace is an error instead of being silently ignored.
    string text{istreambuf_iterator<char>{infile}, istreambuf_iterator<char>{}};

    json doc;
    try {
        doc = json::parse(text);
    }
    catch (const json::parse_error & e) {
        throw GraphFileError{filename, string{"could not be parsed as JSON: "} + e.what(), true};
    }

    if (! doc.is_object())
        throw GraphFileError{filename, "the top level must be a JSON object", true};

    check_keys(doc, "the top-level object", {"format", "version", "directed", "vertices", "edges", "multigraph"}, filename);

    auto format = as_string(require(doc, "format", "the top-level object", filename), "\"format\"", filename);
    if (format != "gss-graph")
        throw GraphFileError{filename, "\"format\" is " + quoted_string(format) + ", but this reader only knows \"gss-graph\"", true};

    const auto & version = require(doc, "version", "the top-level object", filename);
    if (! version.is_number_integer())
        throw GraphFileError{filename, "\"version\" must be an integer", true};
    auto file_version = version.get<long long>();
    if (file_version < 1)
        throw GraphFileError{filename, "\"version\" is " + to_string(file_version) + ", which is not a version", true};
    if (file_version > known_version)
        throw GraphFileError{filename, "this file is version " + to_string(file_version) + ", and this reader understands up to version " + to_string(known_version),
            true};

    // Optional, because its default is the restrictive reading: a file that leaves it
    // out and then gives two edges between one pair is an error, not a different graph.
    bool multigraph = false;
    if (auto it = doc.find("multigraph"); it != doc.end())
        multigraph = as_bool(*it, "\"multigraph\"", filename);

    // Never inferred: both readings of a file without this key parse cleanly and
    // give different answers, so there is no safe default to pick.
    auto directed = as_bool(require(doc, "directed", "the top-level object", filename), "\"directed\"", filename);

    auto vertices = read_vertices(doc, filename);

    map<string, int> by_name;
    size_t vertex_labelled = 0, vertex_costed = 0;
    for (size_t v = 0; v < vertices.size(); ++v) {
        if (vertices[v].name)
            if (! by_name.emplace(*vertices[v].name, int(v)).second)
                throw GraphFileError{filename, "two vertices are both named " + quoted_string(*vertices[v].name), true};
        if (vertices[v].label)
            ++vertex_labelled;
        if (vertices[v].cost)
            ++vertex_costed;
    }

    auto has_vertex_labels = check_all_or_nothing(vertex_labelled, vertices.size(), "vertices", filename);
    auto has_vertex_costs = check_costs_all_or_nothing(vertex_costed, vertices.size(), "vertices", filename);

    auto edges = read_edges(doc, vertices.size(), by_name, directed, multigraph, filename);

    size_t edge_labelled = 0, edge_costed = 0;
    for (const auto & e : edges) {
        if (e.label)
            ++edge_labelled;
        if (e.cost)
            ++edge_costed;
    }
    auto has_edge_labels = check_all_or_nothing(edge_labelled, edges.size(), "edges", filename);
    auto has_edge_costs = check_costs_all_or_nothing(edge_costed, edges.size(), "edges", filename);

    InputGraph result{int(vertices.size()), InputGraphProperties{.has_vertex_labels = has_vertex_labels, .has_edge_labels = has_edge_labels, .directed = directed, .multigraph = multigraph, .has_vertex_costs = has_vertex_costs, .has_edge_costs = has_edge_costs}};

    for (size_t v = 0; v < vertices.size(); ++v) {
        if (vertices[v].name)
            result.set_vertex_name(int(v), *vertices[v].name);
        if (vertices[v].label)
            result.set_vertex_label(int(v), *vertices[v].label);
        if (vertices[v].cost)
            result.set_vertex_cost(int(v), *vertices[v].cost);
    }

    for (const auto & e : edges) {
        auto label = e.label.value_or("");
        if (directed && e.cost)
            result.add_directed_edge(e.from, e.to, label, *e.cost);
        else if (directed)
            result.add_directed_edge(e.from, e.to, label);
        else if (e.cost)
            result.add_edge(e.from, e.to, label, *e.cost);
        else
            result.add_edge(e.from, e.to, label);
    }

    return result;
}

auto write_json_graph(ostream & outfile, const InputGraph & graph) -> void
{
    auto all_named = graph.size() != 0, any_named = false;
    for (int v = 0; v < graph.size(); ++v) {
        if (graph.vertex_has_name(v))
            any_named = true;
        else
            all_named = false;
    }

    // Names are only usable as addresses if every vertex has one.
    auto addressing = all_named ? Addressing::by_name : Addressing::by_index;

    auto endpoint = [&](int v) -> json {
        if (addressing == Addressing::by_name)
            return json(graph.vertex_name(v));
        return json(v);
    };

    outfile << "{\n";
    outfile << "  \"format\": \"gss-graph\",\n";
    outfile << "  \"version\": " << known_version << ",\n";
    outfile << "  \"directed\": " << (graph.directed() ? "true" : "false") << ",\n";
    if (graph.multigraph())
        outfile << "  \"multigraph\": true,\n";

    // With nothing to say about any vertex, the count says it all.
    if (! any_named && ! graph.has_vertex_labels() && ! graph.has_vertex_costs())
        outfile << "  \"vertices\": " << graph.size() << ",\n";
    else if (! graph.has_vertex_labels() && ! graph.has_vertex_costs() && all_named) {
        outfile << "  \"vertices\": [\n";
        for (int v = 0; v < graph.size(); ++v)
            outfile << "    " << json(graph.vertex_name(v)).dump() << (v + 1 < graph.size() ? "," : "") << "\n";
        outfile << "  ],\n";
    }
    else {
        outfile << "  \"vertices\": [\n";
        for (int v = 0; v < graph.size(); ++v) {
            json entry = json::object();
            if (graph.vertex_has_name(v))
                entry["name"] = graph.vertex_name(v);
            if (graph.has_vertex_labels())
                entry["label"] = string{graph.vertex_label(v)};
            if (graph.has_vertex_costs())
                entry["cost"] = graph.vertex_cost(v);
            outfile << "    " << entry.dump() << (v + 1 < graph.size() ? "," : "") << "\n";
        }
        outfile << "  ],\n";
    }

    // Sorted by (from, to, label), so that the output does not depend on the order in
    // which a multigraph's parallel edges were added. An undirected graph holds both
    // directions, and the edge is written once. A cost can only be written in the
    // object form, the array form being frozen.
    vector<tuple<int, int, string, optional<long long>>> edges;
    graph.for_each_edge_and_cost([&](int f, int t, std::string_view l, optional<long long> c) {
        if ((! graph.directed()) && t < f)
            return;
        edges.emplace_back(f, t, string{l}, c);
    });
    sort(edges.begin(), edges.end());

    vector<string> lines;
    for (auto & [f, t, l, c] : edges) {
        json edge;
        if (graph.has_edge_costs()) {
            edge = json::object();
            edge["from"] = endpoint(f);
            edge["to"] = endpoint(t);
            if (graph.has_edge_labels())
                edge["label"] = l;
            edge["cost"] = *c;
        }
        else {
            edge = json::array({endpoint(f), endpoint(t)});
            if (graph.has_edge_labels())
                edge.push_back(l);
        }
        lines.push_back(edge.dump());
    }

    if (lines.empty())
        outfile << "  \"edges\": []\n";
    else {
        outfile << "  \"edges\": [\n";
        for (size_t i = 0; i < lines.size(); ++i)
            outfile << "    " << lines[i] << (i + 1 < lines.size() ? "," : "") << "\n";
        outfile << "  ]\n";
    }

    outfile << "}\n";
}
