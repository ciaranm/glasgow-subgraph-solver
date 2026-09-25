#include <gss/innards/supplemental_graphs.hh>

#include <gss/configuration.hh>
#include <gss/homomorphism.hh>
#include <gss/restarts.hh>

#include <memory>
#include <string>
#include <string_view>
#include <vector>

using std::make_unique;
using std::string_view;
using std::to_string;
using std::vector;

namespace gss::innards
{
    auto build_exact_path_graphs(ProcessedGraphsData & graphs, unsigned size, unsigned & idx,
        unsigned max_graphs, unsigned number_of_exact_path_graphs, bool directed, bool at_most, bool pattern) -> void
    {
        auto & graph_rows = pattern ? graphs.pattern_graph_rows : graphs.target_graph_rows;

        vector<vector<unsigned>> path_counts(size, vector<unsigned>(size, 0));

        // count number of paths from w to v (unless directed, only w >= v, so not v to w)
        for (unsigned v = 0; v < size; ++v) {
            auto nv = graph_rows[v * max_graphs + 0];
            if (at_most)
                nv.set(v);
            for (auto c = nv.find_first(); c != decltype(nv)::npos; c = nv.find_first()) {
                nv.reset(c);
                auto nc = graph_rows[c * max_graphs + 0];
                if (at_most)
                    nc.set(c);
                for (auto w = nc.find_first(); w != decltype(nc)::npos && (directed ? true : w <= v); w = nc.find_first()) {
                    nc.reset(w);
                    ++path_counts[v][w];
                }
            }
        }

        for (unsigned v = 0; v < size; ++v) {
            for (unsigned w = (directed ? 0 : v); w < size; ++w) {
                if (at_most && v == w)
                    graph_rows[v * max_graphs + idx].set(w);
                else {
                    // unless directed, w to v, not v to w, see above
                    unsigned path_count = path_counts[w][v];
                    for (unsigned p = 1; p <= number_of_exact_path_graphs; ++p) {
                        if (path_count >= p) {
                            graph_rows[v * max_graphs + idx + p - 1].set(w);
                            if (! directed)
                                graph_rows[w * max_graphs + idx + p - 1].set(v);
                        }
                    }
                }
            }
        }

        if (pattern)
            for (unsigned p = 1; p <= number_of_exact_path_graphs; ++p)
                graphs.supplemental_graph_names.push_back("exact_path_" + to_string(p));

        idx += number_of_exact_path_graphs;
    }

    auto build_distance3_graphs(ProcessedGraphsData & graphs, unsigned size, unsigned & idx,
        unsigned max_graphs, bool pattern) -> void
    {
        auto & graph_rows = pattern ? graphs.pattern_graph_rows : graphs.target_graph_rows;

        for (unsigned v = 0; v < size; ++v) {
            auto nv = graph_rows[v * max_graphs + 0];
            graph_rows[v * max_graphs + idx] |= nv;
            for (auto c = nv.find_first(); c != decltype(nv)::npos; c = nv.find_first()) {
                nv.reset(c);
                auto nc = graph_rows[c * max_graphs + 0];
                graph_rows[v * max_graphs + idx] |= nc;
                for (auto w = nc.find_first(); w != decltype(nc)::npos; w = nc.find_first()) {
                    nc.reset(w);
                    // v--c--w so v is within distance 3 of w's neighbours
                    graph_rows[v * max_graphs + idx] |= graph_rows[w * max_graphs + 0];
                }
            }
        }

        if (pattern)
            graphs.supplemental_graph_names.push_back("distance3");
        ++idx;
    }

    auto build_k4_graphs(ProcessedGraphsData & graphs, unsigned size, unsigned & idx,
        unsigned max_graphs, bool pattern) -> void
    {
        auto & graph_rows = pattern ? graphs.pattern_graph_rows : graphs.target_graph_rows;

        for (unsigned v = 0; v < size; ++v) {
            auto nv = graph_rows[v * max_graphs + 0];
            for (unsigned w = 0; w < v; ++w) {
                if (nv.test(w)) {
                    // are there two common neighbours with an edge between them?
                    auto common_neighbours = graph_rows[w * max_graphs + 0];
                    common_neighbours &= nv;
                    common_neighbours.reset(v);
                    common_neighbours.reset(w);
                    auto count = common_neighbours.count();
                    if (count >= 2) {
                        bool done = false;
                        auto cn1 = common_neighbours;
                        for (auto x = cn1.find_first(); x != decltype(cn1)::npos && ! done; x = cn1.find_first()) {
                            cn1.reset(x);
                            auto cn2 = common_neighbours;
                            for (auto y = cn2.find_first(); y != decltype(cn2)::npos && ! done; y = cn2.find_first()) {
                                cn2.reset(y);
                                if (v != w && v != x && v != y && w != x && w != y && graph_rows[x * max_graphs + 0].test(y)) {
                                    graph_rows[v * max_graphs + idx].set(w);
                                    graph_rows[w * max_graphs + idx].set(v);
                                    done = true;
                                }
                            }
                        }
                    }
                }
            }
        }

        if (pattern)
            graphs.supplemental_graph_names.push_back("k4");
        ++idx;
    }

    namespace
    {
        // The shape the child search embeds. Its "from" and "to" vertices are pinned by
        // label to the pair being tested, and every other master-graph vertex is labelled
        // "", so every other shape vertex has to carry "" too. A graph file cannot say that,
        // because labels are all or nothing (#87), so here any label other than "from" and
        // "to" means "an interior vertex". Checked up front, because a shape without exactly
        // one of each relates nothing, or everything, and says nothing about why.
        auto shape_for_embedding(const InputGraph & shape) -> InputGraph
        {
            if (shape.has_edge_labels())
                throw UnsupportedConfiguration{"An extra shape graph cannot have edge labels"};

            int froms = 0, tos = 0;
            if (shape.has_vertex_labels())
                for (int v = 0; v < shape.size(); ++v) {
                    froms += shape.vertex_label(v) == "from";
                    tos += shape.vertex_label(v) == "to";
                }
            if (froms != 1 || tos != 1)
                throw UnsupportedConfiguration{"An extra shape graph needs exactly one vertex labelled 'from' and one labelled 'to'"};

            InputGraph result(shape.size(), true, false, shape.directed());
            for (int v = 0; v < shape.size(); ++v) {
                auto label = shape.vertex_label(v);
                result.set_vertex_label(v, (label == "from" || label == "to") ? label : "");
            }
            shape.for_each_edge([&](int f, int t, string_view) {
                if (shape.directed())
                    result.add_directed_edge(f, t, "");
                else
                    result.add_edge(f, t);
            });
            return result;
        }
    }

    auto build_extra_shape(ProcessedGraphsData & graphs, unsigned size, unsigned & idx,
        unsigned max_graphs, const HomomorphismParams & params, InputGraph & shape, bool injective, int count, bool pattern) -> void
    {
        auto & graph_rows = pattern ? graphs.pattern_graph_rows : graphs.target_graph_rows;

        auto child_shape = shape_for_embedding(shape);

        // The shape is embedded into the underlying undirected graph, without loops. A
        // homomorphism of digraphs is also one of their underlying graphs, so this is the
        // same relation on both sides whichever graph is directed. Looking at only the arc
        // v -> w for w < v, as this used to, saw an orientation chosen by vertex numbering,
        // which pattern and target do not share (issue #99, as #97 was for k4).
        InputGraph master_graph(size, true, false);

        for (unsigned v = 0; v < size; ++v) {
            for (unsigned w = 0; w < v; ++w) {
                if (graph_rows[v * max_graphs + 0].test(w) || graph_rows[w * max_graphs + 0].test(v))
                    master_graph.add_edge(v, w);
            }
        }

        auto embeds = [&](unsigned from, unsigned to) -> bool {
            HomomorphismParams child_params;
            child_params.timeout = params.timeout;
            child_params.start_time = params.start_time;
            child_params.induced = false;
            child_params.restarts_schedule = make_unique<NoRestartsSchedule>();
            child_params.clique_detection = false;
            child_params.injectivity = injective ? Injectivity::Injective : Injectivity::NonInjective;
            child_params.no_supplementals = true;

            int seen_count = 0;
            if (count > 1) {
                child_params.count_solutions = true;
                child_params.enumerate_callback = [&](const VertexToVertexMapping &) {
                    return ++seen_count < count;
                };
            }

            master_graph.set_vertex_label(from, "from");
            master_graph.set_vertex_label(to, "to");
            auto child_result = solve_homomorphism_problem(child_shape, master_graph, child_params);
            master_graph.set_vertex_label(from, "");
            master_graph.set_vertex_label(to, "");
            return (count > 1) ? seen_count >= count : ! child_result.mapping.empty();
        };

        // A shape need not look the same from both ends, so a pair is related if the shape
        // embeds in either orientation: the symmetric closure of a preserved relation, so
        // preserved too, and it keeps the supplemental graph symmetric, which is how
        // find_clique() reads it under --cliques-on-supplementals. Testing one orientation
        // and marking both, as this used to, deleted solutions even on undirected graphs.
        for (unsigned v = 0; v < size && ! params.timeout->should_abort(); ++v) {
            for (unsigned w = 0; w < v && ! params.timeout->should_abort(); ++w) {
                if (embeds(v, w) || embeds(w, v)) {
                    graph_rows[v * max_graphs + idx].set(w);
                    graph_rows[w * max_graphs + idx].set(v);
                }
            }
        }

        if (pattern)
            graphs.supplemental_graph_names.push_back("extra_shape");
        ++idx;
    }
}
