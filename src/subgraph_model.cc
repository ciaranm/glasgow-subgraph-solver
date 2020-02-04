/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "subgraph_model.hh"
#include "homomorphism_traits.hh"
#include "configuration.hh"

#include <list>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

using std::list;
using std::map;
using std::max;
using std::pair;
using std::set;
using std::string;
using std::string_view;
using std::to_string;
using std::vector;

namespace
{
    auto calculate_n_shape_graphs(const HomomorphismParams & params) -> unsigned
    {
        return 1 +
            (supports_exact_path_graphs(params) ? params.number_of_exact_path_graphs : 0) +
            (supports_distance3_graphs(params) ? 1 : 0) +
            (supports_k4_graphs(params) ? 1 : 0);
    }
}

struct SubgraphModel::Imp
{
    vector<PatternAdjacencyBitsType> pattern_adjacencies_bits;
    vector<SVOBitset> pattern_graph_rows;
    vector<SVOBitset> target_graph_rows;

    vector<vector<int> > patterns_degrees, targets_degrees;
    int largest_target_degree = 0;
    bool has_less_thans = false;

    vector<int> pattern_vertex_labels, target_vertex_labels, pattern_edge_labels, target_edge_labels;
    vector<int> pattern_loops, target_loops;

    vector<string> pattern_vertex_proof_names, target_vertex_proof_names;
};

SubgraphModel::SubgraphModel(const InputGraph & target, const InputGraph & pattern, const HomomorphismParams & params) :
    _imp(new Imp),
    max_graphs(calculate_n_shape_graphs(params)),
    pattern_size(pattern.size()),
    target_size(target.size())
{
    _imp->patterns_degrees.resize(max_graphs);
    _imp->targets_degrees.resize(max_graphs);

    if (max_graphs > 8 * sizeof(PatternAdjacencyBitsType))
        throw UnsupportedConfiguration{ "Supplemental graphs won't fit in the chosen bitset size" };

    if (pattern.has_edge_labels() && ! params.induced)
        throw UnsupportedConfiguration{ "Currently edge labels only work with --induced" };

    if (params.proof) {
        for (int v = 0 ; v < pattern.size() ; ++v)
            _imp->pattern_vertex_proof_names.push_back(pattern.vertex_name(v));
        for (int v = 0 ; v < target.size() ; ++v)
            _imp->target_vertex_proof_names.push_back(target.vertex_name(v));
    }

    // recode pattern to a bit graph, and strip out loops
    _imp->pattern_graph_rows.resize(pattern_size * max_graphs, SVOBitset(pattern_size, 0));
    _imp->pattern_loops.resize(pattern_size);
    for (unsigned i = 0 ; i < pattern_size ; ++i) {
        for (unsigned j = 0 ; j < pattern_size ; ++j) {
            if (pattern.adjacent(i, j)) {
                if (i == j)
                    _imp->pattern_loops[i] = 1;
                else
                    _imp->pattern_graph_rows[i * max_graphs + 0].set(j);
            }
        }
    }

    // re-encode and store pattern labels
    map<string, int> vertex_labels_map;
    int next_vertex_label = 0;
    if (pattern.has_vertex_labels()) {
        for (unsigned i = 0 ; i < pattern_size ; ++i) {
            if (vertex_labels_map.emplace(pattern.vertex_label(i), next_vertex_label).second)
                ++next_vertex_label;
        }

        _imp->pattern_vertex_labels.resize(pattern_size);
        for (unsigned i = 0 ; i < pattern_size ; ++i)
            _imp->pattern_vertex_labels[i] = vertex_labels_map.find(string{ pattern.vertex_label(i) })->second;
    }

    // re-encode and store edge labels
    map<string, int> edge_labels_map;
    int next_edge_label = 0;
    if (pattern.has_edge_labels()) {
        _imp->pattern_edge_labels.resize(pattern_size * pattern_size);
        for (unsigned i = 0 ; i < pattern_size ; ++i)
            for (unsigned j = 0 ; j < pattern_size ; ++j)
                if (pattern.adjacent(i, j)) {
                    auto r = edge_labels_map.emplace(pattern.edge_label(i, j), next_edge_label);
                    if (r.second)
                        ++next_edge_label;
                    _imp->pattern_edge_labels[i * pattern_size + j] = r.first->second;
                }
    }

    // recode target to a bit graph, and take out loops
    _imp->target_graph_rows.resize(target_size * max_graphs, SVOBitset{ target_size, 0 });
    _imp->target_loops.resize(target_size);
    for (auto e = target.begin_edges(), e_end = target.end_edges() ; e != e_end ; ++e) {
        if (e->first.first == e->first.second)
            _imp->target_loops[e->first.first] = 1;
        else
            _imp->target_graph_rows[e->first.first * max_graphs + 0].set(e->first.second);
    }

    // target vertex labels
    if (pattern.has_vertex_labels()) {
        for (unsigned i = 0 ; i < target_size ; ++i) {
            if (vertex_labels_map.emplace(target.vertex_label(i), next_vertex_label).second)
                ++next_vertex_label;
        }

        _imp->target_vertex_labels.resize(target_size);
        for (unsigned i = 0 ; i < target_size ; ++i)
            _imp->target_vertex_labels[i] = vertex_labels_map.find(string{ target.vertex_label(i) })->second;
    }

    // target edge labels
    if (pattern.has_edge_labels()) {
        _imp->target_edge_labels.resize(target_size * target_size);
        for (auto e = target.begin_edges(), e_end = target.end_edges() ; e != e_end ; ++e) {
            auto r = edge_labels_map.emplace(e->second, next_edge_label);
            if (r.second)
                ++next_edge_label;

            _imp->target_edge_labels[e->first.first * target_size + e->first.second] = r.first->second;
        }
    }

    auto decode = [&] (string_view s) -> int {
        auto n = pattern.vertex_from_name(s);
        if (! n)
            throw UnsupportedConfiguration{ "No vertex named '" + string{ s } + "'" };
        return *n;
    };

    // pattern less than constraints
    if (! params.pattern_less_constraints.empty()) {
        _imp->has_less_thans = true;
        list<pair<unsigned, unsigned> > pattern_less_thans_in_wrong_order;
        for (auto & [ a, b ] : params.pattern_less_constraints) {
            auto a_decoded = decode(a), b_decoded = decode(b);
            pattern_less_thans_in_wrong_order.emplace_back(a_decoded, b_decoded);
        }

        // put them in a convenient order, so we don't need a propagation loop
        while (! pattern_less_thans_in_wrong_order.empty()) {
            bool loop_detect = true;
            set<unsigned> cannot_order_yet;
            for (auto & [ _, b ] : pattern_less_thans_in_wrong_order)
                cannot_order_yet.emplace(b);
            for (auto p = pattern_less_thans_in_wrong_order.begin() ; p != pattern_less_thans_in_wrong_order.end() ; ) {
                if (cannot_order_yet.count(p->first))
                    ++p;
                else {
                    loop_detect = false;
                    pattern_less_thans_in_convenient_order.push_back(*p);
                    pattern_less_thans_in_wrong_order.erase(p++);
                }
            }

            if (loop_detect)
                throw UnsupportedConfiguration{ "Pattern less than constraints form a loop" };
        }
    }
}

SubgraphModel::~SubgraphModel() = default;

auto SubgraphModel::pattern_vertex_for_proof(int v) const -> NamedVertex
{
    if (v < 0 || unsigned(v) >= _imp->pattern_vertex_proof_names.size())
        throw ProofError{ "Oops, there's a bug: v out of range in pattern" };
    return pair{ v, _imp->pattern_vertex_proof_names[v] };
}

auto SubgraphModel::target_vertex_for_proof(int v) const -> NamedVertex
{
    if (v < 0 || unsigned(v) >= _imp->target_vertex_proof_names.size())
        throw ProofError{ "Oops, there's a bug: v out of range in target" };
    return pair{ v, _imp->target_vertex_proof_names[v] };
}

auto SubgraphModel::prepare(const HomomorphismParams & params) -> bool
{
    // pattern and target degrees, for the main graph
    _imp->patterns_degrees.at(0).resize(pattern_size);
    _imp->targets_degrees.at(0).resize(target_size);

    for (unsigned i = 0 ; i < pattern_size ; ++i)
        _imp->patterns_degrees.at(0).at(i) = _imp->pattern_graph_rows[i * max_graphs + 0].count();

    for (unsigned i = 0 ; i < target_size ; ++i)
        _imp->targets_degrees.at(0).at(i) = _imp->target_graph_rows[i * max_graphs + 0].count();

    if (global_degree_is_preserved(params)) {
        vector<pair<int, int> > p_gds, t_gds;
        for (unsigned i = 0 ; i < pattern_size ; ++i)
            p_gds.emplace_back(i, _imp->patterns_degrees.at(0).at(i));
        for (unsigned i = 0 ; i < target_size ; ++i)
            t_gds.emplace_back(i, _imp->targets_degrees.at(0).at(i));

        sort(p_gds.begin(), p_gds.end(), [] (const pair<int, int> & a, const pair<int, int> & b) {
                return a.second > b.second; });
        sort(t_gds.begin(), t_gds.end(), [] (const pair<int, int> & a, const pair<int, int> & b) {
                return a.second > b.second; });

        for (unsigned i = 0 ; i < p_gds.size() ; ++i)
            if (p_gds.at(i).second > t_gds.at(i).second) {
                if (params.proof) {
                    for (unsigned p = 0 ; p <= i ; ++p) {
                        vector<int> n_p;
                        auto np = _imp->pattern_graph_rows[p_gds.at(p).first * max_graphs + 0];
                        for (unsigned j = 0 ; j < pattern_size ; ++j)
                            if (np.test(j))
                                n_p.push_back(j);

                        for (unsigned t = i ; t < t_gds.size() ; ++t) {
                            vector<int> n_t;
                            auto nt = _imp->target_graph_rows[t_gds.at(t).first * max_graphs + 0];
                            for (auto j = nt.find_first() ; j != decltype(nt)::npos ; j = nt.find_first()) {
                                nt.reset(j);
                                n_t.push_back(j);
                            }

                            params.proof->incompatible_by_degrees(0,
                                    pattern_vertex_for_proof(p_gds.at(p).first), n_p,
                                    target_vertex_for_proof(t_gds.at(t).first), n_t);
                        }
                    }

                    vector<int> patterns, targets;
                    for (unsigned p = 0 ; p <= i ; ++p)
                        patterns.push_back(p_gds.at(p).first);
                    for (unsigned t = 0 ; t < i ; ++t)
                        targets.push_back(t_gds.at(t).first);

                    params.proof->emit_hall_set_or_violator(patterns, targets);
                }
                return false;
            }
    }

    unsigned next_pattern_supplemental = 1, next_target_supplemental = 1;
    // build exact path graphs
    if (supports_exact_path_graphs(params)) {
        _build_exact_path_graphs(_imp->pattern_graph_rows, pattern_size, next_pattern_supplemental, params.number_of_exact_path_graphs);
        _build_exact_path_graphs(_imp->target_graph_rows, target_size, next_target_supplemental, params.number_of_exact_path_graphs);

        if (params.proof) {
            for (int g = 1 ; g <= params.number_of_exact_path_graphs ; ++g) {
                for (unsigned p = 0 ; p < pattern_size ; ++p) {
                    for (unsigned q = 0 ; q < pattern_size ; ++q) {
                        if (p == q || ! _imp->pattern_graph_rows[p * max_graphs + g].test(q))
                            continue;

                        auto named_p = pattern_vertex_for_proof(p);
                        auto named_q = pattern_vertex_for_proof(q);

                        auto n_p_q = _imp->pattern_graph_rows[p * max_graphs + 0];
                        n_p_q &= _imp->pattern_graph_rows[q * max_graphs + 0];
                        vector<NamedVertex> between_p_and_q;
                        for (auto v = n_p_q.find_first() ; v != decltype(n_p_q)::npos ; v = n_p_q.find_first()) {
                            n_p_q.reset(v);
                            between_p_and_q.push_back(pattern_vertex_for_proof(v));
                            if (between_p_and_q.size() >= unsigned(g))
                                break;
                        }

                        for (unsigned t = 0 ; t < target_size ; ++t) {
                            auto named_t = target_vertex_for_proof(t);

                            vector<NamedVertex> named_n_t, named_d_n_t;
                            vector<pair<NamedVertex, vector<NamedVertex> > > named_two_away_from_t;
                            auto n_t = _imp->target_graph_rows[t * max_graphs + 0];
                            for (auto w = n_t.find_first() ; w != decltype(n_t)::npos ; w = n_t.find_first()) {
                                n_t.reset(w);
                                named_n_t.push_back(target_vertex_for_proof(w));
                            }

                            auto nd_t = _imp->target_graph_rows[t * max_graphs + g];
                            for (auto w = nd_t.find_first() ; w != decltype(nd_t)::npos ; w = nd_t.find_first()) {
                                nd_t.reset(w);
                                named_d_n_t.push_back(target_vertex_for_proof(w));
                            }

                            auto n2_t = _imp->target_graph_rows[t * max_graphs + 1];
                            for (auto w = n2_t.find_first() ; w != decltype(n2_t)::npos ; w = n2_t.find_first()) {
                                n2_t.reset(w);
                                auto n_t_w = _imp->target_graph_rows[w * max_graphs + 0];
                                n_t_w &= _imp->target_graph_rows[t * max_graphs + 0];
                                vector<NamedVertex> named_n_t_w;
                                for (auto x = n_t_w.find_first() ; x != decltype(n_t_w)::npos ; x = n_t_w.find_first()) {
                                    n_t_w.reset(x);
                                    named_n_t_w.push_back(target_vertex_for_proof(x));
                                }
                                named_two_away_from_t.emplace_back(target_vertex_for_proof(w), named_n_t_w);
                            }

                            params.proof->create_exact_path_graphs(g, named_p, named_q, between_p_and_q,
                                    named_t, named_n_t, named_two_away_from_t, named_d_n_t);
                        }
                    }
                }
            }
        }
    }

    if (supports_distance3_graphs(params)) {
        _build_distance3_graphs(_imp->pattern_graph_rows, pattern_size, next_pattern_supplemental);
        _build_distance3_graphs(_imp->target_graph_rows, target_size, next_target_supplemental);
    }

    if (supports_k4_graphs(params)) {
        _build_k4_graphs(_imp->pattern_graph_rows, pattern_size, next_pattern_supplemental);
        _build_k4_graphs(_imp->target_graph_rows, target_size, next_target_supplemental);
    }

    if (next_pattern_supplemental != max_graphs || next_target_supplemental != max_graphs)
        throw UnsupportedConfiguration{ "something has gone wrong with supplemental graph indexing: " + to_string(next_pattern_supplemental)
            + " " + to_string(next_target_supplemental) + " " + to_string(max_graphs) };

    // pattern and target degrees, for supplemental graphs
    for (unsigned g = 1 ; g < max_graphs ; ++g) {
        _imp->patterns_degrees.at(g).resize(pattern_size);
        _imp->targets_degrees.at(g).resize(target_size);
    }

    for (unsigned g = 1 ; g < max_graphs ; ++g) {
        for (unsigned i = 0 ; i < pattern_size ; ++i)
            _imp->patterns_degrees.at(g).at(i) = _imp->pattern_graph_rows[i * max_graphs + g].count();

        for (unsigned i = 0 ; i < target_size ; ++i)
            _imp->targets_degrees.at(g).at(i) = _imp->target_graph_rows[i * max_graphs + g].count();
    }

    for (unsigned i = 0 ; i < target_size ; ++i)
        _imp->largest_target_degree = max(_imp->largest_target_degree, _imp->targets_degrees[0][i]);

    // pattern adjacencies, compressed
    _imp->pattern_adjacencies_bits.resize(pattern_size * pattern_size);
    for (unsigned g = 0 ; g < max_graphs ; ++g)
        for (unsigned i = 0 ; i < pattern_size ; ++i)
            for (unsigned j = 0 ; j < pattern_size ; ++j)
                if (_imp->pattern_graph_rows[i * max_graphs + g].test(j))
                    _imp->pattern_adjacencies_bits[i * pattern_size + j] |= (1u << g);

    return true;
}

auto SubgraphModel::_build_exact_path_graphs(vector<SVOBitset> & graph_rows, unsigned size, unsigned & idx,
        unsigned number_of_exact_path_graphs) -> void
{
    vector<vector<unsigned> > path_counts(size, vector<unsigned>(size, 0));

    // count number of paths from w to v (only w >= v, so not v to w)
    for (unsigned v = 0 ; v < size ; ++v) {
        auto nv = graph_rows[v * max_graphs + 0];
        for (auto c = nv.find_first() ; c != decltype(nv)::npos ; c = nv.find_first()) {
            nv.reset(c);
            auto nc = graph_rows[c * max_graphs + 0];
            for (auto w = nc.find_first() ; w != decltype(nc)::npos && w <= v ; w = nc.find_first()) {
                nc.reset(w);
                ++path_counts[v][w];
            }
        }
    }

    for (unsigned v = 0 ; v < size ; ++v) {
        for (unsigned w = v ; w < size ; ++w) {
            // w to v, not v to w, see above
            unsigned path_count = path_counts[w][v];
            for (unsigned p = 1 ; p <= number_of_exact_path_graphs ; ++p) {
                if (path_count >= p) {
                    graph_rows[v * max_graphs + idx + p - 1].set(w);
                    graph_rows[w * max_graphs + idx + p - 1].set(v);
                }
            }
        }
    }

    idx += number_of_exact_path_graphs;
}

auto SubgraphModel::_build_distance3_graphs(vector<SVOBitset> & graph_rows, unsigned size, unsigned & idx) -> void
{
    for (unsigned v = 0 ; v < size ; ++v) {
        auto nv = graph_rows[v * max_graphs + 0];
        for (auto c = nv.find_first() ; c != decltype(nv)::npos ; c = nv.find_first()) {
            nv.reset(c);
            auto nc = graph_rows[c * max_graphs + 0];
            for (auto w = nc.find_first() ; w != decltype(nc)::npos ; w = nc.find_first()) {
                nc.reset(w);
                // v--c--w so v is within distance 3 of w's neighbours
                graph_rows[v * max_graphs + idx] |= graph_rows[w * max_graphs + 0];
            }
        }
    }

    ++idx;
}

auto SubgraphModel::_build_k4_graphs(vector<SVOBitset> & graph_rows, unsigned size, unsigned & idx) -> void
{
    for (unsigned v = 0 ; v < size ; ++v) {
        auto nv = graph_rows[v * max_graphs + 0];
        for (unsigned w = 0 ; w < v ; ++w) {
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
                    for (auto x = cn1.find_first() ; x != decltype(cn1)::npos && ! done ; x = cn1.find_first()) {
                        cn1.reset(x);
                        auto cn2 = common_neighbours;
                        for (auto y = cn2.find_first() ; y != decltype(cn2)::npos && ! done ; y = cn2.find_first()) {
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

    ++idx;
}

auto SubgraphModel::pattern_adjacency_bits(int p, int q) const -> PatternAdjacencyBitsType
{
    return _imp->pattern_adjacencies_bits[pattern_size * p + q];
}

auto SubgraphModel::pattern_graph_row(int g, int p) const -> const SVOBitset &
{
    return _imp->pattern_graph_rows[p * max_graphs + g];
}

auto SubgraphModel::target_graph_row(int g, int t) const -> const SVOBitset &
{
    return _imp->target_graph_rows[t * max_graphs + g];
}

auto SubgraphModel::pattern_degree(int g, int p) const -> unsigned
{
    return _imp->patterns_degrees[g][p];
}

auto SubgraphModel::target_degree(int g, int t) const -> unsigned
{
    return _imp->targets_degrees[g][t];
}

auto SubgraphModel::largest_target_degree() const -> unsigned
{
    return _imp->largest_target_degree;
}

auto SubgraphModel::has_vertex_labels() const -> bool
{
    return ! _imp->pattern_vertex_labels.empty();
}

auto SubgraphModel::has_edge_labels() const -> bool
{
    return ! _imp->pattern_edge_labels.empty();
}

auto SubgraphModel::pattern_vertex_label(int p) const -> int
{
    return _imp->pattern_vertex_labels[p];
}

auto SubgraphModel::target_vertex_label(int t) const -> int
{
    return _imp->target_vertex_labels[t];
}

auto SubgraphModel::pattern_edge_label(int p, int q) const -> int
{
    return _imp->pattern_edge_labels[p * pattern_size + q];
}

auto SubgraphModel::target_edge_label(int t, int u) const -> int
{
    return _imp->target_edge_labels[t * target_size + u];
}

auto SubgraphModel::pattern_has_loop(int p) const -> bool
{
    return _imp->pattern_loops[p];
}

auto SubgraphModel::target_has_loop(int t) const -> bool
{
    return _imp->target_loops[t];
}

auto SubgraphModel::has_less_thans() const -> bool
{
    return _imp->has_less_thans;
}

