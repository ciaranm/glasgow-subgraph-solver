/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "homomorphism.hh"
#include "clique.hh"
#include "configuration.hh"
#include "fixed_bit_set.hh"
#include "graph_traits.hh"
#include "homomorphism_traits.hh"
#include "thread_utils.hh"
#include "watches.hh"
#include "proof.hh"

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <functional>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <queue> 
#include <random>
#include <thread>
#include <tuple>
#include <type_traits>
#include <utility>

#include <boost/thread/barrier.hpp>
#include <boost/dynamic_bitset.hpp>

using std::atomic;
using std::conditional_t;
using std::fill;
using std::find_if;
using std::function;
using std::greater;
using std::is_same;
using std::list;
using std::make_optional;
using std::make_unique;
using std::max;
using std::map;
using std::move;
using std::mt19937;
using std::mutex;
using std::numeric_limits;
using std::optional;
using std::pair;
using std::queue;
using std::set;
using std::sort;
using std::stable_sort;
using std::string;
using std::string_view;
using std::swap;
using std::thread;
using std::to_string;
using std::tuple;
using std::uniform_int_distribution;
using std::unique_lock;
using std::unique_ptr;
using std::vector;

using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::steady_clock;
using std::chrono::operator""ms;

using boost::barrier;
using boost::dynamic_bitset;

#include <iostream>
using std::cout;

namespace
{
    auto calculate_n_shape_graphs(const HomomorphismParams & params) -> unsigned
    {
        return 1 +
            (supports_exact_path_graphs(params) ? params.number_of_exact_path_graphs : 0) +
            (supports_common_neighbour_shapes(params) ? params.number_of_common_neighbour_graphs - params.skip_common_neighbour_graphs : 0) +
            (supports_distance3_graphs(params) ? 1 : 0) +
            (supports_k4_graphs(params) ? 1 : 0) +
            (supports_diamond_graphs(params) ? 1 : 0) +
            (params.induced ? 1 : 0);
    }

    template <typename BitSetType_, typename ArrayType_, typename PatternAdjacencyBitsType_>
    struct SubgraphModel
    {
        const unsigned max_graphs;
        unsigned pattern_size, full_pattern_size, target_size;

        vector<PatternAdjacencyBitsType_> pattern_adjacencies_bits;
        vector<dynamic_bitset<> > pattern_graph_rows;
        vector<BitSetType_> target_graph_rows;
        vector<pair<unsigned, unsigned> > pattern_less_thans_in_convenient_order;

        vector<int> pattern_permutation, isolated_vertices;
        vector<vector<int> > patterns_degrees, targets_degrees;
        int largest_target_degree;
        bool has_less_thans;

        vector<int> pattern_vertex_labels, target_vertex_labels, pattern_edge_labels, target_edge_labels;

        vector<pair<int, int> > pattern_dir_degrees, target_dir_degrees;
        vector<pair<bool, bool> > pattern_big_constraints;

        vector<string> pattern_vertex_proof_names, target_vertex_proof_names;

        vector<dynamic_bitset<> > target_graph_reachability, pattern_graph_reachability;
        std::vector<std::pair<bool, std::vector<int> > > target_hyperedges, pattern_hyperedges;

        

        SubgraphModel(const InputGraph & target, const InputGraph & pattern, const HomomorphismParams & params,
                const set<int> & exclude_pattern_vertices) :
            max_graphs(calculate_n_shape_graphs(params)),
            pattern_size(pattern.size()),
            full_pattern_size(pattern.size()),
            target_size(target.size()),
            patterns_degrees(max_graphs),
            targets_degrees(max_graphs),
            largest_target_degree(0),
            has_less_thans(false)
        {
            if (max_graphs > 8 * sizeof(PatternAdjacencyBitsType_))
                throw UnsupportedConfiguration{ "Supplemental graphs won't fit in the chosen bitset size" };

            if (pattern.has_edge_labels() && ! params.induced)
                throw UnsupportedConfiguration{ "Currently edge labels only work with --induced" };

            if (params.proof) {
                for (int v = 0 ; v < pattern.size() ; ++v)
                    pattern_vertex_proof_names.push_back(pattern.vertex_name(v));
                for (int v = 0 ; v < target.size() ; ++v)
                    target_vertex_proof_names.push_back(target.vertex_name(v));
            }

            // strip out isolated vertices in the pattern, and build pattern_permutation
            for (unsigned v = 0 ; v < full_pattern_size ; ++v) {
                if (exclude_pattern_vertices.count(v)) {
                    --pattern_size;
                }
                else if (can_strip_isolated_vertices(params) && 0 == pattern.degree(v)) {
                    isolated_vertices.push_back(v);
                    --pattern_size;
                }
                else
                    pattern_permutation.push_back(v);
            }

            // recode pattern to a bit graph
            pattern_graph_rows.resize(pattern_size * max_graphs, dynamic_bitset<>(pattern_size));
            for (unsigned i = 0 ; i < pattern_size ; ++i)
                for (unsigned j = 0 ; j < pattern_size ; ++j)
                    if (pattern.adjacent(pattern_permutation.at(i), pattern_permutation.at(j)))
                        pattern_graph_rows[i * max_graphs + 0].set(j);

            // re-encode and store pattern labels
            map<string, int> vertex_labels_map;
            int next_vertex_label = 0;
            if (pattern.has_vertex_labels()) {
                for (unsigned i = 0 ; i < pattern_size ; ++i) {
                    if (vertex_labels_map.emplace(pattern.vertex_label(pattern_permutation.at(i)), next_vertex_label).second)
                        ++next_vertex_label;
                }

                pattern_vertex_labels.resize(pattern_size);
                for (unsigned i = 0 ; i < pattern_size ; ++i)
                    pattern_vertex_labels[pattern_permutation[i]] = vertex_labels_map.find(string{ pattern.vertex_label(pattern_permutation[i]) })->second;
            }

            // re-encode and store edge labels
            map<string, int> edge_labels_map;
            int next_edge_label = 0;
            if (pattern.has_edge_labels()) {
                pattern_edge_labels.resize(pattern_size * pattern_size);
                for (unsigned i = 0 ; i < pattern_size ; ++i)
                    for (unsigned j = 0 ; j < pattern_size ; ++j)
                        if (pattern.adjacent(pattern_permutation.at(i), pattern_permutation.at(j))) {
                            auto r = edge_labels_map.emplace(pattern.edge_label(pattern_permutation.at(i), pattern_permutation.at(j)), next_edge_label);
                            if (r.second)
                                ++next_edge_label;
                            pattern_edge_labels[pattern_permutation.at(i) * pattern_size + pattern_permutation.at(j)] = r.first->second;
                        }
            }

            // recode target to a bit graph
            target_graph_rows.resize(target_size * max_graphs, BitSetType_{ target_size, 0 });
            for (auto e = target.begin_edges(), e_end = target.end_edges() ; e != e_end ; ++e)
                target_graph_rows[e->first.first * max_graphs + 0].set(e->first.second);

            // target vertex labels
            if (pattern.has_vertex_labels()) {
                for (unsigned i = 0 ; i < target_size ; ++i) {
                    if (vertex_labels_map.emplace(target.vertex_label(i), next_vertex_label).second)
                        ++next_vertex_label;
                }

                target_vertex_labels.resize(target_size);
                for (unsigned i = 0 ; i < target_size ; ++i)
                    target_vertex_labels[i] = vertex_labels_map.find(string{ target.vertex_label(i) })->second;
            }

            // target edge labels
            if (pattern.has_edge_labels()) {
                target_edge_labels.resize(target_size * target_size);
                for (auto e = target.begin_edges(), e_end = target.end_edges() ; e != e_end ; ++e) {
                    auto r = edge_labels_map.emplace(e->second, next_edge_label);
                    if (r.second)
                        ++next_edge_label;

                    target_edge_labels[e->first.first * target_size + e->first.second] = r.first->second;
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
                has_less_thans = true;
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

            if (params.bigraph) {
                pattern_dir_degrees.resize(pattern_size);
                pattern_big_constraints.resize(pattern_size);
                target_dir_degrees.resize(target_size);

                pattern_graph_reachability.resize(pattern_size, dynamic_bitset<>(pattern_size));
                target_graph_reachability.resize(target_size, dynamic_bitset<>(target_size));
                
                
                target_hyperedges.resize(target.number_of_hyperedges());
                pattern_hyperedges.resize(pattern.number_of_hyperedges());

                for(unsigned int a=0; a<target.number_of_hyperedges();a++) {
                    target_hyperedges[a].second.resize(target_size);
                    target_hyperedges[a] = target.get_hyperedge(a);
                }
                for(unsigned int a=0; a<pattern.number_of_hyperedges();a++) {
                    pattern_hyperedges[a].second.resize(pattern_size);
                    pattern_hyperedges[a] = pattern.get_hyperedge(a);
                }
                


                std::set<int> pattern_unique;
                std::queue<int> pattern_reach;

                std::set<int> target_unique;
                std::queue<int> target_reach;

                // Set in-out degrees of each vertex for bigraph root and site constraints
                for(unsigned int a=0; a<pattern_size;a++) {
                    pattern_dir_degrees[a].first = pattern.in_degree(a);
                    pattern_dir_degrees[a].second = pattern.out_degree(a);
                    pattern_big_constraints[a] = pattern.get_big_constraint(a);

                    pattern_graph_reachability[a][a] = 1;                
                    if(pattern.out_degree(a) == 0) pattern_reach.push(a);
                }

                for(unsigned int a=0; a<target_size;a++) {
                    target_dir_degrees[a].first = target.in_degree(a);
                    target_dir_degrees[a].second = target.out_degree(a);

                    target_graph_reachability[a][a] = 1;
                    if(target.out_degree(a) == 0) target_reach.push(a);            
                }

                
                // Build reachability matrices by bubbling up from sink nodes and performing BFS
                while(!pattern_reach.empty()){
                    int v = pattern_reach.front();
                    pattern_reach.pop();
                    pattern_unique.erase(v);
                    for(unsigned int a=0;a<pattern_size;a++) 
                        if(pattern.adjacent(v,a) && pattern.edge_label(v,a) == "unlabelled") {
                            for(unsigned int b=0;b<pattern_size;b++) pattern_graph_reachability[a][b] |= pattern_graph_reachability[v][b];
                            if(pattern_unique.insert(a).second) pattern_reach.push(a);
                        }
                }
                while(!target_reach.empty()){
                    int v = target_reach.front();
                    target_reach.pop();
                    target_unique.erase(v);
                    for(unsigned int a=0;a<target_size;a++) 
                        if(target.adjacent(v,a) && target.edge_label(v,a) == "unlabelled") {
                            for(unsigned int b=0;b<target_size;b++) target_graph_reachability[a][b] |= target_graph_reachability[v][b];
                            if(target_unique.insert(a).second) target_reach.push(a);
                        }
                }
            }
        }

        auto pattern_vertex_for_proof(int v) const -> NamedVertex
        {
            return pair{ v, pattern_vertex_proof_names[pattern_permutation[v]] };
        }

        auto target_vertex_for_proof(int v) const -> NamedVertex
        {
            return pair{ v, target_vertex_proof_names[v] };
        }


        auto prepare(const HomomorphismParams & params) -> bool
        {
            // pattern and target degrees, for the main graph
            patterns_degrees.at(0).resize(pattern_size);
            targets_degrees.at(0).resize(target_size);

            for (unsigned i = 0 ; i < pattern_size ; ++i)
                patterns_degrees.at(0).at(i) = pattern_graph_rows[i * max_graphs + 0].count();

            for (unsigned i = 0 ; i < target_size ; ++i)
                targets_degrees.at(0).at(i) = target_graph_rows[i * max_graphs + 0].count();

            if (global_degree_is_preserved(params) && ! params.proof) {
                auto pattern_degrees_sorted = patterns_degrees.at(0), target_degrees_sorted = targets_degrees.at(0);
                sort(pattern_degrees_sorted.begin(), pattern_degrees_sorted.end(), greater<int>());
                sort(target_degrees_sorted.begin(), target_degrees_sorted.end(), greater<int>());
                for (unsigned i = 0 ; i < pattern_degrees_sorted.size() ; ++i)
                    if (pattern_degrees_sorted.at(i) > target_degrees_sorted.at(i))
                        return false;
            }

            for (unsigned i = 0 ; i < target_size ; ++i)
                largest_target_degree = max(largest_target_degree, targets_degrees[0][i]);

            unsigned next_pattern_supplemental = 1, next_target_supplemental = 1;
            // build exact path graphs
            if (supports_exact_path_graphs(params)) {
                switch (params.number_of_exact_path_graphs) {
                    case 1:
                        build_exact_path_graphs<1>(pattern_graph_rows, pattern_size, next_pattern_supplemental);
                        build_exact_path_graphs<1>(target_graph_rows, target_size, next_target_supplemental);
                        break;
                    case 2:
                        build_exact_path_graphs<2>(pattern_graph_rows, pattern_size, next_pattern_supplemental);
                        build_exact_path_graphs<2>(target_graph_rows, target_size, next_target_supplemental);
                        break;
                    case 3:
                        build_exact_path_graphs<3>(pattern_graph_rows, pattern_size, next_pattern_supplemental);
                        build_exact_path_graphs<3>(target_graph_rows, target_size, next_target_supplemental);
                        break;
                    case 4:
                        build_exact_path_graphs<4>(pattern_graph_rows, pattern_size, next_pattern_supplemental);
                        build_exact_path_graphs<4>(target_graph_rows, target_size, next_target_supplemental);
                        break;
                    case 5:
                        build_exact_path_graphs<5>(pattern_graph_rows, pattern_size, next_pattern_supplemental);
                        build_exact_path_graphs<5>(target_graph_rows, target_size, next_target_supplemental);
                        break;
                    default:
                        throw UnsupportedConfiguration{ "Unsupported number of exact path graphs" };
                }

                if (supports_common_neighbour_shapes(params)) {
                    build_common_neighbour_graphs(params.number_of_common_neighbour_graphs,
                            params.skip_common_neighbour_graphs, pattern_graph_rows, pattern_size, next_pattern_supplemental);
                    build_common_neighbour_graphs(params.number_of_common_neighbour_graphs,
                            params.skip_common_neighbour_graphs, target_graph_rows, target_size, next_target_supplemental);
                }
            }

            if (supports_distance3_graphs(params)) {
                build_distance3_graphs(pattern_graph_rows, pattern_size, next_pattern_supplemental);
                build_distance3_graphs(target_graph_rows, target_size, next_target_supplemental);
            }

            if (supports_k4_graphs(params)) {
                build_k4_or_diamond_graphs(true, pattern_graph_rows, pattern_size, next_pattern_supplemental);
                build_k4_or_diamond_graphs(true, target_graph_rows, target_size, next_target_supplemental);
            }

            if (supports_diamond_graphs(params)) {
                build_k4_or_diamond_graphs(false, pattern_graph_rows, pattern_size, next_pattern_supplemental);
                build_k4_or_diamond_graphs(false, target_graph_rows, target_size, next_target_supplemental);
            }

            // build complement graphs. these must come last!
            if (params.induced) {
                build_complement_graphs(pattern_graph_rows, pattern_size, next_pattern_supplemental);
                build_complement_graphs(target_graph_rows, target_size, next_target_supplemental);
            }

            if (next_pattern_supplemental != max_graphs || next_target_supplemental != max_graphs)
                throw UnsupportedConfiguration{ "something has gone wrong with supplemental graph indexing: " + to_string(next_pattern_supplemental)
                    + " " + to_string(next_target_supplemental) + " " + to_string(max_graphs) };

            // pattern and target degrees, for supplemental graphs
            for (unsigned g = 1 ; g < max_graphs ; ++g) {
                patterns_degrees.at(g).resize(pattern_size);
                targets_degrees.at(g).resize(target_size);
            }

            for (unsigned g = 1 ; g < max_graphs ; ++g) {
                for (unsigned i = 0 ; i < pattern_size ; ++i)
                    patterns_degrees.at(g).at(i) = pattern_graph_rows[i * max_graphs + g].count();

                for (unsigned i = 0 ; i < target_size ; ++i)
                    targets_degrees.at(g).at(i) = target_graph_rows[i * max_graphs + g].count();
            }

            for (unsigned i = 0 ; i < target_size ; ++i)
                largest_target_degree = max(largest_target_degree, targets_degrees[0][i]);

            // pattern adjacencies, compressed
            pattern_adjacencies_bits.resize(pattern_size * pattern_size);
            for (unsigned g = 0 ; g < max_graphs ; ++g)
                for (unsigned i = 0 ; i < pattern_size ; ++i)
                    for (unsigned j = 0 ; j < pattern_size ; ++j)
                        if (pattern_graph_rows[i * max_graphs + g].test(j))
                            pattern_adjacencies_bits[i * pattern_size + j] |= (1u << g);

            return true;
        }

        template <unsigned number_of_exact_path_graphs_, typename PossiblySomeOtherBitSetType_>
        auto build_exact_path_graphs(vector<PossiblySomeOtherBitSetType_> & graph_rows, unsigned size, unsigned & idx) -> void
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
                    for (unsigned p = 1 ; p <= number_of_exact_path_graphs_ ; ++p) {
                        if (path_count >= p) {
                            graph_rows[v * max_graphs + idx + p - 1].set(w);
                            graph_rows[w * max_graphs + idx + p - 1].set(v);
                        }
                    }
                }
            }

            idx += number_of_exact_path_graphs_;
        }

        template <typename PossiblySomeOtherBitSetType_>
        auto build_distance3_graphs(vector<PossiblySomeOtherBitSetType_> & graph_rows, unsigned size, unsigned & idx) -> void
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

        template <typename PossiblySomeOtherBitSetType_>
        auto build_k4_or_diamond_graphs(bool k4, vector<PossiblySomeOtherBitSetType_> & graph_rows, unsigned size, unsigned & idx) -> void
        {
            for (unsigned v = 0 ; v < size ; ++v) {
                auto nv = graph_rows[v * max_graphs + 0];
                for (unsigned w = 0 ; w < v ; ++w) {
                    // in k4, v -- w, but in a diamond we don't mind
                    if ((! k4) || (nv.test(w))) {
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

        template <typename PossiblySomeOtherBitSetType_>
        auto build_common_neighbour_graphs(
                unsigned number_of_common_neighbour_graphs,
                unsigned skip_common_neighbour_graphs,
                vector<PossiblySomeOtherBitSetType_> & graph_rows, unsigned size, unsigned & idx) -> void
        {
            for (unsigned v = 0 ; v < size ; ++v) {
                auto nv = graph_rows[v * max_graphs + 0];
                for (unsigned w = 0 ; w < v ; ++w) {
                    if (nv.test(w)) {
                        auto common_neighbours = graph_rows[w * max_graphs + 0];
                        common_neighbours &= nv;
                        common_neighbours.reset(v);
                        common_neighbours.reset(w);
                        auto count = common_neighbours.count();

                        for (unsigned p = 1 + skip_common_neighbour_graphs ; p <= number_of_common_neighbour_graphs ; ++p) {
                            if (count >= p) {
                                graph_rows[v * max_graphs + idx + p - 1 - skip_common_neighbour_graphs].set(w);
                                graph_rows[w * max_graphs + idx + p - 1 - skip_common_neighbour_graphs].set(v);
                            }
                        }
                    }
                }
            }

            idx += number_of_common_neighbour_graphs - skip_common_neighbour_graphs;
        }

        template <typename PossiblySomeOtherBitSetType_>
        auto build_complement_graphs(vector<PossiblySomeOtherBitSetType_> & graph_rows, unsigned size, unsigned & idx) -> void
        {
            for (unsigned v = 0 ; v < size ; ++v)
                for (unsigned w = 0 ; w < size ; ++w)
                    if (! graph_rows[v * max_graphs + 0].test(w))
                        graph_rows[v * max_graphs + idx].set(w);

            ++idx;
        }
    };

    enum class SearchResult
    {
        Aborted,
        Unsatisfiable,
        Satisfiable,
        SatisfiableButKeepGoing,
        Restart
    };

    struct Assignment
    {
        unsigned pattern_vertex;
        unsigned target_vertex;

        auto operator== (const Assignment & other) const -> bool
        {
            return pattern_vertex == other.pattern_vertex && target_vertex == other.target_vertex;
        }

        auto operator!= (const Assignment & other) const -> bool
        {
            return ! (*this == other);
        }
    };

    struct AssignmentInformation
    {
        Assignment assignment;
        bool is_decision;
        int discrepancy_count;
        int choice_count;
    };

    struct Assignments
    {
        vector<AssignmentInformation> values;

        bool contains(const Assignment & assignment) const
        {
            // this should not be a linear scan...
            return values.end() != find_if(values.begin(), values.end(), [&] (const auto & a) {
                    return a.assignment == assignment;
                    });
        }
    };

    template <typename BitSetType_>
    struct SubgraphDomain
    {
        unsigned v;
        unsigned count;
        bool fixed = false;
        BitSetType_ values;

        explicit SubgraphDomain(unsigned s) :
            values(s, 0)
        {
        }

        SubgraphDomain(const SubgraphDomain &) = default;
        SubgraphDomain(SubgraphDomain &&) = default;
    };

    template <typename EntryType_>
    struct AssignmentWatchTable
    {
        unsigned target_size;
        vector<EntryType_> data;

        EntryType_ & operator[] (Assignment x)
        {
            return data[target_size * x.pattern_vertex + x.target_vertex];
        }
    };

    template <typename BitSetType_, typename ArrayType_, typename PatternAdjacencyBitsType_>
    struct Searcher
    {
        using Model = SubgraphModel<BitSetType_, ArrayType_, PatternAdjacencyBitsType_>;
        using Domain = SubgraphDomain<BitSetType_>;
        using Domains = vector<Domain>;

        const Model & model;
        const HomomorphismParams & params;

        Watches<Assignment, AssignmentWatchTable> watches;

        mt19937 global_rand;

        Searcher(const Model & m, const HomomorphismParams & p) :
            model(m),
            params(p)
        {
            // set up space for watches
            if (might_have_watches(params)) {
                watches.table.target_size = model.target_size;
                watches.table.data.resize(model.pattern_size * model.target_size);
            }
        }

        auto assignments_as_proof_decisions(const Assignments & assignments) const -> vector<pair<int, int> >
        {
            vector<pair<int, int> > trail;
            for (auto & a : assignments.values)
                if (a.is_decision)
                    trail.emplace_back(model.pattern_permutation[a.assignment.pattern_vertex],
                            a.assignment.target_vertex);
            return trail;
        }

        auto solution_in_proof_form(const Assignments & assignments) const -> vector<pair<NamedVertex, NamedVertex> >
        {
            vector<pair<NamedVertex, NamedVertex> > solution;
            for (auto & a : assignments.values)
                if (solution.end() == find_if(solution.begin(), solution.end(),
                            [&] (const auto & t) { return t.first.first == model.pattern_permutation[a.assignment.pattern_vertex]; }))
                    solution.emplace_back(model.pattern_vertex_for_proof(a.assignment.pattern_vertex), model.target_vertex_for_proof(a.assignment.target_vertex));
            return solution;
        }

        auto expand_to_full_result(const Assignments & assignments, VertexToVertexMapping & mapping) -> void
        {
            for (auto & a : assignments.values)
                mapping.emplace(model.pattern_permutation.at(a.assignment.pattern_vertex), a.assignment.target_vertex);

            // re-add isolated vertices
            int t = 0;
            for (auto & v : model.isolated_vertices) {
                while (mapping.end() != find_if(mapping.begin(), mapping.end(),
                            [&t] (const pair<int, int> & p) { return p.second == t; }))
                        ++t;
                mapping.emplace(v, t);
            }
        }

        auto find_unit_domain(Domains & domains) -> typename Domains::iterator
        {
            return find_if(domains.begin(), domains.end(), [] (Domain & d) {
                    return (! d.fixed) && 1 == d.count;
                    });
        }

        // The max_graphs_ template parameter is so that the for each graph
        // pair loop gets unrolled, which makes an annoyingly large difference
        // to performance. Note that for larger target graphs, half of the
        // total runtime is spent in this function.
        template <int max_graphs_, bool has_edge_labels_>
        auto propagate_adjacency_constraints(Domain & d, const Assignment & current_assignment) -> void
        {
            auto pattern_adjacency_bits = model.pattern_adjacencies_bits[model.pattern_size * current_assignment.pattern_vertex + d.v];

            // for each graph pair...
            for (unsigned g = 0 ; g < max_graphs_ ; ++g) {
                // if we're adjacent...
                if (pattern_adjacency_bits & (1u << g)) {
                    // ...then we can only be mapped to adjacent vertices
                    d.values &= model.target_graph_rows[current_assignment.target_vertex * max_graphs_ + g];
                }
            }

            if constexpr (has_edge_labels_) {
                // if we're adjacent in the original graph, additionally the edge labels need to match up
                if (pattern_adjacency_bits & (1u << 0)) {
                    auto check_d_values = d.values;

                    auto want_forward_label = model.pattern_edge_labels.at(model.pattern_size * d.v + current_assignment.pattern_vertex);
                    auto want_reverse_label = model.pattern_edge_labels.at(model.pattern_size * current_assignment.pattern_vertex + d.v);
                    for (auto c = check_d_values.find_first() ; c != decltype(check_d_values)::npos ; c = check_d_values.find_first()) {
                        check_d_values.reset(c);

                        auto got_forward_label = model.target_edge_labels.at(model.target_size * c + current_assignment.target_vertex);
                        auto got_reverse_label = model.target_edge_labels.at(model.target_size * current_assignment.target_vertex + c);

                        if (got_forward_label != want_forward_label || got_reverse_label != want_reverse_label)
                            d.values.reset(c);
                    }
                }
            }
        }

        auto both_in_the_neighbourhood_of_some_vertex(unsigned v, unsigned w) -> bool
        {
            auto & nv = model.pattern_graph_rows[v * model.max_graphs + 0];
            auto & nw = model.pattern_graph_rows[w * model.max_graphs + 0];
            return ! (nv & nw).empty();
        }

        auto propagate_simple_constraints(Domains & new_domains, const Assignment & current_assignment) -> bool
        {
            // propagate for each remaining domain...
            for (auto & d : new_domains) {
                if (d.fixed)
                    continue;

                // injectivity
                switch (params.injectivity) {
                    case Injectivity::Injective:
                        d.values.reset(current_assignment.target_vertex);
                        break;
                    case Injectivity::LocallyInjective:
                        if (both_in_the_neighbourhood_of_some_vertex(current_assignment.pattern_vertex, d.v))
                            d.values.reset(current_assignment.target_vertex);
                        break;
                    case Injectivity::NonInjective:
                        break;
                }

                // adjacency
                if (model.pattern_edge_labels.empty()) {
                    switch (model.max_graphs) {
                        case 1:  propagate_adjacency_constraints<1, false>(d, current_assignment); break;
                        case 2:  propagate_adjacency_constraints<2, false>(d, current_assignment); break;
                        case 3:  propagate_adjacency_constraints<3, false>(d, current_assignment); break;
                        case 4:  propagate_adjacency_constraints<4, false>(d, current_assignment); break;
                        case 5:  propagate_adjacency_constraints<5, false>(d, current_assignment); break;
                        case 6:  propagate_adjacency_constraints<6, false>(d, current_assignment); break;
                        case 7:  propagate_adjacency_constraints<7, false>(d, current_assignment); break;
                        case 8:  propagate_adjacency_constraints<8, false>(d, current_assignment); break;
                        case 9:  propagate_adjacency_constraints<9, false>(d, current_assignment); break;
                        case 10: propagate_adjacency_constraints<10, false>(d, current_assignment); break;

                        default:
                            throw "you forgot to update the ugly max_graphs hack";
                    }
                }
                else {
                    switch (model.max_graphs) {
                        case 1:  propagate_adjacency_constraints<1, true>(d, current_assignment); break;
                        case 2:  propagate_adjacency_constraints<2, true>(d, current_assignment); break;
                        case 3:  propagate_adjacency_constraints<3, true>(d, current_assignment); break;
                        case 4:  propagate_adjacency_constraints<4, true>(d, current_assignment); break;
                        case 5:  propagate_adjacency_constraints<5, true>(d, current_assignment); break;
                        case 6:  propagate_adjacency_constraints<6, true>(d, current_assignment); break;
                        case 7:  propagate_adjacency_constraints<7, true>(d, current_assignment); break;
                        case 8:  propagate_adjacency_constraints<8, true>(d, current_assignment); break;
                        case 9:  propagate_adjacency_constraints<9, true>(d, current_assignment); break;
                        case 10: propagate_adjacency_constraints<10, true>(d, current_assignment); break;

                        default:
                            throw "you forgot to update the ugly max_graphs hack";
                    }
                }

                // we might have removed values
                d.count = d.values.count();
                if (0 == d.count)
                    return false;
            }

            return true;
        }

        auto propagate_less_thans(Domains & new_domains) -> bool
        {
            ArrayType_ find_domain;
            if constexpr (is_same<ArrayType_, vector<int> >::value)
                find_domain.resize(model.pattern_size);
            fill(find_domain.begin(), find_domain.end(), -1);

            for (unsigned i = 0, i_end = new_domains.size() ; i != i_end ; ++i)
                find_domain[new_domains[i].v] = i;

            for (auto & [ a, b ] : model.pattern_less_thans_in_convenient_order) {
                if (find_domain[a] == -1 || find_domain[b] == -1)
                    continue;
                auto & a_domain = new_domains[find_domain[a]];
                auto & b_domain = new_domains[find_domain[b]];

               // first value of b must be at least one after the first possible value of a
               auto first_a = a_domain.values.find_first();
               if (first_a == decltype(a_domain.values)::npos)
                   return false;
               auto first_allowed_b = first_a + 1;

               if (first_allowed_b >= model.target_size)
                   return false;

               for (auto v = b_domain.values.find_first() ; v != decltype(b_domain.values)::npos ; v = b_domain.values.find_first()) {
                   if (v >= first_allowed_b)
                       break;
                   b_domain.values.reset(v);
               }

               // b might have shrunk (and detect empty before the next bit to make life easier)
               b_domain.count = b_domain.values.count();
               if (0 == b_domain.count)
                   return false;
            }

            for (auto & [ a, b ] : model.pattern_less_thans_in_convenient_order) {
                if (find_domain[a] == -1 || find_domain[b] == -1)
                    continue;
                auto & a_domain = new_domains[find_domain[a]];
                auto & b_domain = new_domains[find_domain[b]];

                // last value of a must be at least one before the last possible value of b
                auto b_values_copy = b_domain.values;
                auto last_b = b_domain.values.find_first();
                for (auto v = last_b ; v != decltype(b_values_copy)::npos ; v = b_values_copy.find_first()) {
                    b_values_copy.reset(v);
                    last_b = v;
                }

                if (last_b == 0)
                    return false;
                auto last_allowed_a = last_b - 1;

                auto a_values_copy = a_domain.values;
                for (auto v = a_values_copy.find_first() ; v != decltype(a_values_copy)::npos ; v = a_values_copy.find_first()) {
                    a_values_copy.reset(v);
                    if (v > last_allowed_a)
                        a_domain.values.reset(v);
                }

                // a might have shrunk
                a_domain.count = a_domain.values.count();
                if (0 == a_domain.count)
                    return false;
            }

            return true;
        }

        auto propagate(Domains & new_domains, Assignments & assignments) -> bool
        {
            // whilst we've got a unit domain...
            for (typename Domains::iterator branch_domain = find_unit_domain(new_domains) ;
                    branch_domain != new_domains.end() ;
                    branch_domain = find_unit_domain(new_domains)) {
                // what are we assigning?
                Assignment current_assignment = { branch_domain->v, unsigned(branch_domain->values.find_first()) };

                // ok, make the assignment
                branch_domain->fixed = true;
                assignments.values.push_back({ current_assignment, false, -1, -1 });

                if (params.proof)
                    params.proof->unit_propagating(
                            model.pattern_vertex_for_proof(current_assignment.pattern_vertex),
                            model.target_vertex_for_proof(current_assignment.target_vertex));

                // propagate watches
                if (might_have_watches(params))
                    watches.propagate(current_assignment,
                            [&] (const Assignment & a) { return ! assignments.contains(a); },
                            [&] (const Assignment & a) {
                                    for (auto & d : new_domains) {
                                        if (d.fixed)
                                            continue;

                                        if (d.v == a.pattern_vertex) {
                                            d.values.reset(a.target_vertex);
                                            break;
                                        }
                                    }
                                });

                // propagate simple all different and adjacency
                if (! propagate_simple_constraints(new_domains, current_assignment))
                    return false;

                // propagate less than
                if (model.has_less_thans && ! propagate_less_thans(new_domains))
                    return false;

                // propagate all different
                if (params.injectivity == Injectivity::Injective)
                    if (! (params.proof ? cheap_all_different<true>(new_domains) : cheap_all_different<false>(new_domains)))
                        return false;
            }

            return true;
        }

        auto find_branch_domain(const Domains & domains) -> const Domain *
        {
            const Domain * result = nullptr;
            for (auto & d : domains)
                if (! d.fixed)
                    if ((! result) ||
                            (d.count < result->count) ||
                            (d.count == result->count && model.patterns_degrees[0][d.v] > model.patterns_degrees[0][result->v]))
                        result = &d;
            return result;
        }

        auto copy_nonfixed_domains_and_make_assignment(
                const Domains & domains,
                unsigned branch_v,
                unsigned f_v) -> Domains
        {
            Domains new_domains;
            new_domains.reserve(domains.size());
            for (auto & d : domains) {
                if (d.fixed)
                    continue;

                new_domains.push_back(d);
                if (d.v == branch_v) {
                    new_domains.back().values.reset();
                    new_domains.back().values.set(f_v);
                    new_domains.back().count = 1;
                }
            }
            return new_domains;
        }

        auto post_nogood(
                const Assignments & assignments)
        {
            if (! might_have_watches(params))
                return;

            Nogood<Assignment> nogood;

            for (auto & a : assignments.values)
                if (a.is_decision)
                    nogood.literals.emplace_back(a.assignment);

            watches.post_nogood(move(nogood));

            if (params.proof)
                params.proof->post_restart_nogood(assignments_as_proof_decisions(assignments));
        }

        auto softmax_shuffle(
                ArrayType_ & branch_v,
                unsigned branch_v_end
                ) -> void
        {
            // repeatedly pick a softmax-biased vertex, move it to the front of branch_v,
            // and then only consider items further to the right in the next iteration.

            // Using floating point here turned out to be way too slow. Fortunately the base
            // of softmax doesn't seem to matter, so we use 2 instead of e, and do everything
            // using bit voodoo.
            auto expish = [largest_target_degree = this->model.largest_target_degree] (int degree) {
                constexpr int sufficient_space_for_adding_up = numeric_limits<long long>::digits - 18;
                auto shift = max<int>(degree - largest_target_degree + sufficient_space_for_adding_up, 0);
                return 1ll << shift;
            };

            long long total = 0;
            for (unsigned v = 0 ; v < branch_v_end ; ++v)
                total += expish(model.targets_degrees[0][branch_v[v]]);

            for (unsigned start = 0 ; start < branch_v_end ; ++start) {
                // pick a random number between 1 and total inclusive
                uniform_int_distribution<long long> dist(1, total);
                long long select_score = dist(global_rand);

                // go over the list until we hit the score
                unsigned select_element = start;
                for ( ; select_element + 1 < branch_v_end ; ++select_element) {
                    select_score -= expish(model.targets_degrees[0][branch_v[select_element]]);
                    if (select_score <= 0)
                        break;
                }

                // move to front
                total -= expish(model.targets_degrees[0][branch_v[select_element]]);
                swap(branch_v[select_element], branch_v[start]);
            }
        }

        auto degree_sort(
                ArrayType_ & branch_v,
                unsigned branch_v_end,
                bool reverse
                ) -> void
        {
            stable_sort(branch_v.begin(), branch_v.begin() + branch_v_end, [&] (int a, int b) -> bool {
                    return (model.targets_degrees[0][a] >= model.targets_degrees[0][b]) ^ reverse;
                    });
        }

        auto bigraph_link_match(
            std::pair<bool, std::vector<int> > pattern_hyperedge,
            std::pair<bool, std::vector<int> > target_hyperedge,
            VertexToVertexMapping mapping) -> bool
        {
            // If pattern is closed and target is open, cannot ever match
            if (target_hyperedge.first && !pattern_hyperedge.first) return false;

            // Closed pattern hyperedge must match exactly with the target hyperedge
            if (!pattern_hyperedge.first) {
                for(unsigned i=0; i<pattern_hyperedge.second.size(); i++)
                    if(pattern_hyperedge.second[i] != target_hyperedge.second[mapping[i]]) return false;
            }
            // Open pattern hyperedges can be subsets of the target hyperedge
            else{          
                for(unsigned i=0; i<pattern_hyperedge.second.size(); i++)
                    if(pattern_hyperedge.second[i] > target_hyperedge.second[mapping[i]]) return false;
            }
            return true;
        }


        auto restarting_search(
                Assignments & assignments,
                const Domains & domains,
                unsigned long long & nodes,
                unsigned long long & propagations,
                unsigned long long & solution_count,
                int depth,
                RestartsSchedule & restarts_schedule) -> SearchResult
        {
            if (params.timeout->should_abort())
                return SearchResult::Aborted;

            ++nodes;

            // find ourselves a domain, or succeed if we're all assigned
            const Domain * branch_domain = find_branch_domain(domains);
            if (! branch_domain) {

                // do stuff here
                // assignments -> contains complete set of everything
                // return SearchResult::Unsatisfiable if transitive closure violated
                //  - or if site/root neighbourhood check is violated
                //  - or if link graph is violated

                if (params.bigraph) {

                    VertexToVertexMapping mapping;
                    expand_to_full_result(assignments, mapping);

                    bool lazy_flag;
                    std::set<int> mapped_targets;

                    // Eliminate closed-pattern/closed-target hyperedges first
                    for(unsigned i=0; i<model.pattern_hyperedges.size(); i++){
                        if(!model.pattern_hyperedges[i].first) {
                            lazy_flag = false;
                            for(unsigned j=0; j<model.target_hyperedges.size(); j++) {
                                if(mapped_targets.find(j) == mapped_targets.end() && 
                                   bigraph_link_match(model.pattern_hyperedges[i], model.target_hyperedges[j], mapping)) {
                                    mapped_targets.insert(j);
                                    lazy_flag = true;             
                                    break;
                                }             
                            }
                            if (!lazy_flag) return SearchResult::Unsatisfiable;
                        }
                    }  
                             
                    // Clean up open pattern hyperedges by matching arbitrarily on remaining target hyperedges    
                    // (by nature of bigraphs, this should be sound)     

                    for(unsigned i=0; i<model.pattern_hyperedges.size(); i++){
                        if(model.pattern_hyperedges[i].first) {
                            lazy_flag = false;
                            for(unsigned j=0; j<model.target_hyperedges.size(); j++) {
                                if(mapped_targets.find(j) == mapped_targets.end() && 
                                   bigraph_link_match(model.pattern_hyperedges[i], model.target_hyperedges[j], mapping)) {
                                   lazy_flag = true;             
                                   break;
                                }             
                            }
                            if (!lazy_flag) return SearchResult::Unsatisfiable;
                        }
                    }

                    //Find transitive closure violations
                    for(unsigned i=0; i<model.pattern_size; i++){
                        if(model.pattern_big_constraints[i].second) { 
                            std::set<int> child_mappings;

                            // Get all corresponding target node's children with a mapping
                            for(unsigned j=0; j<model.pattern_size; j++) 
                                if(i != j &&
                                    model.pattern_graph_rows[i].test(j) && 
                                    model.pattern_graph_reachability[i][j])
                                        child_mappings.insert(mapping[j]);
             

                            // For all target node's children without a mapping, check if it can reach any children of a root node
                            for(unsigned j=0; j<model.target_size; j++)
                                if(mapping[i] != j && 
                                    model.target_graph_rows[mapping[i]].test(j) && 
                                    model.target_graph_reachability[mapping[i]][j] &&
                                    child_mappings.find(j) == child_mappings.end())
                                        for(unsigned k=0; k<model.pattern_size; k++) 
                                            if(model.pattern_big_constraints[k].first && model.target_graph_reachability[j][mapping[k]])
                                                return SearchResult::Unsatisfiable;                                   
                        }
                    }
               }    

                if (params.lackey) {
                    VertexToVertexMapping mapping;
                    expand_to_full_result(assignments, mapping);
                    if (! params.lackey->check_solution(mapping))
                        return SearchResult::Unsatisfiable;
                }

                if (params.proof)
                    params.proof->post_solution(solution_in_proof_form(assignments));

                if (params.count_solutions) {
                    ++solution_count;
                    if (params.enumerate_callback) {
                        VertexToVertexMapping mapping;
                        expand_to_full_result(assignments, mapping);
                        params.enumerate_callback(mapping);
                    }

                    return SearchResult::SatisfiableButKeepGoing;
                }
                else
                    return SearchResult::Satisfiable;
            }

            // pull out the remaining values in this domain for branching
            auto remaining = branch_domain->values;

            ArrayType_ branch_v;
            if constexpr (is_same<ArrayType_, vector<int> >::value)
                branch_v.resize(model.target_size);

            unsigned branch_v_end = 0;
            for (auto f_v = remaining.find_first() ; f_v != decltype(remaining)::npos ; f_v = remaining.find_first()) {
                remaining.reset(f_v);
                branch_v[branch_v_end++] = f_v;
            }

            switch (params.value_ordering_heuristic) {
                case ValueOrdering::Degree:
                    degree_sort(branch_v, branch_v_end, false);
                    break;

                case ValueOrdering::AntiDegree:
                    degree_sort(branch_v, branch_v_end, true);
                    break;

                case ValueOrdering::Biased:
                    softmax_shuffle(branch_v, branch_v_end);
                    break;

                case ValueOrdering::Random:
                    shuffle(branch_v.begin(), branch_v.begin() + branch_v_end, global_rand);
                    break;
            }

            int discrepancy_count = 0;
            bool actually_hit_a_failure = false;

            // for each value remaining...
            for (auto f_v = branch_v.begin(), f_end = branch_v.begin() + branch_v_end ; f_v != f_end ; ++f_v) {
                if (params.proof)
                    params.proof->guessing(depth, model.pattern_vertex_for_proof(branch_domain->v), model.target_vertex_for_proof(*f_v));

                // modified in-place by appending, we can restore by shrinking
                auto assignments_size = assignments.values.size();

                // make the assignment
                assignments.values.push_back({ { branch_domain->v, unsigned(*f_v) }, true, discrepancy_count, int(branch_v_end) });

                // set up new domains
                Domains new_domains = copy_nonfixed_domains_and_make_assignment(domains, branch_domain->v, *f_v);

                // propagate
                ++propagations;
                if (! propagate(new_domains, assignments)) {
                    // failure? restore assignments and go on to the next thing
                    if (params.proof)
                        params.proof->propagation_failure(assignments_as_proof_decisions(assignments), model.pattern_vertex_for_proof(branch_domain->v), model.target_vertex_for_proof(*f_v));

                    assignments.values.resize(assignments_size);
                    actually_hit_a_failure = true;

                    continue;
                }

                if (params.proof)
                    params.proof->start_level(depth + 2);

                // recursive search
                auto search_result = restarting_search(assignments, new_domains, nodes, propagations,
                        solution_count, depth + 1, restarts_schedule);

                switch (search_result) {
                    case SearchResult::Satisfiable:
                        return SearchResult::Satisfiable;

                    case SearchResult::Aborted:
                        return SearchResult::Aborted;

                    case SearchResult::Restart:
                        // restore assignments before posting nogoods, it's easier
                        assignments.values.resize(assignments_size);

                        // post nogoods for everything we've done so far
                        for (auto l = branch_v.begin() ; l != f_v ; ++l) {
                            assignments.values.push_back({ { branch_domain->v, unsigned(*l) }, true, -2, -2 });
                            post_nogood(assignments);
                            assignments.values.pop_back();
                        }

                        return SearchResult::Restart;

                    case SearchResult::SatisfiableButKeepGoing:
                        if (params.proof) {
                            params.proof->back_up_to_level(depth + 1);
                            params.proof->incorrect_guess(assignments_as_proof_decisions(assignments), false);
                        }

                        // restore assignments
                        assignments.values.resize(assignments_size);
                        break;

                    case SearchResult::Unsatisfiable:
                        if (params.proof) {
                            params.proof->back_up_to_level(depth + 1);
                            params.proof->incorrect_guess(assignments_as_proof_decisions(assignments), true);
                        }

                        // restore assignments
                        assignments.values.resize(assignments_size);
                        actually_hit_a_failure = true;

                        break;
                }

                ++discrepancy_count;
            }

            // no values remaining, backtrack, or possibly kick off a restart
            if (params.proof)
                params.proof->out_of_guesses(assignments_as_proof_decisions(assignments));

            if (actually_hit_a_failure)
                restarts_schedule.did_a_backtrack();

            if (restarts_schedule.should_restart()) {
                if (params.proof)
                    params.proof->back_up_to_top();
                post_nogood(assignments);
                return SearchResult::Restart;
            }
            else
                return SearchResult::Unsatisfiable;
        }

        template <bool proof_>
        auto cheap_all_different(Domains & domains) -> bool
        {
            // Pick domains smallest first; ties are broken by smallest .v first.
            // For each count p we have a linked list, whose first member is
            // first[p].  The element following x in one of these lists is next[x].
            // Any domain with a count greater than domains.size() is put
            // int the "count==domains.size()" bucket.
            // The "first" array is sized to be able to hold domains.size()+1
            // elements
            ArrayType_ first, next;

            if constexpr (is_same<ArrayType_, vector<int> >::value)
                first.resize(model.target_size + 1);
            fill(first.begin(), first.begin() + domains.size() + 1, -1);

            if constexpr (is_same<ArrayType_, vector<int> >::value)
                next.resize(model.target_size);
            fill(next.begin(), next.begin() + domains.size(), -1);

            [[ maybe_unused ]] conditional_t<proof_, vector<int>, tuple<> > lhs, hall_lhs, hall_rhs;

            // Iterate backwards, because we insert elements at the head of
            // lists and we want the sort to be stable
            for (int i = int(domains.size()) - 1 ; i >= 0; --i) {
                unsigned count = domains.at(i).count;
                if (count > domains.size())
                    count = domains.size();
                next.at(i) = first.at(count);
                first.at(count) = i;
            }

            // counting all-different
            BitSetType_ domains_so_far{ model.target_size, 0 }, hall{ model.target_size, 0 };
            unsigned neighbours_so_far = 0;

            [[ maybe_unused ]] conditional_t<proof_, unsigned, tuple<> > last_outputted_hall_size{};

            for (unsigned i = 0 ; i <= domains.size() ; ++i) {
                // iterate over linked lists
                int domain_index = first[i];
                while (domain_index != -1) {
                    auto & d = domains.at(domain_index);

                    if constexpr (proof_)
                        lhs.push_back(model.pattern_permutation[d.v]);

                    [[ maybe_unused ]] conditional_t<proof_, unsigned, tuple<> > old_d_values_count;
                    if constexpr (proof_)
                        old_d_values_count = d.values.count();

                    d.values &= ~hall;
                    d.count = d.values.count();

                    if constexpr (proof_)
                        if (last_outputted_hall_size != hall.count() && d.count != old_d_values_count) {
                            last_outputted_hall_size = hall.count();
                            params.proof->emit_hall_set_or_violator(hall_lhs, hall_rhs);
                        }

                    if (0 == d.count)
                        return false;

                    domains_so_far |= d.values;
                    ++neighbours_so_far;

                    unsigned domains_so_far_popcount = domains_so_far.count();

                    if (domains_so_far_popcount < neighbours_so_far) {
                        // hall violator, so we fail (after outputting a proof)
                        if constexpr (proof_) {
                            vector<int> rhs;
                            auto d = domains_so_far;
                            for (auto v = d.find_first() ; v != decltype(d)::npos ; v = d.find_first()) {
                                d.reset(v);
                                rhs.push_back(v);
                            }
                            params.proof->emit_hall_set_or_violator(lhs, rhs);
                        }
                        return false;
                    }
                    else if (domains_so_far_popcount == neighbours_so_far) {
                        // equivalent to hall=domains_so_far
                        hall |= domains_so_far;
                        if constexpr (proof_) {
                            hall_lhs = lhs;
                            hall_rhs.clear();
                            auto d = domains_so_far;
                            for (auto v = d.find_first() ; v != decltype(d)::npos ; v = d.find_first()) {
                                d.reset(v);
                                hall_rhs.push_back(v);
                            }
                        }
                    }
                    domain_index = next[domain_index];
                }
            }

            return true;
        }

        auto save_result(const Assignments & assignments, HomomorphismResult & result) -> void
        {
            expand_to_full_result(assignments, result.mapping);

            string where = "where =";
            for (auto & a : assignments.values)
                where.append(" " + to_string(a.discrepancy_count) + "/" + to_string(a.choice_count));
            result.extra_stats.push_back(where);
        }
    };

    template <typename BitSetType_, typename ArrayType_, typename PatternAdjacencyBitsType_>
    struct BasicSolver
    {
        using Model = SubgraphModel<BitSetType_, ArrayType_, PatternAdjacencyBitsType_>;
        using Domain = SubgraphDomain<BitSetType_>;
        using Domains = vector<Domain>;

        const Model & model;
        const HomomorphismParams & params;

        BasicSolver(const Model & m, const HomomorphismParams & p) :
            model(m),
            params(p)
        {
        }

        auto check_label_compatibility(int p, int t) -> bool
        {
            if (model.pattern_vertex_labels.empty())
                return true;
            else
                return model.pattern_vertex_labels[p] == model.target_vertex_labels[t];
        }

        auto check_loop_compatibility(int p, int t) -> bool
        {
            for (unsigned g = 0 ; g < model.max_graphs ; ++g)
                if (model.pattern_graph_rows[p * model.max_graphs + g].test(p) && ! model.target_graph_rows[t * model.max_graphs + g].test(t))
                    return false;

            return true;
        }

        auto check_bigraph_degree_compatibility(int p, int t) -> bool
        {
            if ((! model.pattern_big_constraints[p].first) && (model.pattern_dir_degrees[p].first != model.target_dir_degrees[t].first))
                return false;
            if (model.pattern_big_constraints[p].first && (model.pattern_dir_degrees[p].first > model.target_dir_degrees[t].first))
                return false;
            if ((! model.pattern_big_constraints[p].second) && (model.pattern_dir_degrees[p].second != model.target_dir_degrees[t].second))
                return false;
            if (model.pattern_big_constraints[p].second && (model.pattern_dir_degrees[p].second > model.target_dir_degrees[t].second))
                return false;
            return true;
        }

        auto check_degree_compatibility(
                int p,
                int t,
                unsigned graphs_to_consider,
                vector<vector<vector<int> > > & patterns_ndss,
                vector<vector<optional<vector<int> > > > & targets_ndss
                ) -> bool
        {
            if (! degree_and_nds_are_preserved(params))
                return true;

            for (unsigned g = 0 ; g < graphs_to_consider ; ++g) {
                if (model.targets_degrees.at(g).at(t) < model.patterns_degrees.at(g).at(p)) {
                    // not ok, degrees differ
                    if (params.proof) {
                        // get the actual neighbours of p and t, in their original terms
                        vector<int> n_p, n_t;

                        auto np = model.pattern_graph_rows[p * model.max_graphs + g];
                        for (auto j = np.find_first() ; j != decltype(np)::npos ; j = np.find_first()) {
                            np.reset(j);
                            n_p.push_back(model.pattern_permutation[j]);
                        }

                        auto nt = model.target_graph_rows[t * model.max_graphs + g];
                        for (auto j = nt.find_first() ; j != decltype(nt)::npos ; j = nt.find_first()) {
                            nt.reset(j);
                            n_t.push_back(j);
                        }

                        params.proof->incompatible_by_degrees(g, model.pattern_vertex_for_proof(p), n_p, model.target_vertex_for_proof(t), n_t);
                    }
                    return false;
                }
                else if (degree_and_nds_are_exact(params, model.full_pattern_size, model.target_size)
                        && model.targets_degrees.at(g).at(t) != model.patterns_degrees.at(g).at(p)) {
                    // not ok, degrees must be exactly the same
                    return false;
                }
            }
            if (params.no_nds)
                return true;

            // full compare of neighbourhood degree sequences
            if (! targets_ndss.at(0).at(t)) {
                for (unsigned g = 0 ; g < graphs_to_consider ; ++g) {
                    targets_ndss.at(g).at(t) = vector<int>{};
                    auto ni = model.target_graph_rows[t * model.max_graphs + g];
                    for (auto j = ni.find_first() ; j != decltype(ni)::npos ; j = ni.find_first()) {
                        ni.reset(j);
                        targets_ndss.at(g).at(t)->push_back(model.targets_degrees.at(g).at(j));
                    }
                    sort(targets_ndss.at(g).at(t)->begin(), targets_ndss.at(g).at(t)->end(), greater<int>());
                }
            }

            for (unsigned g = 0 ; g < graphs_to_consider ; ++g) {
                for (unsigned x = 0 ; x < patterns_ndss.at(g).at(p).size() ; ++x) {
                    if (targets_ndss.at(g).at(t)->at(x) < patterns_ndss.at(g).at(p).at(x)) {
                        if (params.proof)
                            params.proof->incompatible_by_nds(g, p, t);
                        return false;
                    }
                    else if (degree_and_nds_are_exact(params, model.full_pattern_size, model.target_size)
                            && targets_ndss.at(g).at(t)->at(x) != patterns_ndss.at(g).at(p).at(x))
                        return false;
                }
            }

            return true;
        }

        auto initialise_domains(Domains & domains) -> bool
        {
            unsigned graphs_to_consider = model.max_graphs;
            if (params.induced) {
                // when looking at the complement graph, if the largest degree
                // in the pattern is smaller than the smallest degree in the
                // target, then we're going to spend a lot of time doing
                // nothing useful
                auto largest_pattern_c_degree = max_element(model.patterns_degrees[model.max_graphs - 1].begin(), model.patterns_degrees[model.max_graphs - 1].end());
                auto smallest_target_c_degree = min_element(model.targets_degrees[model.max_graphs - 1].begin(), model.targets_degrees[model.max_graphs - 1].end());
                if (largest_pattern_c_degree != model.patterns_degrees[model.max_graphs - 1].end() &&
                        smallest_target_c_degree != model.targets_degrees[model.max_graphs - 1].end() &&
                        *largest_pattern_c_degree < *smallest_target_c_degree)
                    --graphs_to_consider;
            }

            /* pattern and target neighbourhood degree sequences */
            vector<vector<vector<int> > > patterns_ndss(graphs_to_consider);
            vector<vector<optional<vector<int> > > > targets_ndss(graphs_to_consider);

            if (degree_and_nds_are_preserved(params) && ! params.no_nds) {
                for (unsigned g = 0 ; g < graphs_to_consider ; ++g) {
                    patterns_ndss.at(g).resize(model.pattern_size);
                    targets_ndss.at(g).resize(model.target_size);
                }

                for (unsigned g = 0 ; g < graphs_to_consider ; ++g) {
                    for (unsigned i = 0 ; i < model.pattern_size ; ++i) {
                        auto ni = model.pattern_graph_rows[i * model.max_graphs + g];
                        for (auto j = ni.find_first() ; j != decltype(ni)::npos ; j = ni.find_first()) {
                            ni.reset(j);
                            patterns_ndss.at(g).at(i).push_back(model.patterns_degrees.at(g).at(j));
                        }
                        sort(patterns_ndss.at(g).at(i).begin(), patterns_ndss.at(g).at(i).end(), greater<int>());
                    }
                }
            }

            for (unsigned i = 0 ; i < model.pattern_size ; ++i) {
                domains.at(i).v = i;
                domains.at(i).values.reset();

                for (unsigned j = 0 ; j < model.target_size ; ++j) {
                    bool ok = true;

                    if (! check_label_compatibility(i, j))
                        ok = false;
                    else if (! check_loop_compatibility(i, j))
                        ok = false;
                    else if (! check_degree_compatibility(i, j, graphs_to_consider, patterns_ndss, targets_ndss))
                        ok = false;
                    else if (params.bigraph && ! check_bigraph_degree_compatibility(i, j))
                        ok = false;
                    if (ok)
                        domains.at(i).values.set(j);
                }

                domains.at(i).count = domains.at(i).values.count();
                if (0 == domains.at(i).count) {
                    if (params.proof)
                        params.proof->initial_domain_is_empty(i);
                    return false;
                }
            }

            // quick sanity check that we have enough values
            if (is_nonshrinking(params)) {
                BitSetType_ domains_union{ model.target_size, 0 };
                for (auto & d : domains)
                    domains_union |= d.values;

                unsigned domains_union_popcount = domains_union.count();
                if (domains_union_popcount < unsigned(model.pattern_size)) {
                    if (params.proof) {
                        vector<int> hall_lhs, hall_rhs;
                        for (auto & d : domains)
                            hall_lhs.push_back(model.pattern_permutation[d.v]);
                        auto dd = domains_union;
                        for (auto v = dd.find_first() ; v != decltype(dd)::npos ; v = dd.find_first()) {
                            dd.reset(v);
                            hall_rhs.push_back(v);
                        }
                        params.proof->emit_hall_set_or_violator(hall_lhs, hall_rhs);
                    }
                    return false;
                }
            }

            for (auto & d : domains)
                d.count = d.values.count();

            return true;
        }
    };

    template <typename BitSetType_, typename ArrayType_, typename PatternAdjacencyBitsType_>
    struct SequentialSolver :
        BasicSolver<BitSetType_, ArrayType_, PatternAdjacencyBitsType_>
    {
        using BasicSolver<BitSetType_, ArrayType_, PatternAdjacencyBitsType_>::BasicSolver;

        using Model = typename BasicSolver<BitSetType_, ArrayType_, PatternAdjacencyBitsType_>::Model;
        using Domain = typename BasicSolver<BitSetType_, ArrayType_, PatternAdjacencyBitsType_>::Domain;
        using Domains = typename BasicSolver<BitSetType_, ArrayType_, PatternAdjacencyBitsType_>::Domains;

        using BasicSolver<BitSetType_, ArrayType_, PatternAdjacencyBitsType_>::model;
        using BasicSolver<BitSetType_, ArrayType_, PatternAdjacencyBitsType_>::params;

        using BasicSolver<BitSetType_, ArrayType_, PatternAdjacencyBitsType_>::initialise_domains;

        auto solve() -> HomomorphismResult
        {
            HomomorphismResult result;

            // domains
            Domains domains(model.pattern_size, Domain{ model.target_size });
            if (! initialise_domains(domains)) {
                result.complete = true;
                return result;
            }

            // assignments
            Assignments assignments;
            assignments.values.reserve(model.pattern_size);

            // start search timer
            auto search_start_time = steady_clock::now();

            // do the search
            bool done = false;
            unsigned number_of_restarts = 0;

            Searcher<BitSetType_, ArrayType_, PatternAdjacencyBitsType_> searcher(model, params);

            while (! done) {
                ++number_of_restarts;

                // start watching new nogoods
                done = searcher.watches.apply_new_nogoods(
                        [&] (const Assignment & assignment) {
                            for (auto & d : domains)
                                if (d.v == assignment.pattern_vertex) {
                                    d.values.reset(assignment.target_vertex);
                                    d.count = d.values.count();
                                    break;
                                }
                        });

                if (done)
                    break;

                searcher.watches.clear_new_nogoods();

                ++result.propagations;
                if (searcher.propagate(domains, assignments)) {
                    auto assignments_copy = assignments;

                    switch (searcher.restarting_search(assignments_copy, domains, result.nodes, result.propagations,
                                result.solution_count, 0, *params.restarts_schedule)) {
                        case SearchResult::Satisfiable:
                            searcher.save_result(assignments_copy, result);
                            result.complete = true;
                            done = true;
                            break;

                        case SearchResult::SatisfiableButKeepGoing:
                            result.complete = true;
                            done = true;
                            break;

                        case SearchResult::Unsatisfiable:
                            result.complete = true;
                            done = true;
                            break;

                        case SearchResult::Aborted:
                            done = true;
                            break;

                        case SearchResult::Restart:
                            break;
                    }
                }
                else {
                    if (params.proof)
                        params.proof->root_propagation_failed();
                    result.complete = true;
                    done = true;
                }

                params.restarts_schedule->did_a_restart();
            }

            if (params.restarts_schedule->might_restart())
                result.extra_stats.emplace_back("restarts = " + to_string(number_of_restarts));

            result.extra_stats.emplace_back("shape_graphs = " + to_string(model.max_graphs));

            result.extra_stats.emplace_back("search_time = " + to_string(
                        duration_cast<milliseconds>(steady_clock::now() - search_start_time).count()));

            if (might_have_watches(params)) {
                result.extra_stats.emplace_back("nogoods_size = " + to_string(searcher.watches.nogoods.size()));

                map<int, int> nogoods_lengths;
                for (auto & n : searcher.watches.nogoods)
                    nogoods_lengths[n.literals.size()]++;

                string nogoods_lengths_str;
                for (auto & n : nogoods_lengths) {
                    nogoods_lengths_str += " ";
                    nogoods_lengths_str += to_string(n.first) + ":" + to_string(n.second);
                }
                result.extra_stats.emplace_back("nogoods_lengths =" + nogoods_lengths_str);
            }

            if (can_strip_isolated_vertices(params) && ! model.isolated_vertices.empty())
                result.extra_stats.emplace_back("isolated_vertices_removed = " + to_string(model.isolated_vertices.size()));

            return result;
        }
    };

    template <typename BitSetType_, typename ArrayType_, typename PatternAdjacencyBitsType_>
    struct ThreadedSolver :
        BasicSolver<BitSetType_, ArrayType_, PatternAdjacencyBitsType_>
    {
        using Model = typename BasicSolver<BitSetType_, ArrayType_, PatternAdjacencyBitsType_>::Model;
        using Domain = typename BasicSolver<BitSetType_, ArrayType_, PatternAdjacencyBitsType_>::Domain;
        using Domains = typename BasicSolver<BitSetType_, ArrayType_, PatternAdjacencyBitsType_>::Domains;

        using BasicSolver<BitSetType_, ArrayType_, PatternAdjacencyBitsType_>::model;
        using BasicSolver<BitSetType_, ArrayType_, PatternAdjacencyBitsType_>::params;

        using BasicSolver<BitSetType_, ArrayType_, PatternAdjacencyBitsType_>::initialise_domains;

        unsigned n_threads;

        ThreadedSolver(const Model & m, const HomomorphismParams & p, unsigned t) :
            BasicSolver<BitSetType_, ArrayType_, PatternAdjacencyBitsType_>(m, p),
            n_threads(t)
        {
        }

        auto solve() -> HomomorphismResult
        {
            mutex common_result_mutex;
            HomomorphismResult common_result;
            string by_thread_nodes, by_thread_propagations;

            // domains
            Domains common_domains(model.pattern_size, Domain{ model.target_size });
            if (! initialise_domains(common_domains)) {
                common_result.complete = true;
                return common_result;
            }

            // start search timer
            auto search_start_time = steady_clock::now();

            vector<thread> threads;
            threads.reserve(n_threads);

            vector<unique_ptr<Searcher<BitSetType_, ArrayType_, PatternAdjacencyBitsType_> > > searchers{ n_threads };

            barrier wait_for_new_nogoods_barrier{ n_threads }, synced_nogoods_barrier{ n_threads };
            atomic<bool> restart_synchroniser{ false };

            function<auto (unsigned) -> void> work_function = [&searchers, &common_domains, &threads, &work_function,
                        &model = this->model, &params = this->params, n_threads = this->n_threads,
                        &common_result, &common_result_mutex, &by_thread_nodes, &by_thread_propagations,
                        &wait_for_new_nogoods_barrier, &synced_nogoods_barrier, &restart_synchroniser] (unsigned t) -> void
            {
                // do the search
                HomomorphismResult thread_result;

                bool just_the_first_thread = (0 == t) && params.delay_thread_creation;

                searchers[t] = make_unique<Searcher<BitSetType_, ArrayType_, PatternAdjacencyBitsType_> >(model, params);
                if (0 != t)
                    searchers[t]->global_rand.seed(t);

                unsigned number_of_restarts = 0;

                Domains domains = common_domains;

                Assignments thread_assignments;
                thread_assignments.values.reserve(model.pattern_size);

                // each thread needs its own restarts schedule
                unique_ptr<RestartsSchedule> thread_restarts_schedule;
                if (0 == t || ! params.triggered_restarts)
                    thread_restarts_schedule.reset(params.restarts_schedule->clone());
                else
                    thread_restarts_schedule = make_unique<SyncedRestartSchedule>(restart_synchroniser);

                while (true) {
                    ++number_of_restarts;

                    if (! just_the_first_thread) {
                        wait_for_new_nogoods_barrier.wait();

                        for (unsigned u = 0 ; u < n_threads ; ++u)
                            if (t != u)
                                searchers[t]->watches.gather_nogoods_from(searchers[u]->watches);

                        // start watching new nogoods
                        if (searchers[t]->watches.apply_new_nogoods(
                                [&] (const Assignment & assignment) {
                                    for (auto & d : domains)
                                        if (d.v == assignment.pattern_vertex) {
                                            d.values.reset(assignment.target_vertex);
                                            d.count = d.values.count();
                                            break;
                                        }
                                }))
                            break;

                        if (0 == t)
                            restart_synchroniser.store(false);

                        synced_nogoods_barrier.wait();

                        searchers[t]->watches.clear_new_nogoods();
                    }

                    ++thread_result.propagations;
                    if (searchers[t]->propagate(domains, thread_assignments)) {
                        auto assignments_copy = thread_assignments;

                        switch (searchers[t]->restarting_search(assignments_copy, domains, thread_result.nodes, thread_result.propagations,
                                    thread_result.solution_count, 0, *thread_restarts_schedule)) {
                            case SearchResult::Satisfiable:
                                searchers[t]->save_result(assignments_copy, thread_result);
                                thread_result.complete = true;
                                params.timeout->trigger_early_abort();
                                searchers[t]->watches.post_nogood(Nogood<Assignment>{ });
                                break;

                            case SearchResult::SatisfiableButKeepGoing:
                                thread_result.complete = true;
                                params.timeout->trigger_early_abort();
                                searchers[t]->watches.post_nogood(Nogood<Assignment>{ });
                                break;

                            case SearchResult::Unsatisfiable:
                                thread_result.complete = true;
                                params.timeout->trigger_early_abort();
                                searchers[t]->watches.post_nogood(Nogood<Assignment>{ });
                                break;

                            case SearchResult::Aborted:
                                searchers[t]->watches.post_nogood(Nogood<Assignment>{ });
                                break;

                            case SearchResult::Restart:
                                break;
                        }
                    }
                    else {
                        thread_result.complete = true;
                        searchers[t]->watches.post_nogood(Nogood<Assignment>{ });
                        params.timeout->trigger_early_abort();
                    }

                    if (0 == t)
                        restart_synchroniser.store(true);
                    thread_restarts_schedule->did_a_restart();

                    if (params.delay_thread_creation && just_the_first_thread) {
                        if (! thread_result.complete) {
                            just_the_first_thread = false;
                            for (unsigned u = 1 ; u < n_threads ; ++u)
                                threads.emplace_back([&, u] () { work_function(u); });
                        }
                        else
                            break;
                    }
                }

                if (params.delay_thread_creation && 0 == t)
                    for (auto & th : threads)
                        th.join();

                unique_lock<mutex> lock{ common_result_mutex };
                if (! thread_result.mapping.empty())
                    common_result.mapping = move(thread_result.mapping);
                common_result.nodes += thread_result.nodes;
                common_result.propagations += thread_result.propagations;
                common_result.solution_count += thread_result.solution_count;
                common_result.complete = common_result.complete || thread_result.complete;
                for (auto & x : thread_result.extra_stats)
                    common_result.extra_stats.push_back("t" + to_string(t) + "_" + x);

                by_thread_nodes.append(" " + to_string(thread_result.nodes));
                by_thread_propagations.append(" " + to_string(thread_result.propagations));
            };

            if (params.delay_thread_creation)
                work_function(0);
            else {
                for (unsigned u = 0 ; u < n_threads ; ++u)
                    threads.emplace_back([&, u] () { work_function(u); });

                for (auto & th : threads)
                    th.join();
            }

            common_result.extra_stats.emplace_back("by_thread_nodes =" + by_thread_nodes);
            common_result.extra_stats.emplace_back("by_thread_propagations =" + by_thread_propagations);
            common_result.extra_stats.emplace_back("search_time = " + to_string(
                        duration_cast<milliseconds>(steady_clock::now() - search_start_time).count()));

            if (can_strip_isolated_vertices(params) && ! model.isolated_vertices.empty())
                common_result.extra_stats.emplace_back("isolated_vertices_removed = " + to_string(model.isolated_vertices.size()));

            return common_result;
        }
    };

    template <typename BitSetType_, typename ArrayType_, typename PatternAdjacencyBitsType_>
    struct SubgraphRunner
    {
        using Model = SubgraphModel<BitSetType_, ArrayType_, PatternAdjacencyBitsType_>;

        Model model;
        const HomomorphismParams & params;

        SubgraphRunner(const InputGraph & target, const InputGraph & pattern, const HomomorphismParams & p, const set<int> & x) :
            model(target, pattern, p, x),
            params(p)
        {
        }

        auto run() -> HomomorphismResult
        {
            // quick check for size
            if (is_nonshrinking(params) && (model.full_pattern_size > model.target_size)) {
                HomomorphismResult result;
                result.extra_stats.emplace_back("nonshrinking = false");
                if (params.proof)
                    params.proof->failure_due_to_pattern_bigger_than_target();
                return result;
            }

            if (! model.prepare(params)) {
                HomomorphismResult result;
                result.extra_stats.emplace_back("model_consistent = false");
                result.complete = true;
                return result;
            }

            HomomorphismResult result;
            if (1 == params.n_threads) {
                SequentialSolver<BitSetType_, ArrayType_, PatternAdjacencyBitsType_> solver(model, params);
                result = solver.solve();
            }
            else {
                if (! params.restarts_schedule->might_restart())
                    throw UnsupportedConfiguration{ "Threaded search requires restarts" };

                unsigned n_threads = how_many_threads(params.n_threads);
                ThreadedSolver<BitSetType_, ArrayType_, PatternAdjacencyBitsType_> solver(model, params, n_threads);
                result = solver.solve();
            }

            return result;
        }
    };

    // here be witchcraft. figure out the appropriate fixed size data
    // structures for whatever size of graph and combination of parameters we
    // actually have, and run that
    template <
        template <typename, typename, typename> class Algorithm_,
        typename Result_,
        typename Graph_,
        unsigned size_,
        unsigned... other_sizes_,
        typename... Params_>
    auto run_with_appropriate_template_parameters(
            const std::integer_sequence<unsigned, size_, other_sizes_...> &,
            unsigned n_shape_graphs,
            const Graph_ & graph,
            Params_ && ... params) -> Result_
    {
        if (graph.size() < int(size_ * bits_per_word)) {
            if (n_shape_graphs <= 8) {
                Algorithm_<FixedBitSet<size_>, std::array<int, size_ * bits_per_word + 1>, uint8_t> algorithm{ graph, std::forward<Params_>(params)... };
                return algorithm.run();
            }
            else {
                Algorithm_<FixedBitSet<size_>, std::array<int, size_ * bits_per_word + 1>, uint16_t> algorithm{ graph, std::forward<Params_>(params)... };
                return algorithm.run();
            }
        }
        else {
            if constexpr (0 == sizeof...(other_sizes_)) {
                if (n_shape_graphs <= 8) {
                    Algorithm_<boost::dynamic_bitset<>, std::vector<int>, uint8_t> algorithm{ graph, std::forward<Params_>(params)... };
                    return algorithm.run();
                }
                else {
                    Algorithm_<boost::dynamic_bitset<>, std::vector<int>, uint16_t> algorithm{ graph, std::forward<Params_>(params)... };
                    return algorithm.run();
                }
            }
            else
                return run_with_appropriate_template_parameters<Algorithm_, Result_, Graph_>(std::integer_sequence<unsigned, other_sizes_...>{},
                        n_shape_graphs, graph, std::forward<Params_>(params)...);
        }
    }

    using AllGraphSizes = std::integer_sequence<unsigned, 1, 2>;
}

auto solve_homomorphism_problem(const pair<InputGraph, InputGraph> & graphs, const HomomorphismParams & params) -> HomomorphismResult
{
    // start by setting up proof logging, if necessary
    if (params.proof) {
        // proof logging is currently incompatible with a whole load of "extra" features,
        // but can be adapted to support most of them
        if (1 != params.n_threads)
            throw UnsupportedConfiguration{ "Proof logging cannot yet be used with threads" };
        if (! params.no_nds)
            throw UnsupportedConfiguration{ "Proof logging cannot yet be used with neighbourhood degree sequences" };
        if (! params.no_supplementals)
            throw UnsupportedConfiguration{ "Proof logging cannot yet be used with supplemental graphs" };
        if (params.clique_detection)
            throw UnsupportedConfiguration{ "Proof logging cannot yet be used with clique detection" };
        if (params.lackey)
            throw UnsupportedConfiguration{ "Proof logging cannot yet be used with a lackey" };
        if (! params.pattern_less_constraints.empty())
            throw UnsupportedConfiguration{ "Proof logging cannot yet be used with less-constraints" };
        if (params.injectivity != Injectivity::Injective)
            throw UnsupportedConfiguration{ "Proof logging can currently only be used with injectivity" };
        if (params.induced)
            throw UnsupportedConfiguration{ "Proof logging cannot yet be used for induced problems" };
        if (params.minimal_unsat_pattern)
            throw UnsupportedConfiguration{ "Proof logging cannot yet be used for finding a minimal unsat pattern" };
        if (graphs.first.has_vertex_labels() || graphs.first.has_edge_labels())
            throw UnsupportedConfiguration{ "Proof logging cannot yet be used on labelled graphs" };

        // set up our model file, with a set of OPB variables for each CP variable
        for (int n = 0 ; n < graphs.first.size() ; ++n) {
            params.proof->create_cp_variable(n, graphs.second.size(),
                    [&] (int v) { return graphs.first.vertex_name(v); },
                    [&] (int v) { return graphs.second.vertex_name(v); });
        }

        // generate constraints for injectivity
        params.proof->create_injectivity_constraints(graphs.first.size(), graphs.second.size());

        // generate edge constraints, and also handle loops here
        for (int p = 0 ; p < graphs.first.size() ; ++p) {
            for (int t = 0 ; t < graphs.second.size() ; ++t) {
                if (graphs.first.adjacent(p, p) && ! graphs.second.adjacent(t, t))
                    params.proof->create_forbidden_assignment_constraint(p, t);
                else {
                    params.proof->start_adjacency_constraints_for(p, t);
                    // if p can be mapped to t, then each neighbour of p...
                    for (int q = 0 ; q < graphs.first.size() ; ++q)
                        if (q != p && graphs.first.adjacent(p, q)) {
                            // ... must be mapped to a neighbour of t
                            vector<int> permitted;
                            for (int u = 0 ; u < graphs.second.size() ; ++u)
                                if (t != u && graphs.second.adjacent(t, u))
                                    permitted.push_back(u);
                            params.proof->create_adjacency_constraint(p, q, t, permitted);
                        }
                }
            }
        }

        // output the model file
        params.proof->finalise_model();
    }

    // first sanity check: if we're finding an injective mapping, and there
    // aren't enough vertices, fail immediately.
    if (is_nonshrinking(params) && (graphs.first.size() > graphs.second.size())) {
        if (params.proof) {
            params.proof->failure_due_to_pattern_bigger_than_target();
            params.proof->finish_unsat_proof();
        }

        return HomomorphismResult{ };
    }

    // is the pattern a clique? if so, use a clique algorithm instead
    if (can_use_clique(params) && is_simple_clique(graphs.first)) {
        CliqueParams clique_params;
        clique_params.timeout = params.timeout;
        clique_params.start_time = params.start_time;
        clique_params.decide = make_optional(graphs.first.size());
        clique_params.restarts_schedule = make_unique<NoRestartsSchedule>();
        auto clique_result = solve_clique_problem(graphs.second, clique_params);

        // now translate the result back into what we expect
        HomomorphismResult result;
        int v = 0;
        for (auto & m : clique_result.clique) {
            result.mapping.emplace(v++, m);
            // the clique solver can find a bigger clique than we ask for
            if (v >= graphs.first.size())
                break;
        }
        result.nodes = clique_result.nodes;
        result.extra_stats = move(clique_result.extra_stats);
        result.extra_stats.emplace_back("used_clique_solver = true");
        result.complete = clique_result.complete;

        return result;
    }
    else if (params.minimal_unsat_pattern) {
        // if we're finding a minimal unsat pattern, solve the problem
        // repeatedly with shrinking patterns
        HomomorphismResult result = run_with_appropriate_template_parameters<SubgraphRunner, HomomorphismResult>(
                AllGraphSizes(), calculate_n_shape_graphs(params), graphs.second, graphs.first, params, set<int>{});
        if (! result.mapping.empty())
            return result;

        // try to exclude each vertex in turn
        set<int> exclude;
        for (int n = 0 ; n < graphs.first.size() ; ++n) {
            exclude.insert(n);
            result = run_with_appropriate_template_parameters<SubgraphRunner, HomomorphismResult>(
                    AllGraphSizes(), calculate_n_shape_graphs(params), graphs.second, graphs.first, params, exclude);
            if (! result.mapping.empty()) {
                result.mapping.clear();
                exclude.erase(n);
            }
        }

        // build up the minimal unsat pattern
        for (int n = 0 ; n < graphs.first.size() ; ++n) {
            if (! exclude.count(n))
                result.minimal_unsat_pattern.push_back(n);
        }

        return result;
    }
    else {
        // just solve the problem
        auto result = run_with_appropriate_template_parameters<SubgraphRunner, HomomorphismResult>(
                AllGraphSizes(), calculate_n_shape_graphs(params), graphs.second, graphs.first, params, set<int>{});

        if (params.proof && result.mapping.empty())
            params.proof->finish_unsat_proof();

        return result;
    }
}

