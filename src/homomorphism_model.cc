/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "homomorphism_model.hh"
#include "homomorphism_traits.hh"
#include "configuration.hh"

#include <algorithm>
#include <functional>
#include <list>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <string>
#include <utility>
#include <vector>

using std::accumulate;
using std::greater;
using std::list;
using std::map;
using std::max;
using std::optional;
using std::pair;
using std::queue;
using std::set;
using std::string;
using std::string_view;
using std::to_string;
using std::vector;
using std::cout;

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

struct HomomorphismModel::Imp
{
    const HomomorphismParams & params;

    vector<PatternAdjacencyBitsType> pattern_adjacencies_bits;
    vector<SVOBitset> pattern_graph_rows;
    vector<SVOBitset> target_graph_rows, forward_target_graph_rows, reverse_target_graph_rows;

    vector<vector<int> > patterns_degrees, targets_degrees;
    int largest_target_degree = 0;
    bool has_less_thans = false, directed = false;

    vector<int> pattern_vertex_labels, target_vertex_labels, pattern_edge_labels, target_edge_labels;
    vector<int> pattern_loops, target_loops;

    vector<string> pattern_vertex_proof_names, target_vertex_proof_names;

    vector<pair<int, int> > pattern_dir_degrees, target_dir_degrees;
    vector<pair<bool, bool> > pattern_big_constraints;

    vector<SVOBitset> target_graph_reachability, pattern_graph_reachability;
    vector<SVOBitset> pattern_site_reachability, pattern_root_reachability;

    set<int> pattern_anchor_nodes;

    Imp(const HomomorphismParams & p) :
        params(p)
    {
    }
};

HomomorphismModel::HomomorphismModel(const InputGraph & target, const InputGraph & pattern, const HomomorphismParams & params) :
    _imp(new Imp(params)),
    max_graphs(calculate_n_shape_graphs(params)),
    pattern_size(pattern.size()),
    pattern_link_count(pattern.get_no_link_nodes()),
    target_size(target.size()),
    target_link_count(target.get_no_link_nodes())
{
    _imp->patterns_degrees.resize(max_graphs);
    _imp->targets_degrees.resize(max_graphs);

    if (max_graphs > 8 * sizeof(PatternAdjacencyBitsType))
        throw UnsupportedConfiguration{ "Supplemental graphs won't fit in the chosen bitset size" };

    for (int v = 0 ; v < pattern.size() ; ++v) {
        _imp->pattern_vertex_proof_names.push_back(pattern.vertex_name(v));
        if(pattern.vertex_label(v) == "ANCHOR") _imp->pattern_anchor_nodes.insert(v);
    }

    for (int v = 0 ; v < target.size() ; ++v)
        _imp->target_vertex_proof_names.push_back(target.vertex_name(v));

    if (pattern.directed())
        _imp->directed = true;

    // recode pattern to a bit graph, and strip out loops
    _imp->pattern_graph_rows.resize(pattern_size * max_graphs, SVOBitset(pattern_size, 0));
    _imp->pattern_loops.resize(pattern_size);
    for (auto e = pattern.begin_edges(), e_end = pattern.end_edges(); e != e_end; ++e) {
        auto i = e->first.first;
        auto j = e->first.second;
        if (i == j)
            _imp->pattern_loops[i] = 1;
        else
            _imp->pattern_graph_rows[i * max_graphs + 0].set(j);
    }

    // re-encode and store pattern labels
    map<string, int> vertex_labels_map;
    int next_vertex_label = 1;
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
    int next_edge_label = 1;
    if (pattern.has_edge_labels()) {
        _imp->pattern_edge_labels.resize(pattern_size * pattern_size);
        for (unsigned i = 0 ; i < pattern_size ; ++i)
            for (unsigned j = 0 ; j < pattern_size ; ++j)
                if (_imp->pattern_graph_rows[i * max_graphs + 0].test(j)) {
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

    // if directed, do both directions
    if (pattern.directed()) {
        _imp->forward_target_graph_rows.resize(target_size, SVOBitset{ target_size, 0 });
        _imp->reverse_target_graph_rows.resize(target_size, SVOBitset{ target_size, 0 });
        for (auto e = target.begin_edges(), e_end = target.end_edges() ; e != e_end ; ++e) {
            if (e->first.first != e->first.second) {
                _imp->forward_target_graph_rows[e->first.first].set(e->first.second);
                _imp->reverse_target_graph_rows[e->first.second].set(e->first.first);
            }
        }
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

    // Symmetry break bigraph link nodes
    list<pair<unsigned, unsigned> > pattern_less_thans_in_wrong_order;
    if(params.bigraph && ! params.use_bigraph_projection_nogoods) {
        for (unsigned a = 0 ; a != pattern_size-pattern_link_count ; a++) {
            unsigned int prev_no_children = 999999999;
            unsigned int prev_same_clique = 999999999;
            for (unsigned b = pattern_size-pattern_link_count ; b != pattern_size ; b++) {
                if(_imp->pattern_graph_rows[a * max_graphs + 0].test(b)) {
                    if(_imp->pattern_graph_rows[b * max_graphs + 0].count() == 0) {
                        if(prev_no_children != 999999999) {
                            pattern_less_thans_in_wrong_order.push_back(std::make_pair(prev_no_children, b));
                        }
                        prev_no_children = b;
                    }
                    else {
                        if (prev_same_clique != 999999999 && _imp->pattern_graph_rows[prev_same_clique * max_graphs + 0].test(b)) {
                            pattern_less_thans_in_wrong_order.push_back(std::make_pair(prev_same_clique, b));
                        }
                        prev_same_clique = b;
                    }
                } 
            } 
        }
    }

    // pattern less than constraints
    if (! _imp->params.pattern_less_constraints.empty() || ! pattern_less_thans_in_wrong_order.empty()) {
        _imp->has_less_thans = true;
        
        for (auto & [ a, b ] : _imp->params.pattern_less_constraints) {
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

        _imp->pattern_dir_degrees.resize(pattern_size);
        _imp->pattern_big_constraints.resize(pattern_size-pattern_link_count);
        _imp->target_dir_degrees.resize(target_size);
    
        _imp->pattern_graph_reachability.resize(pattern_size-pattern_link_count, SVOBitset{ pattern_size-pattern_link_count, 0 });
        _imp->target_graph_reachability.resize(target_size-target_link_count, SVOBitset{ target_size-target_link_count, 0 });

        int site_size = 0;
        for (int a = 0 ; a != pattern.no_pattern_site_edges() ; ++a)
            if (pattern.get_pattern_site_edge(a).first + 1 > site_size)
                site_size = pattern.get_pattern_site_edge(a).first + 1;
        _imp->pattern_site_reachability.resize(site_size, SVOBitset{ pattern_size-pattern_link_count, 0 });

        for (int a = 0 ; a != pattern.no_pattern_site_edges() ; ++a)
            _imp->pattern_site_reachability[pattern.get_pattern_site_edge(a).first].set(pattern.get_pattern_site_edge(a).second);

        int root_size = 0;
        for (int a = 0 ; a != pattern.no_pattern_root_edges() ; ++a)
            if (pattern.get_pattern_root_edge(a).first + 1 > root_size)
                root_size = pattern.get_pattern_root_edge(a).first + 1;
        _imp->pattern_root_reachability.resize(root_size, SVOBitset{ pattern_size-pattern_link_count, 0 });

        for (int a = 0 ; a != pattern.no_pattern_root_edges() ; ++a)
            _imp->pattern_root_reachability[pattern.get_pattern_root_edge(a).first].set(pattern.get_pattern_root_edge(a).second);

        set<int> pattern_unique;
        queue<int> pattern_reach;

        set<int> target_unique;
        queue<int> target_reach;

        // Set in-out degrees of each vertex for bigraph root and site constraints
        for (unsigned a = 0 ; a < pattern_size ; ++a) {
            _imp->pattern_dir_degrees[a].first = pattern.in_degree(a);
            _imp->pattern_dir_degrees[a].second = pattern.out_degree(a);
        }
        for (unsigned a = 0 ; a < pattern_size-pattern_link_count ; ++a) {
            _imp->pattern_big_constraints[a] = pattern.get_big_constraint(a);
            _imp->pattern_graph_reachability[a].set(a);
            if (pattern.out_degree(a) == 0)
                pattern_reach.push(a);
        }

        for (unsigned a = 0 ; a != target_size ; ++a) {
            _imp->target_dir_degrees[a].first = target.in_degree(a);
            _imp->target_dir_degrees[a].second = target.out_degree(a);
        }

        for (unsigned a = 0 ; a != target_size-target_link_count ; ++a) {
            _imp->target_graph_reachability[a].set(a);
            if (target.out_degree(a) == 0)
                target_reach.push(a);
        }

        // Build reachability matrices by bubbling up from sink nodes and performing BFS
        while (! pattern_reach.empty()) {
            int v = pattern_reach.front();
            pattern_reach.pop();
            pattern_unique.erase(v);
            for (unsigned a = 0 ; a != pattern_size-pattern_link_count ; ++a)
                if (_imp->pattern_graph_rows[a * max_graphs + 0].test(v)) {
                    _imp->pattern_graph_reachability[a] |= _imp->pattern_graph_reachability[v];
                    if (pattern_unique.insert(a).second)
                        pattern_reach.push(a);
                }
        }

        while (! target_reach.empty()) {
            int v = target_reach.front();
            target_reach.pop();
            target_unique.erase(v);
            for (unsigned a = 0 ; a != target_size-target_link_count ; ++a)
                if (_imp->target_graph_rows[a * max_graphs + 0].test(v)) {
                    _imp->target_graph_reachability[a] |= _imp->target_graph_reachability[v];
                    if (target_unique.insert(a).second)
                        target_reach.push(a);
                }
        }
    }
}

HomomorphismModel::~HomomorphismModel() = default;

auto HomomorphismModel::_check_label_compatibility(int p, int t) const -> bool
{
    if (! has_vertex_labels())
        return true;
    else
        return pattern_vertex_label(p) == target_vertex_label(t);
}

auto HomomorphismModel::_check_loop_compatibility(int p, int t) const -> bool
{
    if (pattern_has_loop(p) && ! target_has_loop(t))
        return false;
    else if (_imp->params.induced && (pattern_has_loop(p) != target_has_loop(t)))
        return false;

    return true;
}

auto HomomorphismModel::_check_degree_compatibility(
        int p,
        int t,
        unsigned graphs_to_consider,
        vector<vector<vector<int> > > & patterns_ndss,
        vector<vector<optional<vector<int> > > > & targets_ndss,
        bool do_not_do_nds_yet
        ) const -> bool
{
    if (! degree_and_nds_are_preserved(_imp->params))
        return true;

    for (unsigned g = 0 ; g < graphs_to_consider ; ++g) {
        if (target_degree(g, t) < pattern_degree(g, p)) {
            // not ok, degrees differ
            if (_imp->params.proof) {
                // get the actual neighbours of p and t, in their original terms
                vector<int> n_p, n_t;

                auto np = pattern_graph_row(g, p);
                for (unsigned j = 0 ; j < pattern_size ; ++j)
                    if (np.test(j))
                        n_p.push_back(j);

                auto nt = target_graph_row(g, t);
                for (auto j = nt.find_first() ; j != decltype(nt)::npos ; j = nt.find_first()) {
                    nt.reset(j);
                    n_t.push_back(j);
                }

                _imp->params.proof->incompatible_by_degrees(g, pattern_vertex_for_proof(p), n_p,
                        target_vertex_for_proof(t), n_t);
            }
            return false;
        }
        else if (degree_and_nds_are_exact(_imp->params, pattern_size, target_size)
                && target_degree(g, t) != pattern_degree(g, p)) {
            // not ok, degrees must be exactly the same
            return false;
        }
    }
    if (_imp->params.no_nds || do_not_do_nds_yet)
        return true;

    // full compare of neighbourhood degree sequences
    if (! targets_ndss.at(0).at(t)) {
        for (unsigned g = 0 ; g < graphs_to_consider ; ++g) {
            targets_ndss.at(g).at(t) = vector<int>{};
            auto ni = target_graph_row(g, t);
            for (auto j = ni.find_first() ; j != decltype(ni)::npos ; j = ni.find_first()) {
                ni.reset(j);
                targets_ndss.at(g).at(t)->push_back(target_degree(g, j));
            }
            sort(targets_ndss.at(g).at(t)->begin(), targets_ndss.at(g).at(t)->end(), greater<int>());
        }
    }

    for (unsigned g = 0 ; g < graphs_to_consider ; ++g) {
        for (unsigned x = 0 ; x < patterns_ndss.at(g).at(p).size() ; ++x) {
            if (targets_ndss.at(g).at(t)->at(x) < patterns_ndss.at(g).at(p).at(x)) {
                if (_imp->params.proof) {
                    vector<int> p_subsequence, t_subsequence, t_remaining;

                    // need to know the NDS together with the actual vertices
                    vector<pair<int, int> > p_nds, t_nds;

                    auto np = pattern_graph_row(g, p);
                    for (auto w = np.find_first() ; w != decltype(np)::npos ; w = np.find_first()) {
                        np.reset(w);
                        p_nds.emplace_back(w, pattern_graph_row(g, w).count());
                    }

                    auto nt = target_graph_row(g, t);
                    for (auto w = nt.find_first() ; w != decltype(nt)::npos ; w = nt.find_first()) {
                        nt.reset(w);
                        t_nds.emplace_back(w, target_graph_row(g, w).count());
                    }

                    sort(p_nds.begin(), p_nds.end(), [] (const pair<int, int> & a, const pair<int, int> & b) {
                            return a.second > b.second; });
                    sort(t_nds.begin(), t_nds.end(), [] (const pair<int, int> & a, const pair<int, int> & b) {
                            return a.second > b.second; });

                    for (unsigned y = 0 ; y <= x ; ++y) {
                        p_subsequence.push_back(p_nds[y].first);
                        t_subsequence.push_back(t_nds[y].first);
                    }
                    for (unsigned y = x + 1 ; y < t_nds.size() ; ++y)
                        t_remaining.push_back(t_nds[y].first);

                    _imp->params.proof->incompatible_by_nds(g, pattern_vertex_for_proof(p),
                            target_vertex_for_proof(t), p_subsequence, t_subsequence, t_remaining);
                }
                return false;
            }
            else if (degree_and_nds_are_exact(_imp->params, pattern_size, target_size)
                    && targets_ndss.at(g).at(t)->at(x) != patterns_ndss.at(g).at(p).at(x))
                return false;
        }
    }

    return true;
}

auto HomomorphismModel::_check_open_directed_hyperedges(int p, int t) const -> bool
{
    string name1 = _imp->pattern_vertex_proof_names[p];
    string name2 = _imp->target_vertex_proof_names[t];

    if ((name1.find(":OPX:") == string::npos) || (name2.find(":OPX:") == string::npos))
        return true;

    return name1.substr(name1.length()-1, name1.length()) == name2.substr(name2.length()-1, name2.length());
}

auto HomomorphismModel::_check_bigraph_equality_regions(int p, int t) const -> bool
{
    return ((_imp->pattern_vertex_proof_names[p].find("ROOT") == string::npos) ||
        (_imp->target_vertex_proof_names[t].find("ROOT") == string::npos) ||
        (_imp->pattern_vertex_proof_names[p] == _imp->target_vertex_proof_names[t]));
}

auto HomomorphismModel::_check_bigraph_equality_links(int p, int t) const -> bool
{
    string name1 = _imp->pattern_vertex_proof_names[p];
    string name2 = _imp->target_vertex_proof_names[t];
    
    if ((name1.find(":OPX:") == string::npos) || (name2.find(":OPX:") == string::npos))
        return true;

    name1 = name1.substr(5, name1.length()-1);
    name1 = name1.substr(0, name1.find(':'));
    
    name2 = name2.substr(5, name2.length()-1);
    name2 = name2.substr(0, name2.find(':'));
    
    return name1 == name2;
}

auto HomomorphismModel::_check_bigraph_degree_compatibility(int p, int t) const -> bool
{
    // Cannot match a place node with a link node
    if(p >= pattern_size-pattern_link_count && t < target_size-target_link_count)
        return false;
    if(p < pattern_size-pattern_link_count && t >= target_size-target_link_count)
        return false;

    // Closed hyperedges need same degree to match
    if(_imp->pattern_vertex_proof_names[p].find("C_LINK") != string::npos && _imp->target_vertex_proof_names[t].find("C_LINK") != string::npos)
        if(_imp->pattern_dir_degrees[p].first != _imp->target_dir_degrees[t].first)
            return false;

    if(p >= pattern_size-pattern_link_count && t >= target_size-target_link_count) 
        return true;

    // Place matching degree constraints
    if ((! _imp->pattern_big_constraints[p].first) && (_imp->pattern_dir_degrees[p].first != _imp->target_dir_degrees[t].first))
        return false;
    if (_imp->pattern_big_constraints[p].first && (_imp->pattern_dir_degrees[p].first > _imp->target_dir_degrees[t].first))
        return false;
    if ((! _imp->pattern_big_constraints[p].second) && (_imp->pattern_dir_degrees[p].second != _imp->target_dir_degrees[t].second))
        return false;
    if (_imp->pattern_big_constraints[p].second && (_imp->pattern_dir_degrees[p].second > _imp->target_dir_degrees[t].second))
        return false;
    return true;
}

auto HomomorphismModel::initialise_domains(vector<HomomorphismDomain> & domains) const -> bool
{
    unsigned graphs_to_consider = max_graphs;

    /* pattern and target neighbourhood degree sequences */
    vector<vector<vector<int> > > patterns_ndss(graphs_to_consider);
    vector<vector<optional<vector<int> > > > targets_ndss(graphs_to_consider);

    if (degree_and_nds_are_preserved(_imp->params) && ! _imp->params.no_nds) {
        for (unsigned g = 0 ; g < graphs_to_consider ; ++g) {
            patterns_ndss.at(g).resize(pattern_size);
            targets_ndss.at(g).resize(target_size);
        }

        for (unsigned g = 0 ; g < graphs_to_consider ; ++g) {
            for (unsigned i = 0 ; i < pattern_size ; ++i) {
                auto ni = pattern_graph_row(g, i);
                for (auto j = ni.find_first() ; j != decltype(ni)::npos ; j = ni.find_first()) {
                    ni.reset(j);
                    patterns_ndss.at(g).at(i).push_back(pattern_degree(g, j));
                }
                sort(patterns_ndss.at(g).at(i).begin(), patterns_ndss.at(g).at(i).end(), greater<int>());
            }
        }
    }

    for (unsigned i = 0 ; i < pattern_size ; ++i) {
        domains.at(i).v = i;
        domains.at(i).values.reset();

        for (unsigned j = 0 ; j < target_size ; ++j) {
            bool ok = true;

            if (! _check_label_compatibility(i, j))
                ok = false;
            else if (! _check_loop_compatibility(i, j))
                ok = false;
            else if (! _check_degree_compatibility(i, j, graphs_to_consider, patterns_ndss, targets_ndss, _imp->params.proof.get()))
                ok = false;
            else if (_imp->params.bigraph && ! _check_bigraph_degree_compatibility(i, j))
                ok = false;
            else if (_imp->params.bigraph && _imp->params.equality_check && ! _check_bigraph_equality_regions(i, j))
                ok = false;
            else if (_imp->params.bigraph && _imp->params.equality_check && ! _check_bigraph_equality_links(i, j))
                ok = false;
            else if (_imp->params.bigraph && _imp->params.directed && ! _check_open_directed_hyperedges(i, j))
                ok = false;
            if (ok)
                domains.at(i).values.set(j);
        }

        domains.at(i).count = domains.at(i).values.count();
        if (0 == domains.at(i).count)
            return false;
    }

    // for proof logging, we need degree information before we can output nds proofs
    if (_imp->params.proof && degree_and_nds_are_preserved(_imp->params) && ! _imp->params.no_nds) {
        for (unsigned i = 0 ; i < pattern_size ; ++i) {
            for (unsigned j = 0 ; j < target_size ; ++j) {
                if (domains.at(i).values.test(j) &&
                        ! _check_degree_compatibility(i, j, graphs_to_consider, patterns_ndss, targets_ndss, false)) {
                    domains.at(i).values.reset(j);
                    if (0 == --domains.at(i).count)
                        return false;
                }
            }
        }
    }

    // quick sanity check that we have enough values
    if (is_nonshrinking(_imp->params)) {
        SVOBitset domains_union{ target_size, 0 };
        for (auto & d : domains)
            domains_union |= d.values;

        unsigned domains_union_popcount = domains_union.count();
        if (domains_union_popcount < unsigned(pattern_size)) {
            if (_imp->params.proof) {
                vector<int> hall_lhs, hall_rhs;
                for (auto & d : domains)
                    hall_lhs.push_back(d.v);
                auto dd = domains_union;
                for (auto v = dd.find_first() ; v != decltype(dd)::npos ; v = dd.find_first()) {
                    dd.reset(v);
                    hall_rhs.push_back(v);
                }
                _imp->params.proof->emit_hall_set_or_violator(hall_lhs, hall_rhs);
            }
            return false;
        }
    }

    for (auto & d : domains) {
        d.count = d.values.count();
        if (0 == d.count && _imp->params.proof) {
            _imp->params.proof->initial_domain_is_empty(d.v);
            return false;
        }
    }

    return true;
}

auto HomomorphismModel::pattern_vertex_for_proof(int v) const -> NamedVertex
{
    if (v < 0 || unsigned(v) >= _imp->pattern_vertex_proof_names.size())
        throw ProofError{ "Oops, there's a bug: v out of range in pattern" };
    return pair{ v, _imp->pattern_vertex_proof_names[v] };
}

auto HomomorphismModel::target_vertex_for_proof(int v) const -> NamedVertex
{
    if (v < 0 || unsigned(v) >= _imp->target_vertex_proof_names.size())
        throw ProofError{ "Oops, there's a bug: v out of range in target" };
    return pair{ v, _imp->target_vertex_proof_names[v] };
}

auto HomomorphismModel::prepare() -> bool
{
    if (is_nonshrinking(_imp->params) && (pattern_size > target_size))
        return false;

    // pattern and target degrees, for the main graph
    _imp->patterns_degrees.at(0).resize(pattern_size);
    _imp->targets_degrees.at(0).resize(target_size);

    for (unsigned i = 0 ; i < pattern_size ; ++i)
        _imp->patterns_degrees.at(0).at(i) = _imp->pattern_graph_rows[i * max_graphs + 0].count();

    for (unsigned i = 0 ; i < target_size ; ++i)
        _imp->targets_degrees.at(0).at(i) = _imp->target_graph_rows[i * max_graphs + 0].count();

    if (global_degree_is_preserved(_imp->params)) {
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
                if (_imp->params.proof) {
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

                            _imp->params.proof->incompatible_by_degrees(0,
                                    pattern_vertex_for_proof(p_gds.at(p).first), n_p,
                                    target_vertex_for_proof(t_gds.at(t).first), n_t);
                        }
                    }

                    vector<int> patterns, targets;
                    for (unsigned p = 0 ; p <= i ; ++p)
                        patterns.push_back(p_gds.at(p).first);
                    for (unsigned t = 0 ; t < i ; ++t)
                        targets.push_back(t_gds.at(t).first);

                    _imp->params.proof->emit_hall_set_or_violator(patterns, targets);
                }
                return false;
            }
    }

    unsigned next_pattern_supplemental = 1, next_target_supplemental = 1;
    // build exact path graphs
    if (supports_exact_path_graphs(_imp->params)) {
        _build_exact_path_graphs(_imp->pattern_graph_rows, pattern_size, next_pattern_supplemental, _imp->params.number_of_exact_path_graphs, _imp->directed);
        _build_exact_path_graphs(_imp->target_graph_rows, target_size, next_target_supplemental, _imp->params.number_of_exact_path_graphs, _imp->directed);

        if (_imp->params.proof) {
            for (int g = 1 ; g <= _imp->params.number_of_exact_path_graphs ; ++g) {
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

                            _imp->params.proof->create_exact_path_graphs(g, named_p, named_q, between_p_and_q,
                                    named_t, named_n_t, named_two_away_from_t, named_d_n_t);
                        }
                    }
                }
            }
        }
    }

    if (supports_distance3_graphs(_imp->params)) {
        _build_distance3_graphs(_imp->pattern_graph_rows, pattern_size, next_pattern_supplemental);
        _build_distance3_graphs(_imp->target_graph_rows, target_size, next_target_supplemental);
    }

    if (supports_k4_graphs(_imp->params)) {
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

auto HomomorphismModel::_build_exact_path_graphs(vector<SVOBitset> & graph_rows, unsigned size, unsigned & idx,
        unsigned number_of_exact_path_graphs, bool directed) -> void
{
    vector<vector<unsigned> > path_counts(size, vector<unsigned>(size, 0));

    // count number of paths from w to v (unless directed, only w >= v, so not v to w)
    for (unsigned v = 0 ; v < size ; ++v) {
        auto nv = graph_rows[v * max_graphs + 0];
        for (auto c = nv.find_first() ; c != decltype(nv)::npos ; c = nv.find_first()) {
            nv.reset(c);
            auto nc = graph_rows[c * max_graphs + 0];
            for (auto w = nc.find_first() ; w != decltype(nc)::npos && (directed ? true : w <= v) ; w = nc.find_first()) {
                nc.reset(w);
                ++path_counts[v][w];
            }
        }
    }

    for (unsigned v = 0 ; v < size ; ++v) {
        for (unsigned w = (directed ? 0 : v) ; w < size ; ++w) {
            // nuless directed, w to v, not v to w, see above
            unsigned path_count = path_counts[w][v];
            for (unsigned p = 1 ; p <= number_of_exact_path_graphs ; ++p) {
                if (path_count >= p) {
                    graph_rows[v * max_graphs + idx + p - 1].set(w);
                    if (! directed)
                        graph_rows[w * max_graphs + idx + p - 1].set(v);
                }
            }
        }
    }

    idx += number_of_exact_path_graphs;
}

auto HomomorphismModel::_build_distance3_graphs(vector<SVOBitset> & graph_rows, unsigned size, unsigned & idx) -> void
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

auto HomomorphismModel::_build_k4_graphs(vector<SVOBitset> & graph_rows, unsigned size, unsigned & idx) -> void
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

auto HomomorphismModel::pattern_adjacency_bits(int p, int q) const -> PatternAdjacencyBitsType
{
    return _imp->pattern_adjacencies_bits[pattern_size * p + q];
}

auto HomomorphismModel::pattern_graph_row(int g, int p) const -> const SVOBitset &
{
    return _imp->pattern_graph_rows[p * max_graphs + g];
}

auto HomomorphismModel::target_graph_row(int g, int t) const -> const SVOBitset &
{
    return _imp->target_graph_rows[t * max_graphs + g];
}

auto HomomorphismModel::forward_target_graph_row(int t) const -> const SVOBitset &
{
    return _imp->forward_target_graph_rows[t];
}

auto HomomorphismModel::reverse_target_graph_row(int t) const -> const SVOBitset &
{
    return _imp->reverse_target_graph_rows[t];
}

auto HomomorphismModel::pattern_degree(int g, int p) const -> unsigned
{
    return _imp->patterns_degrees[g][p];
}

auto HomomorphismModel::target_degree(int g, int t) const -> unsigned
{
    return _imp->targets_degrees[g][t];
}

auto HomomorphismModel::largest_target_degree() const -> unsigned
{
    return _imp->largest_target_degree;
}

auto HomomorphismModel::has_vertex_labels() const -> bool
{
    return ! _imp->pattern_vertex_labels.empty();
}

auto HomomorphismModel::has_edge_labels() const -> bool
{
    return ! _imp->pattern_edge_labels.empty();
}

auto HomomorphismModel::pattern_vertex_label(int p) const -> int
{
    return _imp->pattern_vertex_labels[p];
}

auto HomomorphismModel::target_vertex_label(int t) const -> int
{
    return _imp->target_vertex_labels[t];
}

auto HomomorphismModel::pattern_edge_label(int p, int q) const -> int
{
    return _imp->pattern_edge_labels[p * pattern_size + q];
}

auto HomomorphismModel::target_edge_label(int t, int u) const -> int
{
    return _imp->target_edge_labels[t * target_size + u];
}

auto HomomorphismModel::pattern_has_loop(int p) const -> bool
{
    return _imp->pattern_loops[p];
}

auto HomomorphismModel::target_has_loop(int t) const -> bool
{
    return _imp->target_loops[t];
}

auto HomomorphismModel::has_less_thans() const -> bool
{
    return _imp->has_less_thans;
}

auto HomomorphismModel::directed() const -> bool
{
    return _imp->directed;
}

auto HomomorphismModel::is_pattern_anchor(int p) const -> bool
{
    return _imp->pattern_anchor_nodes.find(p) != _imp->pattern_anchor_nodes.end();
}

auto HomomorphismModel::check_extra_bigraph_constraints(const VertexToVertexMapping & mapping) const -> bool
{
    // Find transitive closure violations
    for (unsigned i = 0 ; i != pattern_size-pattern_link_count ; ++i) {
        if (_imp->pattern_big_constraints[i].second) {
            set<int> child_mappings;

            // Get all corresponding target node's children with a mapping
            for (unsigned j = 0 ; j != pattern_size-pattern_link_count ; ++j)
                if (i != j &&
                        _imp->pattern_graph_rows[i * max_graphs].test(j) &&
                        _imp->pattern_graph_reachability[i].test(j))
                    child_mappings.insert(mapping.find(j)->second);

            // Check if pattern node only points to sites that are shared
            bool sites_okay = false;
            vector<int> site_mappings;
            for (unsigned j = 0 ; j != _imp->pattern_site_reachability.size() ; ++j) {
                if (_imp->pattern_site_reachability[j].test(i) && _imp->pattern_site_reachability[j].count() == 1) {
                    sites_okay = true;
                    break;
                }
                if (_imp->pattern_site_reachability[j].test(i))
                    site_mappings.push_back(j);
            }

            // For all target node's children without a mapping, check if it can reach any children of a root node
            for (unsigned j = 0 ; j != target_size-target_link_count ; ++j) {
                if (mapping.find(i)->second != int(j) &&
                        _imp->target_graph_rows[mapping.find(i)->second * max_graphs].test(j) &&
                        _imp->target_graph_reachability[mapping.find(i)->second].test(j) &&
                        child_mappings.find(j) == child_mappings.end()) {
                    for (unsigned k = 0 ; k != pattern_size-pattern_link_count ; ++k)
                        if (_imp->pattern_big_constraints[k].first && _imp->target_graph_reachability[j].test(mapping.find(k)->second))
                            return false;

                    // If only points to shared sites, check that all unmapped children have an adjacency superset of a site
                    if (! sites_okay) {
                        bool site_sat = true;
                        for (unsigned k = 0 ; k != site_mappings.size() ; ++k) {
                            site_sat = true;
                            for (unsigned l = 0 ; l != pattern_size-pattern_link_count ; ++l) {
                                if (l != i && _imp->pattern_site_reachability[site_mappings[k]].test(l) &&
                                        ! (_imp->target_graph_rows[mapping.find(l)->second * max_graphs].test(j) &&
                                            _imp->target_graph_reachability[mapping.find(l)->second].test(j))) {
                                    site_sat = false;
                                    break;
                                }
                            }
                            if (site_sat)
                                break;
                        }

                        if (! site_sat)
                            return false;
                    }
                }
            }
        }

        // Check if pattern node only is pointed to by roots that are shared
        if (_imp->pattern_big_constraints[i].first) {
            vector<int> root_mappings;
            bool roots_okay = false;
            for (unsigned j = 0 ; j != _imp->pattern_root_reachability.size() ; ++j) {
                if (_imp->pattern_root_reachability[j].test(i) && _imp->pattern_root_reachability[j].count() == 1) {
                    roots_okay = true;
                    break;
                }
                if (_imp->pattern_root_reachability[j].test(i))
                    root_mappings.push_back(j);
            }

            if (! roots_okay) {
                set<int> parent_mappings;
                for (unsigned j = 0 ; j != pattern_size-pattern_link_count ; ++j)
                    if (i != j &&
                            _imp->pattern_graph_rows[i * max_graphs].test(j) &&
                            _imp->pattern_graph_reachability[j].test(i))
                        parent_mappings.insert(mapping.find(j)->second);

                // If only points to shared roots, check that all unmapped parents have an adjacency superset of a root
                for (unsigned j = 0 ; j != target_size-target_link_count ; ++j) {
                    if (mapping.find(i)->second != int(j) &&
                            _imp->target_graph_rows[j * max_graphs].test(mapping.find(i)->second) &&
                            _imp->target_graph_reachability[j].test(mapping.find(i)->second) &&
                            ! parent_mappings.count(j)) {
                        bool root_sat = true;
                        for (unsigned k = 0 ; k != root_mappings.size() ; ++k) {
                            root_sat = true;
                            for (unsigned l = 0 ; l < pattern_size-pattern_link_count ; ++l){
                                if (l != i && _imp->pattern_root_reachability[root_mappings[k]].test(l) &&
                                        ! (_imp->target_graph_rows[j * max_graphs].test(mapping.find(l)->second) &&
                                            _imp->target_graph_reachability[j].test(mapping.find(l)->second))) {
                                    root_sat = false;
                                    break;
                                }
                            }
                            if (root_sat)
                                break;
                        }
                        if (! root_sat)
                            return false;
                    }
                }
            }
        }
    }

    return true;
}

