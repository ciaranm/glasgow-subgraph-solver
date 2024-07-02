#include <gss/clique.hh>
#include <gss/configuration.hh>
#include <gss/innards/homomorphism_model.hh>
#include <gss/innards/homomorphism_traits.hh>

#include <chrono>
#include <functional>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace gss;
using namespace gss::innards;

using std::greater;
using std::list;
using std::make_optional;
using std::make_shared;
using std::make_unique;
using std::map;
using std::max;
using std::nullopt;
using std::optional;
using std::pair;
using std::set;
using std::shared_ptr;
using std::string;
using std::string_view;
using std::stringstream;
using std::to_string;
using std::vector;

using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::chrono::steady_clock;

namespace
{
    auto calculate_n_shape_graphs(const HomomorphismParams & params) -> unsigned
    {
        return 1 +
            (supports_exact_path_graphs(params) ? params.number_of_exact_path_graphs : 0) +
            (supports_distance2_graphs(params) ? 1 : 0) +
            (supports_distance3_graphs(params) ? 1 : 0) +
            (supports_k4_graphs(params) ? 1 : 0) +
            params.extra_shapes.size();
    }

    auto find_clique(
        const shared_ptr<Timeout> & timeout,
        unsigned size,
        const vector<SVOBitset> & rows,
        unsigned g,
        unsigned max_graphs,
        unsigned v,
        optional<int> largest_if_target,
        vector<int> & best_knowns,
        list<string> & build_times,
        list<string> & solve_times,
        list<string> & find_nodes,
        list<string> & prove_nodes) -> int
    {
        if (largest_if_target && (best_knowns[v] >= *largest_if_target))
            return best_knowns[v];

        auto start_time = steady_clock::now();

        CliqueParams params;
        params.timeout = timeout;
        params.start_time = steady_clock::now();
        params.decide = nullopt;
        params.restarts_schedule = make_unique<NoRestartsSchedule>();
        if (largest_if_target)
            params.stop_after_finding = *largest_if_target - 1;

        vector<int> include(size, -1), invinclude(size, 0);
        int count = 0;
        for (unsigned w = 0; w < size; ++w)
            if (w != v && rows[w * max_graphs + g].test(v)) {
                include[w] = count;
                invinclude[count] = w;
                ++count;
            }

        InputGraph gv(count, false, false);
        for (unsigned f = 0; f < size; ++f)
            if (include[f] != -1)
                for (unsigned t = 0; t < size; ++t)
                    if (f != t && include[t] != -1 && rows[f * max_graphs + g].test(t))
                        gv.add_edge(include[f], include[t]);

        build_times.push_back(to_string(duration_cast<milliseconds>(steady_clock::now() - start_time).count()));

        start_time = steady_clock::now();

        auto result = solve_clique_problem(gv, params);

        solve_times.push_back(to_string(duration_cast<milliseconds>(steady_clock::now() - start_time).count()));
        find_nodes.push_back(to_string(result.find_nodes));
        prove_nodes.push_back(to_string(result.prove_nodes));

        best_knowns[v] = max<int>(best_knowns[v], result.clique.size() + 1);
        for (auto & w : result.clique)
            best_knowns[invinclude[w]] = max<int>(best_knowns[invinclude[w]], result.clique.size() + 1);

        return result.clique.size() + 1;
    }
}

struct HomomorphismModel::Imp
{
    const HomomorphismParams & params;
    shared_ptr<Proof> proof;

    vector<PatternAdjacencyBitsType> pattern_adjacencies_bits;
    vector<SVOBitset> pattern_graph_rows;
    vector<SVOBitset> target_graph_rows, forward_target_graph_rows, reverse_target_graph_rows;

    vector<vector<int>> patterns_degrees, targets_degrees;
    int largest_target_degree = 0;
    bool has_less_thans = false, has_occur_less_thans = false, directed = false;

    vector<int> pattern_vertex_labels, target_vertex_labels, pattern_edge_labels, target_edge_labels;
    vector<int> pattern_loops, target_loops;

    vector<string> pattern_vertex_proof_names, target_vertex_proof_names;

    mutable bool has_pattern_cliques_sizes = false;
    mutable vector<vector<int>> pattern_cliques_sizes, target_cliques_sizes, pattern_cliques_best_knowns, target_cliques_best_knowns;
    mutable vector<int> largest_pattern_clique;

    unsigned max_graphs_for_clique_size_constraints = 0;
    mutable list<string> pattern_cliques_build_times, pattern_cliques_solve_times, pattern_cliques_solve_find_nodes, pattern_cliques_solve_prove_nodes;
    mutable list<string> target_cliques_build_times, target_cliques_solve_times, target_cliques_solve_find_nodes, target_cliques_solve_prove_nodes;

    mutable list<string> supplemental_graph_names;

    Imp(const HomomorphismParams & p, const std::shared_ptr<Proof> & r) :
        params(p),
        proof(r)
    {
    }
};

HomomorphismModel::HomomorphismModel(const InputGraph & target, const InputGraph & pattern, const HomomorphismParams & params,
    const std::shared_ptr<Proof> & proof) :
    _imp(new Imp(params, proof)),
    max_graphs(calculate_n_shape_graphs(params)),
    pattern_size(pattern.size()),
    target_size(target.size())
{
    if (_imp->params.clique_size_constraints)
        _imp->max_graphs_for_clique_size_constraints = (_imp->params.clique_size_constraints_on_supplementals ? max_graphs : 1);

    _imp->patterns_degrees.resize(max_graphs);
    _imp->targets_degrees.resize(max_graphs);

    if (max_graphs > 8 * sizeof(PatternAdjacencyBitsType))
        throw UnsupportedConfiguration{"Supplemental graphs won't fit in the chosen bitset size"};

    if (_imp->proof) {
        for (int v = 0; v < pattern.size(); ++v)
            _imp->pattern_vertex_proof_names.push_back(pattern.vertex_name(v));
        for (int v = 0; v < target.size(); ++v)
            _imp->target_vertex_proof_names.push_back(target.vertex_name(v));
    }

    if (pattern.directed())
        _imp->directed = true;

    // recode pattern to a bit graph, and strip out loops
    _imp->pattern_graph_rows.resize(pattern_size * max_graphs, SVOBitset(pattern_size, 0));
    _imp->pattern_loops.resize(pattern_size);
    for (unsigned i = 0; i < pattern_size; ++i) {
        for (unsigned j = 0; j < pattern_size; ++j) {
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
    int next_vertex_label = 1;
    if (pattern.has_vertex_labels()) {
        for (unsigned i = 0; i < pattern_size; ++i) {
            if (vertex_labels_map.emplace(pattern.vertex_label(i), next_vertex_label).second)
                ++next_vertex_label;
        }

        _imp->pattern_vertex_labels.resize(pattern_size);
        for (unsigned i = 0; i < pattern_size; ++i)
            _imp->pattern_vertex_labels[i] = vertex_labels_map.find(string{pattern.vertex_label(i)})->second;
    }

    // re-encode and store edge labels
    map<string, int> edge_labels_map;
    int next_edge_label = 1;
    if (pattern.has_edge_labels()) {
        _imp->pattern_edge_labels.resize(pattern_size * pattern_size);
        for (unsigned i = 0; i < pattern_size; ++i)
            for (unsigned j = 0; j < pattern_size; ++j)
                if (pattern.adjacent(i, j)) {
                    auto r = edge_labels_map.emplace(pattern.edge_label(i, j), next_edge_label);
                    if (r.second)
                        ++next_edge_label;
                    _imp->pattern_edge_labels[i * pattern_size + j] = r.first->second;
                }
    }

    // recode target to a bit graph, and take out loops
    _imp->target_graph_rows.resize(target_size * max_graphs, SVOBitset{target_size, 0});
    _imp->target_loops.resize(target_size);
    target.for_each_edge([&](int f, int t, string_view) {
        if (f == t)
            _imp->target_loops[f] = 1;
        else
            _imp->target_graph_rows[f * max_graphs + 0].set(t);
    });

    // if directed, do both directions
    if (pattern.directed()) {
        _imp->forward_target_graph_rows.resize(target_size, SVOBitset{target_size, 0});
        _imp->reverse_target_graph_rows.resize(target_size, SVOBitset{target_size, 0});
        target.for_each_edge([&](int f, int t, string_view l) {
            if (f != t && l != "unlabelled") {
                _imp->forward_target_graph_rows[f].set(t);
                _imp->reverse_target_graph_rows[t].set(f);
            }
        });
    }

    // target vertex labels
    if (pattern.has_vertex_labels()) {
        for (unsigned i = 0; i < target_size; ++i) {
            if (vertex_labels_map.emplace(target.vertex_label(i), next_vertex_label).second)
                ++next_vertex_label;
        }

        _imp->target_vertex_labels.resize(target_size);
        for (unsigned i = 0; i < target_size; ++i)
            _imp->target_vertex_labels[i] = vertex_labels_map.find(string{target.vertex_label(i)})->second;
    }

    // target edge labels
    if (pattern.has_edge_labels()) {
        _imp->target_edge_labels.resize(target_size * target_size);
        target.for_each_edge([&](int f, int t, string_view l) {
            auto r = edge_labels_map.emplace(l, next_edge_label);
            if (r.second)
                ++next_edge_label;

            _imp->target_edge_labels[f * target_size + t] = r.first->second;
        });
    }

    auto decode = [&](const InputGraph & g, string_view s) -> int {
        auto n = g.vertex_from_name(s);
        if (! n)
            throw UnsupportedConfiguration{"No vertex named '" + string{s} + "'"};
        return *n;
    };

    // pattern less than constraints
    if (! _imp->params.pattern_less_constraints.empty()) {
        _imp->has_less_thans = true;
        list<pair<unsigned, unsigned>> pattern_less_thans_in_wrong_order;
        for (auto & [a, b] : _imp->params.pattern_less_constraints) {
            auto a_decoded = decode(pattern, a), b_decoded = decode(pattern, b);
            pattern_less_thans_in_wrong_order.emplace_back(a_decoded, b_decoded);
        }

        // put them in a convenient order, so we don't need a propagation loop
        while (! pattern_less_thans_in_wrong_order.empty()) {
            bool loop_detect = true;
            set<unsigned> cannot_order_yet;
            for (auto & [_, b] : pattern_less_thans_in_wrong_order)
                cannot_order_yet.emplace(b);
            for (auto p = pattern_less_thans_in_wrong_order.begin(); p != pattern_less_thans_in_wrong_order.end();) {
                if (cannot_order_yet.count(p->first))
                    ++p;
                else {
                    loop_detect = false;
                    pattern_less_thans_in_convenient_order.push_back(*p);
                    pattern_less_thans_in_wrong_order.erase(p++);
                }
            }

            if (loop_detect)
                throw UnsupportedConfiguration{"Pattern less than constraints form a loop"};
        }
    }

    // target less than constraints
    if (! _imp->params.target_occur_less_constraints.empty()) {
        _imp->has_occur_less_thans = true;
        list<pair<unsigned, unsigned>> target_occur_less_thans_in_wrong_order;
        for (auto & [a, b] : _imp->params.target_occur_less_constraints) {
            auto a_decoded = decode(target, a), b_decoded = decode(target, b);
            target_occur_less_thans_in_wrong_order.emplace_back(a_decoded, b_decoded);
        }

        // put them in a convenient order, so we don't need a propagation loop
        while (! target_occur_less_thans_in_wrong_order.empty()) {
            bool loop_detect = true;
            set<unsigned> cannot_order_yet;
            for (auto & [_, b] : target_occur_less_thans_in_wrong_order)
                cannot_order_yet.emplace(b);
            for (auto t = target_occur_less_thans_in_wrong_order.begin(); t != target_occur_less_thans_in_wrong_order.end();) {
                if (cannot_order_yet.count(t->first))
                    ++t;
                else {
                    loop_detect = false;
                    target_occur_less_thans_in_convenient_order.push_back(*t);
                    target_occur_less_thans_in_wrong_order.erase(t++);
                }
            }

            if (loop_detect)
                throw UnsupportedConfiguration{"Target less than constraints form a loop"};
        }
    }

    // make space for clique constraints
    if (params.clique_size_constraints) {
        for (unsigned g = 0; g < _imp->max_graphs_for_clique_size_constraints; ++g) {
            _imp->pattern_cliques_sizes.push_back(vector<int>(pattern.size(), 0));
            _imp->target_cliques_sizes.push_back(vector<int>(target.size(), 0));
            _imp->pattern_cliques_best_knowns.push_back(vector<int>(pattern.size(), 0));
            _imp->target_cliques_best_knowns.push_back(vector<int>(target.size(), 0));
        }
        _imp->largest_pattern_clique.resize(_imp->max_graphs_for_clique_size_constraints);
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
    if (pattern_has_loop(p) && ! target_has_loop(t)) {
        if (_imp->proof)
            _imp->proof->incompatible_by_loops(pattern_vertex_for_proof(p), target_vertex_for_proof(t));
        return false;
    }
    else if (_imp->params.induced && (pattern_has_loop(p) != target_has_loop(t)))
        return false;

    return true;
}

auto HomomorphismModel::_build_pattern_clique_sizes() const -> void
{
    for (unsigned g = 0; g < _imp->max_graphs_for_clique_size_constraints; ++g) {
        for (unsigned v = 0; v < pattern_size; ++v) {
            auto c = find_clique(_imp->params.timeout, pattern_size, _imp->pattern_graph_rows, g, max_graphs, v, nullopt,
                _imp->pattern_cliques_best_knowns[g], _imp->pattern_cliques_build_times, _imp->pattern_cliques_solve_times,
                _imp->pattern_cliques_solve_find_nodes, _imp->pattern_cliques_solve_prove_nodes);
            _imp->pattern_cliques_sizes[g][v] = c;
            _imp->largest_pattern_clique[g] = max(_imp->largest_pattern_clique[g], c);
        }
        _imp->has_pattern_cliques_sizes = true;
    }
}

auto HomomorphismModel::_build_target_clique_size(int v) const -> void
{
    if (0 == _imp->target_cliques_sizes[0][v])
        for (unsigned g = 0; g < _imp->max_graphs_for_clique_size_constraints; ++g) {
            _imp->target_cliques_sizes[g][v] = find_clique(_imp->params.timeout, target_size, _imp->target_graph_rows, g, max_graphs, v,
                _imp->largest_pattern_clique[g], _imp->target_cliques_best_knowns[0], _imp->target_cliques_build_times,
                _imp->target_cliques_solve_times, _imp->target_cliques_solve_find_nodes, _imp->target_cliques_solve_prove_nodes);
        }
}

auto HomomorphismModel::_check_clique_compatibility(int p, int t) const -> bool
{
    if (! _imp->params.clique_size_constraints)
        return true;

    if (! _imp->has_pattern_cliques_sizes)
        _build_pattern_clique_sizes();

    _build_target_clique_size(t);

    for (unsigned g = 0; g < _imp->max_graphs_for_clique_size_constraints; ++g) {
        if (_imp->pattern_cliques_sizes[g][p] > _imp->target_cliques_sizes[g][t]) {
            if (_imp->proof)
                _prove_no_clique(g, p, t);
            return false;
        }
    }

    return true;
}

auto HomomorphismModel::_prove_no_clique(
    unsigned g,
    int p,
    int tt) const -> void
{
    vector<NamedVertex> p_clique;
    map<int, NamedVertex> t_clique_neighbourhood;
    unsigned decide_size;

    {
        vector<int> include(pattern_size, -1), invinclude(pattern_size, 0);
        int count = 0;
        for (int w = 0; w < int(pattern_size); ++w)
            if (w != p && _imp->pattern_graph_rows[w * max_graphs + g].test(p)) {
                include[w] = count;
                invinclude[count] = w;
                ++count;
            }

        InputGraph gv(count, false, false);
        for (unsigned f = 0; f < pattern_size; ++f)
            if (include[f] != -1)
                for (unsigned t = 0; t < pattern_size; ++t)
                    if (f != t && include[t] != -1 && _imp->pattern_graph_rows[f * max_graphs + g].test(t))
                        gv.add_edge(include[f], include[t]);

        CliqueParams params;
        params.timeout = _imp->params.timeout;
        params.start_time = steady_clock::now();
        params.restarts_schedule = make_unique<NoRestartsSchedule>();
        auto result = solve_clique_problem(gv, params);
        for (auto & v : result.clique)
            p_clique.push_back(pattern_vertex_for_proof(invinclude[v]));
        decide_size = result.clique.size();
    }

    {
        vector<int> include(target_size, -1), invinclude(target_size, 0);
        int count = 0;
        for (int w = 0; w < int(target_size); ++w)
            if (w != tt && _imp->target_graph_rows[w * max_graphs + g].test(tt)) {
                t_clique_neighbourhood.emplace(count, target_vertex_for_proof(w));
                include[w] = count;
                invinclude[count] = w;
                ++count;
            }

        _imp->proof->prepare_hom_clique_proof(pattern_vertex_for_proof(p), target_vertex_for_proof(tt), decide_size);

        InputGraph gv(count, false, false);
        for (unsigned f = 0; f < target_size; ++f)
            if (include[f] != -1)
                for (unsigned t = 0; t < target_size; ++t) {
                    if (f != t && include[t] != -1) {
                        if (_imp->target_graph_rows[f * max_graphs + g].test(t))
                            gv.add_edge(include[f], include[t]);
                        else if (f < t)
                            _imp->proof->add_hom_clique_non_edge(
                                pattern_vertex_for_proof(p), target_vertex_for_proof(tt),
                                p_clique, target_vertex_for_proof(f), target_vertex_for_proof(t));
                    }
                }

        _imp->proof->start_hom_clique_proof(pattern_vertex_for_proof(p), move(p_clique), target_vertex_for_proof(tt), move(t_clique_neighbourhood));

        CliqueParams params;
        params.timeout = _imp->params.timeout;
        params.start_time = steady_clock::now();
        params.decide = make_optional(decide_size);
        params.restarts_schedule = make_unique<NoRestartsSchedule>();
        params.extend_proof = _imp->proof;
        params.proof_is_for_hom = true;

        auto result = solve_clique_problem(gv, params);
        if (result.complete && ! result.clique.empty())
            throw ProofError{"Oops, found a clique that shound't exist"};
        _imp->proof->finish_hom_clique_proof(pattern_vertex_for_proof(p), target_vertex_for_proof(tt), decide_size);
    }
}

auto HomomorphismModel::_check_degree_compatibility(
    int p,
    int t,
    unsigned graphs_to_consider,
    vector<vector<vector<int>>> & patterns_ndss,
    vector<vector<optional<vector<int>>>> & targets_ndss,
    bool do_not_do_nds_yet) const -> bool
{
    if (! degree_and_nds_are_preserved(_imp->params))
        return true;

    for (unsigned g = 0; g < graphs_to_consider; ++g) {
        if (target_degree(g, t) < pattern_degree(g, p)) {
            // not ok, degrees differ
            if (_imp->proof) {
                // get the actual neighbours of p and t, in their original terms
                vector<int> n_p, n_t;

                auto np = pattern_graph_row(g, p);
                for (unsigned j = 0; j < pattern_size; ++j)
                    if (np.test(j))
                        n_p.push_back(j);

                auto nt = target_graph_row(g, t);
                for (auto j = nt.find_first(); j != decltype(nt)::npos; j = nt.find_first()) {
                    nt.reset(j);
                    n_t.push_back(j);
                }

                _imp->proof->incompatible_by_degrees(g, pattern_vertex_for_proof(p), n_p,
                    target_vertex_for_proof(t), n_t);
            }
            return false;
        }
        else if (degree_and_nds_are_exact(_imp->params, pattern_size, target_size) && target_degree(g, t) != pattern_degree(g, p)) {
            // not ok, degrees must be exactly the same
            return false;
        }
    }
    if (_imp->params.no_nds || do_not_do_nds_yet)
        return true;

    // full compare of neighbourhood degree sequences
    if (! targets_ndss.at(0).at(t)) {
        for (unsigned g = 0; g < graphs_to_consider; ++g) {
            targets_ndss.at(g).at(t) = vector<int>{};
            auto ni = target_graph_row(g, t);
            for (auto j = ni.find_first(); j != decltype(ni)::npos; j = ni.find_first()) {
                ni.reset(j);
                targets_ndss.at(g).at(t)->push_back(target_degree(g, j));
            }
            sort(targets_ndss.at(g).at(t)->begin(), targets_ndss.at(g).at(t)->end(), greater<int>());
        }
    }

    for (unsigned g = 0; g < graphs_to_consider; ++g) {
        for (unsigned x = 0; x < patterns_ndss.at(g).at(p).size(); ++x) {
            if (targets_ndss.at(g).at(t)->at(x) < patterns_ndss.at(g).at(p).at(x)) {
                if (_imp->proof) {
                    vector<int> p_subsequence, t_subsequence, t_remaining;

                    // need to know the NDS together with the actual vertices
                    vector<pair<int, int>> p_nds, t_nds;

                    auto np = pattern_graph_row(g, p);
                    for (auto w = np.find_first(); w != decltype(np)::npos; w = np.find_first()) {
                        np.reset(w);
                        p_nds.emplace_back(w, pattern_graph_row(g, w).count());
                    }

                    auto nt = target_graph_row(g, t);
                    for (auto w = nt.find_first(); w != decltype(nt)::npos; w = nt.find_first()) {
                        nt.reset(w);
                        t_nds.emplace_back(w, target_graph_row(g, w).count());
                    }

                    sort(p_nds.begin(), p_nds.end(), [](const pair<int, int> & a, const pair<int, int> & b) { return a.second > b.second; });
                    sort(t_nds.begin(), t_nds.end(), [](const pair<int, int> & a, const pair<int, int> & b) { return a.second > b.second; });

                    for (unsigned y = 0; y <= x; ++y) {
                        p_subsequence.push_back(p_nds[y].first);
                        t_subsequence.push_back(t_nds[y].first);
                    }
                    for (unsigned y = x + 1; y < t_nds.size(); ++y)
                        t_remaining.push_back(t_nds[y].first);

                    _imp->proof->incompatible_by_nds(g, pattern_vertex_for_proof(p),
                        target_vertex_for_proof(t), p_subsequence, t_subsequence, t_remaining);
                }
                return false;
            }
            else if (degree_and_nds_are_exact(_imp->params, pattern_size, target_size) && targets_ndss.at(g).at(t)->at(x) != patterns_ndss.at(g).at(p).at(x))
                return false;
        }
    }

    return true;
}

auto HomomorphismModel::initialise_domains(vector<HomomorphismDomain> & domains) const -> bool
{
    unsigned max_graphs_for_degree_things = (_imp->params.injectivity == Injectivity::LocallyInjective ? 1 : max_graphs);

    /* pattern and target neighbourhood degree sequences */
    vector<vector<vector<int>>> patterns_ndss(max_graphs_for_degree_things);
    vector<vector<optional<vector<int>>>> targets_ndss(max_graphs_for_degree_things);

    if (degree_and_nds_are_preserved(_imp->params) && ! _imp->params.no_nds) {
        for (unsigned g = 0; g < max_graphs_for_degree_things; ++g) {
            patterns_ndss.at(g).resize(pattern_size);
            targets_ndss.at(g).resize(target_size);
        }

        for (unsigned g = 0; g < max_graphs_for_degree_things; ++g) {
            for (unsigned i = 0; i < pattern_size; ++i) {
                auto ni = pattern_graph_row(g, i);
                for (auto j = ni.find_first(); j != decltype(ni)::npos; j = ni.find_first()) {
                    ni.reset(j);
                    patterns_ndss.at(g).at(i).push_back(pattern_degree(g, j));
                }
                sort(patterns_ndss.at(g).at(i).begin(), patterns_ndss.at(g).at(i).end(), greater<int>());
            }
        }
    }

    for (unsigned i = 0; i < pattern_size; ++i) {
        domains.at(i).v = i;
        domains.at(i).values.reset();

        for (unsigned j = 0; j < target_size; ++j) {
            bool ok = true;

            if (! _check_label_compatibility(i, j))
                ok = false;
            else if (! _check_loop_compatibility(i, j))
                ok = false;
            else if (! _check_degree_compatibility(i, j, max_graphs_for_degree_things, patterns_ndss, targets_ndss, _imp->proof.get()))
                ok = false;
            else if (! _check_clique_compatibility(i, j))
                ok = false;

            if (ok)
                domains.at(i).values.set(j);
        }

        domains.at(i).count = domains.at(i).values.count();
        if (0 == domains.at(i).count) {
            if (_imp->proof)
                _imp->proof->initial_domain_is_empty(domains.at(i).v, "compatibility stage");
            return false;
        }
    }

    // for proof logging, we need degree information before we can output nds proofs
    if (_imp->proof && degree_and_nds_are_preserved(_imp->params) && ! _imp->params.no_nds) {
        for (unsigned i = 0; i < pattern_size; ++i) {
            for (unsigned j = 0; j < target_size; ++j) {
                if (domains.at(i).values.test(j) &&
                    ! _check_degree_compatibility(i, j, max_graphs_for_degree_things, patterns_ndss, targets_ndss, false)) {
                    domains.at(i).values.reset(j);
                    if (0 == --domains.at(i).count) {
                        if (_imp->proof)
                            _imp->proof->initial_domain_is_empty(domains.at(i).v, "nds stage");
                        return false;
                    }
                }
            }
        }
    }

    // quick sanity check that we have enough values
    if (is_nonshrinking(_imp->params)) {
        SVOBitset domains_union{target_size, 0};
        for (auto & d : domains)
            domains_union |= d.values;

        unsigned domains_union_popcount = domains_union.count();
        if (domains_union_popcount < unsigned(pattern_size)) {
            if (_imp->proof) {
                vector<NamedVertex> hall_lhs, hall_rhs;
                for (auto & d : domains)
                    hall_lhs.push_back(pattern_vertex_for_proof(d.v));
                auto dd = domains_union;
                for (auto v = dd.find_first(); v != decltype(dd)::npos; v = dd.find_first()) {
                    dd.reset(v);
                    hall_rhs.push_back(target_vertex_for_proof(v));
                }
                _imp->proof->emit_hall_set_or_violator(hall_lhs, hall_rhs);
            }
            return false;
        }
    }

    for (auto & d : domains) {
        d.count = d.values.count();
        if (0 == d.count && _imp->proof) {
            _imp->proof->initial_domain_is_empty(d.v, "post-initialisation stage");
            return false;
        }
    }

    if (_imp->params.lackey) {
        // If we're dealing with a model from Essence, it's possible some values will
        // be completely eliminated from the upper and lower bounds of certain domains,
        // in which case they won't show up as being deleted during propagation.
        bool wipeout = false;
        if (! _imp->params.lackey->reduce_initial_bounds([&](int p, int t) -> void {
                for (auto & d : domains)
                    if (d.v == unsigned(p)) {
                        if (d.values.test(t)) {
                            d.values.reset(t);
                            if (0 == --d.count)
                                wipeout = true;
                        }
                        break;
                    }
            }) ||
            wipeout) {
            return false;
        }
    }

    return true;
}

auto HomomorphismModel::pattern_vertex_for_proof(int v) const -> NamedVertex
{
    if (v < 0 || unsigned(v) >= _imp->pattern_vertex_proof_names.size())
        throw ProofError{"Oops, there's a bug: v out of range in pattern"};
    return pair{v, _imp->pattern_vertex_proof_names[v]};
}

auto HomomorphismModel::target_vertex_for_proof(int v) const -> NamedVertex
{
    if (v < 0 || unsigned(v) >= _imp->target_vertex_proof_names.size())
        throw ProofError{"Oops, there's a bug: v out of range in target"};
    return pair{v, _imp->target_vertex_proof_names[v]};
}

auto HomomorphismModel::prepare() -> bool
{
    if (is_nonshrinking(_imp->params) && (pattern_size > target_size))
        return false;

    _imp->supplemental_graph_names.push_back("original");

    // pattern and target degrees, for the main graph
    _imp->patterns_degrees.at(0).resize(pattern_size);
    _imp->targets_degrees.at(0).resize(target_size);

    for (unsigned i = 0; i < pattern_size; ++i)
        _imp->patterns_degrees.at(0).at(i) = _imp->pattern_graph_rows[i * max_graphs + 0].count();

    for (unsigned i = 0; i < target_size; ++i)
        _imp->targets_degrees.at(0).at(i) = _imp->target_graph_rows[i * max_graphs + 0].count();

    if (global_degree_is_preserved(_imp->params)) {
        vector<pair<int, int>> p_gds, t_gds;
        for (unsigned i = 0; i < pattern_size; ++i)
            p_gds.emplace_back(i, _imp->patterns_degrees.at(0).at(i));
        for (unsigned i = 0; i < target_size; ++i)
            t_gds.emplace_back(i, _imp->targets_degrees.at(0).at(i));

        sort(p_gds.begin(), p_gds.end(), [](const pair<int, int> & a, const pair<int, int> & b) { return a.second > b.second; });
        sort(t_gds.begin(), t_gds.end(), [](const pair<int, int> & a, const pair<int, int> & b) { return a.second > b.second; });

        for (unsigned i = 0; i < p_gds.size(); ++i)
            if (p_gds.at(i).second > t_gds.at(i).second) {
                if (_imp->proof) {
                    for (unsigned p = 0; p <= i; ++p) {
                        vector<int> n_p;
                        auto np = _imp->pattern_graph_rows[p_gds.at(p).first * max_graphs + 0];
                        for (unsigned j = 0; j < pattern_size; ++j)
                            if (np.test(j))
                                n_p.push_back(j);

                        for (unsigned t = i; t < t_gds.size(); ++t) {
                            vector<int> n_t;
                            auto nt = _imp->target_graph_rows[t_gds.at(t).first * max_graphs + 0];
                            for (auto j = nt.find_first(); j != decltype(nt)::npos; j = nt.find_first()) {
                                nt.reset(j);
                                n_t.push_back(j);
                            }

                            _imp->proof->incompatible_by_degrees(0,
                                pattern_vertex_for_proof(p_gds.at(p).first), n_p,
                                target_vertex_for_proof(t_gds.at(t).first), n_t);
                        }
                    }

                    vector<NamedVertex> patterns, targets;
                    for (unsigned p = 0; p <= i; ++p)
                        patterns.push_back(pattern_vertex_for_proof(p_gds.at(p).first));
                    for (unsigned t = 0; t < i; ++t)
                        targets.push_back(target_vertex_for_proof(t_gds.at(t).first));

                    _imp->proof->emit_hall_set_or_violator(patterns, targets);
                }
                return false;
            }
    }

    unsigned next_pattern_supplemental = 1, next_target_supplemental = 1;
    // build exact path graphs
    if (supports_exact_path_graphs(_imp->params)) {
        _build_exact_path_graphs(_imp->pattern_graph_rows, pattern_size, next_pattern_supplemental, _imp->params.number_of_exact_path_graphs, _imp->directed, false, true);
        _build_exact_path_graphs(_imp->target_graph_rows, target_size, next_target_supplemental, _imp->params.number_of_exact_path_graphs, _imp->directed, false, false);

        if (_imp->proof) {
            for (int g = 1; g <= _imp->params.number_of_exact_path_graphs; ++g) {
                for (unsigned p = 0; p < pattern_size; ++p) {
                    for (unsigned q = 0; q < pattern_size; ++q) {
                        if (p == q || ! _imp->pattern_graph_rows[p * max_graphs + g].test(q))
                            continue;

                        auto named_p = pattern_vertex_for_proof(p);
                        auto named_q = pattern_vertex_for_proof(q);

                        auto n_p_q = _imp->pattern_graph_rows[p * max_graphs + 0];
                        n_p_q &= _imp->pattern_graph_rows[q * max_graphs + 0];
                        vector<NamedVertex> between_p_and_q;
                        for (auto v = n_p_q.find_first(); v != decltype(n_p_q)::npos; v = n_p_q.find_first()) {
                            n_p_q.reset(v);
                            between_p_and_q.push_back(pattern_vertex_for_proof(v));
                            if (between_p_and_q.size() >= unsigned(g))
                                break;
                        }

                        for (unsigned t = 0; t < target_size; ++t) {
                            auto named_t = target_vertex_for_proof(t);

                            vector<NamedVertex> named_n_t, named_d_n_t;
                            vector<pair<NamedVertex, vector<NamedVertex>>> named_two_away_from_t;
                            auto n_t = _imp->target_graph_rows[t * max_graphs + 0];
                            for (auto w = n_t.find_first(); w != decltype(n_t)::npos; w = n_t.find_first()) {
                                n_t.reset(w);
                                named_n_t.push_back(target_vertex_for_proof(w));
                            }

                            auto nd_t = _imp->target_graph_rows[t * max_graphs + g];
                            for (auto w = nd_t.find_first(); w != decltype(nd_t)::npos; w = nd_t.find_first()) {
                                nd_t.reset(w);
                                named_d_n_t.push_back(target_vertex_for_proof(w));
                            }

                            auto n2_t = _imp->target_graph_rows[t * max_graphs + 1];
                            for (auto w = n2_t.find_first(); w != decltype(n2_t)::npos; w = n2_t.find_first()) {
                                n2_t.reset(w);
                                auto n_t_w = _imp->target_graph_rows[w * max_graphs + 0];
                                n_t_w &= _imp->target_graph_rows[t * max_graphs + 0];
                                vector<NamedVertex> named_n_t_w;
                                for (auto x = n_t_w.find_first(); x != decltype(n_t_w)::npos; x = n_t_w.find_first()) {
                                    n_t_w.reset(x);
                                    named_n_t_w.push_back(target_vertex_for_proof(x));
                                }
                                named_two_away_from_t.emplace_back(target_vertex_for_proof(w), named_n_t_w);
                            }

                            _imp->proof->create_exact_path_graphs(g, named_p, named_q, between_p_and_q,
                                named_t, named_n_t, named_two_away_from_t, named_d_n_t);
                        }
                    }
                }
            }
        }
    }

    if (supports_distance2_graphs(_imp->params)) {
        _build_exact_path_graphs(_imp->pattern_graph_rows, pattern_size, next_pattern_supplemental, 1, _imp->directed, true, true);
        _build_exact_path_graphs(_imp->target_graph_rows, target_size, next_target_supplemental, 1, _imp->directed, true, false);
    }

    if (supports_distance3_graphs(_imp->params)) {
        _build_distance3_graphs(_imp->pattern_graph_rows, pattern_size, next_pattern_supplemental, true);
        _build_distance3_graphs(_imp->target_graph_rows, target_size, next_target_supplemental, false);

        if (_imp->proof) {
            for (unsigned p = 0; p < pattern_size; ++p) {
                for (unsigned q = 0; q < pattern_size; ++q) {
                    auto named_p = pattern_vertex_for_proof(p);
                    auto named_q = pattern_vertex_for_proof(q);

                    // only do this if they're actually adjacent
                    if (p == q || ! _imp->pattern_graph_rows[p * max_graphs + next_pattern_supplemental - 1].test(q))
                        continue;

                    bool actually_adjacent = false;
                    optional<NamedVertex> path_from_p_to_q_1 = nullopt, path_from_p_to_q_2 = nullopt;

                    auto n_p = _imp->pattern_graph_rows[p * max_graphs + 0];

                    // are they actually distance 1 apart?
                    if (n_p.test(q))
                        actually_adjacent = true;
                    else {
                        auto n_q = _imp->pattern_graph_rows[q * max_graphs + 0];

                        auto n_p_q = n_p;
                        n_p_q &= n_q;
                        n_p_q.reset(p);
                        n_p_q.reset(q);

                        if (n_p_q.any()) {
                            // they're actually distance 2 apart
                            path_from_p_to_q_1 = pattern_vertex_for_proof(n_p_q.find_first());
                        }
                        else {
                            // find a path of length 3
                            n_p.reset(p);
                            n_p.reset(q);
                            for (auto v = n_p.find_first(); v != decltype(n_p)::npos && ! path_from_p_to_q_1; v = n_p.find_first()) {
                                n_p.reset(v);
                                auto n_v = _imp->pattern_graph_rows[v * max_graphs + 0];
                                n_v.reset(v);
                                n_v.reset(p);
                                n_v.reset(q);
                                for (auto w = n_v.find_first(); w != decltype(n_v)::npos && ! path_from_p_to_q_1; w = n_v.find_first()) {
                                    n_v.reset(w);
                                    if (_imp->pattern_graph_rows[w * max_graphs + 0].test(q)) {
                                        path_from_p_to_q_1 = pattern_vertex_for_proof(v);
                                        path_from_p_to_q_2 = pattern_vertex_for_proof(w);
                                    }
                                }
                            }
                        }

                        if (! path_from_p_to_q_1)
                            throw ProofError{"Oops, there's a bug: missing path from " + named_p.second + " to " + named_q.second};
                    }

                    for (unsigned t = 0; t < target_size; ++t) {
                        auto named_t = target_vertex_for_proof(t);

                        vector<NamedVertex> d1_from_t, d2_from_t, d3_from_t;
                        set<NamedVertex> d2_from_t_set, d3_from_t_set;
                        auto n_t = _imp->target_graph_rows[t * max_graphs + 0];
                        n_t.set(t);
                        for (auto v = n_t.find_first(); v != decltype(n_t)::npos; v = n_t.find_first()) {
                            n_t.reset(v);
                            d1_from_t.push_back(target_vertex_for_proof(v));
                            auto n_v = _imp->target_graph_rows[v * max_graphs + 0];
                            n_v.set(v);
                            for (auto w = n_v.find_first(); w != decltype(n_v)::npos; w = n_v.find_first()) {
                                n_v.reset(w);
                                d2_from_t_set.insert(target_vertex_for_proof(w));
                                auto n_w = _imp->target_graph_rows[w * max_graphs + 0];
                                n_w.set(w);
                                for (auto x = n_w.find_first(); x != decltype(n_w)::npos; x = n_w.find_first()) {
                                    n_w.reset(x);
                                    d3_from_t_set.insert(target_vertex_for_proof(x));
                                }
                            }
                        }

                        d2_from_t.assign(d2_from_t_set.begin(), d2_from_t_set.end());
                        d3_from_t.assign(d3_from_t_set.begin(), d3_from_t_set.end());

                        if (actually_adjacent)
                            _imp->proof->create_distance3_graphs_but_actually_distance_1(next_pattern_supplemental - 1, named_p, named_q, named_t, d3_from_t);
                        else if (path_from_p_to_q_2)
                            _imp->proof->create_distance3_graphs(next_pattern_supplemental - 1, named_p, named_q, *path_from_p_to_q_1,
                                *path_from_p_to_q_2, named_t, d1_from_t, d2_from_t, d3_from_t);
                        else
                            _imp->proof->create_distance3_graphs_but_actually_distance_2(next_pattern_supplemental - 1, named_p, named_q, *path_from_p_to_q_1,
                                named_t, d1_from_t, d2_from_t, d3_from_t);
                    }
                }
            }
        }
    }

    if (supports_k4_graphs(_imp->params)) {
        _build_k4_graphs(_imp->pattern_graph_rows, pattern_size, next_pattern_supplemental, true);
        _build_k4_graphs(_imp->target_graph_rows, target_size, next_target_supplemental, false);
    }

    for (auto & [shape, injective, count] : _imp->params.extra_shapes) {
        _build_extra_shape(_imp->pattern_graph_rows, pattern_size, next_pattern_supplemental, *shape, injective, count, true);
        _build_extra_shape(_imp->target_graph_rows, target_size, next_target_supplemental, *shape, injective, count, false);

        if (_imp->proof) {
            for (unsigned p = 0; p < pattern_size; ++p) {
                for (unsigned q = 0; q < pattern_size; ++q) {
                    auto named_p = pattern_vertex_for_proof(p);
                    auto named_q = pattern_vertex_for_proof(q);

                    // only do this if they're actually adjacent
                    if (! _imp->pattern_graph_rows[p * max_graphs + next_pattern_supplemental - 1].test(q))
                        continue;

                    for (unsigned t = 0; t < target_size; ++t) {
                        auto named_t = target_vertex_for_proof(t);
                        vector<NamedVertex> named_n_t;
                        auto n_t = _imp->target_graph_rows[t * max_graphs + next_pattern_supplemental - 1];
                        for (auto v = n_t.find_first(); v != decltype(n_t)::npos; v = n_t.find_first()) {
                            n_t.reset(v);
                            named_n_t.push_back(target_vertex_for_proof(v));
                        }
                        _imp->proof->hack_in_shape_graph(next_pattern_supplemental - 1, named_p, named_q, named_t, named_n_t);
                    }
                }
            }
        }
    }

    if (next_pattern_supplemental != max_graphs || next_target_supplemental != max_graphs ||
            next_pattern_supplemental != _imp->supplemental_graph_names.size())
        throw UnsupportedConfiguration{"something has gone wrong with supplemental graph indexing: " + to_string(next_pattern_supplemental) + " " + to_string(next_target_supplemental) + " " + to_string(max_graphs) + " "
        + to_string(_imp->supplemental_graph_names.size())};

    // pattern and target degrees, for supplemental graphs
    for (unsigned g = 1; g < max_graphs; ++g) {
        _imp->patterns_degrees.at(g).resize(pattern_size);
        _imp->targets_degrees.at(g).resize(target_size);
    }

    for (unsigned g = 1; g < max_graphs; ++g) {
        for (unsigned i = 0; i < pattern_size; ++i)
            _imp->patterns_degrees.at(g).at(i) = _imp->pattern_graph_rows[i * max_graphs + g].count();

        for (unsigned i = 0; i < target_size; ++i)
            _imp->targets_degrees.at(g).at(i) = _imp->target_graph_rows[i * max_graphs + g].count();
    }

    for (unsigned i = 0; i < target_size; ++i)
        _imp->largest_target_degree = max(_imp->largest_target_degree, _imp->targets_degrees[0][i]);

    // re-add loops
    for (unsigned i = 0; i < pattern_size; ++i)
        if (_imp->pattern_loops[i])
            _imp->pattern_graph_rows[i * max_graphs + 0].set(i);

    for (unsigned i = 0; i < target_size; ++i)
        if (_imp->target_loops[i])
            _imp->target_graph_rows[i * max_graphs + 0].set(i);

    // pattern adjacencies, compressed
    _imp->pattern_adjacencies_bits.resize(pattern_size * pattern_size);
    for (unsigned g = 0; g < max_graphs; ++g)
        for (unsigned i = 0; i < pattern_size; ++i)
            for (unsigned j = 0; j < pattern_size; ++j)
                if (_imp->pattern_graph_rows[i * max_graphs + g].test(j))
                    _imp->pattern_adjacencies_bits[i * pattern_size + j] |= (1u << g);

    return true;
}

auto HomomorphismModel::_build_exact_path_graphs(vector<SVOBitset> & graph_rows, unsigned size, unsigned & idx,
    unsigned number_of_exact_path_graphs, bool directed, bool at_most, bool pattern) -> void
{
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
            _imp->supplemental_graph_names.push_back("exact_path_" + to_string(p));

    idx += number_of_exact_path_graphs;
}

auto HomomorphismModel::_build_distance3_graphs(vector<SVOBitset> & graph_rows, unsigned size, unsigned & idx,
        bool pattern) -> void
{
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
        _imp->supplemental_graph_names.push_back("distance3");
    ++idx;
}

auto HomomorphismModel::_build_k4_graphs(vector<SVOBitset> & graph_rows, unsigned size, unsigned & idx, bool pattern) -> void
{
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
        _imp->supplemental_graph_names.push_back("k4");
    ++idx;
}

auto HomomorphismModel::_build_extra_shape(vector<SVOBitset> & graph_rows, unsigned size, unsigned & idx, InputGraph & shape,
    bool injective, int count, bool pattern) -> void
{
    InputGraph master_graph(size, true, false);

    for (unsigned v = 0; v < size; ++v) {
        auto nv = graph_rows[v * max_graphs + 0];
        for (unsigned w = 0; w < v; ++w) {
            if (nv.test(w))
                master_graph.add_edge(v, w);
        }
    }

    for (unsigned v = 0; v < size && ! _imp->params.timeout->should_abort(); ++v) {
        for (unsigned w = 0; w < v && ! _imp->params.timeout->should_abort(); ++w) {
            HomomorphismParams child_params;
            child_params.timeout = _imp->params.timeout;
            child_params.start_time = _imp->params.start_time;
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

            master_graph.set_vertex_label(v, "from");
            master_graph.set_vertex_label(w, "to");
            auto child_result = solve_homomorphism_problem(shape, master_graph, child_params);
            if ((count > 1) ? seen_count >= count : ! child_result.mapping.empty()) {
                graph_rows[v * max_graphs + idx].set(w);
                graph_rows[w * max_graphs + idx].set(v);
            }
            master_graph.set_vertex_label(v, "");
            master_graph.set_vertex_label(w, "");
        }
    }

    if (pattern)
        _imp->supplemental_graph_names.push_back("extra_shape");
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

auto HomomorphismModel::has_occur_less_thans() const -> bool
{
    return _imp->has_occur_less_thans;
}

auto HomomorphismModel::directed() const -> bool
{
    return _imp->directed;
}

auto HomomorphismModel::add_extra_stats(list<string> & x) const -> void
{
    auto join = [](string_view t, auto & l) -> string {
        stringstream s;
        s << t;
        for (auto & i : l)
            s << " " << i;
        return s.str();
    };

    if (! _imp->pattern_cliques_sizes.empty()) {
        x.emplace_back(join("pattern_cliques_build_times =", _imp->pattern_cliques_build_times));
        x.emplace_back(join("pattern_cliques_solve_times =", _imp->pattern_cliques_solve_times));
        x.emplace_back(join("pattern_cliques_solve_find_nodes =", _imp->pattern_cliques_solve_find_nodes));
        x.emplace_back(join("pattern_cliques_solve_prove_nodes =", _imp->pattern_cliques_solve_prove_nodes));

        x.emplace_back(join("target_cliques_build_times =", _imp->target_cliques_build_times));
        x.emplace_back(join("target_cliques_solve_times =", _imp->target_cliques_solve_times));
        x.emplace_back(join("target_cliques_solve_find_nodes =", _imp->target_cliques_solve_find_nodes));
        x.emplace_back(join("target_cliques_solve_prove_nodes =", _imp->target_cliques_solve_prove_nodes));
    }

    x.emplace_back(join("supplemental_graph_names =", _imp->supplemental_graph_names));
}
