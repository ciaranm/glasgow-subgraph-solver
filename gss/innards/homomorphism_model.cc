#include <gss/clique.hh>
#include <gss/configuration.hh>
#include <gss/innards/clique_size_constraints.hh>
#include <gss/innards/homomorphism_model.hh>
#include <gss/innards/homomorphism_proofs.hh>
#include <gss/innards/homomorphism_traits.hh>
#include <gss/innards/processed_graphs_data.hh>
#include <gss/innards/supplemental_graphs.hh>

#include <chrono>
#include <functional>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
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
    // One entry of the shape-graph plan: a supplemental graph (or run of graphs) to
    // build, in build order. The plan is the single source of truth for which
    // supplemental graphs the model has -- max_graphs is derived from it (see
    // number_of_shape_graphs) and prepare() builds from it -- so the count, the build
    // order, and the per-graph applicability all live in one place, rather than in a
    // separate count function and a parallel sequence of build blocks that have to be
    // kept in step (dev_docs/preprocessor-refactor.md, S1).
    struct ShapeGraphSpec
    {
        enum class Kind
        {
            ExactPath,
            Distance2,
            Distance3,
            K4,
            ExtraShape
        };

        Kind kind;
        unsigned slot_count;
        const std::tuple<std::unique_ptr<InputGraph>, bool, int> * extra_shape = nullptr;
    };

    auto make_shape_graph_plan(const HomomorphismParams & params, bool has_loops) -> vector<ShapeGraphSpec>
    {
        vector<ShapeGraphSpec> plan;
        if (supports_exact_path_graphs(params, has_loops))
            plan.push_back({ShapeGraphSpec::Kind::ExactPath, unsigned(params.number_of_exact_path_graphs)});
        if (supports_distance2_graphs(params, has_loops))
            plan.push_back({ShapeGraphSpec::Kind::Distance2, 1});
        if (supports_distance3_graphs(params))
            plan.push_back({ShapeGraphSpec::Kind::Distance3, 1});
        if (supports_k4_graphs(params, has_loops))
            plan.push_back({ShapeGraphSpec::Kind::K4, 1});
        for (auto & shape : params.extra_shapes)
            plan.push_back({ShapeGraphSpec::Kind::ExtraShape, 1, &shape});
        return plan;
    }

    // The original graph plus every slot the plan allocates: this is the bitset stride
    // (max_graphs), computed once at construction so the graph rows can be sized.
    auto number_of_shape_graphs(const HomomorphismParams & params, bool has_loops) -> unsigned
    {
        unsigned n = 1;
        for (const auto & spec : make_shape_graph_plan(params, has_loops))
            n += spec.slot_count;
        return n;
    }

}

struct HomomorphismModel::Imp
{
    const HomomorphismParams & params;
    shared_ptr<Proof> proof;

    // The recoded graphs and everything derived from them (see processed_graphs_data.hh).
    ProcessedGraphsData graphs;

    bool has_less_thans = false, has_occur_less_thans = false;

    // The solver-proofs middle layer (owns vertex naming + hom-specific derivations).
    // Non-owning: it is created and owned by the solve pipeline (so it can emit the model
    // before the model is built) and outlives this object. Null when proof logging is off.
    HomomorphismProofs * proofs = nullptr;

    // Clique-size filtering caches (see clique_size_constraints.hh); mutable because the
    // sizes are computed lazily during the const domain-compatibility checks.
    mutable CliqueSizeData clique_data;

    Imp(const HomomorphismParams & p, const std::shared_ptr<Proof> & r, HomomorphismProofs * pr) :
        params(p),
        proof(r),
        proofs(pr)
    {
    }
};

HomomorphismModel::HomomorphismModel(const InputGraph & target, const InputGraph & pattern, const HomomorphismParams & params,
    const std::shared_ptr<Proof> & proof, HomomorphismProofs * proofs) :
    _imp(make_unique<Imp>(params, proof, proofs)),
    max_graphs(number_of_shape_graphs(params, pattern.loopy() || target.loopy())),
    pattern_size(pattern.size()),
    target_size(target.size())
{
    _imp->graphs.has_loops = pattern.loopy() || target.loopy();

    _imp->graphs.patterns_degrees.resize(max_graphs);
    _imp->graphs.targets_degrees.resize(max_graphs);

    if (max_graphs > 8 * sizeof(PatternAdjacencyBitsType))
        throw UnsupportedConfiguration{"Supplemental graphs won't fit in the chosen bitset size"};

    if (pattern.directed())
        _imp->graphs.directed = true;

    // recode pattern to a bit graph, and strip out loops
    _imp->graphs.pattern_graph_rows.resize(pattern_size * max_graphs, SVOBitset(pattern_size, 0));
    _imp->graphs.pattern_loops.resize(pattern_size);
    for (unsigned i = 0; i < pattern_size; ++i) {
        for (unsigned j = 0; j < pattern_size; ++j) {
            if (pattern.adjacent(i, j)) {
                if (i == j)
                    _imp->graphs.pattern_loops[i] = 1;
                else
                    _imp->graphs.pattern_graph_rows[i * max_graphs + 0].set(j);
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

        _imp->graphs.pattern_vertex_labels.resize(pattern_size);
        for (unsigned i = 0; i < pattern_size; ++i)
            _imp->graphs.pattern_vertex_labels[i] = vertex_labels_map.find(string{pattern.vertex_label(i)})->second;
    }

    // re-encode and store edge labels
    map<string, int> edge_labels_map;
    int next_edge_label = 1;
    if (pattern.has_edge_labels()) {
        _imp->graphs.pattern_edge_labels.resize(pattern_size * pattern_size);
        for (unsigned i = 0; i < pattern_size; ++i)
            for (unsigned j = 0; j < pattern_size; ++j)
                if (pattern.adjacent(i, j)) {
                    auto r = edge_labels_map.emplace(pattern.edge_label(i, j), next_edge_label);
                    if (r.second)
                        ++next_edge_label;
                    _imp->graphs.pattern_edge_labels[i * pattern_size + j] = r.first->second;
                }
    }

    // recode target to a bit graph, and take out loops
    _imp->graphs.target_graph_rows.resize(target_size * max_graphs, SVOBitset{target_size, 0});
    _imp->graphs.target_loops.resize(target_size);
    target.for_each_edge([&](int f, int t, string_view) {
        if (f == t)
            _imp->graphs.target_loops[f] = 1;
        else
            _imp->graphs.target_graph_rows[f * max_graphs + 0].set(t);
    });

    // if directed, do both directions
    if (pattern.directed()) {
        _imp->graphs.forward_target_graph_rows.resize(target_size, SVOBitset{target_size, 0});
        _imp->graphs.reverse_target_graph_rows.resize(target_size, SVOBitset{target_size, 0});
        target.for_each_edge([&](int f, int t, string_view l) {
            if (f != t && l != "unlabelled") {
                _imp->graphs.forward_target_graph_rows[f].set(t);
                _imp->graphs.reverse_target_graph_rows[t].set(f);
            }
        });
    }

    // target vertex labels
    if (pattern.has_vertex_labels()) {
        for (unsigned i = 0; i < target_size; ++i) {
            if (vertex_labels_map.emplace(target.vertex_label(i), next_vertex_label).second)
                ++next_vertex_label;
        }

        _imp->graphs.target_vertex_labels.resize(target_size);
        for (unsigned i = 0; i < target_size; ++i)
            _imp->graphs.target_vertex_labels[i] = vertex_labels_map.find(string{target.vertex_label(i)})->second;
    }

    // target edge labels
    if (pattern.has_edge_labels()) {
        _imp->graphs.target_edge_labels.resize(target_size * target_size);
        target.for_each_edge([&](int f, int t, string_view l) {
            auto r = edge_labels_map.emplace(l, next_edge_label);
            if (r.second)
                ++next_edge_label;

            _imp->graphs.target_edge_labels[f * target_size + t] = r.first->second;
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

    // set up the clique-size filtering caches
    init_clique_size_data(_imp->clique_data, params, max_graphs, pattern.size(), target.size());
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

auto HomomorphismModel::_check_clique_compatibility(int p, int t) const -> bool
{
    return check_clique_compatibility(_imp->clique_data, _imp->graphs, max_graphs, pattern_size, target_size,
        _imp->params, _imp->proofs, p, t);
}

auto HomomorphismModel::_check_degree_compatibility(
    int p,
    int t,
    unsigned graphs_to_consider,
    vector<vector<vector<int>>> & patterns_ndss,
    vector<vector<optional<vector<int>>>> & targets_ndss,
    bool do_not_do_nds_yet) const -> bool
{
    if (! degree_and_nds_are_preserved(_imp->params, _imp->graphs.has_loops))
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

    if (degree_and_nds_are_preserved(_imp->params, _imp->graphs.has_loops) && ! _imp->params.no_nds) {
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
    if (_imp->proof && degree_and_nds_are_preserved(_imp->params, _imp->graphs.has_loops) && ! _imp->params.no_nds) {
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

    return true;
}

auto HomomorphismModel::pattern_vertex_for_proof(int v) const -> NamedVertex
{
    return _imp->proofs->pattern_vertex(v);
}

auto HomomorphismModel::target_vertex_for_proof(int v) const -> NamedVertex
{
    return _imp->proofs->target_vertex(v);
}

auto HomomorphismModel::prepare() -> bool
{
    if (is_nonshrinking(_imp->params) && (pattern_size > target_size))
        return false;

    _imp->graphs.supplemental_graph_names.push_back("original");

    // Emit every loop-based incompatibility up front, as a unit clause. The target graph
    // rows have self-loops stripped, so the adjacency, supplemental-graph and degree
    // derivations all carry stray "maps to a loopy vertex" terms; having ~x_p_t available
    // as a unit lets unit propagation (and the search-failure rups) discharge them. This
    // also covers the induced loop-mismatch case, which otherwise prunes p->t silently
    // without a proof line (issue #56).
    if (_imp->proof) {
        for (unsigned p = 0; p < pattern_size; ++p)
            for (unsigned t = 0; t < target_size; ++t)
                if ((pattern_has_loop(p) && ! target_has_loop(t)) ||
                    (_imp->params.induced && (pattern_has_loop(p) != target_has_loop(t))))
                    _imp->proof->incompatible_by_loops(pattern_vertex_for_proof(p), target_vertex_for_proof(t));
    }

    // pattern and target degrees, for the main graph
    _imp->graphs.patterns_degrees.at(0).resize(pattern_size);
    _imp->graphs.targets_degrees.at(0).resize(target_size);

    for (unsigned i = 0; i < pattern_size; ++i)
        _imp->graphs.patterns_degrees.at(0).at(i) = _imp->graphs.pattern_graph_rows[i * max_graphs + 0].count();

    for (unsigned i = 0; i < target_size; ++i)
        _imp->graphs.targets_degrees.at(0).at(i) = _imp->graphs.target_graph_rows[i * max_graphs + 0].count();

    if (global_degree_is_preserved(_imp->params)) {
        vector<pair<int, int>> p_gds, t_gds;
        for (unsigned i = 0; i < pattern_size; ++i)
            p_gds.emplace_back(i, _imp->graphs.patterns_degrees.at(0).at(i));
        for (unsigned i = 0; i < target_size; ++i)
            t_gds.emplace_back(i, _imp->graphs.targets_degrees.at(0).at(i));

        sort(p_gds.begin(), p_gds.end(), [](const pair<int, int> & a, const pair<int, int> & b) { return a.second > b.second; });
        sort(t_gds.begin(), t_gds.end(), [](const pair<int, int> & a, const pair<int, int> & b) { return a.second > b.second; });

        for (unsigned i = 0; i < p_gds.size(); ++i)
            if (p_gds.at(i).second > t_gds.at(i).second) {
                if (_imp->proof) {
                    for (unsigned p = 0; p <= i; ++p) {
                        vector<int> n_p;
                        auto np = _imp->graphs.pattern_graph_rows[p_gds.at(p).first * max_graphs + 0];
                        for (unsigned j = 0; j < pattern_size; ++j)
                            if (np.test(j))
                                n_p.push_back(j);

                        for (unsigned t = i; t < t_gds.size(); ++t) {
                            vector<int> n_t;
                            auto nt = _imp->graphs.target_graph_rows[t_gds.at(t).first * max_graphs + 0];
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

    // Build every supplemental graph the plan registers, in plan order, each into the
    // next free slot(s), then (when proving) derive it through the solver-proofs layer.
    // The plan also fixes max_graphs, so the bump counters land exactly on max_graphs at
    // the end (checked below).
    for (const auto & spec : make_shape_graph_plan(_imp->params, _imp->graphs.has_loops)) {
        switch (spec.kind) {
        case ShapeGraphSpec::Kind::ExactPath:
            build_exact_path_graphs(_imp->graphs, pattern_size, next_pattern_supplemental, max_graphs, _imp->params.number_of_exact_path_graphs, _imp->graphs.directed, false, true);
            build_exact_path_graphs(_imp->graphs, target_size, next_target_supplemental, max_graphs, _imp->params.number_of_exact_path_graphs, _imp->graphs.directed, false, false);
            if (_imp->proof) {
                // exact-path graph g (1..number_of_exact_path_graphs) was built into slot
                // base + g - 1; pair each index with its slot for the derivation.
                unsigned base = next_pattern_supplemental - _imp->params.number_of_exact_path_graphs;
                vector<pair<int, unsigned>> exact_path_index_and_slot;
                for (int g = 1; g <= _imp->params.number_of_exact_path_graphs; ++g)
                    exact_path_index_and_slot.emplace_back(g, base + g - 1);
                _imp->proofs->prove_exact_path_graphs(_imp->graphs, max_graphs, exact_path_index_and_slot, base);
            }
            break;

        case ShapeGraphSpec::Kind::Distance2:
            build_exact_path_graphs(_imp->graphs, pattern_size, next_pattern_supplemental, max_graphs, 1, _imp->graphs.directed, true, true);
            build_exact_path_graphs(_imp->graphs, target_size, next_target_supplemental, max_graphs, 1, _imp->graphs.directed, true, false);
            break;

        case ShapeGraphSpec::Kind::Distance3:
            build_distance3_graphs(_imp->graphs, pattern_size, next_pattern_supplemental, max_graphs, true);
            build_distance3_graphs(_imp->graphs, target_size, next_target_supplemental, max_graphs, false);
            if (_imp->proof)
                _imp->proofs->prove_distance3_graphs(_imp->graphs, max_graphs, next_pattern_supplemental - 1);
            break;

        case ShapeGraphSpec::Kind::K4:
            build_k4_graphs(_imp->graphs, pattern_size, next_pattern_supplemental, max_graphs, true);
            build_k4_graphs(_imp->graphs, target_size, next_target_supplemental, max_graphs, false);
            break;

        case ShapeGraphSpec::Kind::ExtraShape: {
            auto & [shape, injective, count] = *spec.extra_shape;
            build_extra_shape(_imp->graphs, pattern_size, next_pattern_supplemental, max_graphs, _imp->params, *shape, injective, count, true);
            build_extra_shape(_imp->graphs, target_size, next_target_supplemental, max_graphs, _imp->params, *shape, injective, count, false);
            if (_imp->proof)
                _imp->proofs->prove_extra_shape(_imp->graphs, max_graphs, next_pattern_supplemental - 1);
        } break;
        }
    }

    // Defensive invariant: every plan slot has now been built, so the bump counters must
    // land on max_graphs (derived from the same plan) and match the recorded names.
    if (next_pattern_supplemental != max_graphs || next_target_supplemental != max_graphs ||
            next_pattern_supplemental != _imp->graphs.supplemental_graph_names.size())
        throw UnsupportedConfiguration{"something has gone wrong with supplemental graph indexing: " + to_string(next_pattern_supplemental) + " " + to_string(next_target_supplemental) + " " + to_string(max_graphs) + " "
        + to_string(_imp->graphs.supplemental_graph_names.size())};

    // pattern and target degrees, for supplemental graphs
    for (unsigned g = 1; g < max_graphs; ++g) {
        _imp->graphs.patterns_degrees.at(g).resize(pattern_size);
        _imp->graphs.targets_degrees.at(g).resize(target_size);
    }

    for (unsigned g = 1; g < max_graphs; ++g) {
        for (unsigned i = 0; i < pattern_size; ++i)
            _imp->graphs.patterns_degrees.at(g).at(i) = _imp->graphs.pattern_graph_rows[i * max_graphs + g].count();

        for (unsigned i = 0; i < target_size; ++i)
            _imp->graphs.targets_degrees.at(g).at(i) = _imp->graphs.target_graph_rows[i * max_graphs + g].count();
    }

    for (unsigned i = 0; i < target_size; ++i)
        _imp->graphs.largest_target_degree = max(_imp->graphs.largest_target_degree, _imp->graphs.targets_degrees[0][i]);

    // re-add loops
    for (unsigned i = 0; i < pattern_size; ++i)
        if (_imp->graphs.pattern_loops[i])
            _imp->graphs.pattern_graph_rows[i * max_graphs + 0].set(i);

    for (unsigned i = 0; i < target_size; ++i)
        if (_imp->graphs.target_loops[i])
            _imp->graphs.target_graph_rows[i * max_graphs + 0].set(i);

    // pattern adjacencies, compressed
    _imp->graphs.pattern_adjacencies_bits.resize(pattern_size * pattern_size);
    for (unsigned g = 0; g < max_graphs; ++g)
        for (unsigned i = 0; i < pattern_size; ++i)
            for (unsigned j = 0; j < pattern_size; ++j)
                if (_imp->graphs.pattern_graph_rows[i * max_graphs + g].test(j))
                    _imp->graphs.pattern_adjacencies_bits[i * pattern_size + j] |= (1u << g);

    return true;
}

auto HomomorphismModel::pattern_adjacency_bits(int p, int q) const -> PatternAdjacencyBitsType
{
    return _imp->graphs.pattern_adjacencies_bits[pattern_size * p + q];
}

auto HomomorphismModel::pattern_graph_row(int g, int p) const -> const SVOBitset &
{
    return _imp->graphs.pattern_graph_rows[p * max_graphs + g];
}

auto HomomorphismModel::target_graph_row(int g, int t) const -> const SVOBitset &
{
    return _imp->graphs.target_graph_rows[t * max_graphs + g];
}

auto HomomorphismModel::forward_target_graph_row(int t) const -> const SVOBitset &
{
    return _imp->graphs.forward_target_graph_rows[t];
}

auto HomomorphismModel::reverse_target_graph_row(int t) const -> const SVOBitset &
{
    return _imp->graphs.reverse_target_graph_rows[t];
}

auto HomomorphismModel::pattern_degree(int g, int p) const -> unsigned
{
    return _imp->graphs.patterns_degrees[g][p];
}

auto HomomorphismModel::target_degree(int g, int t) const -> unsigned
{
    return _imp->graphs.targets_degrees[g][t];
}

auto HomomorphismModel::largest_target_degree() const -> unsigned
{
    return _imp->graphs.largest_target_degree;
}

auto HomomorphismModel::has_vertex_labels() const -> bool
{
    return ! _imp->graphs.pattern_vertex_labels.empty();
}

auto HomomorphismModel::has_edge_labels() const -> bool
{
    return ! _imp->graphs.pattern_edge_labels.empty();
}

auto HomomorphismModel::pattern_vertex_label(int p) const -> int
{
    return _imp->graphs.pattern_vertex_labels[p];
}

auto HomomorphismModel::target_vertex_label(int t) const -> int
{
    return _imp->graphs.target_vertex_labels[t];
}

auto HomomorphismModel::pattern_edge_label(int p, int q) const -> int
{
    return _imp->graphs.pattern_edge_labels[p * pattern_size + q];
}

auto HomomorphismModel::target_edge_label(int t, int u) const -> int
{
    return _imp->graphs.target_edge_labels[t * target_size + u];
}

auto HomomorphismModel::pattern_has_loop(int p) const -> bool
{
    return _imp->graphs.pattern_loops[p];
}

auto HomomorphismModel::target_has_loop(int t) const -> bool
{
    return _imp->graphs.target_loops[t];
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
    return _imp->graphs.directed;
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

    if (! _imp->clique_data.pattern_cliques_sizes.empty()) {
        x.emplace_back(join("pattern_cliques_build_times =", _imp->clique_data.pattern_cliques_build_times));
        x.emplace_back(join("pattern_cliques_solve_times =", _imp->clique_data.pattern_cliques_solve_times));
        x.emplace_back(join("pattern_cliques_solve_find_nodes =", _imp->clique_data.pattern_cliques_solve_find_nodes));
        x.emplace_back(join("pattern_cliques_solve_prove_nodes =", _imp->clique_data.pattern_cliques_solve_prove_nodes));

        x.emplace_back(join("target_cliques_build_times =", _imp->clique_data.target_cliques_build_times));
        x.emplace_back(join("target_cliques_solve_times =", _imp->clique_data.target_cliques_solve_times));
        x.emplace_back(join("target_cliques_solve_find_nodes =", _imp->clique_data.target_cliques_solve_find_nodes));
        x.emplace_back(join("target_cliques_solve_prove_nodes =", _imp->clique_data.target_cliques_solve_prove_nodes));
    }

    x.emplace_back(join("supplemental_graph_names =", _imp->graphs.supplemental_graph_names));
}
