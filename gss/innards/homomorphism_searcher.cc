#include <gss/innards/cheap_all_different.hh>
#include <gss/innards/homomorphism_searcher.hh>

#include <optional>
#include <tuple>
#include <type_traits>

using namespace gss;
using namespace gss::innards;

using std::conditional_t;
using std::make_optional;
using std::max;
using std::move;
using std::mt19937;
using std::nullopt;
using std::numeric_limits;
using std::optional;
using std::pair;
using std::string;
using std::swap;
using std::to_string;
using std::tuple;
using std::uniform_int_distribution;
using std::vector;

HomomorphismSearcher::HomomorphismSearcher(const HomomorphismModel & m, const HomomorphismParams & p,
    const DuplicateSolutionFilterer & d, const std::shared_ptr<Proof> & f) :
    model(m),
    params(p),
    _duplicate_solution_filterer(d),
    _phase_size(0),
    proof(f)
{
    if (might_have_watches(params)) {
        watches.table.target_size = model.target_size;
        watches.table.data.resize(model.pattern_size * model.target_size);
    }

    for (unsigned v = 0 ; v < m.pattern_size ; ++v)
        _branch_scores.push_back(m.pattern_degree(0, v));

    if (params.phase_saving)
        _phases.resize(model.pattern_size);
}

auto HomomorphismSearcher::assignments_as_proof_decisions(const HomomorphismAssignments & assignments) const -> vector<pair<int, int>>
{
    vector<pair<int, int>> trail;
    for (auto & a : assignments.values)
        if (a.is_decision)
            trail.emplace_back(a.assignment.pattern_vertex, a.assignment.target_vertex);
    return trail;
}

auto HomomorphismSearcher::solution_in_proof_form(const HomomorphismAssignments & assignments) const -> vector<pair<NamedVertex, NamedVertex>>
{
    vector<pair<NamedVertex, NamedVertex>> solution;
    for (auto & a : assignments.values)
        if (solution.end() == find_if(solution.begin(), solution.end(), [&](const auto & t) { return unsigned(t.first.first) == a.assignment.pattern_vertex; }))
            solution.emplace_back(
                model.pattern_vertex_for_proof(a.assignment.pattern_vertex),
                model.target_vertex_for_proof(a.assignment.target_vertex));
    return solution;
}

auto HomomorphismSearcher::expand_to_full_result(const HomomorphismAssignments & assignments, VertexToVertexMapping & mapping) -> void
{
    for (auto & a : assignments.values)
        mapping.emplace(a.assignment.pattern_vertex, a.assignment.target_vertex);
}

auto HomomorphismSearcher::save_result(const HomomorphismAssignments & assignments, HomomorphismResult & result) -> void
{
    expand_to_full_result(assignments, result.mapping);

    string where = "where =";
    for (auto & a : assignments.values)
        where.append(" " + to_string(a.discrepancy_count) + "/" + to_string(a.choice_count));
    result.extra_stats.push_back(where);
}

auto HomomorphismSearcher::restarting_search(
    HomomorphismAssignments & assignments,
    const Domains & domains,
    unsigned long long & nodes,
    unsigned long long & propagations,
    loooong & solution_count,
    int depth,
    RestartsSchedule & restarts_schedule) -> SearchResult
{
    if (proof && proof->super_extra_verbose()) {
        vector<pair<NamedVertex, vector<NamedVertex>>> proof_domains;
        for (auto & d : domains) {
            proof_domains.push_back(pair{model.pattern_vertex_for_proof(d.v), vector<NamedVertex>{}});
            auto values = d.values;
            for (auto v = values.find_first(); v != decltype(values)::npos; v = values.find_first()) {
                values.reset(v);
                proof_domains.back().second.push_back(model.target_vertex_for_proof(v));
            }
        }
        proof->show_domains("entering depth " + to_string(depth), proof_domains);
    }

    if (params.timeout->should_abort())
        return SearchResult::Aborted;

    ++nodes;

    // find ourselves a domain, or succeed if we're all assigned
    const HomomorphismDomain * branch_domain = find_branch_domain(domains);
    if (! branch_domain) {
        if (params.lackey) {
            VertexToVertexMapping mapping;
            expand_to_full_result(assignments, mapping);
            if (! params.lackey->check_solution(mapping, false, params.count_solutions, {})) {
                switch (params.propagate_using_lackey) {
                case PropagateUsingLackey::RootAndBackjump:
                    return SearchResult::UnsatisfiableAndBackjumpUsingLackey;
                case PropagateUsingLackey::Never:
                case PropagateUsingLackey::Root:
                case PropagateUsingLackey::Always:
                    return SearchResult::Unsatisfiable;
                }
            }
        }

        if (proof)
            proof->post_solution(solution_in_proof_form(assignments));

        if (params.count_solutions) {
            // we could be finding duplicate solutions, in threaded search
            if (_duplicate_solution_filterer(assignments)) {
                ++solution_count;
                if (params.enumerate_callback) {
                    VertexToVertexMapping mapping;
                    expand_to_full_result(assignments, mapping);
                    if (! params.enumerate_callback(mapping))
                        return SearchResult::Satisfiable;
                }
            }

            return SearchResult::SatisfiableButKeepGoing;
        }
        else
            return SearchResult::Satisfiable;
    }

    // pull out the remaining values in this domain for branching
    auto remaining = branch_domain->values;

    vector<int> branch_v(model.target_size);

    unsigned branch_v_start = 0, branch_v_end = 0;
    optional<unsigned> save = params.phase_saving ? _phases.at(branch_domain->v) : nullopt;
    for (auto f_v = remaining.find_first(); f_v != decltype(remaining)::npos; f_v = remaining.find_first()) {
        remaining.reset(f_v);
        branch_v[branch_v_end] = f_v;
        if (save && *save == f_v) {
            swap(branch_v[branch_v_start], branch_v[branch_v_end]);
            ++branch_v_start;
        }
        ++branch_v_end;
    }

    switch (params.value_ordering_heuristic) {
    case ValueOrdering::None:
        break;

    case ValueOrdering::Degree:
        degree_sort(branch_v, branch_v_start, branch_v_end, false);
        break;

    case ValueOrdering::AntiDegree:
        degree_sort(branch_v, branch_v_start, branch_v_end, true);
        break;

    case ValueOrdering::Biased:
        softmax_shuffle(branch_v, branch_v_start, branch_v_end);
        break;

    case ValueOrdering::Random:
        shuffle(branch_v.begin() + branch_v_start, branch_v.begin() + branch_v_end, global_rand);
        break;
    }

    int discrepancy_count = 0;
    bool actually_hit_a_failure = false;

    // override whether we use the lackey for propagation, in case we are inside a backjump
    bool use_lackey_for_propagation = false;

    // for each value remaining...
    for (auto f_v = branch_v.begin(), f_end = branch_v.begin() + branch_v_end; f_v != f_end; ++f_v) {
        if (proof)
            proof->guessing(depth, model.pattern_vertex_for_proof(branch_domain->v), model.target_vertex_for_proof(*f_v));

        // modified in-place by appending, we can restore by shrinking
        auto assignments_size = assignments.values.size();

        // make the assignment
        assignments.values.push_back({{branch_domain->v, unsigned(*f_v)}, true, discrepancy_count, int(branch_v_end)});

        // phase saving?
        if (params.phase_saving) {
            auto how_many = count_if(assignments.values.begin(), assignments.values.end(), [&](const auto & a) { return ! a.is_decision; });
            if (how_many > _phase_size) {
                _phase_size = how_many;
                fill(_phases.begin(), _phases.end(), nullopt);
                for (const auto & a : assignments.values)
                    if (a.is_decision)
                        _phases.at(a.assignment.pattern_vertex) = a.assignment.target_vertex;
            }
        }

        // set up new domains
        Domains new_domains = copy_nonfixed_domains_and_make_assignment(domains, branch_domain->v, *f_v);

        // propagate
        ++propagations;
        optional<FailingVertex> failing_vertex;
        if (nullopt != (failing_vertex = propagate(false, new_domains, assignments, use_lackey_for_propagation || (params.propagate_using_lackey == PropagateUsingLackey::Always)))) {
            // failure? restore assignments and go on to the next thing
            if (proof)
                proof->propagation_failure(assignments_as_proof_decisions(assignments), model.pattern_vertex_for_proof(branch_domain->v), model.target_vertex_for_proof(*f_v));

            assignments.values.resize(assignments_size);
            actually_hit_a_failure = true;

            if (params.variable_ordering_heuristic == VariableOrdering::DomOverWDeg || params.variable_ordering_heuristic == VariableOrdering::DomThenWDeg) {
                if (failing_vertex->pattern_vertex_1 && assignments.values.end() == find_if(
                            assignments.values.begin(), assignments.values.end(), [&] (const auto & a) { return a.assignment.pattern_vertex == *failing_vertex->pattern_vertex_1; })) {
                    auto & v = _branch_scores.at(*failing_vertex->pattern_vertex_1);
                    if (++v > (1 << 24)) {
                        for (auto & w : _branch_scores)
                            w /= 2;
                    }
                }
                if (failing_vertex->pattern_vertex_2 && assignments.values.end() == find_if(
                            assignments.values.begin(), assignments.values.end(), [&] (const auto & a) { return a.assignment.pattern_vertex == *failing_vertex->pattern_vertex_2; })) {
                    auto & v = _branch_scores.at(*failing_vertex->pattern_vertex_2);
                    if (++v > (1 << 24)) {
                        for (auto & w : _branch_scores)
                            w /= 2;
                    }
                }
            }

            continue;
        }

        if (proof)
            proof->start_level(depth + 2);

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
            for (auto l = branch_v.begin(); l != f_v; ++l) {
                assignments.values.push_back({{branch_domain->v, unsigned(*l)}, true, -2, -2});
                post_nogood(assignments);
                assignments.values.pop_back();
            }

            return SearchResult::Restart;

        case SearchResult::SatisfiableButKeepGoing:
            if (proof) {
                proof->back_up_to_level(depth + 1);
                proof->incorrect_guess(assignments_as_proof_decisions(assignments), false);
                proof->forget_level(depth + 2);
            }

            // restore assignments
            assignments.values.resize(assignments_size);
            break;

        case SearchResult::UnsatisfiableAndBackjumpUsingLackey:
            use_lackey_for_propagation = true;
            [[std::fallthrough]];

        case SearchResult::Unsatisfiable:
            if (proof) {
                proof->back_up_to_level(depth + 1);
                proof->incorrect_guess(assignments_as_proof_decisions(assignments), true);
                proof->forget_level(depth + 2);
            }

            // restore assignments
            assignments.values.resize(assignments_size);
            actually_hit_a_failure = true;
            break;
        }

        ++discrepancy_count;
    }

    // no values remaining, backtrack, or possibly kick off a restart
    if (proof)
        proof->out_of_guesses(assignments_as_proof_decisions(assignments));

    if (actually_hit_a_failure)
        restarts_schedule.did_a_backtrack();

    if (restarts_schedule.should_restart()) {
        if (proof)
            proof->back_up_to_top();
        post_nogood(assignments);
        return SearchResult::Restart;
    }
    else
        return use_lackey_for_propagation ? SearchResult::UnsatisfiableAndBackjumpUsingLackey : SearchResult::Unsatisfiable;
}

auto HomomorphismSearcher::degree_sort(
    vector<int> & branch_v,
    unsigned branch_v_start,
    unsigned branch_v_end,
    bool reverse) -> void
{
    stable_sort(branch_v.begin() + branch_v_start, branch_v.begin() + branch_v_end, [&](int a, int b) -> bool {
        return reverse
            ? model.target_degree(0, a) < model.target_degree(0, b)
            : -model.target_degree(0, a) < -model.target_degree(0, b);
    });
}

auto HomomorphismSearcher::softmax_shuffle(
    vector<int> & branch_v,
    unsigned branch_v_start,
    unsigned branch_v_end) -> void
{
    // repeatedly pick a softmax-biased vertex, move it to the front of branch_v,
    // and then only consider items further to the right in the next iteration.

    // Using floating point here turned out to be way too slow. Fortunately the base
    // of softmax doesn't seem to matter, so we use 2 instead of e, and do everything
    // using bit voodoo.
    auto expish = [largest_target_degree = this->model.largest_target_degree()](int degree) {
        constexpr int sufficient_space_for_adding_up = numeric_limits<long long>::digits - 18;
        auto shift = max<int>(degree - largest_target_degree + sufficient_space_for_adding_up, 0);
        return 1ll << shift;
    };

    long long total = 0;
    for (unsigned v = 0; v < branch_v_end; ++v)
        total += expish(model.target_degree(0, branch_v[v]));

    for (unsigned start = branch_v_start; start < branch_v_end; ++start) {
        // pick a random number between 1 and total inclusive
        uniform_int_distribution<long long> dist(1, total);
        long long select_score = dist(global_rand);

        // go over the list until we hit the score
        unsigned select_element = start;
        for (; select_element + 1 < branch_v_end; ++select_element) {
            select_score -= expish(model.target_degree(0, branch_v[select_element]));
            if (select_score <= 0)
                break;
        }

        // move to front
        total -= expish(model.target_degree(0, branch_v[select_element]));
        swap(branch_v[select_element], branch_v[start]);
    }
}

auto HomomorphismSearcher::post_nogood(const HomomorphismAssignments & assignments) -> void
{
    if (! might_have_watches(params))
        return;

    Nogood<HomomorphismAssignment> nogood;

    for (auto & a : assignments.values)
        if (a.is_decision)
            nogood.literals.emplace_back(a.assignment);

    watches.post_nogood(move(nogood));

    if (proof)
        proof->post_restart_nogood(assignments_as_proof_decisions(assignments));
}

auto HomomorphismSearcher::copy_nonfixed_domains_and_make_assignment(
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

auto HomomorphismSearcher::find_branch_domain(const Domains & domains) -> const HomomorphismDomain *
{
    const HomomorphismDomain * result = nullptr;

    switch (params.variable_ordering_heuristic) {
    case VariableOrdering::DomThenDeg:
    case VariableOrdering::DomThenWDeg:
        for (auto & d : domains)
            if (! d.fixed)
                if ((! result) ||
                    (d.count < result->count) ||
                    (d.count == result->count && branch_score(d.v) > branch_score(result->v)))
                    result = &d;
        break;

    case VariableOrdering::DomOverDeg:
    case VariableOrdering::DomOverWDeg:
        for (auto & d : domains)
            if (! d.fixed)
                if ((! result) ||
                    (d.count < result->count) ||
                    (d.count / (1.0 + branch_score(d.v)) < (result->count / (1.0 + branch_score(result->v)))))
                    result = &d;
        break;
    }

    return result;
}

template <bool directed_, bool has_edge_labels_, bool induced_, bool verbose_proofs_>
auto HomomorphismSearcher::propagate_adjacency_constraints(HomomorphismDomain & d, const HomomorphismAssignment & current_assignment) -> void
{
    const auto & graph_pairs_to_consider = model.pattern_adjacency_bits(current_assignment.pattern_vertex, d.v);

    [[maybe_unused]] conditional_t<verbose_proofs_, SVOBitset, tuple<>> before;
    if constexpr (verbose_proofs_) {
        before = d.values;
    }

    if constexpr (! directed_) {
        // for the original graph pair, if we're adjacent...
        if (graph_pairs_to_consider & (1u << 0)) {
            // ...then we can only be mapped to adjacent vertices
            d.values &= model.target_graph_row(0, current_assignment.target_vertex);
        }
        else {
            if constexpr (induced_) {
                // ...otherwise we can only be mapped to adjacent vertices
                d.values.intersect_with_complement(model.target_graph_row(0, current_assignment.target_vertex));
            }
        }
    }
    else {
        // both forward and reverse edges to consider
        if (graph_pairs_to_consider & (1u << 0)) {
            // ...then we can only be mapped to adjacent vertices
            d.values &= model.forward_target_graph_row(current_assignment.target_vertex);
        }
        else {
            if constexpr (induced_) {
                // ...otherwise we can only be mapped to adjacent vertices
                d.values.intersect_with_complement(model.forward_target_graph_row(current_assignment.target_vertex));
            }
        }

        const auto & reverse_edge_graph_pairs_to_consider = model.pattern_adjacency_bits(d.v, current_assignment.pattern_vertex);

        if (reverse_edge_graph_pairs_to_consider & (1u << 0)) {
            // ...then we can only be mapped to adjacent vertices
            d.values &= model.reverse_target_graph_row(current_assignment.target_vertex);
        }
        else {
            if constexpr (induced_) {
                // ...otherwise we can only be mapped to adjacent vertices
                d.values.intersect_with_complement(model.reverse_target_graph_row(current_assignment.target_vertex));
            }
        }
    }

    if constexpr (verbose_proofs_) {
        if (before.count() != d.values.count())
            proof->propagated(model.pattern_vertex_for_proof(current_assignment.pattern_vertex), model.target_vertex_for_proof(current_assignment.target_vertex),
                0, before.count() - d.values.count(), model.pattern_vertex_for_proof(d.v));
        before = d.values;
    }

    // and for each remaining graph pair...
    for (unsigned g = 1; g < model.max_graphs; ++g) {
        // if we're adjacent...
        if (graph_pairs_to_consider & (1u << g)) {
            // ...then we can only be mapped to adjacent vertices
            d.values &= model.target_graph_row(g, current_assignment.target_vertex);
        }

        if constexpr (verbose_proofs_) {
            if (before.count() != d.values.count())
                proof->propagated(model.pattern_vertex_for_proof(current_assignment.pattern_vertex), model.target_vertex_for_proof(current_assignment.target_vertex),
                    g, before.count() - d.values.count(), model.pattern_vertex_for_proof(d.v));
            before = d.values;
        }
    }

    if constexpr (has_edge_labels_) {
        // if we're adjacent in the original graph, additionally the edge labels need to match up
        if (graph_pairs_to_consider & (1u << 0)) {
            auto check_d_values = d.values;

            auto want_forward_label = model.pattern_edge_label(current_assignment.pattern_vertex, d.v);
            for (auto c = check_d_values.find_first(); c != decltype(check_d_values)::npos; c = check_d_values.find_first()) {
                check_d_values.reset(c);

                auto got_forward_label = model.target_edge_label(current_assignment.target_vertex, c);
                if (got_forward_label != want_forward_label)
                    d.values.reset(c);
            }
        }

        const auto & reverse_edge_graph_pairs_to_consider = model.pattern_adjacency_bits(d.v, current_assignment.pattern_vertex);
        if (reverse_edge_graph_pairs_to_consider & (1u << 0)) {
            auto check_d_values = d.values;

            auto want_reverse_label = model.pattern_edge_label(d.v, current_assignment.pattern_vertex);
            for (auto c = check_d_values.find_first(); c != decltype(check_d_values)::npos; c = check_d_values.find_first()) {
                check_d_values.reset(c);

                auto got_reverse_label = model.target_edge_label(c, current_assignment.target_vertex);
                if (got_reverse_label != want_reverse_label)
                    d.values.reset(c);
            }
        }
    }
}

auto HomomorphismSearcher::both_in_the_neighbourhood_of_some_vertex(unsigned v, unsigned w) -> bool
{
    auto i = model.pattern_graph_row(0, v);
    i &= model.pattern_graph_row(0, w);
    return i.any();
}

auto HomomorphismSearcher::propagate_simple_constraints(Domains & new_domains, const HomomorphismAssignment & current_assignment) -> optional<FailingVertex>
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
        if (! model.has_edge_labels()) {
            if (params.induced) {
                if (model.directed()) {
                    if ((! proof) || (! proof->super_extra_verbose()))
                        propagate_adjacency_constraints<true, false, true, false>(d, current_assignment);
                    else
                        propagate_adjacency_constraints<true, false, true, true>(d, current_assignment);
                }
                else {
                    if ((! proof) || (! proof->super_extra_verbose()))
                        propagate_adjacency_constraints<false, false, true, false>(d, current_assignment);
                    else
                        propagate_adjacency_constraints<false, false, true, true>(d, current_assignment);
                }
            }
            else {
                if (model.directed()) {
                    if ((! proof) || (! proof->super_extra_verbose()))
                        propagate_adjacency_constraints<true, false, false, false>(d, current_assignment);
                    else
                        propagate_adjacency_constraints<true, false, false, true>(d, current_assignment);
                }
                else {
                    if ((! proof) || (! proof->super_extra_verbose()))
                        propagate_adjacency_constraints<false, false, false, false>(d, current_assignment);
                    else
                        propagate_adjacency_constraints<false, false, false, true>(d, current_assignment);
                }
            }
        }
        else {
            // edge labels are always directed
            if (params.induced) {
                if ((! proof) || (! proof->super_extra_verbose()))
                    propagate_adjacency_constraints<true, true, true, false>(d, current_assignment);
                else
                    propagate_adjacency_constraints<true, true, true, true>(d, current_assignment);
            }
            else {
                if ((! proof) || (! proof->super_extra_verbose()))
                    propagate_adjacency_constraints<true, true, false, false>(d, current_assignment);
                else
                    propagate_adjacency_constraints<true, true, false, true>(d, current_assignment);
            }
        }

        // we might have removed values
        d.count = d.values.count();
        if (0 == d.count)
            return FailingVertex{current_assignment.pattern_vertex, d.v};
    }

    return nullopt;
}

auto HomomorphismSearcher::propagate_less_thans(Domains & new_domains) -> optional<FailingVertex>
{
    vector<int> find_domain(model.pattern_size, -1);

    for (unsigned i = 0, i_end = new_domains.size(); i != i_end; ++i)
        find_domain[new_domains[i].v] = i;

    for (auto & [a, b] : model.pattern_less_thans_in_convenient_order) {
        if (find_domain[a] == -1 || find_domain[b] == -1)
            continue;
        auto & a_domain = new_domains[find_domain[a]];
        auto & b_domain = new_domains[find_domain[b]];

        // first value of b must be at least one after the first possible value of a
        auto first_a = a_domain.values.find_first();
        if (first_a == decltype(a_domain.values)::npos)
            return FailingVertex{a_domain.v, b_domain.v};
        auto first_allowed_b = first_a + 1;

        if (first_allowed_b >= model.target_size)
            return FailingVertex{a_domain.v, b_domain.v};

        for (auto v = b_domain.values.find_first(); v != decltype(b_domain.values)::npos; v = b_domain.values.find_first()) {
            if (v >= first_allowed_b)
                break;
            b_domain.values.reset(v);
        }

        // b might have shrunk (and detect empty before the next bit to make life easier)
        b_domain.count = b_domain.values.count();
        if (0 == b_domain.count)
            return FailingVertex{a_domain.v, b_domain.v};
    }

    for (auto & [a, b] : model.pattern_less_thans_in_convenient_order) {
        if (find_domain[a] == -1 || find_domain[b] == -1)
            continue;
        auto & a_domain = new_domains[find_domain[a]];
        auto & b_domain = new_domains[find_domain[b]];

        // last value of a must be at least one before the last possible value of b
        auto b_values_copy = b_domain.values;
        auto last_b = b_domain.values.find_first();
        for (auto v = last_b; v != decltype(b_values_copy)::npos; v = b_values_copy.find_first()) {
            b_values_copy.reset(v);
            last_b = v;
        }

        if (last_b == 0)
            return FailingVertex{a_domain.v, b_domain.v};
        auto last_allowed_a = last_b - 1;

        auto a_values_copy = a_domain.values;
        for (auto v = a_values_copy.find_first(); v != decltype(a_values_copy)::npos; v = a_values_copy.find_first()) {
            a_values_copy.reset(v);
            if (v > last_allowed_a)
                a_domain.values.reset(v);
        }

        // a might have shrunk
        a_domain.count = a_domain.values.count();
        if (0 == a_domain.count)
            return FailingVertex{a_domain.v, b_domain.v};
    }

    return nullopt;
}

auto HomomorphismSearcher::propagate_occur_less_thans(
    const optional<HomomorphismAssignment> & current_assignment,
    const HomomorphismAssignments & assignments,
    Domains & new_domains) -> optional<FailingVertex>
{
    vector<optional<SVOBitset>> occurs(model.target_size);

    auto build_occurs = [&](int p) -> void {
        if (occurs[p])
            return;

        occurs[p] = make_optional<SVOBitset>(model.pattern_size, 0);
        for (auto & d : new_domains)
            if (d.values.test(p))
                occurs[p]->set(d.v);
    };

    for (auto & [a, b] : model.target_occur_less_thans_in_convenient_order) {
        build_occurs(a);
        build_occurs(b);
    }

    for (auto & a : assignments.values)
        if (occurs[a.assignment.target_vertex])
            occurs[a.assignment.target_vertex]->set(a.assignment.pattern_vertex);

    // propagate lower bounds
    for (auto & [a, b] : model.target_occur_less_thans_in_convenient_order) {
        auto first_a = occurs[a]->find_first();
        if (first_a == SVOBitset::npos) {
            // no occurrence of value a, value b cannot be used either
            occurs[b]->reset();
            for (auto & d : new_domains)
                if (d.values.test(b)) {
                    d.values.reset(b);
                    if (0 == --d.count)
                        return FailingVertex{d.v};
                }
        }
        else {
            // value a first occurs in variable x, value b cannot be used in a variable lower than x
            for (auto & d : new_domains) {
                if (d.v < first_a && d.values.test(b)) {
                    occurs[b]->reset(d.v);
                    d.values.reset(b);
                    if (0 == --d.count)
                        return FailingVertex{d.v};
                }
            }
        }
    }

    // propagate other way: if value b must occur (because it has been assigned) then
    // value a must go before
    if (current_assignment) {
        for (auto & [a, b] : model.target_occur_less_thans_in_convenient_order) {
            if (b != current_assignment->target_vertex)
                continue;

            bool saw_an_a = false;
            for (auto & d : new_domains) {
                if (d.v < current_assignment->pattern_vertex) {
                    // it's before
                    if (d.values.test(a))
                        saw_an_a = true;
                }
                else if (d.v > current_assignment->pattern_vertex) {
                    // comes after, can't use a
                    if (d.values.test(a)) {
                        occurs[a]->reset(d.v);
                        d.values.reset(a);
                        if (0 == --d.count)
                            return FailingVertex{d.v};
                    }
                }
            }

            for (auto & d : assignments.values)
                if (d.assignment.pattern_vertex < current_assignment->pattern_vertex && a == d.assignment.target_vertex)
                    saw_an_a = true;

            if (! saw_an_a)
                return FailingVertex{current_assignment->pattern_vertex};
        }
    }

    return nullopt;
}

auto HomomorphismSearcher::propagate(bool initial, Domains & new_domains, HomomorphismAssignments & assignments, bool propagate_using_lackey) -> optional<FailingVertex>
{
    optional<FailingVertex> wipeout = nullopt;

    // nogoods might be watching things in initial assignments. this is possibly not the
    // best place to put this...
    if (initial && might_have_watches(params)) {
        for (auto & a : assignments.values) {
            HomomorphismAssignment current_assignment = {a.assignment.pattern_vertex, a.assignment.target_vertex};
            watches.propagate(
                current_assignment,
                [&](const HomomorphismAssignment & a) { return ! assignments.contains(a); },
                [&](const HomomorphismAssignment & a) {
                    for (auto & d : new_domains) {
                        if (d.v == a.pattern_vertex) {
                            if (d.values.test(a.target_vertex)) {
                                d.values.reset(a.target_vertex);
                                if (0 == --d.count)
                                    wipeout = FailingVertex{NoIdentifiableCause{}};
                            }
                            break;
                        }
                    }
                });

            if (wipeout)
                return wipeout;
        }
    }

    auto find_unit_domain = [&]() {
        return find_if(new_domains.begin(), new_domains.end(), [](HomomorphismDomain & d) {
            return (! d.fixed) && 1 == d.count;
        });
    };

    bool done_globals_at_least_once = false;

    // whilst we've got a unit domain...
    for (typename Domains::iterator branch_domain = find_unit_domain();
         branch_domain != new_domains.end() || ! done_globals_at_least_once;
         branch_domain = find_unit_domain()) {
        optional<HomomorphismAssignment> current_assignment;
        if (branch_domain != new_domains.end()) {
            // what are we assigning?
            current_assignment = HomomorphismAssignment{branch_domain->v, unsigned(branch_domain->values.find_first())};

            // ok, make the assignment
            branch_domain->fixed = true;
            assignments.values.push_back({*current_assignment, false, -1, -1});

            if (proof)
                proof->unit_propagating(
                    model.pattern_vertex_for_proof(current_assignment->pattern_vertex),
                    model.target_vertex_for_proof(current_assignment->target_vertex));

            // propagate watches
            if (might_have_watches(params)) {
                watches.propagate(
                    *current_assignment,
                    [&](const HomomorphismAssignment & a) { return ! assignments.contains(a); },
                    [&](const HomomorphismAssignment & a) {
                        for (auto & d : new_domains) {
                            if (d.fixed)
                                continue;

                            if (d.v == a.pattern_vertex) {
                                if (d.values.test(a.target_vertex)) {
                                    d.values.reset(a.target_vertex);
                                    if (0 == --d.count)
                                        wipeout = FailingVertex{NoIdentifiableCause{}};
                                }
                                break;
                            }
                        }
                    });

                if (wipeout)
                    return wipeout;
            }

            // propagate simple all different and adjacency
            wipeout = propagate_simple_constraints(new_domains, *current_assignment);
            if (wipeout)
                return wipeout;
        }

        // propagate less thans
        if (model.has_less_thans()) {
            wipeout = propagate_less_thans(new_domains);
            if (wipeout)
                return wipeout;
        }

        if (model.has_occur_less_thans()) {
            wipeout = propagate_occur_less_thans(current_assignment, assignments, new_domains);
            if (wipeout)
                return wipeout;
        }

        // propagate all different
        if (params.injectivity == Injectivity::Injective) {
            wipeout = cheap_all_different(model.target_size, new_domains, proof, &model);
            if (wipeout)
                return wipeout;
        }
        done_globals_at_least_once = true;
    }

    int dcount = 0;
    if (params.lackey && (propagate_using_lackey || params.send_partials_to_lackey)) {
        VertexToVertexMapping mapping;
        expand_to_full_result(assignments, mapping);

        if (! propagate_using_lackey) {
            wipeout = params.lackey->check_solution(mapping, true, false, Lackey::DeletionFunction{});
            if (wipeout)
                return wipeout;
        }
        else {
            vector<int> find_domain(model.pattern_size, -1);
            for (unsigned i = 0; i < new_domains.size(); ++i)
                find_domain[new_domains[i].v] = i;

            auto deletion = [&](int p, int t) -> bool {
                if (! wipeout) {
                    if (int d = find_domain[p]; d != -1) {
                        if (new_domains[d].values.test(t)) {
                            ++dcount;
                            new_domains[d].values.reset(t);
                            if (0 == --new_domains[d].count)
                                wipeout = FailingVertex{NoIdentifiableCause{}};
                            return true;
                        }
                    }
                }
                return false;
            };

            wipeout = params.lackey->check_solution(mapping, true, false, deletion);

            if (wipeout)
                return wipeout;
        }
    }

    return wipeout;
}

auto HomomorphismSearcher::set_seed(int t) -> void
{
    global_rand.seed(t);
}

auto HomomorphismSearcher::branch_score(unsigned branch_v) -> int
{
    return _branch_scores.at(branch_v);
}
