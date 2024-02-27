#include <gss/innards/cheap_all_different.hh>
#include <gss/innards/homomorphism_searcher.hh>

#include <optional>
#include <tuple>
#include <type_traits>
#include <boost/range/adaptor/reversed.hpp>

using namespace gss;
using namespace gss::innards;

using std::conditional_t;
using std::make_optional;
using std::max;
using std::move;
using std::mt19937;
using std::numeric_limits;
using std::optional;
using std::pair;
using std::string;
using std::string_view;
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
    proof(f)
{
    if (might_have_watches(params)) {
        watches.table.target_size = model.target_size;
        watches.table.data.resize(model.pattern_size * model.target_size);
    }
    if (model.has_less_thans()) {
        symmetric_value_displacement.resize(model.target_size);
    }
    if (model.do_dynamic_occur_less_thans()) {
        std::vector<innards::SVOBitset> adjacency_matrix;
        for (size_t i = 0; i < model.target_size; i++) {
            adjacency_matrix.emplace_back(model.target_graph_row(0,i));
        }
        initialise_dynamic_structure(t_rschreier, adjacency_matrix);
    }
    if (model.do_dynamic_less_thans()) {
        std::vector<innards::SVOBitset> adjacency_matrix;
        for (size_t i = 0; i < model.pattern_size; i++) {
            adjacency_matrix.emplace_back(model.pattern_graph_row(0,i));
        }
        initialise_dynamic_structure(p_rschreier, adjacency_matrix);
    }
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
    Domains & domains,
    unsigned long long & nodes,
    unsigned long long & propagations,
    loooong & solution_count,
    int depth,
    RestartsSchedule & restarts_schedule,
    vector<int> & pattern_orbit_base,
    vector<int> & target_orbit_base) -> SearchResult
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

    unsigned branch_v_end = 0;
    for (auto f_v = remaining.find_first(); f_v != decltype(remaining)::npos; f_v = remaining.find_first()) {
        remaining.reset(f_v);
        branch_v[branch_v_end++] = f_v;
    }

    switch (params.value_ordering_heuristic) {
    case ValueOrdering::None:
        break;

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

    // override whether we use the lackey for propagation, in case we are inside a backjump
    bool use_lackey_for_propagation = false;

    vector<uint8_t> skip_target_vertex_due_to_orbits(model.target_size, 0);

    dejavu::groups::orbit pattern_orbit_partition{model.pattern_size};
    bool this_vertex_has_pattern_orbit = false;
    if (model.has_pattern_orbits) {
        if (pattern_orbit_base != model.pattern_orbit_base) {
            model.pattern_orbit_base = pattern_orbit_base;
            model.pattern_orbits_schreier->set_base(model.pattern_orbit_base);
        }
        pattern_orbit_base.push_back(branch_domain->v);
        model.pattern_orbit_base.push_back(branch_domain->v);
        model.pattern_orbits_schreier->get_stabilizer_orbit(pattern_orbit_base.size(), pattern_orbit_partition);
        if (pattern_orbit_partition.orbit_size(branch_domain->v) == 1) {
            pattern_orbit_base.pop_back();
            model.pattern_orbit_base.pop_back();
        }
        else
            this_vertex_has_pattern_orbit = true;
    }

    // for each value remaining...
    for (auto f_v = branch_v.begin(), f_end = branch_v.begin() + branch_v_end; f_v != f_end; ++f_v) {
        dejavu::groups::orbit target_orbit_partition{model.target_size};
        bool this_vertex_has_target_orbit = false;

        if (! skip_target_vertex_due_to_orbits.at(*f_v)) {
            if (proof)
                proof->guessing(depth, model.pattern_vertex_for_proof(branch_domain->v), model.target_vertex_for_proof(*f_v));

            // modified in-place by appending, we can restore by shrinking
            auto assignments_size = assignments.values.size();

            // make the assignment
            assignments.values.push_back({{branch_domain->v, unsigned(*f_v)}, true, discrepancy_count, int(branch_v_end)});

            // set up new domains
            Domains new_domains = copy_nonfixed_domains_and_make_assignment(domains, branch_domain->v, *f_v);

            if (model.has_target_orbits) {
                if (target_orbit_base != model.target_orbit_base) {
                    model.target_orbit_base = target_orbit_base;
                    model.target_orbits_schreier->set_base(model.target_orbit_base);
                }
                target_orbit_base.push_back(*f_v);
                model.target_orbit_base.push_back(*f_v);
                model.target_orbits_schreier->get_stabilizer_orbit(target_orbit_base.size(), target_orbit_partition);
                if (target_orbit_partition.orbit_size(*f_v) == 1) {
                    target_orbit_base.pop_back();
                    model.target_orbit_base.pop_back();
                }
                else
                    this_vertex_has_target_orbit = true;
            }

            // propagate
            ++propagations;
            if (! propagate(false, new_domains, assignments, use_lackey_for_propagation || (params.propagate_using_lackey == PropagateUsingLackey::Always))) {
                // failure? restore assignments and go on to the next thing
                if (proof)
                    proof->propagation_failure(assignments_as_proof_decisions(assignments), model.pattern_vertex_for_proof(branch_domain->v), model.target_vertex_for_proof(*f_v));

                assignments.values.resize(assignments_size);
                actually_hit_a_failure = true;
            }
            else {
                if (proof)
                    proof->start_level(depth + 2);

                // recursive search
                auto pattern_orbit_base_copy = pattern_orbit_base, target_orbit_base_copy = target_orbit_base;
                auto search_result = restarting_search(assignments, new_domains, nodes, propagations,
                    solution_count, depth + 1, restarts_schedule, pattern_orbit_base_copy, target_orbit_base_copy);

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
        }

        if (this_vertex_has_pattern_orbit) {
            for (auto & d : domains) {
                if (d.v != branch_domain->v && pattern_orbit_partition.are_in_same_orbit(d.v, branch_domain->v)) {
                    if (d.values.test(*f_v)) {
                        d.values.reset(*f_v);
                        if (0 == --d.count)
                            return SearchResult::Unsatisfiable;
                    }
                }
            }
        }

        if (this_vertex_has_target_orbit) {
            for (auto f_w = next(f_v); f_w != f_end; ++f_w)
                if (target_orbit_partition.are_in_same_orbit(*f_v, *f_w))
                    skip_target_vertex_due_to_orbits.at(*f_w) = 1;
        }
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
    unsigned branch_v_end,
    bool reverse) -> void
{
    stable_sort(branch_v.begin(), branch_v.begin() + branch_v_end, [&](int a, int b) -> bool {
        return reverse
            ? model.target_degree(0, a) < model.target_degree(0, b)
            : -model.target_degree(0, a) < -model.target_degree(0, b);
    });
}

auto HomomorphismSearcher::softmax_shuffle(
    vector<int> & branch_v,
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

    for (unsigned start = 0; start < branch_v_end; ++start) {
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
    for (auto & d : domains)
        if (! d.fixed)
            if ((! result) ||
                (d.count < result->count) ||
                (d.count == result->count && model.pattern_degree(0, d.v) > model.pattern_degree(0, result->v)))
                result = &d;
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

auto HomomorphismSearcher::propagate_simple_constraints(Domains & new_domains, const HomomorphismAssignment & current_assignment) -> bool
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
            return false;
    }

    return true;
}

auto HomomorphismSearcher::propagate_less_thans(Domains & new_domains) -> bool {
    return propagate_less_thans(new_domains, model.pattern_less_thans_in_convenient_order);
}

auto HomomorphismSearcher::propagate_less_thans(Domains & new_domains, const std::vector<std::pair<unsigned int, unsigned int>> & constraints) -> bool
{
    vector<int> find_domain(model.pattern_size, -1);

    for (unsigned i = 0, i_end = new_domains.size(); i != i_end; ++i)
        find_domain[new_domains[i].v] = i;

    for (auto & [a, b] : constraints) {
        if (find_domain[a] == -1 || find_domain[b] == -1)
            continue;
        auto & a_domain = new_domains[find_domain[a]];
        auto & b_domain = new_domains[find_domain[b]];

        // first value of b must be at least one after the first possible value of a
        auto first_a = a_domain.values.find_first();
        if (first_a == decltype(a_domain.values)::npos)
            return false;
        if (model.do_dynamic_occur_less_thans()) {
            for (unsigned int i = 0; i < model.target_size; i++) {
                if (a_domain.values.test(i) && (symmetric_value_displacement[i] + i < symmetric_value_displacement[first_a] + first_a)) {
                    first_a = i;
                }
            }
        }
        
        auto first_allowed_b = symmetric_value_displacement[first_a] + first_a + 1;

        if (first_allowed_b >= model.target_size)
            return false;

        if (model.do_dynamic_occur_less_thans()) {
            auto b_domain_values_copy = b_domain.values;
            for (auto v = b_domain_values_copy.find_first(); v != decltype(b_domain_values_copy)::npos; v = b_domain_values_copy.find_first()) {
                if (b_domain.values.test(v) && (symmetric_value_displacement[v] + v) < first_allowed_b) {
                    b_domain.values.reset(v);
                }
                b_domain_values_copy.reset(v);
            }
            // for (unsigned int v = 0; v < model.target_size; v++) {      // This version is probably slightly less efficient
            //     if (b_domain.values.test(v) && (symmetric_value_displacement[v] + v) < first_allowed_b) {
            //         b_domain.values.reset(v);
            //     }
            // }
        }
        else {
            for (auto v = b_domain.values.find_first(); v != decltype(b_domain.values)::npos; v = b_domain.values.find_first()) {
                if (v >= first_allowed_b)
                    break;
                b_domain.values.reset(v);
            }
        }

        // b might have shrunk (and detect empty before the next bit to make life easier)
        b_domain.count = b_domain.values.count();
        if (0 == b_domain.count)
            return false;
    }

    for (auto & [a, b] : boost::adaptors::reverse(constraints)) {
        if (find_domain[a] == -1 || find_domain[b] == -1)
            continue;
        auto & a_domain = new_domains[find_domain[a]];
        auto & b_domain = new_domains[find_domain[b]];

        // last value of a must be at least one before the last possible value of b
        if (model.do_dynamic_occur_less_thans()) {
            auto b_values_copy = b_domain.values;
            auto last_b = b_domain.values.find_first();
            for (auto v = last_b; v != decltype(b_values_copy)::npos; v = b_values_copy.find_first()) {
                b_values_copy.reset(v);
                if ((symmetric_value_displacement[last_b] + last_b) < (symmetric_value_displacement[v] + v))
                    last_b = v;
            }

            if (last_b + symmetric_value_displacement[last_b] == 0) // last_b is the first value in the ordering
                return false;
            auto last_allowed_a = symmetric_value_displacement[last_b] + last_b - 1;   // last a must be at least one before last b in the ordering

            auto a_values_copy = a_domain.values;
            for (auto v = a_values_copy.find_first(); v != decltype(a_values_copy)::npos; v = a_values_copy.find_first()) {
                a_values_copy.reset(v);
                if (v + symmetric_value_displacement[v] > last_allowed_a)
                    a_domain.values.reset(v);
            }
        }
        else {
            auto b_values_copy = b_domain.values;
            auto last_b = b_domain.values.find_first();
            for (auto v = last_b; v != decltype(b_values_copy)::npos; v = b_values_copy.find_first()) {
                b_values_copy.reset(v);
                last_b = v;
            }

            if (last_b == 0)
                return false;
            auto last_allowed_a = last_b - 1;

            auto a_values_copy = a_domain.values;
            for (auto v = a_values_copy.find_first(); v != decltype(a_values_copy)::npos; v = a_values_copy.find_first()) {
                a_values_copy.reset(v);
                if (v > last_allowed_a)
                    a_domain.values.reset(v);
            }
        }

        // a might have shrunk
        a_domain.count = a_domain.values.count();
        if (0 == a_domain.count)
            return false;
    }

    return true;
}

auto HomomorphismSearcher::propagate_occur_less_thans(
    const optional<HomomorphismAssignment> & current_assignment,
    const HomomorphismAssignments & assignments,
    Domains & new_domains) -> bool
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

    // std::cout << current_assignment->pattern_vertex << " -> " << current_assignment->target_vertex << "\n";
    // std::cout << "[";
    // for (auto & d : new_domains) {
    //     std::cout << d.v << " : [";
    //     for (int i = 0; i < model.target_size; i++) {
    //         std::cout << d.values.test(i) << ",";
    //     }
    //     std::cout << "]\n";
    // }
    // std::cout << "]\n";

    // for (auto j = 0; j < occurs.size(); j++) {
    //     if (occurs[j]) {
    //         std::cout << j << " " << occurs[j]->count() << "\n";
    //     }
        
    // }

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
                        return false;
                }
        }
        else {
            // value a first occurs in variable x, value b cannot be used in a variable lower than x
            for (auto & d : new_domains) {
                if (d.v < first_a && d.values.test(b)) {
                    occurs[b]->reset(d.v);
                    d.values.reset(b);
                    if (0 == --d.count)
                        return false;
                }
            }
        }
    }

    // propagate other way: if value b must occur (because it has been assigned) then
    // value a must go before
    if (current_assignment) {
        for (auto & [a, b] : model.target_occur_less_thans_in_convenient_order) {
            if (b != current_assignment->target_vertex) {
                // TODO sanity check more efficiently
                if (occurs[a]->count() == 0) {  // a has not been and will not be assigned
                    for (auto & d: assignments.values) {
                        if (b == d.assignment.target_vertex) return false;  // b has been assigned previously
                    }
                }
                continue;
            }

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
                            return false;
                    }
                }
            }

            for (auto & d : assignments.values)
                if (d.assignment.pattern_vertex < current_assignment->pattern_vertex && a == d.assignment.target_vertex)
                    saw_an_a = true;

            if (! saw_an_a)
                return false;
        }
    }

    return true;
}

auto HomomorphismSearcher::propagate_dynamic_occur_less_thans(
    const optional<HomomorphismAssignment> & current_assignment,
    const HomomorphismAssignments & assignments,
    Domains & new_domains) -> bool
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

    for (auto & [a, b] : useful_target_constraints) {
        build_occurs(a);
        build_occurs(b);
    }

    // std::cout << "[";
    // for (auto d: new_domains) {
    //     std::cout << "[";
    //     for (int i = 0; i < 4; i++) {
    //         std::cout << d.values.test(i) << ",";
    //     }
    //     std::cout << "]\n";
    // }
    // std::cout << "]\n";

    for (auto & a : assignments.values)
        if (occurs[a.assignment.target_vertex])
            occurs[a.assignment.target_vertex]->set(a.assignment.pattern_vertex);

    // propagate lower bounds
    for (auto & [a, b] : useful_target_constraints) {
        auto first_a = occurs[a]->find_first();
        if (first_a == SVOBitset::npos) {
            // no occurrence of value a, value b cannot be used either
            occurs[b]->reset();
            for (auto & d : new_domains)
                if (d.values.test(b)) {
                    d.values.reset(b);
                    if (0 == --d.count)
                        return false;
                }
        }
        else {
            // value a first occurs in variable x, value b cannot be used in a variable lower than x
            for (auto & d : new_domains) {
                if (d.v < first_a && d.values.test(b)) {        // TODO this line might get funky with pattern symmetries
                    occurs[b]->reset(d.v);
                    d.values.reset(b);
                    if (0 == --d.count)
                        return false;
                }
            }
        }
    }

    // propagate other way: if value b must occur (because it has been assigned) then
    // value a must go before
    if (current_assignment) {
        // for (auto & [a, b] : model.target_occur_less_thans_in_convenient_order) {
        for (auto & [a, b] : useful_target_constraints) {
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
                            return false;
                    }
                }
            }

            for (auto & d : assignments.values)
                if (d.assignment.pattern_vertex < current_assignment->pattern_vertex && a == d.assignment.target_vertex)
                    saw_an_a = true;

            if (! saw_an_a)
                return false;
        }
    }

    return true;
}

// Returns true if new useful constraints were added, false otherwise;
auto HomomorphismSearcher::make_useful_target_constraints(
    const std::optional<HomomorphismAssignment> & current_assignment,
    std::vector<std::pair<unsigned int, unsigned int>> & useful_constraints
) -> bool
{
    if(!current_assignment) return false;
    unsigned int t = current_assignment->target_vertex;

    if (std::find(target_base.begin(), target_base.end(), t) == target_base.end()) {
        target_base.push_back(t);       // Add this vertex as a new base point

        unsigned int size_before = useful_constraints.size();

        innards::dynamic_order_constraints(model.target_size, target_base, t_rschreier, useful_constraints);    // Compute constraints at the new base point

        if (model.has_less_thans() && model.do_dynamic_occur_less_thans()) {        // We need to record value order if we're also doing pattern symmetry
            for (unsigned int i = size_before; i < useful_constraints.size(); i++) {       // for each a<b just added
                std::pair<unsigned int, unsigned int> con = useful_constraints[i];
                if ((con.first + symmetric_value_displacement[con.first]) > (con.second + symmetric_value_displacement[con.second])) {    // if b is currently before a in the value ordering
                    int count = 0;
                    for (unsigned int i = 0; i < model.target_size; i++) {   // we have to check every vertex since we only know offset, not explicit order
                        if ((i + symmetric_value_displacement[i]) >= (con.second + symmetric_value_displacement[con.second]) &&
                            (i + symmetric_value_displacement[i]) < (con.first + symmetric_value_displacement[con.first])) {    // if b <= i < a in the current ordering
                            symmetric_value_displacement[i]++;      // Shift i up one index
                            count++;
                        }
                    }
                    symmetric_value_displacement[con.first] -= count;       // Shift a left by `count` indices
                }
            }
        }

        return (useful_constraints.size() - size_before) > 0;       // Return true if new constraints were added

    }

    return false;
}

auto HomomorphismSearcher::make_useful_pattern_constraints(
    const std::optional<HomomorphismAssignment> & current_assignment,
    std::vector<std::pair<unsigned int, unsigned int>> & useful_constraints
) -> bool
{
    if(!current_assignment) return false;
    unsigned int p = current_assignment->pattern_vertex;

    if (std::find(pattern_base.begin(), pattern_base.end(), p) == pattern_base.end()) {
        pattern_base.push_back(p);      // Add this vertex as a new base point

        int size_before = useful_constraints.size();

        innards::dynamic_order_constraints(model.pattern_size, pattern_base, p_rschreier, useful_constraints);      // Compute new constraints at new base point

        return (useful_constraints.size() - size_before) > 0;       // Return true if new constraints were added

    }

    return false;
}

auto HomomorphismSearcher::propagate(bool initial, Domains & new_domains, HomomorphismAssignments & assignments, bool propagate_using_lackey) -> bool
{
    // nogoods might be watching things in initial assignments. this is possibly not the
    // best place to put this...
    if (initial && might_have_watches(params)) {
        for (auto & a : assignments.values) {
            HomomorphismAssignment current_assignment = {a.assignment.pattern_vertex, a.assignment.target_vertex};
            bool wipeout = false;
            watches.propagate(
                current_assignment,
                [&](const HomomorphismAssignment & a) { return ! assignments.contains(a); },
                [&](const HomomorphismAssignment & a) {
                    for (auto & d : new_domains) {
                        if (d.v == a.pattern_vertex) {
                            if (d.values.test(a.target_vertex)) {
                                d.values.reset(a.target_vertex);
                                if (0 == --d.count)
                                    wipeout = true;
                            }
                            break;
                        }
                    }
                });

            if (wipeout)
                return false;
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
                bool wipeout = false;
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
                                        wipeout = true;
                                }
                                break;
                            }
                        }
                    });

                if (wipeout)
                    return false;
            }

            // propagate simple all different and adjacency
            if (! propagate_simple_constraints(new_domains, *current_assignment))
                return false;
        }

        // propagate less thans
        if (model.has_less_thans() && !model.do_dynamic_less_thans()) {
            if (!propagate_less_thans(new_domains)) {
                return false;
            }
        }
        if (model.has_occur_less_thans()) {
            if (model.do_dynamic_occur_less_thans()) {
                make_useful_target_constraints(current_assignment, useful_target_constraints);
                if (!propagate_dynamic_occur_less_thans(current_assignment, assignments, new_domains))
                    return false;
            }
            else if (! propagate_occur_less_thans(current_assignment, assignments, new_domains)) 
                return false;
        }
        if (model.has_less_thans() && model.do_dynamic_less_thans()) {
                make_useful_pattern_constraints(current_assignment, useful_pattern_constraints);
                if (!propagate_less_thans(new_domains, useful_pattern_constraints)) return false;
            }

        // propagate all different
        if (params.injectivity == Injectivity::Injective)
            if (! cheap_all_different(model.target_size, new_domains, proof, &model))
                return false;
        done_globals_at_least_once = true;
    }

    int dcount = 0;
    if (params.lackey && (propagate_using_lackey || params.send_partials_to_lackey)) {
        VertexToVertexMapping mapping;
        expand_to_full_result(assignments, mapping);

        if (! propagate_using_lackey) {
            if (! params.lackey->check_solution(mapping, true, false, Lackey::DeletionFunction{}))
                return false;
        }
        else {
            bool wipeout = false;
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
                                wipeout = true;
                            return true;
                        }
                    }
                }
                return false;
            };

            if (! params.lackey->check_solution(mapping, true, false, deletion) || wipeout)
                return false;
        }
    }

    return true;
}

auto HomomorphismSearcher::set_seed(int t) -> void
{
    global_rand.seed(t);
}

auto HomomorphismSearcher::print_pattern_constraints() -> void {
    std::cout << "pattern constraints = ";
    for (auto & [a,b] : useful_pattern_constraints) {
        std::cout << a << "<" << b << " ";
    }
    std::cout << "\n";
}

auto HomomorphismSearcher::print_target_constraints() -> void {
    std::cout << "target constraints = ";
    for (auto & [a,b] : useful_target_constraints) {
        std::cout << a << "<" << b << " ";
    }
    std::cout << "\n";
}
