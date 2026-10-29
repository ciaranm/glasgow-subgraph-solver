#include <gss/innards/cheap_all_different.hh>
#include <gss/innards/homomorphism_proofs.hh>
#include <gss/innards/homomorphism_searcher.hh>

#include <algorithm>
#include <cmath>
#include <optional>
#include <set>
#include <stdexcept>
#include <unordered_map>

using namespace gss;
using namespace gss::innards;

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
using std::uniform_int_distribution;
using std::vector;

namespace
{
    struct TallyKeyHash
    {
        auto operator()(const vector<int> & key) const -> size_t
        {
            size_t result = key.size();
            for (auto k : key)
                result = (result * 1000003u) ^ size_t(unsigned(k));
            return result;
        }
    };

    using TallyCounts = std::unordered_map<vector<int>, loooong, TallyKeyHash>;

    // A split's subproblem, as a cache key: see split_into_components.
    using CacheKey = vector<unsigned long long>;

    struct CacheKeyHash
    {
        auto operator()(const CacheKey & key) const -> size_t
        {
            size_t result = key.size();
            for (auto k : key)
                result = (result * 1000003u) ^ size_t(k ^ (k >> 32));
            return result;
        }
    };

    auto append_to_key(CacheKey & key, const SVOBitset & bits) -> void
    {
        key.insert(key.end(), bits.words(), bits.words() + bits.number_of_words());
    }
}

namespace gss::innards
{
    // What a split contributes to the enclosing tally (or, at the top, to the count): keyed by
    // the images of layout, which are the enclosing tally's vertices among those unassigned at
    // the split. Empty if the split has no solutions.
    struct SplitResult
    {
        vector<unsigned> layout;
        TallyCounts combined;
    };

    // A split's result, and a component's tally, depend only on the domains of the vertices
    // involved and on which images are recorded, not on how the search got there, so each is
    // cached on exactly those.
    struct SplitCaches
    {
        std::unordered_map<CacheKey, SplitResult, CacheKeyHash> splits;
        std::unordered_map<CacheKey, TallyCounts, CacheKeyHash> components;
        unsigned long long split_hits = 0, component_hits = 0;

        // Roughly what the caches hold, and how many times they have been emptied for
        // reaching HomomorphismParams::split_cache_limit.
        unsigned long long bytes = 0, clears = 0;

        // Make room for an entry of about this many bytes: emptying the caches loses only
        // the work in them, since every entry can be worked out again.
        auto make_room(unsigned long long entry_bytes, unsigned long long limit) -> void
        {
            if (bytes + entry_bytes > limit && bytes > 0) {
                splits.clear();
                components.clear();
                bytes = 0;
                ++clears;
            }
            bytes += entry_bytes;
        }
    };

    // Roughly the memory a cache entry takes, counting the hash table's nodes and the keys.
    auto approximate_bytes(const CacheKey & key, const TallyCounts & counts, std::size_t width) -> unsigned long long
    {
        return 64 + 8 * key.size() + counts.size() * (96 + 4 * width);
    }

    // The solutions of one component of a split, keyed by the images of the vertices that the
    // split still needs to see: those that could collide with another component, and those
    // the enclosing tally needs, if this split is itself inside a component. An image that
    // cannot matter to any of them is recorded as -1, so solutions that differ only there
    // share an entry.
    struct ComponentTally
    {
        // ascending
        vector<unsigned> vertices;

        // for each of vertices, the images that are recorded rather than collapsed to -1
        vector<SVOBitset> relevant;

        TallyCounts counts;

        auto index_of(unsigned v) const -> optional<unsigned>
        {
            auto i = std::lower_bound(vertices.begin(), vertices.end(), v);
            if (i == vertices.end() || *i != v)
                return nullopt;
            return unsigned(i - vertices.begin());
        }

        auto image(unsigned i, int t) const -> int
        {
            return (t != -1 && relevant[i].test(t)) ? t : -1;
        }
    };
}

HomomorphismSearcher::HomomorphismSearcher(const HomomorphismModel & m, const HomomorphismParams & p,
    const DuplicateSolutionFilterer & d, const std::shared_ptr<Proof> & f,
    Watches<HomomorphismAssignment, HomomorphismAssignmentWatchTable> & w) :
    model(m),
    params(p),
    _duplicate_solution_filterer(d),
    proof(f),
    _record_filter_activations(p.record_filter_activations),
    _verbose_proof_comments(f && f->super_extra_verbose()),
    _track_removals(p.record_filter_activations || (f && f->super_extra_verbose())),
    _filter_activations(m.max_graphs),
    watches(w)
{
    if (might_have_watches(params)) {
        watches.table.target_size = model.target_size;
        watches.table.data.resize(model.pattern_size * model.target_size);
    }
}

HomomorphismSearcher::~HomomorphismSearcher() = default;

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

auto HomomorphismSearcher::add_extra_stats(std::list<std::string> & stats) const -> void
{
    if (_record_filter_activations)
        _filter_activations.add_search_stats(stats);
    if (params.split_components) {
        stats.push_back("splits = " + to_string(_splits));
        stats.push_back("nested_splits = " + to_string(_nested_splits));
        stats.push_back("split_cache_hits = " + to_string(_split_caches ? _split_caches->split_hits : 0));
        stats.push_back("split_cache_entries = " + to_string(_split_caches ? _split_caches->splits.size() : 0));
        stats.push_back("component_cache_hits = " + to_string(_split_caches ? _split_caches->component_hits : 0));
        stats.push_back("component_cache_entries = " + to_string(_split_caches ? _split_caches->components.size() : 0));
        stats.push_back("split_cache_clears = " + to_string(_split_caches ? _split_caches->clears : 0));
    }
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
        vector<pair<int, vector<int>>> proof_domains;
        for (auto & d : domains) {
            proof_domains.push_back(pair{int(d.v), vector<int>{}});
            auto values = d.values;
            for (auto v = values.find_first(); v != decltype(values)::npos; v = values.find_first()) {
                values.reset(v);
                proof_domains.back().second.push_back(int(v));
            }
        }
        model.proofs()->show_domains("entering depth " + to_string(depth), proof_domains);
    }

    if (params.timeout->should_abort())
        return SearchResult::Aborted;

    ++nodes;

    // find ourselves a domain, or succeed if we're all assigned
    const HomomorphismDomain * branch_domain = find_branch_domain(domains);
    if (! branch_domain) {
        if (_tally) {
            vector<int> key(_tally->vertices.size(), -1);
            for (auto & a : assignments.values)
                if (auto i = _tally->index_of(a.assignment.pattern_vertex))
                    key[*i] = _tally->image(*i, int(a.assignment.target_vertex));
            ++_tally->counts[key];
            return SearchResult::SatisfiableButKeepGoing;
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

    if (params.split_components)
        if (auto result = split_into_components(assignments, domains, nodes, propagations, solution_count, depth, restarts_schedule))
            return *result;

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

    // for each value remaining...
    for (auto f_v = branch_v.begin(), f_end = branch_v.begin() + branch_v_end; f_v != f_end; ++f_v) {
        if (proof) {
            if (auto hp = model.proofs()) {
                hp->guessing(depth, int(branch_domain->v), int(*f_v));
                // branching on branch_v->*f_v forward-checks using its supplemental adjacency; emit
                // any still-pending supplemental for this antecedent before the propagation.
                hp->materialise_adjacency_for(int(branch_domain->v), int(*f_v));
            }
        }

        // modified in-place by appending, we can restore by shrinking
        auto assignments_size = assignments.values.size();

        // make the assignment
        assignments.values.push_back({{branch_domain->v, unsigned(*f_v)}, true, discrepancy_count, int(branch_v_end)});

        // set up new domains
        Domains new_domains = copy_nonfixed_domains_and_make_assignment(domains, branch_domain->v, *f_v);

        // propagate
        ++propagations;
        if (! propagate(false, new_domains, assignments)) {
            // failure? restore assignments and go on to the next thing
            if (proof)
                model.proofs()->propagation_failure(assignments_as_proof_decisions(assignments), int(branch_domain->v), int(*f_v));

            assignments.values.resize(assignments_size);
            actually_hit_a_failure = true;

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
        return SearchResult::Unsatisfiable;
}

auto HomomorphismSearcher::split_into_components(
    HomomorphismAssignments & assignments,
    const Domains & domains,
    unsigned long long & nodes,
    unsigned long long & propagations,
    loooong & solution_count,
    int depth,
    RestartsSchedule & restarts_schedule) -> optional<SearchResult>
{
    // Every constraint between an assigned and an unassigned vertex has already been
    // propagated into the domains, so what is left is the constraints among the unassigned
    // vertices. Those inside a connected component are the component's own business. Between
    // two components there are no edges, and since induced mappings are refused, the only
    // constraint left is (local) injectivity: two vertices may not share an image. So the
    // solutions below this node are the combinations of one solution per component in which
    // no two such vertices collide, and each component can be searched on its own, keeping
    // only enough of each solution to tell which combinations collide.
    vector<const HomomorphismDomain *> unassigned;
    SVOBitset unassigned_mask(model.pattern_size, 0);
    for (auto & d : domains)
        if (! d.fixed) {
            unassigned.push_back(&d);
            unassigned_mask.set(d.v);
        }

    if (unassigned.size() < 2)
        return nullopt;

    // edges in either direction join a component
    vector<int> component(model.pattern_size, -1);
    int n_components = 0;
    SVOBitset unvisited = unassigned_mask;
    vector<unsigned> todo;
    for (auto * d : unassigned) {
        if (! unvisited.test(d->v))
            continue;
        unvisited.reset(d->v);
        component[d->v] = n_components;
        todo.push_back(d->v);
        while (! todo.empty()) {
            auto v = todo.back();
            todo.pop_back();
            auto row = model.pattern_graph_row(0, v);
            row |= model.pattern_in_neighbour_row(v);
            row &= unvisited;
            unvisited.intersect_with_complement(row);
            for (auto w = row.find_first(); w != decltype(row)::npos; w = row.find_first()) {
                row.reset(w);
                component[w] = n_components;
                todo.push_back(w);
            }
        }
        ++n_components;
    }

    if (n_components < 2)
        return nullopt;

    auto may_collide = [&](unsigned v, unsigned w) -> bool {
        switch (params.injectivity) {
        case Injectivity::Injective: return true;
        case Injectivity::LocallyInjective: return both_in_the_neighbourhood_of_some_vertex(v, w);
        case Injectivity::NonInjective: return false;
        }
        return true;
    };

    // If this split is itself inside a component of an enclosing split, its solutions go to
    // that component's tally, which needs the images of some of these vertices too.
    ComponentTally * const enclosing = _tally;
    auto enclosing_index = [&](unsigned v) -> optional<unsigned> {
        return enclosing ? enclosing->index_of(v) : nullopt;
    };

    // The pairs of vertices in different components that could take the same image, by index
    // into unassigned, and for each vertex in such a pair, the images that a partner could
    // also take. Every table the combining builds is keyed by some of the recorded images, so
    // the number of possible signatures bounds them all, and a split with too many is not
    // made. Each vertex in a pair at least doubles that number, so the scan can stop early.
    const double log_limit = std::log2(double(*params.split_components));
    vector<pair<unsigned, unsigned>> collisions;
    vector<optional<SVOBitset>> shared(unassigned.size());
    unsigned vertices_in_collisions = 0;
    for (unsigned i = 0; i < unassigned.size(); ++i)
        for (unsigned j = i + 1; j < unassigned.size(); ++j)
            if (component[unassigned[i]->v] != component[unassigned[j]->v] && may_collide(unassigned[i]->v, unassigned[j]->v)) {
                auto both = unassigned[i]->values;
                both &= unassigned[j]->values;
                if (both.any()) {
                    collisions.emplace_back(i, j);
                    for (auto k : {i, j}) {
                        if (! shared[k]) {
                            if (++vertices_in_collisions > log_limit)
                                return nullopt;
                            shared[k] = SVOBitset(model.target_size, 0);
                        }
                        *shared[k] |= both;
                    }
                }
            }

    // A recorded image is one of the relevant ones, or -1.
    double log_signatures = 0.0;
    for (unsigned i = 0; i < unassigned.size(); ++i) {
        auto e = enclosing_index(unassigned[i]->v);
        if (! e && ! shared[i])
            continue;
        SVOBitset relevant(model.target_size, 0);
        if (e)
            relevant |= enclosing->relevant[*e];
        if (shared[i])
            relevant |= *shared[i];
        log_signatures += std::log2(double(relevant.count()) + 1.0);
    }
    if (log_signatures > log_limit)
        return nullopt;

    ++_splits;
    if (_tally)
        ++_nested_splits;

    // Hand a result to the enclosing tally, or at the top, to the count.
    auto deliver = [&](const SplitResult & result) -> SearchResult {
        if (result.combined.empty())
            return SearchResult::Unsatisfiable;

        if (! enclosing) {
            // layout is empty, and so combined has the one entry
            for (auto & [_, count] : result.combined)
                solution_count += count;
        }
        else {
            // The enclosing tally's vertices are either assigned already, or in layout.
            vector<int> key(enclosing->vertices.size(), -1);
            for (auto & a : assignments.values)
                if (auto e = enclosing->index_of(a.assignment.pattern_vertex))
                    key[*e] = enclosing->image(*e, int(a.assignment.target_vertex));
            vector<unsigned> layout_to_enclosing;
            for (auto v : result.layout)
                layout_to_enclosing.push_back(*enclosing->index_of(v));
            for (auto & [combined_key, count] : result.combined) {
                for (unsigned k = 0; k < result.layout.size(); ++k)
                    key[layout_to_enclosing[k]] = enclosing->image(layout_to_enclosing[k], combined_key[k]);
                enclosing->counts[key] += count;
            }
        }

        return SearchResult::SatisfiableButKeepGoing;
    };

    // The result depends only on the unassigned vertices' domains, and on which of their
    // images the enclosing tally records. Everything else about this split, from the
    // components to the order they are combined in, follows from those.
    vector<const HomomorphismDomain *> by_vertex = unassigned;
    std::sort(by_vertex.begin(), by_vertex.end(), [](auto * a, auto * b) { return a->v < b->v; });
    CacheKey split_key;
    split_key.push_back(by_vertex.size());
    for (auto * d : by_vertex) {
        split_key.push_back(d->v);
        append_to_key(split_key, d->values);
    }
    for (auto * d : by_vertex)
        if (auto e = enclosing_index(d->v)) {
            split_key.push_back(d->v);
            append_to_key(split_key, enclosing->relevant[*e]);
        }

    if (! _split_caches)
        _split_caches = std::make_unique<SplitCaches>();
    if (auto cached = _split_caches->splits.find(split_key); cached != _split_caches->splits.end()) {
        ++_split_caches->split_hits;
        return deliver(cached->second);
    }

    // Search the components that look smallest first, so that a component with no
    // solutions is found as early as possible.
    vector<double> size_estimate(n_components, 0.0);
    for (auto * d : unassigned)
        size_estimate[component[d->v]] += std::log2(double(d->count));
    vector<int> order(n_components);
    for (int c = 0; c < n_components; ++c)
        order[c] = c;
    std::stable_sort(order.begin(), order.end(), [&](int a, int b) { return size_estimate[a] < size_estimate[b]; });
    vector<int> position(n_components);
    for (int s = 0; s < n_components; ++s)
        position[order[s]] = s;

    // For each unassigned vertex, the last position in the order of a component holding a
    // vertex it could collide with.
    vector<int> last_partner(unassigned.size(), -1);
    std::set<pair<unsigned, unsigned>> colliding_vertices;
    for (auto [i, j] : collisions) {
        last_partner[i] = max(last_partner[i], position[component[unassigned[j]->v]]);
        last_partner[j] = max(last_partner[j], position[component[unassigned[i]->v]]);
        colliding_vertices.emplace(unassigned[i]->v, unassigned[j]->v);
        colliding_vertices.emplace(unassigned[j]->v, unassigned[i]->v);
    }

    vector<int> index_in_unassigned(model.pattern_size, -1);
    for (unsigned i = 0; i < unassigned.size(); ++i)
        index_in_unassigned[unassigned[i]->v] = i;

    // Is v's image still needed once the components up to position s have been combined?
    auto needed_after = [&](unsigned v, int s) -> bool {
        return enclosing_index(v) || last_partner[index_in_unassigned[v]] > s;
    };

    SplitResult result;

    // Count each component on its own, into its own tally, and combine the tallies into
    // result. Returns false if aborted, and leaves result empty if there are no solutions.
    auto count_and_combine = [&]() -> bool {
        vector<ComponentTally> tallies(n_components);
        for (int s = 0; s < n_components; ++s) {
            int c = order[s];
            auto & tally = tallies[s];
            Domains component_domains;
            for (auto * d : by_vertex)
                if (component[d->v] == c) {
                    component_domains.push_back(*d);
                    if (enclosing_index(d->v) || shared[index_in_unassigned[d->v]])
                        tally.vertices.push_back(d->v);
                }
            for (auto v : tally.vertices) {
                SVOBitset relevant(model.target_size, 0);
                if (auto e = enclosing_index(v))
                    relevant |= enclosing->relevant[*e];
                if (auto & mine = shared[index_in_unassigned[v]])
                    relevant |= *mine;
                tally.relevant.push_back(relevant);
            }

            CacheKey component_key;
            component_key.push_back(component_domains.size());
            for (auto & d : component_domains) {
                component_key.push_back(d.v);
                append_to_key(component_key, d.values);
            }
            for (unsigned i = 0; i < tally.vertices.size(); ++i) {
                component_key.push_back(tally.vertices[i]);
                append_to_key(component_key, tally.relevant[i]);
            }

            if (auto cached = _split_caches->components.find(component_key); cached != _split_caches->components.end()) {
                ++_split_caches->component_hits;
                tally.counts = cached->second;
            }
            else {
                _tally = &tally;
                auto search_result = restarting_search(assignments, component_domains, nodes, propagations, solution_count, depth, restarts_schedule);
                _tally = enclosing;

                switch (search_result) {
                case SearchResult::Aborted:
                    return false;
                case SearchResult::Unsatisfiable:
                case SearchResult::SatisfiableButKeepGoing:
                    break;
                case SearchResult::Satisfiable:
                case SearchResult::Restart:
                    // refused in solve_homomorphism_problem: no enumerate callback, no restarts
                    throw std::logic_error{"split component search did not run to completion"};
                }

                _split_caches->make_room(approximate_bytes(component_key, tally.counts, tally.vertices.size()), params.split_cache_limit);
                _split_caches->components.emplace(move(component_key), tally.counts);
            }

            if (tally.counts.empty())
                return true;
        }

        // Combine the tallies one component at a time, keeping in the key only the images
        // that a later component could collide with, or that the enclosing tally needs.
        vector<unsigned> layout;
        TallyCounts combined;
        combined.emplace(vector<int>{}, loooong{1});
        for (int s = 0; s < n_components; ++s) {
            auto & tally = tallies[s];

            vector<pair<unsigned, unsigned>> compare;
            for (unsigned a = 0; a < layout.size(); ++a)
                for (unsigned b = 0; b < tally.vertices.size(); ++b)
                    if (colliding_vertices.contains({layout[a], tally.vertices[b]}))
                        compare.emplace_back(a, b);

            // where each entry of the new key comes from: the combined key, or the tally's
            vector<unsigned> new_layout;
            vector<pair<bool, unsigned>> take;
            for (unsigned a = 0; a < layout.size(); ++a)
                if (needed_after(layout[a], s)) {
                    new_layout.push_back(layout[a]);
                    take.emplace_back(false, a);
                }
            for (unsigned b = 0; b < tally.vertices.size(); ++b)
                if (needed_after(tally.vertices[b], s)) {
                    new_layout.push_back(tally.vertices[b]);
                    take.emplace_back(true, b);
                }

            TallyCounts next;
            vector<int> key(take.size());
            unsigned long long pairs_considered = 0;
            for (auto & [left_key, left_count] : combined)
                for (auto & [right_key, right_count] : tally.counts) {
                    // the signatures bound this, but the bound can still be large
                    if (0 == (++pairs_considered % 65536) && params.timeout->should_abort())
                        return false;

                    bool compatible = true;
                    for (auto [a, b] : compare)
                        if (left_key[a] != -1 && left_key[a] == right_key[b]) {
                            compatible = false;
                            break;
                        }
                    if (! compatible)
                        continue;

                    for (unsigned k = 0; k < take.size(); ++k)
                        key[k] = take[k].first ? right_key[take[k].second] : left_key[take[k].second];
                    next[key].add_product(left_count, right_count);
                }

            if (next.empty())
                return true;

            combined = move(next);
            layout = move(new_layout);
        }

        result.layout = move(layout);
        result.combined = move(combined);
        return true;
    };

    if (! count_and_combine())
        return SearchResult::Aborted;

    _split_caches->make_room(approximate_bytes(split_key, result.combined, result.layout.size()), params.split_cache_limit);
    auto & stored = _split_caches->splits.emplace(move(split_key), move(result)).first->second;
    return deliver(stored);
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

template <bool directed_, bool has_edge_labels_, bool induced_, bool track_removals_>
auto HomomorphismSearcher::propagate_adjacency_constraints(HomomorphismDomain & d, const HomomorphismAssignment & current_assignment) -> void
{
    const auto & graph_pairs_to_consider = model.pattern_adjacency_bits(current_assignment.pattern_vertex, d.v);

    // Lazy supplemental proofs: every value this forward-check removes from dom(d.v) is an
    // adjacency-based elimination whose branch-consequence a later Hall confinement may cite,
    // so materialise the supplemental for (d.v, removed value) below. Snapshot dom(d.v) now;
    // skipped once nothing is left to defer (or with proofs off).
    const bool lazy_proof = proof && model.proofs() && model.proofs()->has_pending_supplementals();
    SVOBitset lazy_proof_before;
    if (lazy_proof)
        lazy_proof_before = d.values;

    // Attributing a removal to the graph pair that made it costs a popcount per graph pair,
    // which is why it is a template parameter rather than a runtime flag: measured on a
    // 2-million-node unsatisfiable instance, a perfectly predicted branch here still costs
    // 3.5%. Two things want the attribution -- verbose proof comments and filter-activation
    // recording -- and neither is on in a normal solve, so they share the one instantiation
    // (see _track_removals).
    [[maybe_unused]] unsigned long long remaining = 0;
    if constexpr (track_removals_)
        remaining = d.values.count();

    [[maybe_unused]] const auto track = [&](unsigned g, unsigned long long & activation_counter) {
        if constexpr (track_removals_) {
            auto now = d.values.count();
            if (now != remaining) {
                if (_verbose_proof_comments)
                    model.proofs()->propagated(int(current_assignment.pattern_vertex), int(current_assignment.target_vertex),
                        int(g), remaining - now, int(d.v));
                if (_record_filter_activations)
                    activation_counter += remaining - now;
                remaining = now;
            }
        }
    };

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

    track(0, _filter_activations.search[0]);

    // and for each remaining graph pair...
    for (unsigned g = 1; g < model.max_graphs; ++g) {
        // if we're adjacent...
        if (graph_pairs_to_consider & (1u << g)) {
            // ...then we can only be mapped to adjacent vertices
            d.values &= model.target_graph_row(g, current_assignment.target_vertex);
        }

        track(g, _filter_activations.search[g]);
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

    if constexpr (has_edge_labels_ && track_removals_) {
        // no proof comment for this one: the edge-label check has never had one, and the
        // proof's propagated() line is indexed by graph pair, which this is not
        if (_record_filter_activations)
            _filter_activations.search_edge_labels += remaining - d.values.count();
    }

    if (lazy_proof) {
        auto removed = lazy_proof_before;
        removed.intersect_with_complement(d.values); // before & ~after = removed values
        for (auto c = removed.find_first(); c != decltype(removed)::npos; c = removed.find_first()) {
            removed.reset(c);
            // pattern_vertex_for_proof / target_vertex_for_proof are the identity on the index,
            // so (d.v, c) is already the proof-space (p, t) antecedent.
            model.proofs()->materialise_adjacency_for(int(d.v), int(c));
        }
    }
}

auto HomomorphismSearcher::both_in_the_neighbourhood_of_some_vertex(unsigned v, unsigned w) -> bool
{
    // Is there a pattern vertex u with u -> v and u -> w? That is what local injectivity
    // asks about, and it is what the verifier and the proof encoding both mean by it. On an
    // undirected pattern the in-neighbourhood is just the row, so this is unchanged there;
    // on a directed one, intersecting the rows would ask the unrelated question of whether v
    // and w share a successor (issue #96).
    auto i = model.pattern_in_neighbour_row(v);
    i &= model.pattern_in_neighbour_row(w);
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
                    if (! _track_removals)
                        propagate_adjacency_constraints<true, false, true, false>(d, current_assignment);
                    else
                        propagate_adjacency_constraints<true, false, true, true>(d, current_assignment);
                }
                else {
                    if (! _track_removals)
                        propagate_adjacency_constraints<false, false, true, false>(d, current_assignment);
                    else
                        propagate_adjacency_constraints<false, false, true, true>(d, current_assignment);
                }
            }
            else {
                if (model.directed()) {
                    if (! _track_removals)
                        propagate_adjacency_constraints<true, false, false, false>(d, current_assignment);
                    else
                        propagate_adjacency_constraints<true, false, false, true>(d, current_assignment);
                }
                else {
                    if (! _track_removals)
                        propagate_adjacency_constraints<false, false, false, false>(d, current_assignment);
                    else
                        propagate_adjacency_constraints<false, false, false, true>(d, current_assignment);
                }
            }
        }
        else {
            // edge labels are always directed
            if (params.induced) {
                if (! _track_removals)
                    propagate_adjacency_constraints<true, true, true, false>(d, current_assignment);
                else
                    propagate_adjacency_constraints<true, true, true, true>(d, current_assignment);
            }
            else {
                if (! _track_removals)
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

auto HomomorphismSearcher::propagate_less_thans(Domains & new_domains) -> bool
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
            return false;
        auto first_allowed_b = first_a + 1;

        if (first_allowed_b >= model.target_size)
            return false;

        for (auto v = b_domain.values.find_first(); v != decltype(b_domain.values)::npos; v = b_domain.values.find_first()) {
            if (v >= first_allowed_b)
                break;
            b_domain.values.reset(v);
        }

        // b might have shrunk (and detect empty before the next bit to make life easier)
        b_domain.count = b_domain.values.count();
        if (0 == b_domain.count)
            return false;
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
            return false;
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

auto HomomorphismSearcher::propagate(bool initial, Domains & new_domains, HomomorphismAssignments & assignments) -> bool
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

            if (proof) {
                if (auto hp = model.proofs()) {
                    hp->unit_propagating(int(current_assignment->pattern_vertex), int(current_assignment->target_vertex));
                    // a forced var->val also forward-checks using its supplemental adjacency.
                    hp->materialise_adjacency_for(int(current_assignment->pattern_vertex), int(current_assignment->target_vertex));
                }
            }

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
        if (model.has_less_thans() && ! propagate_less_thans(new_domains))
            return false;
        if (model.has_occur_less_thans() && ! propagate_occur_less_thans(current_assignment, assignments, new_domains))
            return false;

        // propagate all different
        if (params.injectivity == Injectivity::Injective)
            if (! cheap_all_different(model.target_size, new_domains, proof, &model))
                return false;
        done_globals_at_least_once = true;
    }

    return true;
}

auto HomomorphismSearcher::set_seed(int t) -> void
{
    global_rand.seed(t);
}
