#include <gss/configuration.hh>
#include <gss/counting_statement.hh>
#include <gss/homomorphism.hh>
#include <gss/innards/proof.hh>
#include <gss/loooong.hh>
#include <gss/sip_decomposer.hh>

#include <algorithm>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace gss;
using namespace gss::innards;

using std::map;
using std::optional;
using std::pair;
using std::set;
using std::sort;
using std::string;
using std::to_string;
using std::vector;

namespace
{
    auto find_removable_isolated_pattern_vertices(
        const InputGraph & pattern,
        const HomomorphismParams & params,
        set<int> & isolated_pattern_vertices) -> void
    {
        // the solution count fix-up below is a falling factorial, which is only right if
        // the isolated vertices must land on distinct target vertices not used by the rest
        if (params.induced || params.injectivity != Injectivity::Injective || (! params.pattern_less_constraints.empty()) || (! params.target_occur_less_constraints.empty()) || pattern.has_vertex_labels() || (params.count_solutions && params.enumerate_callback) || pattern.directed())
            return;

        for (int i = 0; i < pattern.size(); ++i)
            if (0 == pattern.degree(i))
                isolated_pattern_vertices.emplace(i);
    }

    template <typename I_>
    auto n_choose_k(I_ n, I_ k) -> I_
    {
        I_ r = 1;
        for (I_ d = 1; d <= k; ++d) {
            I_ numerator = n - k + d;
            I_ denominator = d;
            I_ cf = gss::gcd(r, denominator);
            r /= cf;
            denominator /= cf;
            numerator /= denominator;
            r *= numerator;
        }
        return r;
    }

    template <typename I_>
    auto factorial(I_ n) -> I_
    {
        I_ r = 1;
        for (I_ d = 2; d <= n; ++d)
            r *= d;
        return r;
    }

    auto falling_factorial(long n, long k) -> loooong
    {
        loooong r = 1;
        for (long i = 0; i < k; ++i)
            r *= loooong{n - i};
        return r;
    }

    // A graph small enough to want a canonical form for: vertex labels (or none),
    // and labelled undirected edges, keyed with the smaller end first.
    struct SmallGraph
    {
        int size = 0;
        vector<string> labels;
        map<pair<int, int>, string> edges;
    };

    auto encode(const SmallGraph & g, const vector<int> & order) -> string
    {
        vector<int> position(g.size);
        for (int i = 0; i < g.size; ++i)
            position[order[i]] = i;

        string result = to_string(g.size) + ";";
        if (! g.labels.empty())
            for (auto v : order)
                result += to_string(g.labels[v].size()) + ":" + g.labels[v];
        vector<std::tuple<int, int, const string *>> edges;
        for (auto & [e, label] : g.edges) {
            auto [a, b] = std::minmax(position[e.first], position[e.second]);
            edges.emplace_back(a, b, &label);
        }
        sort(edges.begin(), edges.end(), [](const auto & x, const auto & y) {
            return std::tie(std::get<0>(x), std::get<1>(x), *std::get<2>(x)) < std::tie(std::get<0>(y), std::get<1>(y), *std::get<2>(y));
        });
        for (auto & [a, b, label] : edges)
            result += ";" + to_string(a) + "," + to_string(b) + ":" + to_string(label->size()) + ":" + *label;
        return result;
    }

    // A string that is equal for two graphs if and only if they are isomorphic,
    // computed by colour refinement and then trying every ordering that respects the
    // colours. When that is too many orderings, we give up on isomorphism and use the
    // graph as it is, which is still correct, just less good as a cache key.
    auto canonical_key(const SmallGraph & g) -> string
    {
        vector<vector<pair<int, const string *>>> neighbours(g.size);
        for (auto & [e, label] : g.edges) {
            neighbours[e.first].emplace_back(e.second, &label);
            if (e.first != e.second)
                neighbours[e.second].emplace_back(e.first, &label);
        }

        // colours are numbered by sorting their descriptions, never by vertex number,
        // so that isomorphic graphs get the same colours
        auto renumber = [&](const vector<string> & descriptions) {
            vector<string> sorted = descriptions;
            sort(sorted.begin(), sorted.end());
            sorted.erase(std::unique(sorted.begin(), sorted.end()), sorted.end());
            vector<int> colours(g.size);
            for (int v = 0; v < g.size; ++v)
                colours[v] = std::lower_bound(sorted.begin(), sorted.end(), descriptions[v]) - sorted.begin();
            return pair{colours, int(sorted.size())};
        };

        vector<string> descriptions(g.size);
        for (int v = 0; v < g.size; ++v)
            descriptions[v] = g.labels.empty() ? string{} : g.labels[v];
        auto [colours, number_of_colours] = renumber(descriptions);

        while (true) {
            for (int v = 0; v < g.size; ++v) {
                vector<string> seen;
                for (auto & [w, label] : neighbours[v])
                    seen.push_back(to_string(colours[w]) + ":" + *label);
                sort(seen.begin(), seen.end());
                descriptions[v] = to_string(colours[v]);
                for (auto & d : seen)
                    descriptions[v] += ";" + to_string(d.size()) + ":" + d;
            }
            auto [new_colours, new_number_of_colours] = renumber(descriptions);
            colours = new_colours;
            if (new_number_of_colours == number_of_colours)
                break;
            number_of_colours = new_number_of_colours;
        }

        vector<vector<int>> classes(number_of_colours);
        for (int v = 0; v < g.size; ++v)
            classes[colours[v]].push_back(v);

        constexpr long most_orderings = 5040;
        long orderings = 1;
        for (auto & c : classes)
            for (long i = 2; i <= long(c.size()); ++i)
                if ((orderings *= i) > most_orderings) {
                    vector<int> identity(g.size);
                    for (int v = 0; v < g.size; ++v)
                        identity[v] = v;
                    return "exact;" + encode(g, identity);
                }

        optional<string> best;
        vector<int> order;
        auto try_orderings = [&](auto & self, unsigned c) -> void {
            if (c == classes.size()) {
                auto e = encode(g, order);
                if ((! best) || e < *best)
                    best = e;
                return;
            }
            auto members = classes[c];
            do {
                order.insert(order.end(), members.begin(), members.end());
                self(self, c + 1);
                order.resize(order.size() - members.size());
            } while (std::next_permutation(members.begin(), members.end()));
        };
        try_orderings(try_orderings, 0);

        return "canonical;" + *best;
    }

    // H / sigma, for a partition sigma of the pattern vertices vertices[0 .. n - 1],
    // written as the block number of each, numbered in order of first occurrence. Gluing
    // vertices with different labels, or two edges with different labels into one, is
    // impossible for an embedding, and gives nothing.
    struct Quotient
    {
        int size = 0;
        map<pair<int, int>, string> edges;
        vector<string> labels;
        vector<vector<int>> components;
    };

    auto build_quotient(const InputGraph & pattern, const vector<int> & vertices, const vector<int> & sigma) -> optional<Quotient>
    {
        Quotient q;
        q.size = *std::max_element(sigma.begin(), sigma.end()) + 1;

        if (pattern.has_vertex_labels()) {
            vector<optional<string>> labels(q.size);
            for (unsigned v = 0; v < vertices.size(); ++v) {
                auto label = string{pattern.vertex_label(vertices[v])};
                if (! labels[sigma[v]])
                    labels[sigma[v]] = label;
                else if (*labels[sigma[v]] != label)
                    return std::nullopt;
            }
            for (auto & l : labels)
                q.labels.push_back(*l);
        }

        // two pattern edges can be glued into one quotient edge, in which case the
        // target edge must carry both labels, which it cannot unless they agree
        for (unsigned v = 0; v < vertices.size(); ++v)
            for (unsigned w = v; w < vertices.size(); ++w)
                if (pattern.adjacent(vertices[v], vertices[w])) {
                    auto label = pattern.has_edge_labels() ? string{pattern.edge_label(vertices[v], vertices[w])} : string{};
                    auto key = std::minmax(sigma[v], sigma[w]);
                    auto [where, inserted] = q.edges.emplace(key, label);
                    if ((! inserted) && where->second != label)
                        return std::nullopt;
                }

        vector<int> parent(q.size);
        for (int b = 0; b < q.size; ++b)
            parent[b] = b;
        auto find = [&](int b) {
            while (parent[b] != b)
                b = parent[b] = parent[parent[b]];
            return b;
        };
        for (auto & [e, _] : q.edges)
            parent[find(e.first)] = find(e.second);

        map<int, int> component_number;
        for (int b = 0; b < q.size; ++b) {
            auto [where, inserted] = component_number.emplace(find(b), q.components.size());
            if (inserted)
                q.components.emplace_back();
            q.components[where->second].push_back(b);
        }

        return q;
    }

    // How many target vertices carry each label, with every vertex of an unlabelled
    // target under the empty label.
    auto count_target_labels(const InputGraph & pattern, const InputGraph & target) -> map<string, long>
    {
        map<string, long> result;
        for (int t = 0; t < target.size(); ++t)
            ++result[pattern.has_vertex_labels() ? string{target.vertex_label(t)} : string{}];
        return result;
    }

    // An injective embedding of H / sigma needs a different target vertex for each block,
    // carrying that block's label, so there is none if some label has more blocks than
    // target vertices. This is what stops a pattern bigger than its target from
    // exploring every way of gluing it, all of which are empty.
    auto more_blocks_than_targets(const Quotient & q, const map<string, long> & target_with_label) -> bool
    {
        map<string, long> blocks_with_label;
        for (int b = 0; b < q.size; ++b) {
            auto label = q.labels.empty() ? string{} : q.labels[b];
            auto t = target_with_label.find(label);
            if (++blocks_with_label[label] > (t == target_with_label.end() ? 0 : t->second))
                return true;
        }
        return false;
    }

    // The connected component of a quotient made up of the given blocks, renumbered
    // from zero in the order given.
    auto component_graph(const InputGraph & pattern, const Quotient & q, const vector<int> & blocks) -> SmallGraph
    {
        map<int, int> index_of_block;
        for (unsigned i = 0; i < blocks.size(); ++i)
            index_of_block.emplace(blocks[i], i);

        SmallGraph component;
        component.size = blocks.size();
        if (pattern.has_vertex_labels())
            for (auto b : blocks)
                component.labels.push_back(q.labels[b]);
        for (auto & [e, label] : q.edges)
            if (auto a = index_of_block.find(e.first); a != index_of_block.end())
                component.edges.emplace(std::minmax(a->second, index_of_block.at(e.second)), label);
        return component;
    }

    auto to_input_graph(const InputGraph & pattern, const SmallGraph & component) -> InputGraph
    {
        InputGraph graph(component.size, pattern.has_vertex_labels(), pattern.has_edge_labels());
        for (int i = 0; i < component.size && pattern.has_vertex_labels(); ++i)
            graph.set_vertex_label(i, component.labels[i]);
        for (auto & [e, label] : component.edges)
            graph.add_edge(e.first, e.second, label);
        return graph;
    }

    // Calls f(tau) for each coarsening tau != sigma of sigma that glues each block of
    // H / sigma only to blocks from other components: block b goes into an existing group
    // that has nothing from its component yet, or starts a new one. Groups are opened in
    // block order, so tau is already numbered in order of first occurrence. Stops early if
    // f returns false.
    template <typename F_>
    auto for_each_coarsening(const Quotient & q, const vector<int> & sigma, const F_ & f) -> void
    {
        vector<int> component_of_block(q.size);
        for (unsigned c = 0; c < q.components.size(); ++c)
            for (auto b : q.components[c])
                component_of_block[b] = c;

        vector<int> group_of_block(q.size);
        vector<vector<bool>> group_has_component;
        bool stopped = false;
        auto enumerate = [&](auto & self, int b) -> void {
            if (stopped)
                return;
            if (b == q.size) {
                if (int(group_has_component.size()) == q.size)
                    return;
                vector<int> tau(sigma.size());
                for (unsigned v = 0; v < sigma.size(); ++v)
                    tau[v] = group_of_block[sigma[v]];
                if (! f(tau))
                    stopped = true;
                return;
            }

            for (unsigned g = 0; g < group_has_component.size(); ++g)
                if (! group_has_component[g][component_of_block[b]]) {
                    group_of_block[b] = g;
                    group_has_component[g][component_of_block[b]] = true;
                    self(self, b + 1);
                    group_has_component[g][component_of_block[b]] = false;
                }

            group_of_block[b] = group_has_component.size();
            group_has_component.emplace_back(q.components.size(), false);
            group_has_component.back()[component_of_block[b]] = true;
            self(self, b + 1);
            group_has_component.pop_back();
        };
        enumerate(enumerate, 0);
    }

    // Counting by components: what we need to be able to say that the number of
    // embeddings of H is a function of the numbers of embeddings of its connected
    // components, and of the graphs obtained by gluing them together.
    auto can_count_by_components(const InputGraph & pattern, const HomomorphismParams & params) -> bool
    {
        return params.count_solutions && (! params.enumerate_callback) && params.injectivity == Injectivity::Injective && (! params.induced) && params.pattern_less_constraints.empty() && params.target_occur_less_constraints.empty() && (! pattern.directed());
    }

    // Counts injective embeddings of a disconnected pattern into the target, by
    // counting each connected component separately and correcting for the
    // embeddings in which components overlap by inclusion-exclusion. This is the
    // c = 0 case of the algorithm in Section 3.6 of Jessica Ryan's PhD thesis
    // (Glasgow, 2023), with the solver doing the per-component counts.
    //
    // The vertices being decomposed are numbered 0 .. n - 1 here, and are
    // pattern vertices vertices[0 .. n - 1]. A partition of them is written as the
    // block number of each vertex, with blocks numbered in order of first
    // occurrence, so that equal partitions are equal vectors. A partition is
    // valid if no block holds two vertices of the same connected component of
    // the pattern: an embedding is injective on each component, so these are the
    // only ways in which embeddings of the components can overlap. For a valid
    // partition sigma, the quotient H / sigma glues together the vertices in each
    // block; it has no new loops, because adjacent vertices share a component.
    //
    // Write none(sigma) for the number of injective embeddings of H / sigma, and
    // product(sigma) for the product, over the connected components of H / sigma,
    // of the number of injective embeddings of each. A tuple of embeddings counted
    // by product(sigma) is injective on each component of H / sigma, so it is an
    // injective embedding of H / tau for exactly one tau that coarsens sigma by
    // gluing together blocks from distinct components of H / sigma. Hence
    //
    //     none(sigma) = product(sigma) - sum over such tau != sigma of none(tau),
    //
    // which is Lemma 3.8 of the thesis, and what we want is none of the discrete
    // partition. Every term is non-negative, so if product(sigma) is zero then so
    // is none(tau) for every such tau, and we need not look at them.
    class ComponentCounter
    {
    private:
        const InputGraph & pattern;
        const InputGraph & target;
        const HomomorphismParams & params;

        vector<int> vertices;
        int number_of_components = 0;

        map<string, long> target_with_label;
        map<string, loooong> component_counts;
        map<string, loooong> none_counts;

        unsigned long long solver_calls = 0, coarsenings_considered = 0;

        auto count_component(const string & key, const SmallGraph & component) -> loooong
        {
            if (auto c = component_counts.find(key); c != component_counts.end())
                return c->second;

            auto graph = to_input_graph(pattern, component);
            ++solver_calls;
            auto result = solve_homomorphism_problem(graph, target, params);
            nodes += result.nodes;
            propagations += result.propagations;
            if (! result.complete)
                aborted = true;

            return component_counts.emplace(key, result.solution_count).first->second;
        }

        auto none(const vector<int> & sigma) -> loooong
        {
            auto q = build_quotient(pattern, vertices, sigma);
            if ((! q) || more_blocks_than_targets(*q, target_with_label))
                return 0;

            // none(sigma) depends only on H / sigma up to isomorphism, and a disjoint
            // union is determined up to isomorphism by its components' isomorphism
            // classes, so this is the memo key, and also how we find components that
            // have already been counted
            vector<SmallGraph> components;
            vector<string> keys;
            for (auto & c : q->components) {
                components.push_back(component_graph(pattern, *q, c));
                keys.push_back(canonical_key(components.back()));
            }
            auto sorted_keys = keys;
            sort(sorted_keys.begin(), sorted_keys.end());
            string quotient_key;
            for (auto & k : sorted_keys)
                quotient_key += to_string(k.size()) + ":" + k;

            if (auto c = none_counts.find(quotient_key); c != none_counts.end())
                return c->second;

            loooong total = 1;
            for (unsigned c = 0; c < components.size(); ++c) {
                total *= count_component(keys[c], components[c]);
                if (aborted)
                    return 0;
                if (total == 0)
                    return none_counts.emplace(quotient_key, 0).first->second;
            }

            for_each_coarsening(*q, sigma, [&](const vector<int> & tau) {
                ++coarsenings_considered;
                total -= none(tau);
                return ! aborted;
            });

            if (aborted)
                return 0;
            return none_counts.emplace(quotient_key, total).first->second;
        }

    public:
        unsigned long long nodes = 0, propagations = 0;
        bool aborted = false;

        ComponentCounter(const InputGraph & p, const InputGraph & t, const HomomorphismParams & r, const vector<int> & v) :
            pattern(p),
            target(t),
            params(r),
            vertices(v),
            target_with_label(count_target_labels(p, t))
        {
        }

        auto count() -> loooong
        {
            if (vertices.empty())
                return 1;

            vector<int> discrete(vertices.size());
            for (unsigned v = 0; v < vertices.size(); ++v)
                discrete[v] = v;

            auto q = build_quotient(pattern, vertices, discrete);
            number_of_components = q ? q->components.size() : 0;
            return none(discrete);
        }

        auto stats(std::list<string> & extra_stats) const -> void
        {
            extra_stats.emplace_back("decomposition_components = " + to_string(number_of_components));
            extra_stats.emplace_back("decomposition_solver_calls = " + to_string(solver_calls));
            extra_stats.emplace_back("decomposition_quotients = " + to_string(none_counts.size()));
            extra_stats.emplace_back("decomposition_coarsenings_considered = " + to_string(coarsenings_considered));
        }
    };

    auto count_by_components(const InputGraph & pattern, const InputGraph & target, const HomomorphismParams & params) -> HomomorphismResult
    {
        // Isolated vertices are not worth a solver call each, nor the inclusion-exclusion
        // that gluing them on would cost: every other pattern vertex with label L uses
        // exactly one target vertex labelled L, so the isolated ones labelled L can go on
        // any of the rest. (A looped vertex has a nonzero degree, so is not isolated.)
        vector<int> vertices;
        map<string, long> isolated_with_label, others_with_label;
        auto target_with_label = count_target_labels(pattern, target);
        for (int v = 0; v < pattern.size(); ++v) {
            auto label = pattern.has_vertex_labels() ? string{pattern.vertex_label(v)} : string{};
            if (0 == pattern.degree(v))
                ++isolated_with_label[label];
            else {
                ++others_with_label[label];
                vertices.push_back(v);
            }
        }

        HomomorphismResult result;
        ComponentCounter counter(pattern, target, params, vertices);
        result.solution_count = counter.count();
        for (auto & [label, k] : isolated_with_label)
            result.solution_count *= falling_factorial(target_with_label[label] - others_with_label[label], k);

        result.nodes = counter.nodes;
        result.propagations = counter.propagations;
        result.complete = ! counter.aborted;
        if (counter.aborted)
            result.solution_count = loooong{0};

        long isolated = 0;
        for (auto & [_, k] : isolated_with_label)
            isolated += k;
        result.extra_stats.emplace_back("isolated_pattern_vertices = " + to_string(isolated));
        counter.stats(result.extra_stats);
        return result;
    }

    // The parameters for a subcall that enumerates, as far as they can be copied: the
    // extra shapes cannot be, and are only filtering, so leaving them out changes nothing
    // but speed.
    auto enumerating_params(const HomomorphismParams & params) -> HomomorphismParams
    {
        HomomorphismParams result;
        result.timeout = params.timeout;
        result.start_time = params.start_time;
        result.induced = params.induced;
        result.injectivity = params.injectivity;
        result.count_solutions = true;
        result.value_ordering_heuristic = params.value_ordering_heuristic;
        result.restarts_schedule.reset(params.restarts_schedule->clone());
        result.nogood_size_limit = params.nogood_size_limit;
        result.n_threads = params.n_threads;
        result.delay_thread_creation = params.delay_thread_creation;
        result.triggered_restarts = params.triggered_restarts;
        result.clique_detection = params.clique_detection;
        result.distance3 = params.distance3;
        result.k4 = params.k4;
        result.no_supplementals = params.no_supplementals;
        result.number_of_exact_path_graphs = params.number_of_exact_path_graphs;
        result.clique_size_constraints = params.clique_size_constraints;
        result.clique_size_constraints_on_supplementals = params.clique_size_constraints_on_supplementals;
        result.no_nds = params.no_nds;
        result.staged = params.staged;
        result.staged_first_round_backtracks = params.staged_first_round_backtracks;
        return result;
    }

    // Builds the counting-language statement for the same recursion as ComponentCounter
    // (see dev_docs/decomposition-counting.md, "Compiling Ryan's algorithm"):
    //
    //     Exact_sigma = Prod_sigma minus (the sum of Exact_tau over coarsenings tau),
    //     Prod_sigma  = the product of Leaf_K over the components K of H / sigma,
    //     Leaf_K      = the sum of one cube per embedding of K,
    //
    // where a cube gives every pattern vertex in a block the block's image. Nothing is
    // shared up to isomorphism, because a node is about particular pattern vertices, and
    // isolated vertices are components like any other, so this does far more work than
    // ComponentCounter: it is for certifying a count, not for finding one quickly.
    //
    // Exact_tau is empty whenever some component of H / tau has no embedding, or H / tau
    // has more vertices with some label than the target, and is left out of the sums that
    // would mention it; Exact of the discrete partition is the root, so is always built.
    class StatementBuilder
    {
    private:
        const InputGraph & pattern;
        const InputGraph & target;
        const HomomorphismParams & params;
        CountingStatement & statement;

        vector<int> vertices;
        map<string, long> target_with_label;
        map<vector<vector<int>>, pair<int, bool>> leaves;
        map<vector<int>, optional<int>> exacts;

        auto leaf(const vector<vector<int>> & members, const SmallGraph & component) -> pair<int, bool>
        {
            if (auto l = leaves.find(members); l != leaves.end())
                return l->second;

            auto graph = to_input_graph(pattern, component);
            auto sub_params = enumerating_params(params);
            vector<int> cubes;
            sub_params.enumerate_callback = [&](const VertexToVertexMapping & mapping) -> bool {
                vector<pair<int, int>> literals;
                for (auto & [i, t] : mapping)
                    for (auto & v : members[i])
                        literals.emplace_back(v, t);
                sort(literals.begin(), literals.end());
                cubes.push_back(statement.add_cube(std::move(literals)));
                return true;
            };

            ++solver_calls;
            auto result = solve_homomorphism_problem(graph, target, sub_params);
            nodes += result.nodes;
            propagations += result.propagations;
            if (! result.complete)
                aborted = true;

            pair<int, bool> node;
            if (cubes.empty()) {
                vector<int> scope;
                for (auto & m : members)
                    scope.insert(scope.end(), m.begin(), m.end());
                sort(scope.begin(), scope.end());
                node = pair{statement.add_false(std::move(scope)), true};
            }
            else
                node = pair{statement.add_sum(std::move(cubes)), false};

            return leaves.emplace(members, node).first->second;
        }

        auto exact(const vector<int> & sigma, bool is_root) -> optional<int>
        {
            if (auto e = exacts.find(sigma); e != exacts.end())
                return e->second;

            auto q = build_quotient(pattern, vertices, sigma);
            if (! q)
                return exacts.emplace(sigma, std::nullopt).first->second;

            if (more_blocks_than_targets(*q, target_with_label)) {
                optional<int> result;
                if (is_root)
                    result = statement.add_false(vertices);
                return exacts.emplace(sigma, result).first->second;
            }

            vector<int> children;
            bool any_empty = false;
            for (auto & blocks : q->components) {
                vector<vector<int>> members(blocks.size());
                for (unsigned v = 0; v < sigma.size(); ++v)
                    if (auto b = std::find(blocks.begin(), blocks.end(), sigma[v]); b != blocks.end())
                        members[b - blocks.begin()].push_back(vertices[v]);

                auto [node, empty] = leaf(members, component_graph(pattern, *q, blocks));
                if (aborted)
                    return std::nullopt;
                children.push_back(node);
                any_empty = any_empty || empty;
            }

            if (any_empty && ! is_root)
                return exacts.emplace(sigma, std::nullopt).first->second;

            auto product = children.size() == 1 ? children.front() : statement.add_product(std::move(children));
            if (any_empty)
                return exacts.emplace(sigma, product).first->second;

            vector<int> overlaps;
            for_each_coarsening(*q, sigma, [&](const vector<int> & tau) {
                if (auto e = exact(tau, false))
                    overlaps.push_back(*e);
                return ! aborted;
            });
            if (aborted)
                return std::nullopt;

            optional<int> result = product;
            if (! overlaps.empty())
                result = statement.add_difference(product, overlaps.size() == 1 ? overlaps.front() : statement.add_sum(std::move(overlaps)));
            return exacts.emplace(sigma, result).first->second;
        }

    public:
        unsigned long long nodes = 0, propagations = 0, solver_calls = 0;
        bool aborted = false;

        StatementBuilder(const InputGraph & p, const InputGraph & t, const HomomorphismParams & r, CountingStatement & s) :
            pattern(p),
            target(t),
            params(r),
            statement(s),
            target_with_label(count_target_labels(p, t))
        {
        }

        auto build() -> void
        {
            vector<string> target_names;
            for (int t = 0; t < target.size(); ++t)
                target_names.push_back(target.vertex_name(t));
            auto domain = statement.add_domain(std::move(target_names));
            for (int v = 0; v < pattern.size(); ++v) {
                statement.add_variable(pattern.vertex_name(v), domain);
                vertices.push_back(v);
            }

            if (vertices.empty()) {
                statement.root = statement.add_cube({});
                return;
            }

            vector<int> discrete(vertices.size());
            for (unsigned v = 0; v < vertices.size(); ++v)
                discrete[v] = v;
            if (auto root = exact(discrete, true))
                statement.root = *root;
        }

        auto stats(std::list<string> & extra_stats) const -> void
        {
            extra_stats.emplace_back("decomposition_solver_calls = " + to_string(solver_calls));
            extra_stats.emplace_back("decomposition_leaves = " + to_string(leaves.size()));
            extra_stats.emplace_back("decomposition_quotients = " + to_string(exacts.size()));
            extra_stats.emplace_back("statement_nodes = " + to_string(statement.nodes.size()));
        }
    };

    auto build_statement(const InputGraph & pattern, const InputGraph & target, const HomomorphismParams & params,
        CountingStatement & statement) -> HomomorphismResult
    {
        HomomorphismResult result;
        StatementBuilder builder(pattern, target, params, statement);
        builder.build();

        result.nodes = builder.nodes;
        result.propagations = builder.propagations;
        result.complete = ! builder.aborted;
        if (result.complete) {
            auto count = count_counting_statement(statement);
            result.solution_count = count.count;
            result.extra_stats.emplace_back("statement_obligations = " + to_string(count.obligations.size()));
        }
        builder.stats(result.extra_stats);
        return result;
    }
}

auto gss::solve_sip_by_decomposition(const InputGraph & pattern, const InputGraph & target,
    const HomomorphismParams & params, CountingStatement * statement) -> HomomorphismResult
{
    // each piece would write its own proof over the top of the last one, and none of
    // them would be a proof about the whole instance
    if (params.proof_options)
        throw UnsupportedConfiguration{"Proof logging cannot be used with decomposition"};

    if (statement) {
        if (! can_count_by_components(pattern, params))
            throw UnsupportedConfiguration{"A counting statement needs injective, non-induced, undirected counting, without enumeration or symmetry constraints"};
        return build_statement(pattern, target, params, *statement);
    }

    if (can_count_by_components(pattern, params))
        return count_by_components(pattern, target, params);

    set<int> isolated_pattern_vertices;
    find_removable_isolated_pattern_vertices(pattern, params, isolated_pattern_vertices);

    if (! isolated_pattern_vertices.empty()) {
        InputGraph reduced_pattern(pattern.size() - isolated_pattern_vertices.size(), pattern.has_vertex_labels(),
            pattern.has_edge_labels(), pattern.directed());

        vector<int> original_to_reduced, reduced_to_original;

        for (int i = 0; i < pattern.size(); ++i) {
            if (isolated_pattern_vertices.count(i))
                original_to_reduced.push_back(-1);
            else {
                original_to_reduced.push_back(reduced_to_original.size());
                reduced_to_original.push_back(i);
            }
        }

        for (int i = 0; i < pattern.size(); ++i) {
            auto r_i = original_to_reduced.at(i);
            if (r_i == -1)
                continue;

            for (int j = 0; j < pattern.size(); ++j) {
                auto r_j = original_to_reduced.at(j);
                if (r_j == -1)
                    continue;

                // Mirror the edge as it is in the pattern: dropping its label would
                // silently turn a labelled pattern into one demanding empty labels,
                // and add_directed_edge() would mark an undirected reduced pattern
                // directed. (find_removable_isolated_pattern_vertices() declines to
                // decompose a directed pattern, so only the second arm runs today.)
                if (pattern.adjacent(i, j)) {
                    if (pattern.directed())
                        reduced_pattern.add_directed_edge(r_i, r_j, pattern.edge_label(i, j));
                    else
                        reduced_pattern.add_edge(r_i, r_j, pattern.edge_label(i, j));
                }
            }
        }
        auto result = solve_homomorphism_problem(reduced_pattern, target, params);

        result.extra_stats.emplace_back("isolated_pattern_vertices = " + to_string(isolated_pattern_vertices.size()));

        // fix up the result to be in the non-decomposed form
        if (! result.mapping.empty()) {
            // need to re-add isolated vertices, arbitrarily. what's available?
            set<int> t_avail;
            for (int i = 0; i < target.size(); ++i)
                t_avail.emplace(i);

            // convert back to original indices, and remove used vertices
            auto sub_mapping = result.mapping;
            result.mapping.clear();
            for (auto & [p, t] : sub_mapping) {
                t_avail.erase(t);
                result.mapping.emplace(reduced_to_original.at(p), t);
            }

            // add in isolated vertices
            for (auto & p : isolated_pattern_vertices) {
                result.mapping.emplace(p, *t_avail.begin());
                t_avail.erase(t_avail.begin());
            }
        }

        // fix up the solution count
        if (params.count_solutions && result.solution_count > 0) {
            loooong unmapped_target_vertices = target.size() - reduced_pattern.size();
            loooong solution_multiplier = n_choose_k<loooong>(unmapped_target_vertices, isolated_pattern_vertices.size());
            loooong isolated_symmetry_multiplier = factorial<loooong>(isolated_pattern_vertices.size());
            result.solution_count *= solution_multiplier * isolated_symmetry_multiplier;
        }

        return result;
    }
    else
        return solve_homomorphism_problem(pattern, target, params);
}
