#include <gss/configuration.hh>
#include <gss/counting_statement.hh>
#include <gss/homomorphism.hh>
#include <gss/innards/decomposition_proof.hh>
#include <gss/innards/homomorphism_model.hh>
#include <gss/innards/proof.hh>
#include <gss/loooong.hh>
#include <gss/sip_decomposer.hh>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <ostream>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

using namespace gss;
using namespace gss::innards;

using std::make_unique;
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
        // the list of target vertices each vertex may go to, as a list number (see
        // ComponentCounter), or empty if there are no lists
        vector<int> lists;
        // only for searching under proof, where they must be the pattern's names
        vector<string> names;
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
        if (! g.lists.empty()) {
            result += "|";
            for (auto v : order)
                result += to_string(g.lists[v]) + ",";
        }
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

    // A key that is equal for two graphs if and only if they are isomorphic, and the
    // ordering of the vertices that produced it: order[p] is the vertex in position p,
    // so for two graphs with equal keys, the vertices in the same position correspond.
    struct CanonicalForm
    {
        string key;
        vector<int> order;
    };

    // Computed by colour refinement and then trying every ordering that respects the
    // colours. When that is too many orderings, we give up on isomorphism and use the
    // graph as it is, which is still correct, just less good as a cache key.
    auto canonical_form(const SmallGraph & g) -> CanonicalForm
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

        // a vertex's list is part of its colour, like its label: an isomorphism must
        // preserve both
        vector<string> descriptions(g.size);
        for (int v = 0; v < g.size; ++v) {
            descriptions[v] = g.labels.empty() ? string{} : g.labels[v];
            if (! g.lists.empty())
                descriptions[v] += "|" + to_string(g.lists[v]);
        }
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
                    return CanonicalForm{"exact;" + encode(g, identity), identity};
                }

        optional<string> best;
        vector<int> best_order, order;
        auto try_orderings = [&](auto & self, unsigned c) -> void {
            if (c == classes.size()) {
                auto e = encode(g, order);
                if ((! best) || e < *best) {
                    best = e;
                    best_order = order;
                }
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

        return CanonicalForm{"canonical;" + *best, best_order};
    }

    // The vertices of each connected component of the pattern restricted to vertices,
    // which must contain every neighbour of each of its members.
    auto pattern_components(const InputGraph & pattern, const vector<int> & vertices) -> vector<vector<int>>
    {
        vector<int> parent(pattern.size());
        for (int v = 0; v < pattern.size(); ++v)
            parent[v] = v;
        auto find = [&](int v) {
            while (parent[v] != v)
                v = parent[v] = parent[parent[v]];
            return v;
        };
        for (auto v : vertices)
            for (auto w : vertices)
                if (pattern.adjacent(v, w))
                    parent[find(v)] = find(w);

        map<int, vector<int>> components;
        for (auto v : vertices)
            components[find(v)].push_back(v);
        vector<vector<int>> result;
        for (auto & [_, c] : components)
            result.push_back(c);
        return result;
    }

    // How many partitions of the vertices of these components put no two vertices of
    // one component in the same block: the coarsenings that the inclusion-exclusion can
    // visit. Adding a component of size s to a partition with b blocks puts j of its
    // vertices into distinct existing blocks and opens s - j new ones.
    auto number_of_valid_partitions(const vector<vector<int>> & components) -> double
    {
        vector<double> with_blocks{1};
        for (auto & component : components) {
            int s = component.size();
            vector<double> next(with_blocks.size() + s, 0);
            for (int b = 0; b < int(with_blocks.size()); ++b) {
                double ways = 1;
                for (int j = 0; j <= s && j <= b; ++j) {
                    next[b + s - j] += with_blocks[b] * ways;
                    // from C(s, j) b! / (b - j)! to C(s, j + 1) b! / (b - j - 1)!
                    ways *= double(s - j) / (j + 1) * (b - j);
                }
            }
            with_blocks = std::move(next);
        }
        double result = 0;
        for (auto w : with_blocks)
            result += w;
        return result;
    }

    // log(exp(a) + exp(b)), for adding up estimates that are too big for a double
    auto log_add(double a, double b) -> double
    {
        auto hi = std::max(a, b), lo = std::min(a, b);
        if (std::isinf(hi) || std::isinf(lo))
            return hi;
        return hi + std::log1p(std::exp(lo - hi));
    }

    // The hubs are the target vertices that pattern vertices are guessed onto, rather than
    // found by search: the c target vertices of highest degree, ties going to the lower
    // numbered vertex. In Ryan's terms, what is left has maximum degree Delta. Any choice
    // of hubs gives the right count; it only affects how long that takes. If hubs is
    // given, it is c. Otherwise c is chosen to minimise a rough estimate of the work,
    //
    //     (number of guesses) x P + (embeddings of the components and the glued one).
    //
    // Each guess runs its own inclusion-exclusion, over at most P coarsenings, the
    // number of valid partitions of the m vertices. The guesses put any j <= c of the m
    // vertices onto distinct hubs, which is at most the sum over j of C(m, j) c! / (c - j)!
    // ways. That only grows with c, so once it alone costs more than the best so far, no
    // larger c can do better.
    //
    // The embeddings are those that enumeration visits: of each component, and of the
    // biggest glued component, which glues one vertex of greatest degree from each
    // component together. A graph with e edges and degrees g(v) has about
    //
    //     product over v of F(g(v)) / F(1)^e,  where F(j) = sum over non-hubs t of d(t)_j,
    //
    // embeddings into a random graph with the degrees d(t) of what is left once the hubs
    // are gone, writing d_j for the falling factorial. That is exact for a star, which is
    // what makes a heavy-tailed degree distribution worth guessing on and a regular one
    // not, and it charges little for a dense pattern in a sparse target, where search
    // visits few embeddings, and where guessing gains little.
    auto choose_hubs(const InputGraph & pattern, const InputGraph & target, const vector<int> & vertices, optional<int> hubs) -> vector<int>
    {
        int n = target.size();
        vector<int> by_degree(n);
        for (int t = 0; t < n; ++t)
            by_degree[t] = t;
        std::stable_sort(by_degree.begin(), by_degree.end(), [&](int a, int b) { return target.degree(a) > target.degree(b); });

        // the lists are made from undirected adjacency
        if (hubs && *hubs > 0 && target.directed())
            throw UnsupportedConfiguration{"Decomposition cannot yet use hubs with a directed target"};
        if (hubs && (*hubs < 0 || *hubs > n))
            throw UnsupportedConfiguration{"The number of hubs must be between 0 and the number of target vertices"};

        int c = 0;
        if (hubs)
            c = *hubs;
        else if (! (target.directed() || vertices.empty())) {
            auto components = pattern_components(pattern, vertices);
            double m = vertices.size();
            double log_partitions = std::log(number_of_valid_partitions(components));

            // the degrees of each component's vertices, never counting a loop, its
            // number of edges, and the same for the glued component
            struct Shape
            {
                vector<int> degrees;
                long edges = 0;
            };
            vector<Shape> shapes;
            Shape glued;
            glued.degrees.push_back(0);
            for (auto & component : components) {
                Shape shape;
                for (auto v : component) {
                    int g = 0;
                    for (auto w : component)
                        if (v != w && pattern.adjacent(v, w))
                            ++g;
                    shape.degrees.push_back(g);
                    shape.edges += g;
                }
                shape.edges /= 2;
                auto greatest = std::max_element(shape.degrees.begin(), shape.degrees.end());
                glued.degrees.front() += *greatest;
                for (auto i = shape.degrees.begin(); i != shape.degrees.end(); ++i)
                    if (i != greatest)
                        glued.degrees.push_back(*i);
                glued.edges += shape.edges;
                shapes.push_back(std::move(shape));
            }
            if (shapes.size() > 1)
                shapes.push_back(std::move(glued));

            set<int> wanted{1};
            for (auto & shape : shapes)
                wanted.insert(shape.degrees.begin(), shape.degrees.end());

            vector<vector<int>> neighbours(n);
            target.for_each_edge([&](int f, int t, std::string_view) {
                if (f != t)
                    neighbours[f].push_back(t);
            });
            vector<int> degree(n);
            vector<double> log_factorial(1, 0);
            for (int t = 0; t < n; ++t) {
                degree[t] = neighbours[t].size();
                while (int(log_factorial.size()) <= degree[t])
                    log_factorial.push_back(log_factorial.back() + std::log(double(log_factorial.size())));
            }
            vector<bool> removed(n, false);

            // log F(j) for each wanted j, from scratch, because subtracting from F(j) as
            // hubs go would lose all precision for large j
            auto log_enumeration = [&]() {
                map<int, double> log_f;
                for (auto j : wanted) {
                    double most = -std::numeric_limits<double>::infinity();
                    for (int t = 0; t < n; ++t)
                        if ((! removed[t]) && degree[t] >= j)
                            most = std::max(most, log_factorial[degree[t]] - log_factorial[degree[t] - j]);
                    double sum = 0;
                    if (! std::isinf(most))
                        for (int t = 0; t < n; ++t)
                            if ((! removed[t]) && degree[t] >= j)
                                sum += std::exp(log_factorial[degree[t]] - log_factorial[degree[t] - j] - most);
                    log_f[j] = std::isinf(most) ? most : most + std::log(sum);
                }

                double result = -std::numeric_limits<double>::infinity();
                for (auto & shape : shapes) {
                    double estimate = -double(shape.edges) * log_f[1];
                    for (auto g : shape.degrees)
                        if (std::isinf(log_f[g])) {
                            estimate = -std::numeric_limits<double>::infinity();
                            break;
                        }
                        else
                            estimate += log_f[g];
                    result = log_add(result, estimate);
                }
                return result;
            };

            double best = log_add(log_partitions, log_enumeration());
            for (int h = 1; h < n; ++h) {
                auto hub = by_degree[h - 1];
                removed[hub] = true;
                for (auto u : neighbours[hub])
                    if (! removed[u])
                        --degree[u];

                double log_guesses = 0, log_term = 0;
                for (int j = 1; j <= h && j <= m; ++j) {
                    log_term += std::log((m - j + 1) / j * (h - j + 1));
                    log_guesses = log_add(log_guesses, log_term);
                }
                if (log_guesses + log_partitions >= best)
                    break;
                auto cost = log_add(log_guesses + log_partitions, log_enumeration());
                if (cost < best) {
                    best = cost;
                    c = h;
                }
            }
        }

        vector<int> result(by_degree.begin(), by_degree.begin() + c);
        sort(result.begin(), result.end());
        return result;
    }

    // The plain solver's checks before search look at the whole instance, so they can
    // rule it out where every component would fit on its own: by comparing degree
    // sequences, for instance. Decomposing never sees the whole instance, so ask first.
    auto rejected_before_search(const InputGraph & pattern, const InputGraph & target, const HomomorphismParams & params) -> bool
    {
        innards::HomomorphismModel model(target, pattern, params, nullptr, nullptr);
        return ! model.prepare();
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
    // algorithm in Section 3.6 of Jessica Ryan's PhD thesis (Glasgow, 2023), with the
    // solver doing the per-component counts.
    //
    // First, some pattern vertices are guessed onto hubs (see choose_hubs), in every
    // possible way: injectively, preserving labels, loops and the edges between them.
    // Every embedding puts exactly one set of pattern vertices onto the hubs, so the
    // count is the sum over the guesses. In each guess, the other pattern vertices go
    // to target vertices that are not hubs, and each gets a list of the ones it may go
    // to: those adjacent to the hubs that its neighbours were guessed onto. With no
    // hubs, there is one guess, which puts nothing on a hub, and there are no lists.
    //
    // Then, for each guess, the vertices being decomposed are the ones left over,
    // numbered 0 .. n - 1 here, and are pattern vertices vertices[0 .. n - 1]. A
    // partition of them is written as the block number of each vertex, with blocks
    // numbered in order of first occurrence, so that equal partitions are equal
    // vectors. A partition is valid if no block holds two vertices of the same
    // connected component of what is left of the pattern: an embedding is injective on
    // each component, so these are the only ways in which embeddings of the components
    // can overlap. For a valid partition sigma, the quotient H / sigma glues together
    // the vertices in each block; it has no new loops, because adjacent vertices share
    // a component. A block's list is the intersection of its members' lists.
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
    // is none(tau) for every such tau, and we need not look at them. All of this
    // holds with lists, which only say which embeddings exist.
    class ComponentCounter
    {
    private:
        const InputGraph & pattern;
        const InputGraph & target;
        const HomomorphismParams & params;

        // With any hubs, lists are numbered: list_values[l] is the target vertices of
        // list l, in order. Isomorphic graphs whose vertices have the same lists can
        // share a count, so equal lists must get equal numbers.
        vector<int> hubs;
        vector<bool> is_hub;
        vector<vector<int>> list_values;
        map<vector<int>, int> list_numbers;
        map<pair<int, int>, int> intersections;
        map<pair<int, string>, int> hub_neighbourhoods;
        map<int, int> non_hubs_with_degree_at_least;

        // degrees and neighbours, never counting a loop
        vector<vector<int>> pattern_neighbours;
        vector<int> target_degree, non_hub_degree;

        // the vertices being decomposed, for the guess at hand, and their lists, which
        // are empty if there are no hubs
        vector<int> vertices;
        vector<int> vertex_lists;
        int number_of_components = 0;

        // how many target vertices carry each label, with every vertex of an unlabelled
        // target under the empty label, and how many of those are not hubs
        map<string, long> target_with_label, non_hubs_with_label;

        map<string, loooong> component_counts;
        map<string, loooong> none_counts;

        unsigned long long solver_calls = 0, coarsenings_considered = 0, guesses = 0;

        auto list_number(vector<int> values) -> int
        {
            auto [where, inserted] = list_numbers.emplace(std::move(values), list_values.size());
            if (inserted)
                list_values.push_back(where->first);
            return where->second;
        }

        auto intersection(int a, int b) -> int
        {
            if (a == b)
                return a;
            pair<int, int> key = std::minmax(a, b);
            if (auto i = intersections.find(key); i != intersections.end())
                return i->second;
            vector<int> both;
            std::set_intersection(list_values[a].begin(), list_values[a].end(), list_values[b].begin(), list_values[b].end(), std::back_inserter(both));
            auto result = list_number(std::move(both));
            intersections.emplace(key, result);
            return result;
        }

        auto pattern_edge_label(int v, int w) const -> string
        {
            return pattern.has_edge_labels() ? string{pattern.edge_label(v, w)} : string{};
        }

        auto target_edge_label(int t, int u) const -> string
        {
            return pattern.has_edge_labels() ? string{target.edge_label(t, u)} : string{};
        }

        // The non-hubs with at least d neighbours that are also non-hubs.
        auto with_non_hub_degree_at_least(int d) -> int
        {
            if (auto n = non_hubs_with_degree_at_least.find(d); n != non_hubs_with_degree_at_least.end())
                return n->second;
            vector<int> values;
            for (int t = 0; t < target.size(); ++t)
                if ((! is_hub[t]) && non_hub_degree[t] >= d)
                    values.push_back(t);
            auto result = list_number(std::move(values));
            non_hubs_with_degree_at_least.emplace(d, result);
            return result;
        }

        // The non-hubs that a pattern vertex can go to if it has a neighbour on hub h,
        // joined by an edge with this label.
        auto hub_neighbourhood(int h, const string & label) -> int
        {
            if (auto n = hub_neighbourhoods.find({h, label}); n != hub_neighbourhoods.end())
                return n->second;
            vector<int> values;
            for (int t = 0; t < target.size(); ++t)
                if ((! is_hub[t]) && target.adjacent(t, h) && target_edge_label(t, h) == label)
                    values.push_back(t);
            auto result = list_number(std::move(values));
            hub_neighbourhoods.emplace(pair{h, label}, result);
            return result;
        }

        // Can pattern vertex v go on hub h, given the guess so far? Its neighbours need
        // distinct neighbours of h, because embeddings are injective.
        auto fits_on_hub(int v, int h, const vector<pair<int, int>> & cube) const -> bool
        {
            if (pattern.has_vertex_labels() && pattern.vertex_label(v) != target.vertex_label(h))
                return false;
            if (int(pattern_neighbours[v].size()) > target_degree[h])
                return false;
            if (pattern.adjacent(v, v) && ! (target.adjacent(h, h) && pattern_edge_label(v, v) == target_edge_label(h, h)))
                return false;
            for (auto & [w, g] : cube)
                if (pattern.adjacent(v, w) && ! (target.adjacent(h, g) && pattern_edge_label(v, w) == target_edge_label(h, g)))
                    return false;
            return true;
        }

        // Calls f(cube, rest, lists) for each guess of which candidates go on which hubs,
        // where cube is the guess, as (pattern vertex, hub) in order, rest is the other
        // candidates, and lists is their lists. A vertex not on a hub must go to a non-hub
        // adjacent to the hubs of its neighbours that are, and its other neighbours need
        // distinct non-hub neighbours of where it goes. A guess that leaves some vertex an
        // empty list has no embeddings, and is skipped. The candidates must be in order,
        // and must contain every neighbour of each of their members.
        template <typename F_>
        auto for_each_guess(const vector<int> & candidates, const F_ & f) -> void
        {
            if (hubs.empty()) {
                ++guesses;
                f(vector<pair<int, int>>{}, candidates, vector<int>{});
                return;
            }

            vector<bool> on_hub(pattern.size(), false), used(target.size(), false);
            vector<pair<int, int>> cube;
            auto guess = [&](auto & self, unsigned i) -> void {
                if (out_of_time())
                    return;
                if (i == candidates.size()) {
                    vector<int> rest, lists;
                    for (auto v : candidates)
                        if (! on_hub[v]) {
                            int off_hubs = 0;
                            for (auto w : pattern_neighbours[v])
                                if (! on_hub[w])
                                    ++off_hubs;
                            int list = with_non_hub_degree_at_least(off_hubs);
                            for (auto & [s, h] : cube)
                                if (pattern.adjacent(v, s))
                                    list = intersection(list, hub_neighbourhood(h, pattern_edge_label(v, s)));
                            if (list_values[list].empty())
                                return;
                            rest.push_back(v);
                            lists.push_back(list);
                        }
                    ++guesses;
                    f(cube, rest, lists);
                    return;
                }

                auto v = candidates[i];
                self(self, i + 1);
                for (auto h : hubs)
                    if ((! used[h]) && fits_on_hub(v, h, cube)) {
                        on_hub[v] = used[h] = true;
                        cube.emplace_back(v, h);
                        self(self, i + 1);
                        cube.pop_back();
                        on_hub[v] = used[h] = false;
                    }
            };
            guess(guess, 0);
        }

        // The guesses and the coarsenings can be numerous enough to outlast the timeout
        // without any solver call noticing, so they check it themselves.
        auto out_of_time() -> bool
        {
            if (params.timeout && params.timeout->should_abort())
                aborted = true;
            return aborted;
        }

        auto set_vertices(const vector<int> & v, const vector<int> & lists) -> void
        {
            vertices = v;
            vertex_lists = lists;
        }

        auto discrete() const -> vector<int>
        {
            vector<int> result(vertices.size());
            for (unsigned v = 0; v < vertices.size(); ++v)
                result[v] = v;
            return result;
        }

        // The pattern vertices in each block of sigma, in order, with the block's list, or
        // -1 if there are no lists.
        using Blocks = vector<pair<vector<int>, int>>;

        auto blocks_of(const vector<int> & sigma) -> Blocks
        {
            Blocks blocks(*std::max_element(sigma.begin(), sigma.end()) + 1, pair{vector<int>{}, -1});
            for (unsigned v = 0; v < vertices.size(); ++v) {
                auto & [members, list] = blocks[sigma[v]];
                members.push_back(vertices[v]);
                if (! vertex_lists.empty())
                    list = (list == -1) ? vertex_lists[v] : intersection(list, vertex_lists[v]);
            }
            return blocks;
        }

        struct Quotient
        {
            int size = 0;
            map<pair<int, int>, string> edges;
            vector<string> labels;
            vector<int> lists;
            vector<vector<int>> components;
        };

        auto build_quotient(const vector<int> & sigma) -> optional<Quotient>
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

            // a block whose members' lists have nothing in common cannot go anywhere
            if (! vertex_lists.empty())
                for (auto & [_, list] : blocks_of(sigma)) {
                    if (list_values[list].empty())
                        return std::nullopt;
                    q.lists.push_back(list);
                }

            // two pattern edges can be glued into one quotient edge, in which case the
            // target edge must carry both labels, which it cannot unless they agree
            for (unsigned v = 0; v < vertices.size(); ++v)
                for (unsigned w = v; w < vertices.size(); ++w)
                    if (pattern.adjacent(vertices[v], vertices[w])) {
                        auto label = pattern_edge_label(vertices[v], vertices[w]);
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

        // An injective embedding of H / sigma needs a different target vertex for each
        // block, carrying that block's label, so there is none if some label has more
        // blocks than target vertices. Without this, a pattern bigger than its target
        // explores every way of gluing itself together, all of them empty. Blocks never
        // go on hubs, so only the other target vertices are available.
        auto more_blocks_than_targets(const Quotient & q, const map<string, long> & available) const -> bool
        {
            map<string, long> blocks_with_label;
            for (int b = 0; b < q.size; ++b) {
                auto label = q.labels.empty() ? string{} : q.labels[b];
                auto t = available.find(label);
                if (++blocks_with_label[label] > (t == available.end() ? 0 : t->second))
                    return true;
            }
            return false;
        }

        auto component_graph(const Quotient & q, const vector<int> & blocks) const -> SmallGraph
        {
            map<int, int> index_of_block;
            for (unsigned i = 0; i < blocks.size(); ++i)
                index_of_block.emplace(blocks[i], i);

            SmallGraph component;
            component.size = blocks.size();
            if (pattern.has_vertex_labels())
                for (auto b : blocks)
                    component.labels.push_back(q.labels[b]);
            if (! q.lists.empty())
                for (auto b : blocks)
                    component.lists.push_back(q.lists[b]);
            for (auto & [e, label] : q.edges)
                if (auto a = index_of_block.find(e.first); a != index_of_block.end())
                    component.edges.emplace(std::minmax(a->second, index_of_block.at(e.second)), label);
            return component;
        }

        // Each vertex of the component may only go to the vertices on its list.
        auto restriction(const SmallGraph & component) const -> std::function<auto(int, int)->bool>
        {
            auto allowed = std::make_shared<vector<vector<bool>>>(component.size, vector<bool>(target.size(), false));
            for (int v = 0; v < component.size; ++v)
                for (auto t : list_values[component.lists[v]])
                    (*allowed)[v][t] = true;
            return [allowed](int p, int t) { return (*allowed)[p][t]; };
        }

        auto count_component(const string & key, const SmallGraph & component) -> loooong
        {
            if (auto c = component_counts.find(key); c != component_counts.end())
                return c->second;

            ++solver_calls;
            HomomorphismResult result;
            if (component.lists.empty())
                result = solve_homomorphism_problem(as_input_graph(component), target, params);
            else {
                auto restricted = enumeration_params();
                restricted.restrict_initial_domains = restriction(component);
                result = solve_homomorphism_problem(as_input_graph(component), target, restricted);
            }
            nodes += result.nodes;
            propagations += result.propagations;
            if (! result.complete)
                aborted = true;

            return component_counts.emplace(key, result.solution_count).first->second;
        }

        auto as_input_graph(const SmallGraph & component) const -> InputGraph
        {
            InputGraph graph(component.size, pattern.has_vertex_labels(), pattern.has_edge_labels());
            for (int i = 0; i < component.size && ! component.names.empty(); ++i)
                graph.set_vertex_name(i, component.names[i]);
            for (int i = 0; i < component.size && pattern.has_vertex_labels(); ++i)
                graph.set_vertex_label(i, component.labels[i]);
            for (auto & [e, label] : component.edges)
                graph.add_edge(e.first, e.second, label);
            return graph;
        }

        // Calls f(tau) for each coarsening tau != sigma of sigma that glues each block
        // to blocks from other components only: block b goes into an existing group that
        // has nothing from its component yet, or starts a new one. Groups are opened in
        // block order, so tau is already numbered in order of first occurrence.
        template <typename F_>
        auto for_each_coarsening(const Quotient & q, const vector<int> & sigma, const F_ & f) -> void
        {
            vector<int> component_of_block(q.size);
            for (unsigned c = 0; c < q.components.size(); ++c)
                for (auto b : q.components[c])
                    component_of_block[b] = c;

            vector<int> group_of_block(q.size);
            vector<vector<bool>> group_has_component;
            auto enumerate = [&](auto & self, int b) -> void {
                if (out_of_time())
                    return;
                if (b == q.size) {
                    if (int(group_has_component.size()) == q.size)
                        return;
                    ++coarsenings_considered;
                    vector<int> tau(sigma.size());
                    for (unsigned v = 0; v < sigma.size(); ++v)
                        tau[v] = group_of_block[sigma[v]];
                    f(tau);
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

        auto none(const vector<int> & sigma) -> loooong
        {
            auto q = build_quotient(sigma);
            if ((! q) || more_blocks_than_targets(*q, non_hubs_with_label))
                return 0;

            // none(sigma) depends only on H / sigma up to isomorphism, and a disjoint
            // union is determined up to isomorphism by its components' isomorphism
            // classes, so this is the memo key, and also how we find components that
            // have already been counted. Lists are part of the keys, so this holds across
            // guesses too.
            vector<SmallGraph> components;
            vector<string> keys;
            for (auto & c : q->components) {
                components.push_back(component_graph(*q, c));
                keys.push_back(canonical_form(components.back()).key);
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

            for_each_coarsening(*q, sigma, [&](const vector<int> & tau) { total -= none(tau); });

            if (aborted)
                return 0;
            return none_counts.emplace(quotient_key, total).first->second;
        }

        // Writing a statement (see dev_docs/counting-statements.md) instead of just
        // counting. The recursion is the same, but every quotient gets its own node,
        // over the actual pattern vertices, so isomorphic quotients cannot share a memo
        // entry: the memo is keyed on the blocks themselves, with their lists, which
        // also lets guesses share nodes where their lists agree. Isomorphic components
        // still share one enumeration, which is renamed into place.
        CountingStatementWriter * statement = nullptr;
        map<string, vector<vector<int>>> embeddings_by_key;
        map<Blocks, int> leaf_nodes;
        map<Blocks, optional<pair<int, loooong>>> exact_nodes;

        // What a proof needs to know about the statement's structure, by node.
        struct ExactNode
        {
            int product;
            vector<int> leaves;
            optional<int> overlaps;
            vector<int> children;
        };
        map<int, ExactNode> exact_structure;
        map<int, vector<int>> leaf_cubes;
        set<int> glued_leaves;
        innards::DecompositionProof * decomposition_proof = nullptr;

        auto enumeration_params() const -> HomomorphismParams
        {
            HomomorphismParams result;
            result.timeout = params.timeout;
            result.start_time = params.start_time;
            result.restarts_schedule = make_unique<NoRestartsSchedule>();
            result.count_solutions = true;
            result.value_ordering_heuristic = params.value_ordering_heuristic;
            result.nogood_size_limit = params.nogood_size_limit;
            result.no_supplementals = params.no_supplementals;
            result.no_nds = params.no_nds;
            result.number_of_exact_path_graphs = params.number_of_exact_path_graphs;
            result.distance3 = params.distance3;
            result.k4 = params.k4;
            result.clique_size_constraints = params.clique_size_constraints;
            result.clique_size_constraints_on_supplementals = params.clique_size_constraints_on_supplementals;
            return result;
        }

        // Every injective embedding of the component, each as the target vertex of each
        // position in its canonical form, so that an isomorphic component can reuse them.
        auto embeddings(const SmallGraph & component, const CanonicalForm & form) -> const vector<vector<int>> &
        {
            if (auto e = embeddings_by_key.find(form.key); e != embeddings_by_key.end())
                return e->second;

            vector<vector<int>> found;
            auto child_params = enumeration_params();
            if (! component.lists.empty())
                child_params.restrict_initial_domains = restriction(component);
            child_params.enumerate_callback = [&](const VertexToVertexMapping & mapping) -> bool {
                vector<int> embedding(form.order.size());
                for (unsigned p = 0; p < form.order.size(); ++p)
                    embedding[p] = mapping.at(form.order[p]);
                found.push_back(std::move(embedding));
                return true;
            };

            ++solver_calls;
            auto result = solve_homomorphism_problem(as_input_graph(component), target, child_params);
            nodes += result.nodes;
            propagations += result.propagations;
            if (! result.complete)
                aborted = true;

            return embeddings_by_key.emplace(form.key, std::move(found)).first->second;
        }

        // The leaf for one component of a quotient: an or of one cube per embedding,
        // with every pattern vertex in a block taking that block's image.
        auto leaf(const Blocks & members, const vector<vector<int>> & found, const CanonicalForm & form) -> int
        {
            if (auto l = leaf_nodes.find(members); l != leaf_nodes.end())
                return l->second;

            vector<int> position(form.order.size());
            for (unsigned p = 0; p < form.order.size(); ++p)
                position[form.order[p]] = p;

            vector<int> cubes;
            for (auto & embedding : found) {
                vector<pair<int, int>> assignment;
                for (unsigned i = 0; i < members.size(); ++i)
                    for (auto v : members[i].first)
                        assignment.emplace_back(v, embedding[position[i]]);
                sort(assignment.begin(), assignment.end());
                cubes.push_back(statement->add_cube(assignment));
            }
            auto node = statement->add_or(cubes);
            leaf_cubes.emplace(node, cubes);
            if (std::any_of(members.begin(), members.end(), [](const auto & m) { return m.first.size() > 1; }))
                glued_leaves.insert(node);
            return leaf_nodes.emplace(members, node).first->second;
        }

        // Does H / sigma have any embeddings, as far as exact can tell without writing
        // anything? Exactly when this holds, exact(sigma) gives a node.
        auto might_embed(const vector<int> & sigma) -> bool
        {
            auto q = build_quotient(sigma);
            if ((! q) || more_blocks_than_targets(*q, non_hubs_with_label))
                return false;
            for (auto & c : q->components) {
                auto component = component_graph(*q, c);
                if (embeddings(component, canonical_form(component)).empty() || aborted)
                    return false;
            }
            return true;
        }

        // The node for the embeddings of H / sigma, and its count, or nothing if
        // product(sigma) is zero, in which case no node is written at all.
        auto exact(const vector<int> & sigma) -> optional<pair<int, loooong>>
        {
            auto blocks = blocks_of(sigma);
            sort(blocks.begin(), blocks.end());
            if (auto e = exact_nodes.find(blocks); e != exact_nodes.end())
                return e->second;

            // a pigeonholed quotient has no node, which leaves it out of its parents' sums:
            // nothing is lost, because it has no embeddings to count
            auto q = build_quotient(sigma);
            if ((! q) || more_blocks_than_targets(*q, non_hubs_with_label))
                return exact_nodes.emplace(blocks, std::nullopt).first->second;

            // enumerate every component before writing any of them, so that nothing is
            // written for a quotient with an unembeddable component
            vector<SmallGraph> components;
            vector<CanonicalForm> forms;
            vector<const vector<vector<int>> *> found;
            for (auto & c : q->components) {
                components.push_back(component_graph(*q, c));
                forms.push_back(canonical_form(components.back()));
                found.push_back(&embeddings(components.back(), forms.back()));
                if (aborted)
                    return std::nullopt;
                if (found.back()->empty())
                    return exact_nodes.emplace(blocks, std::nullopt).first->second;
            }

            auto members_of_block = blocks_of(sigma);
            vector<int> leaves;
            loooong product = 1;
            for (unsigned c = 0; c < q->components.size(); ++c) {
                Blocks members;
                for (auto b : q->components[c])
                    members.push_back(members_of_block[b]);
                leaves.push_back(leaf(members, *found[c], forms[c]));
                product *= loooong{static_cast<unsigned long>(found[c]->size())};
            }
            int product_node = leaves.size() == 1 ? leaves.front() : statement->add_and(leaves);

            vector<int> overlaps;
            loooong overlap = 0;
            for_each_coarsening(*q, sigma, [&](const vector<int> & tau) {
                if (auto e = exact(tau)) {
                    overlaps.push_back(e->first);
                    overlap += e->second;
                }
            });
            if (aborted)
                return std::nullopt;

            optional<int> overlaps_node;
            if (! overlaps.empty())
                overlaps_node = statement->add_or(overlaps);
            int node = overlaps_node ? statement->add_minus(product_node, *overlaps_node) : product_node;
            exact_structure.emplace(node, ExactNode{product_node, leaves, overlaps_node, overlaps});
            return exact_nodes.emplace(blocks, pair{node, product - overlap}).first->second;
        }

    public:
        unsigned long long nodes = 0, propagations = 0;
        bool aborted = false;

        // whether to start with the plain solver's checks on the whole instance, and
        // whether they ruled it out
        bool check_whole_instance = true, rejected = false;

        auto rejected_before_search() -> bool
        {
            rejected = check_whole_instance && ::rejected_before_search(pattern, target, params);
            return rejected;
        }

        ComponentCounter(const InputGraph & p, const InputGraph & t, const HomomorphismParams & r, const vector<int> & h) :
            pattern(p),
            target(t),
            params(r),
            hubs(h),
            is_hub(t.size(), false)
        {
            for (auto v : hubs)
                is_hub[v] = true;

            pattern_neighbours.resize(pattern.size());
            for (int v = 0; v < pattern.size(); ++v)
                for (int w = 0; w < pattern.size(); ++w)
                    if (v != w && pattern.adjacent(v, w))
                        pattern_neighbours[v].push_back(w);
            target_degree.resize(target.size());
            non_hub_degree.resize(target.size());
            target.for_each_edge([&](int f, int t, std::string_view) {
                if (f != t) {
                    ++target_degree[f];
                    if (! is_hub[t])
                        ++non_hub_degree[f];
                }
            });
            for (int w = 0; w < target.size(); ++w) {
                auto label = pattern.has_vertex_labels() ? string{target.vertex_label(w)} : string{};
                ++target_with_label[label];
                if (! is_hub[w])
                    ++non_hubs_with_label[label];
            }
        }

        // The number of embeddings of the pattern restricted to these vertices, which
        // must be in order, and must contain every neighbour of each of their members.
        auto count(const vector<int> & candidates) -> loooong
        {
            if (candidates.empty())
                return 1;

            number_of_components = pattern_components(pattern, candidates).size();
            loooong total = 0;
            for_each_guess(candidates, [&](const vector<pair<int, int>> &, const vector<int> & rest, const vector<int> & lists) {
                if (rest.empty()) {
                    total += 1;
                    return;
                }
                set_vertices(rest, lists);
                total += none(discrete());
            });
            return total;
        }

        // Writes a statement whose count is the number of embeddings, and returns that
        // count, unless the pattern cannot fit by pigeonhole or no guess has any
        // embeddings at all. Then nothing is written, because a statement cannot say
        // zero that way, and no_statement_because says why.
        //
        // The root is the or of a node for each guess that has embeddings, or just that
        // node if there is only one. A guess's node is the and of a cube for the vertices
        // on hubs and an exact node for the rest, or just whichever of them is not empty.
        // Every embedding puts one set of vertices on hubs, so the guesses are disjoint.
        auto write_statement(std::ostream & out, const optional<ProofOptions> & proof_options) -> optional<loooong>
        {
            vector<int> all(pattern.size());
            for (int v = 0; v < pattern.size(); ++v)
                all[v] = v;

            if (proof_options) {
                if (! hubs.empty())
                    throw UnsupportedConfiguration{"Proof logging cannot yet be used with decomposition hubs"};
                return write_statement_with_proof(out, *proof_options, all);
            }

            if (all.empty()) {
                CountingStatementWriter writer{out, {}};
                writer.finish(writer.add_cube({}));
                statement_nodes = writer.number_of_nodes();
                return loooong{1};
            }

            set_vertices(all, {});
            auto q = build_quotient(discrete());
            number_of_components = q->components.size();
            if (more_blocks_than_targets(*q, target_with_label)) {
                no_statement_because = "the pattern has more vertices of some label than the target";
                return std::nullopt;
            }
            // not yet under proof, which would need the plain solver's proof of it
            if (rejected_before_search()) {
                no_statement_because = "the checks before search rule out every embedding";
                return std::nullopt;
            }

            // find the guesses with embeddings before writing anything
            struct Guess
            {
                vector<pair<int, int>> cube;
                vector<int> rest, lists;
            };
            vector<Guess> viable;
            for_each_guess(all, [&](const vector<pair<int, int>> & cube, const vector<int> & rest, const vector<int> & lists) {
                set_vertices(rest, lists);
                if (rest.empty() || might_embed(discrete()))
                    viable.push_back(Guess{cube, rest, lists});
            });
            if (aborted)
                return std::nullopt;
            if (viable.empty()) {
                no_statement_because = hubs.empty() ? "some component has no embeddings" : "no way of putting pattern vertices on hubs leaves the rest any embeddings";
                return std::nullopt;
            }

            CountingStatementWriter writer{out, vector<int>(pattern.size(), target.size())};
            statement = &writer;
            vector<int> guess_nodes;
            loooong total = 0;
            for (auto & g : viable) {
                optional<int> rest_node;
                if (! g.rest.empty()) {
                    set_vertices(g.rest, g.lists);
                    auto e = exact(discrete());
                    if (aborted)
                        return std::nullopt;
                    if (! e)
                        throw UnsupportedConfiguration{"a guess that might embed has no node, which cannot happen"};
                    rest_node = e->first;
                    total += e->second;
                }
                else
                    total += 1;

                if (g.cube.empty())
                    guess_nodes.push_back(*rest_node);
                else {
                    auto cube = writer.add_cube(g.cube);
                    guess_nodes.push_back(rest_node ? writer.add_and({cube, *rest_node}) : cube);
                }
            }
            statement = nullptr;
            writer.finish(guess_nodes.size() == 1 ? guess_nodes.front() : writer.add_or(guess_nodes));
            statement_nodes = writer.number_of_nodes();
            return total;
        }

        // As write_statement, but also writing a proof that the solver's encoding of the
        // instance is equienumerable with the statement's counting formula (phases 4 and 5
        // of decomposition-counting.md), or that it has no solutions when there is no
        // statement. Only without hubs, for now.
        auto write_statement_with_proof(std::ostream & out, const ProofOptions & proof_options, const vector<int> & all) -> optional<loooong>
        {
            innards::DecompositionProof proof(proof_options, pattern, target);
            decomposition_proof = &proof;

            set_vertices(all, {});
            auto discrete_partition = discrete();
            auto q = build_quotient(discrete_partition);
            number_of_components = q->components.size();
            ++guesses;
            if (more_blocks_than_targets(*q, target_with_label)) {
                no_statement_because = "the pattern has more vertices of some label than the target";
                proof.finish_pattern_bigger_than_target();
                return std::nullopt;
            }

            CountingStatementWriter writer{out, vector<int>(pattern.size(), target.size())};
            writer.set_observer([&](int i, const CountingNode & node) { proof.define_gate(i, node); });
            statement = &writer;

            // Each component of the pattern is searched under proof, which concludes that
            // its leaf holds. Isomorphic components cannot share these searches, because
            // each proof is about its own variables, but they still fill the enumeration
            // cache for the glued components later on.
            // for now, without the preprocessing whose derivations might not survive a
            // nested solve
            auto search_params = [&] {
                auto result = enumeration_params();
                result.clique_detection = false;
                result.no_supplementals = true;
                result.no_nds = true;
                return result;
            };
            innards::DecompositionProof::Structure structure;
            for (auto & c : q->components) {
                auto component = component_graph(*q, c);
                vector<int> component_vertices;
                for (auto b : c) {
                    component_vertices.push_back(vertices[b]);
                    component.names.push_back(pattern.vertex_name(vertices[b]));
                }

                ++solver_calls;
                auto searched = proof.search_component(as_input_graph(component), component_vertices, search_params(), writer);
                if (! searched.complete) {
                    aborted = true;
                    return std::nullopt;
                }
                if (! searched.leaf) {
                    no_statement_because = "some component has no embeddings";
                    proof.finish_no_solutions();
                    return std::nullopt;
                }

                Blocks members;
                for (auto v : component_vertices)
                    members.emplace_back(vector{v}, -1);
                leaf_nodes.emplace(members, *searched.leaf);
                leaf_cubes.emplace(*searched.leaf, searched.cubes);
                structure.components.push_back({component_vertices, *searched.leaf, searched.cubes, searched.embeddings});

                auto form = canonical_form(component);
                if (! embeddings_by_key.contains(form.key)) {
                    vector<vector<int>> by_position;
                    for (auto & e : searched.embeddings) {
                        vector<int> p(form.order.size());
                        for (unsigned i = 0; i < form.order.size(); ++i)
                            p[i] = e[form.order[i]];
                        by_position.push_back(std::move(p));
                    }
                    embeddings_by_key.emplace(form.key, std::move(by_position));
                }
            }

            auto root = exact(discrete_partition);
            statement = nullptr;
            if (aborted || ! root)
                return std::nullopt;

            // The root: every component's leaf holds, so the product does, and every
            // overlap has a glued leaf, every cube of which sends two pattern vertices to
            // one target vertex, so no overlap holds.
            auto & r = exact_structure.at(root->first);
            if (r.leaves.size() > 1)
                proof.derive(counting_gate_name(r.product));
            if (r.overlaps) {
                for (auto & c : r.children) {
                    auto & e = exact_structure.at(c);
                    auto glued = std::find_if(e.leaves.begin(), e.leaves.end(), [&](int l) { return glued_leaves.contains(l); });
                    if (glued == e.leaves.end())
                        throw UnsupportedConfiguration{"an overlap has no glued leaf, which cannot happen"};
                    proof.refute_glued_leaf(*glued, leaf_cubes.at(*glued));
                    proof.derive("~" + counting_gate_name(c));
                }
                proof.derive("~" + counting_gate_name(*r.overlaps));
            }

            // Then the root implies the encoding, which is what lets the encoding's
            // constraints be deleted, and that needs each overlap's partition.
            structure.root = root->first;
            structure.root_line = proof.derive(counting_gate_name(root->first));
            structure.product = r.product;
            structure.overlaps = r.overlaps;
            set<int> children(r.children.begin(), r.children.end());
            for (auto & [blocks, e] : exact_nodes)
                if (e && children.contains(e->first)) {
                    vector<vector<int>> partition;
                    for (auto & [members, _] : blocks)
                        partition.push_back(members);
                    structure.overlaps_by_partition.emplace(partition, e->first);
                }
            proof.finish(structure);

            writer.finish(root->first);
            statement_nodes = writer.number_of_nodes();
            return root->second;
        }

        int statement_nodes = 0;
        string no_statement_because;

        auto stats(std::list<string> & extra_stats) const -> void
        {
            if (rejected)
                extra_stats.emplace_back("decomposition_rejected_before_search = true");
            extra_stats.emplace_back("decomposition_components = " + to_string(number_of_components));
            extra_stats.emplace_back("decomposition_hubs = " + to_string(hubs.size()));
            extra_stats.emplace_back("decomposition_guesses = " + to_string(guesses));
            extra_stats.emplace_back("decomposition_solver_calls = " + to_string(solver_calls));
            extra_stats.emplace_back("decomposition_quotients = " + to_string(statement_nodes ? exact_nodes.size() : none_counts.size()));
            extra_stats.emplace_back("decomposition_coarsenings_considered = " + to_string(coarsenings_considered));
            if (statement_nodes)
                extra_stats.emplace_back("decomposition_statement_nodes = " + to_string(statement_nodes));
        }
    };

    auto count_by_components(const InputGraph & pattern, const InputGraph & target, const HomomorphismParams & params,
        const DecompositionOptions & options) -> HomomorphismResult
    {
        // Isolated vertices are not worth a solver call each, nor the inclusion-exclusion
        // that gluing them on would cost: every other pattern vertex with label L uses
        // exactly one target vertex labelled L, so the isolated ones labelled L can go on
        // any of the rest. (A looped vertex has a nonzero degree, so is not isolated.)
        vector<int> vertices;
        map<string, long> isolated_with_label, others_with_label, target_with_label;
        for (int v = 0; v < pattern.size(); ++v) {
            auto label = pattern.has_vertex_labels() ? string{pattern.vertex_label(v)} : string{};
            if (0 == pattern.degree(v))
                ++isolated_with_label[label];
            else {
                ++others_with_label[label];
                vertices.push_back(v);
            }
        }
        for (int t = 0; t < target.size(); ++t)
            ++target_with_label[pattern.has_vertex_labels() ? string{target.vertex_label(t)} : string{}];

        HomomorphismResult result;
        ComponentCounter counter(pattern, target, params, choose_hubs(pattern, target, vertices, options.hubs));
        counter.check_whole_instance = options.check_whole_instance;
        if (! counter.rejected_before_search()) {
            result.solution_count = counter.count(vertices);
            for (auto & [label, k] : isolated_with_label)
                result.solution_count *= falling_factorial(target_with_label[label] - others_with_label[label], k);
        }

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

    // The same count, but stating the argument as it goes. Isolated vertices are
    // components like any other here, because the statement language has no way of
    // saying "a falling factorial". Under proof there are no hubs, for now.
    auto write_statement_by_components(const InputGraph & pattern, const InputGraph & target, const HomomorphismParams & params,
        std::ostream & statement, const DecompositionOptions & options) -> HomomorphismResult
    {
        auto hubs = options.hubs;
        vector<int> vertices(pattern.size());
        for (int v = 0; v < pattern.size(); ++v)
            vertices[v] = v;

        if (params.proof_options && ! hubs)
            hubs = 0;

        HomomorphismResult result;
        ComponentCounter counter(pattern, target, params, choose_hubs(pattern, target, vertices, hubs));
        counter.check_whole_instance = options.check_whole_instance;
        auto count = counter.write_statement(statement, params.proof_options);
        result.nodes = counter.nodes;
        result.propagations = counter.propagations;
        result.complete = ! counter.aborted;
        result.solution_count = (count && ! counter.aborted) ? *count : loooong{0};
        if (! count && ! counter.aborted)
            result.extra_stats.emplace_back("decomposition_statement = none, because " + counter.no_statement_because);
        counter.stats(result.extra_stats);
        return result;
    }
}

auto gss::solve_sip_by_decomposition(const InputGraph & pattern, const InputGraph & target,
    const HomomorphismParams & params, const DecompositionOptions & options) -> HomomorphismResult
{
    auto statement = options.statement;

    // Without a statement, each piece would write its own proof over the top of the last
    // one, and none of them would be a proof about the whole instance. With one, the
    // proof is about the statement (see dev_docs/decomposition-counting.md).
    if (params.proof_options && ! statement)
        throw UnsupportedConfiguration{"Proof logging can only be used with decomposition when writing a statement"};
    if (params.proof_options && (pattern.has_vertex_labels() || pattern.has_edge_labels()))
        throw UnsupportedConfiguration{"Proof logging cannot yet be used on labelled graphs"};

    if (statement) {
        if (! can_count_by_components(pattern, params))
            throw UnsupportedConfiguration{"A decomposition statement can only be written when counting injective, non-induced embeddings of an undirected pattern, with no symmetry constraints"};
        return write_statement_by_components(pattern, target, params, *statement, options);
    }

    if (can_count_by_components(pattern, params))
        return count_by_components(pattern, target, params, options);

    if (options.hubs)
        throw UnsupportedConfiguration{"Decomposition can only use hubs when counting injective, non-induced embeddings of an undirected pattern, with no symmetry constraints"};

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
