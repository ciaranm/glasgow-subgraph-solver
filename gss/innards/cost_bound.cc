#include <gss/configuration.hh>
#include <gss/innards/cost_bound.hh>

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace gss;
using namespace gss::innards;

using std::map;
using std::max;
using std::min;
using std::pair;
using std::to_string;
using std::vector;

namespace
{
    // A cost no mapping can reach, and the arithmetic that keeps it where it is. The
    // constructor checks that real sums stay far below it, so anything at or above half
    // of it can only have come from an infinity.
    constexpr long long inf = std::numeric_limits<long long>::max() / 4;

    auto is_inf(long long a) -> bool
    {
        return a >= inf / 2;
    }

    auto add(long long a, long long b) -> long long
    {
        if (is_inf(a) || is_inf(b))
            return inf;
        return a + b;
    }

    auto floor_half(long long a) -> long long
    {
        return a >= 0 ? a / 2 : -((-a + 1) / 2);
    }

    // Rectangular minimum-cost assignment of rows to distinct columns (rows <= columns),
    // by the Hungarian algorithm with potentials, in exact integers. Returns inf if every
    // assignment uses an infinite entry. On success the potentials are dual feasible:
    // row_potential[i] + column_potential[j] <= cost[i][j] wherever that is finite, and
    // every column potential is at most 0. Note that a column left unmatched may have a
    // negative potential, so the assignment's cost plus a reduced cost is *not* a bound
    // on forcing an entry; see the filtering in propagate() for one that is.
    auto assign(const vector<vector<long long>> & cost, unsigned columns,
        vector<long long> & row_potential, vector<long long> & column_potential) -> long long
    {
        unsigned rows = cost.size();
        if (rows > columns)
            return inf;

        vector<long long> u(rows + 1, 0), v(columns + 1, 0);
        vector<unsigned> p(columns + 1, 0), way(columns + 1, 0);
        for (unsigned i = 1; i <= rows; ++i) {
            p[0] = i;
            unsigned j0 = 0;
            vector<long long> minv(columns + 1, inf);
            vector<char> used(columns + 1, false);
            do {
                used[j0] = true;
                unsigned i0 = p[j0], j1 = 0;
                long long delta = inf;
                for (unsigned j = 1; j <= columns; ++j)
                    if (! used[j]) {
                        long long c = cost[i0 - 1][j - 1];
                        if (! is_inf(c)) {
                            long long cur = c - u[i0] - v[j];
                            if (cur < minv[j]) {
                                minv[j] = cur;
                                way[j] = j0;
                            }
                        }
                        if (minv[j] < delta) {
                            delta = minv[j];
                            j1 = j;
                        }
                    }

                // Nothing finite left to reach: no assignment avoids an infinite entry.
                if (is_inf(delta))
                    return inf;

                for (unsigned j = 0; j <= columns; ++j)
                    if (used[j]) {
                        u[p[j]] += delta;
                        v[j] -= delta;
                    }
                    else if (! is_inf(minv[j]))
                        minv[j] -= delta;
                j0 = j1;
            } while (p[j0] != 0);

            do {
                unsigned j1 = way[j0];
                p[j0] = p[j1];
                j0 = j1;
            } while (j0);
        }

        long long result = 0;
        for (unsigned j = 1; j <= columns; ++j)
            if (p[j])
                result = add(result, cost[p[j] - 1][j - 1]);

        row_potential.assign(u.begin() + 1, u.end());
        column_potential.assign(v.begin() + 1, v.end());
        return result;
    }

    // One pair of unassigned original pattern vertices joined by at least one
    // edge-vertex, and the cost of each pair of their candidate images.
    struct PairTerm
    {
        unsigned row_a, row_b;
        vector<long long> cost; // candidates of a by candidates of b
        vector<long long> message_a, message_b;
    };
}

CostBound::CostBound(const CostData & data, unsigned pattern_size, unsigned target_size, int dual_sweeps, bool pruning) :
    _data(data),
    _pattern_size(pattern_size),
    _target_size(target_size),
    _dual_sweeps(dual_sweeps),
    _pruning(pruning)
{
    // Every sum the bound forms is at most the pattern size times the largest cost, with
    // the dual messages bounded by the same, so this keeps all of them far from inf.
    long long largest = 0;
    for (auto c : _data.target_costs) {
        if (c == std::numeric_limits<long long>::min())
            throw UnsupportedConfiguration{"Target costs are too large to sum safely"};
        largest = max(largest, std::llabs(c));
    }
    if (largest != 0 && (inf / 64) / largest < (long long)(pattern_size) + 1)
        throw UnsupportedConfiguration{"Target costs are too large to sum safely over a mapping of this pattern: the largest is " + to_string(largest)};

    _scores.resize(_data.pattern_original_size);
}

auto CostBound::infinity() -> long long
{
    return inf;
}

auto CostBound::is_original_pattern_vertex(unsigned p) const -> bool
{
    return int(p) < _data.pattern_original_size;
}

auto CostBound::pattern_original_size() const -> unsigned
{
    return _data.pattern_original_size;
}

auto CostBound::cost_of(const vector<int> & assigned) const -> long long
{
    long long result = 0;
    for (auto t : assigned)
        result += _data.target_costs.at(t);
    return result;
}

auto CostBound::score(unsigned p, unsigned t) const -> long long
{
    // Only original vertices are scored. Search can still branch on an edge-vertex, when
    // a pattern edge without a label has several parallel target edges to choose from,
    // and then there is nothing to say about which to try first.
    if (int(p) >= _data.pattern_original_size || _scores[p].empty() || int(t) >= _data.target_original_size)
        return 0;
    return _scores[p][t];
}

auto CostBound::add_extra_stats(std::list<std::string> & stats) const -> void
{
    stats.emplace_back("cost_bound_calls = " + to_string(_calls));
    stats.emplace_back("cost_bound_removals = " + to_string(_removals));
    stats.emplace_back("cost_bound_time = " + to_string(std::chrono::duration_cast<std::chrono::milliseconds>(_time).count()));
}

auto CostBound::propagate(const vector<int> & assigned, vector<HomomorphismDomain> & domains,
    long long upper_bound, bool & changed) -> bool
{
    if (! _pruning)
        return true;

    ++_calls;
    auto start = std::chrono::steady_clock::now();
    struct AddTime
    {
        std::chrono::nanoseconds & total;
        std::chrono::steady_clock::time_point start;
        ~AddTime()
        {
            total += std::chrono::steady_clock::now() - start;
        }
    } add_time{_time, start};

    return propagate_timed(assigned, domains, upper_bound, changed);
}

auto CostBound::propagate_timed(const vector<int> & assigned, vector<HomomorphismDomain> & domains,
    long long upper_bound, bool & changed) -> bool
{
    const int pattern_original = _data.pattern_original_size;
    const int target_original = _data.target_original_size;

    vector<HomomorphismDomain *> domain_of(_pattern_size, nullptr);
    for (auto & d : domains)
        if (! d.fixed)
            domain_of[d.v] = &d;

    // The cost of what is already decided.
    long long fixed_cost = 0;
    for (unsigned p = 0; p < _pattern_size; ++p)
        if (assigned[p] != -1)
            fixed_cost = add(fixed_cost, _data.target_costs[assigned[p]]);

    // Rows: the unassigned original pattern vertices, each with its candidate images,
    // and for each candidate the unary cost it carries so far.
    vector<unsigned> rows;
    vector<int> row_of(pattern_original, -1);
    vector<vector<unsigned>> candidates;
    vector<vector<int>> position; // per row, per original target vertex, index in candidates or -1
    vector<vector<long long>> unary;
    for (int p = 0; p < pattern_original; ++p) {
        if (assigned[p] != -1 || ! domain_of[p])
            continue;
        row_of[p] = rows.size();
        rows.push_back(p);
        candidates.emplace_back();
        position.emplace_back(target_original, -1);
        unary.emplace_back();
        domain_of[p]->values.for_each([&](unsigned t) {
            // Labels keep an original vertex off edge-vertices, but say so rather than trust it.
            if (int(t) >= target_original)
                return;
            position.back()[t] = candidates.back().size();
            candidates.back().push_back(t);
            unary.back().push_back(_data.target_costs[t]);
        });
    }

    // Edge-vertices: fold each into a constant, a unary cost, or a pairwise cost,
    // according to how many of its endpoints are assigned.
    map<pair<unsigned, unsigned>, unsigned> pair_index;
    vector<PairTerm> pairs;

    auto target_endpoints = [&](unsigned t) -> const EdgeVertexEndpoints & {
        return _data.target_edge_vertices[t - target_original];
    };

    for (unsigned k = 0; k < _data.pattern_edge_vertices.size(); ++k) {
        unsigned e = pattern_original + k;
        if (assigned[e] != -1 || ! domain_of[e])
            continue;

        auto [from, to] = _data.pattern_edge_vertices[k];
        int f_from = assigned[from], f_to = assigned[to];

        // The cheapest image left for an edge-vertex whose endpoints are both decided.
        if (f_from != -1 && f_to != -1) {
            long long best = inf;
            domain_of[e]->values.for_each([&](unsigned t) { best = min(best, _data.target_costs[t]); });
            fixed_cost = add(fixed_cost, best);
            continue;
        }

        // One endpoint decided (or a loop with its one endpoint not): a cost on the
        // other endpoint's candidates.
        if (f_from != -1 || f_to != -1 || from == to) {
            unsigned free_vertex = (f_from == -1) ? from : to;
            int decided = (f_from == -1) ? f_to : f_from;
            bool free_is_from = (free_vertex == unsigned(from));
            if (row_of[free_vertex] == -1)
                continue;
            unsigned r = row_of[free_vertex];

            vector<long long> best(candidates[r].size(), inf);
            domain_of[e]->values.for_each([&](unsigned t) {
                auto & te = target_endpoints(t);
                auto consider = [&](int free_image, int other_image) {
                    if (from == to) {
                        if (free_image != other_image)
                            return;
                    }
                    else if (other_image != decided)
                        return;
                    if (free_image >= target_original)
                        return;
                    int i = position[r][free_image];
                    if (i != -1)
                        best[i] = min(best[i], _data.target_costs[t]);
                };
                if (_data.directed)
                    free_is_from ? consider(te.from, te.to) : consider(te.to, te.from);
                else {
                    consider(te.from, te.to);
                    consider(te.to, te.from);
                }
            });

            for (unsigned i = 0; i < candidates[r].size(); ++i)
                unary[r][i] = add(unary[r][i], best[i]);
            continue;
        }

        // Neither endpoint decided: a cost on the pair of their images.
        unsigned a = min<unsigned>(from, to), b = max<unsigned>(from, to);
        bool a_is_from = (a == unsigned(from));
        if (row_of[a] == -1 || row_of[b] == -1)
            continue;
        auto [it, fresh] = pair_index.emplace(pair{a, b}, pairs.size());
        if (fresh) {
            PairTerm term;
            term.row_a = row_of[a];
            term.row_b = row_of[b];
            term.cost.assign(candidates[term.row_a].size() * candidates[term.row_b].size(), 0);
            term.message_a.assign(candidates[term.row_a].size(), 0);
            term.message_b.assign(candidates[term.row_b].size(), 0);
            pairs.push_back(std::move(term));
        }
        auto & term = pairs[it->second];
        unsigned width = candidates[term.row_b].size();

        vector<long long> this_edge(term.cost.size(), inf);
        domain_of[e]->values.for_each([&](unsigned t) {
            auto & te = target_endpoints(t);
            auto consider = [&](int image_a, int image_b) {
                if (image_a >= target_original || image_b >= target_original)
                    return;
                int i = position[term.row_a][image_a], j = position[term.row_b][image_b];
                if (i != -1 && j != -1)
                    this_edge[i * width + j] = min(this_edge[i * width + j], _data.target_costs[t]);
            };
            if (_data.directed)
                a_is_from ? consider(te.from, te.to) : consider(te.to, te.from);
            else {
                consider(te.from, te.to);
                consider(te.to, te.from);
            }
        });

        for (unsigned i = 0; i < term.cost.size(); ++i)
            term.cost[i] = add(term.cost[i], this_edge[i]);
    }

    // Only an infinite part can conclude anything yet: costs may be negative, so the
    // part decided so far reaching upper_bound says nothing about the whole.
    if (is_inf(fixed_cost))
        return false;

    // Dual ascent. theta is each row's unary cost plus the messages the pairs have sent
    // it; each pair's residual cost is its cost minus the two messages it sent.
    auto & theta = unary;
    for (int sweep = 0; sweep < _dual_sweeps && ! pairs.empty(); ++sweep) {
        for (auto & term : pairs) {
            auto & ta = theta[term.row_a];
            auto & tb = theta[term.row_b];
            unsigned na = ta.size(), nb = tb.size();

            for (unsigned i = 0; i < na; ++i)
                if (! is_inf(ta[i]))
                    ta[i] -= term.message_a[i];
            for (unsigned j = 0; j < nb; ++j)
                if (! is_inf(tb[j]))
                    tb[j] -= term.message_b[j];

            vector<long long> best_a(na, inf), best_b(nb, inf);
            for (unsigned i = 0; i < na; ++i)
                for (unsigned j = 0; j < nb; ++j) {
                    long long c = term.cost[i * nb + j];
                    if (is_inf(c))
                        continue;
                    best_a[i] = min(best_a[i], add(c, tb[j]));
                    best_b[j] = min(best_b[j], add(c, ta[i]));
                }

            // Rounded down, so that the messages are integers and the residuals exact.
            for (unsigned i = 0; i < na; ++i) {
                if (is_inf(ta[i]))
                    term.message_a[i] = 0;
                else if (is_inf(best_a[i])) {
                    term.message_a[i] = 0;
                    ta[i] = inf;
                }
                else
                    term.message_a[i] = floor_half(best_a[i] - ta[i]);
            }
            for (unsigned j = 0; j < nb; ++j) {
                if (is_inf(tb[j]))
                    term.message_b[j] = 0;
                else if (is_inf(best_b[j])) {
                    term.message_b[j] = 0;
                    tb[j] = inf;
                }
                else
                    term.message_b[j] = floor_half(best_b[j] - tb[j]);
            }

            for (unsigned i = 0; i < na; ++i)
                if (! is_inf(ta[i]))
                    ta[i] += term.message_a[i];
            for (unsigned j = 0; j < nb; ++j)
                if (! is_inf(tb[j]))
                    tb[j] += term.message_b[j];
        }
    }

    // What is left in each pair after the messages, at its cheapest.
    long long bound = fixed_cost;
    for (auto & term : pairs) {
        auto & ta = theta[term.row_a];
        auto & tb = theta[term.row_b];
        unsigned nb = tb.size();
        long long residual = inf;
        for (unsigned i = 0; i < ta.size(); ++i) {
            if (is_inf(ta[i]))
                continue;
            for (unsigned j = 0; j < nb; ++j) {
                long long c = term.cost[i * nb + j];
                if (is_inf(c) || is_inf(tb[j]))
                    continue;
                residual = min(residual, c - term.message_a[i] - term.message_b[j]);
            }
        }
        bound = add(bound, residual);
    }

    if (is_inf(bound))
        return false;

    // Combine the reparametrised unary costs by a minimum-cost assignment over the
    // original target vertices that are still candidates somewhere.
    vector<int> column_of(target_original, -1);
    vector<unsigned> columns;
    for (auto & c : candidates)
        for (auto t : c)
            if (column_of[t] == -1) {
                column_of[t] = columns.size();
                columns.push_back(t);
            }

    vector<vector<long long>> matrix(rows.size(), vector<long long>(columns.size(), inf));
    for (unsigned r = 0; r < rows.size(); ++r)
        for (unsigned i = 0; i < candidates[r].size(); ++i)
            matrix[r][column_of[candidates[r][i]]] = theta[r][i];

    vector<long long> row_potential, column_potential;
    long long assignment = rows.empty() ? 0 : assign(matrix, columns.size(), row_potential, column_potential);
    long long before_assignment = bound;
    bound = add(bound, assignment);
    if (is_inf(bound) || bound >= upper_bound)
        return false;

    // What forcing row r onto column j costs at least. Any assignment s using it costs
    // sum_k (u_k + v_s(k) + reduced_k) >= sum_k u_k + reduced_rj + sum_k v_s(k), and the
    // column potentials are at most 0, so the last sum is at least the most negative
    // total of k potentials that includes v_j. Summing only the matched columns'
    // potentials instead, as the assignment's own cost does, would over-prune.
    long long row_total = 0;
    for (auto u : row_potential)
        row_total += u;
    vector<long long> sorted_columns(column_potential);
    std::sort(sorted_columns.begin(), sorted_columns.end());
    unsigned k = rows.size();
    long long most_negative_k = 0, most_negative_k_minus_one = 0;
    for (unsigned i = 0; i < k && i < sorted_columns.size(); ++i) {
        most_negative_k += sorted_columns[i];
        if (i + 1 < k)
            most_negative_k_minus_one += sorted_columns[i];
    }
    long long kth_most_negative = (k >= 1 && k <= sorted_columns.size()) ? sorted_columns[k - 1] : 0;
    auto forced_bound = [&](unsigned r, unsigned t, long long th) -> long long {
        long long v = column_potential[column_of[t]];
        long long columns_part = (v <= kth_most_negative) ? most_negative_k : most_negative_k_minus_one + v;
        long long reduced = th - row_potential[r] - v;
        return add(before_assignment, row_total + columns_part + reduced);
    };

    // Remove every value that cannot be in a mapping cheaper than upper_bound, and
    // every value no edge-vertex supports.
    for (unsigned r = 0; r < rows.size(); ++r) {
        auto & d = *domain_of[rows[r]];
        auto & scores = _scores[rows[r]];
        scores.assign(target_original, inf);
        for (unsigned i = 0; i < candidates[r].size(); ++i) {
            unsigned t = candidates[r][i];
            long long th = theta[r][i];
            bool remove = is_inf(th) || forced_bound(r, t, th) >= upper_bound;
            if (remove) {
                d.values.reset(t);
                --d.count;
                ++_removals;
                changed = true;
            }
            else
                scores[t] = th;
        }
        if (0 == d.count)
            return false;
    }

    return true;
}
