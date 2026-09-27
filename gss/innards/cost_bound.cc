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

    // Minimum-cost assignment of rows to distinct columns (rows <= columns), by the
    // Hungarian algorithm on the square problem padded with zero-cost rows, in exact
    // integers. Returns false if every assignment uses an infinite entry, and then
    // hall_rows and hall_columns are a Hall violator: rows whose finite entries all lie
    // in fewer columns than there are rows. Otherwise it gives an exact dual of the
    // rectangular problem, row potentials alpha of any sign and column potentials
    // beta >= 0, with alpha[i] - beta[j] <= cost[i][j] wherever that is finite: the dual
    // value, sum alpha - sum beta, is then a lower bound on any assignment, and
    // cost[i][j] - alpha[i] + beta[j] >= 0 a lower bound on how much forcing i onto j
    // adds to it.
    //
    // (The unpadded algorithm gives potentials whose unmatched columns may be negative,
    // and then the dual value falls short of the assignment's cost, and cost plus a
    // reduced cost is not a valid forced bound. Padding is what makes it exact.)
    auto assign(const vector<vector<long long>> & cost, unsigned columns,
        vector<long long> & alpha, vector<long long> & beta,
        vector<unsigned> & hall_rows, vector<unsigned> & hall_columns) -> bool
    {
        unsigned rows = cost.size(), n = columns;
        if (rows > n) {
            for (unsigned i = 0; i < rows; ++i)
                hall_rows.push_back(i);
            for (unsigned j = 0; j < n; ++j)
                hall_columns.push_back(j);
            return false;
        }

        auto entry = [&](unsigned i, unsigned j) -> long long {
            return i < rows ? cost[i][j] : 0;
        };

        vector<long long> u(n + 1, 0), v(n + 1, 0);
        vector<unsigned> p(n + 1, 0), way(n + 1, 0);
        for (unsigned i = 1; i <= n; ++i) {
            p[0] = i;
            unsigned j0 = 0;
            vector<long long> minv(n + 1, inf);
            vector<char> used(n + 1, false);
            do {
                used[j0] = true;
                unsigned i0 = p[j0], j1 = 0;
                long long delta = inf;
                for (unsigned j = 1; j <= n; ++j)
                    if (! used[j]) {
                        long long c = entry(i0 - 1, j - 1);
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

                // Nothing finite left to reach. The real rows come first and a padding
                // row reaches everything, so this is a real row, and the rows of the
                // alternating tree only reach the columns in it, each already matched
                // to one of them.
                if (is_inf(delta)) {
                    for (unsigned j = 0; j <= n; ++j)
                        if (used[j]) {
                            hall_rows.push_back(p[j] - 1);
                            if (j != 0)
                                hall_columns.push_back(j - 1);
                        }
                    return false;
                }

                for (unsigned j = 0; j <= n; ++j)
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

        // Shift to the rectangular problem's dual: alpha = u + s and beta = s - v, which
        // keeps alpha - beta = u + v. beta >= 0 needs v <= s: with padding rows, every
        // padding potential w has w + v <= 0, so s = -(the largest) will do; without,
        // s is the largest v.
        long long s;
        if (rows < n) {
            long long largest = std::numeric_limits<long long>::min();
            for (unsigned i = rows + 1; i <= n; ++i)
                largest = max(largest, u[i]);
            s = -largest;
        }
        else {
            s = std::numeric_limits<long long>::min();
            for (unsigned j = 1; j <= n; ++j)
                s = max(s, v[j]);
        }

        alpha.assign(rows, 0);
        beta.assign(n, 0);
        for (unsigned i = 0; i < rows; ++i)
            alpha[i] = u[i + 1] + s;
        for (unsigned j = 0; j < n; ++j)
            beta[j] = s - v[j + 1];
        return true;
    }

    // A pair of adjacent original pattern vertices, and the cost of each pair of their
    // candidate images: the sum over the edge-vertices between them of the cheapest image
    // consistent with those, or inf if one has none.
    struct PairTerm
    {
        unsigned a, b; // pattern vertices, a < b, also their rows
        vector<long long> cost; // candidates of a by candidates of b
        vector<long long> message_a, message_b;
        long long residual = 0;
        bool folded = false; // one side has a single candidate, so it is not a dual variable
    };
}

CostBound::CostBound(const CostData & data, unsigned pattern_size, unsigned target_size, int dual_sweeps) :
    _data(data),
    _pattern_size(pattern_size),
    _target_size(target_size),
    _dual_sweeps(dual_sweeps)
{
    // A reparametrised unary cost is at most the pattern size times the largest cost (it
    // gathers messages from pairs, each at most a pair's cost), and the assignment step's
    // potentials at most the number of columns times that, so this keeps every sum the
    // bound forms far from inf.
    long long largest = 0;
    for (auto c : _data.target_costs) {
        if (c == std::numeric_limits<long long>::min())
            throw UnsupportedConfiguration{"Target costs are too large to sum safely"};
        largest = max(largest, std::llabs(c));
    }
    if (largest != 0 && ((inf / 64) / largest) / ((long long)(pattern_size) + 1) < (long long)(target_size) + 1)
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

auto CostBound::want_certificates() -> void
{
    _want_certificates = true;
}

auto CostBound::certificate() const -> const CostBoundCertificate &
{
    return _certificate;
}

auto CostBound::propagate_timed(const vector<int> & assigned, vector<HomomorphismDomain> & domains,
    long long upper_bound, bool & changed) -> bool
{
    const int pattern_original = _data.pattern_original_size;
    const int target_original = _data.target_original_size;

    _certificate = CostBoundCertificate{};

    vector<HomomorphismDomain *> domain_of(_pattern_size, nullptr);
    for (auto & d : domains)
        if (! d.fixed)
            domain_of[d.v] = &d;

    // The values a pattern vertex might still take: its value, if it has one, else its
    // domain.
    auto values_of = [&](unsigned p, auto && f) {
        if (assigned[p] != -1)
            f(unsigned(assigned[p]));
        else if (domain_of[p])
            domain_of[p]->values.for_each(f);
    };

    // Rows: every original pattern vertex, with its candidate images and the unary cost
    // each carries so far.
    vector<vector<unsigned>> candidates(pattern_original);
    vector<vector<int>> position(pattern_original, vector<int>(target_original, -1));
    vector<vector<long long>> theta(pattern_original);
    for (int p = 0; p < pattern_original; ++p) {
        values_of(p, [&](unsigned t) {
            // Labels keep an original vertex off edge-vertices, but say so rather than trust it.
            if (int(t) >= target_original)
                return;
            position[p][t] = candidates[p].size();
            candidates[p].push_back(t);
            theta[p].push_back(_data.target_costs[t]);
        });
        if (candidates[p].empty())
            return false;
    }

    auto target_endpoints = [&](unsigned t) -> const EdgeVertexEndpoints & {
        return _data.target_edge_vertices[t - target_original];
    };

    // Edge-vertices: a loop is a unary cost on its vertex, and anything else a pairwise
    // cost on its endpoints.
    //
    // For proofs, these costs must be exactly the model's: P(x, y) must be the cost of
    // z(a, b, x, y) for every pair of live candidates, so an edge-vertex's domain must
    // still hold every label-compatible target edge-vertex whose endpoints are both still
    // candidates. That is true because nothing prunes an edge-vertex's domain otherwise
    // under proof logging (injectivity on edge-vertices, the degree filters and the
    // supplemental graphs are all off then). Something that did would leave a positive
    // coefficient on a z that propagation cannot falsify, and the proof would fail to
    // verify rather than be wrong.
    map<pair<unsigned, unsigned>, unsigned> pair_index;
    vector<PairTerm> pairs;
    for (unsigned k = 0; k < _data.pattern_edge_vertices.size(); ++k) {
        unsigned e = pattern_original + k;
        auto [from, to] = _data.pattern_edge_vertices[k];

        if (from == to) {
            vector<long long> best(candidates[from].size(), inf);
            values_of(e, [&](unsigned t) {
                auto & te = target_endpoints(t);
                if (te.from == te.to && te.from < target_original && position[from][te.from] != -1)
                    best[position[from][te.from]] = min(best[position[from][te.from]], _data.target_costs[t]);
            });
            for (unsigned i = 0; i < candidates[from].size(); ++i)
                theta[from][i] = add(theta[from][i], best[i]);
            continue;
        }

        unsigned a = min<unsigned>(from, to), b = max<unsigned>(from, to);
        bool a_is_from = (a == unsigned(from));
        auto [it, fresh] = pair_index.emplace(pair{a, b}, pairs.size());
        if (fresh) {
            PairTerm term;
            term.a = a;
            term.b = b;
            term.cost.assign(candidates[a].size() * candidates[b].size(), 0);
            term.message_a.assign(candidates[a].size(), 0);
            term.message_b.assign(candidates[b].size(), 0);
            pairs.push_back(std::move(term));
        }
        auto & term = pairs[it->second];
        unsigned width = candidates[b].size();

        vector<long long> this_edge(term.cost.size(), inf);
        values_of(e, [&](unsigned t) {
            auto & te = target_endpoints(t);
            auto consider = [&](int image_a, int image_b) {
                if (image_a >= target_original || image_b >= target_original)
                    return;
                int i = position[a][image_a], j = position[b][image_b];
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

    // A pair with one candidate on a side folds exactly: all of its cost goes to the other
    // side (or to the residual, if both sides have one).
    for (auto & term : pairs) {
        unsigned na = candidates[term.a].size(), nb = candidates[term.b].size();
        if (na == 1 && nb == 1) {
            term.folded = true;
            term.residual = term.cost[0];
        }
        else if (na == 1) {
            term.folded = true;
            for (unsigned j = 0; j < nb; ++j) {
                term.message_b[j] = is_inf(term.cost[j]) ? 0 : term.cost[j];
                theta[term.b][j] = is_inf(term.cost[j]) ? inf : add(theta[term.b][j], term.cost[j]);
            }
        }
        else if (nb == 1) {
            term.folded = true;
            for (unsigned i = 0; i < na; ++i) {
                term.message_a[i] = is_inf(term.cost[i]) ? 0 : term.cost[i];
                theta[term.a][i] = is_inf(term.cost[i]) ? inf : add(theta[term.a][i], term.cost[i]);
            }
        }
    }

    // Dual ascent over the rest. theta is each row's unary cost plus the messages the
    // pairs have sent it; each pair's residual cost is its cost minus the two messages
    // it sent.
    for (int sweep = 0; sweep < _dual_sweeps; ++sweep) {
        for (auto & term : pairs) {
            if (term.folded)
                continue;
            auto & ta = theta[term.a];
            auto & tb = theta[term.b];
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

    // What is left in each pair after the messages, at its cheapest over the values still
    // alive on both sides.
    long long bound = 0;
    for (auto & term : pairs) {
        if (! term.folded) {
            auto & ta = theta[term.a];
            auto & tb = theta[term.b];
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
            term.residual = residual;
        }
        bound = add(bound, term.residual);
    }

    // Nothing supports any pair of images for some pair: propagation over the linking
    // equalities finds that for itself.
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

    vector<vector<long long>> matrix(pattern_original, vector<long long>(columns.size(), inf));
    for (int p = 0; p < pattern_original; ++p)
        for (unsigned i = 0; i < candidates[p].size(); ++i)
            matrix[p][column_of[candidates[p][i]]] = theta[p][i];

    vector<long long> alpha, beta;
    vector<unsigned> hall_rows, hall_columns;
    if (! assign(matrix, columns.size(), alpha, beta, hall_rows, hall_columns)) {
        if (_want_certificates) {
            _certificate.kind = CostBoundCertificate::Kind::HallViolator;
            _certificate.hall_rows = hall_rows;
            for (auto j : hall_columns)
                _certificate.hall_columns.push_back(columns[j]);
        }
        return false;
    }

    for (auto a : alpha)
        bound += a;
    for (auto b : beta)
        bound -= b;

    auto record_certificate = [&]() {
        if (! _want_certificates)
            return;
        _certificate.kind = CostBoundCertificate::Kind::Bound;
        _certificate.bound = bound;

        // A pair's residual rides on the exactly-one of its first vertex, and on the
        // linking equalities of that vertex's candidates: sum over x of
        // (sum_y z(a, b, x, y) - x(a, x)) plus sum_x x(a, x) = 1 is sum z = 1.
        vector<long long> exactly_one(alpha);
        for (auto & term : pairs)
            if (term.residual != 0)
                exactly_one[term.a] += term.residual;
        for (int p = 0; p < pattern_original; ++p)
            if (exactly_one[p] != 0)
                _certificate.exactly_one.emplace_back(p, exactly_one[p]);
        for (unsigned j = 0; j < columns.size(); ++j)
            if (beta[j] != 0)
                _certificate.at_most_one.emplace_back(columns[j], beta[j]);

        for (auto & term : pairs) {
            for (unsigned i = 0; i < candidates[term.a].size(); ++i) {
                long long m = (is_inf(theta[term.a][i]) ? 0 : term.message_a[i]) + term.residual;
                if (m != 0)
                    _certificate.links.push_back({term.a, term.b, true, candidates[term.a][i], m});
            }
            for (unsigned j = 0; j < candidates[term.b].size(); ++j) {
                long long m = is_inf(theta[term.b][j]) ? 0 : term.message_b[j];
                if (m != 0)
                    _certificate.links.push_back({term.a, term.b, false, candidates[term.b][j], m});
            }
        }
    };

    if (bound >= upper_bound) {
        record_certificate();
        return false;
    }

    // Remove every value that cannot be in a mapping cheaper than upper_bound, and
    // every value nothing supports.
    bool any_by_bound = false;
    for (int p = 0; p < pattern_original; ++p) {
        auto & scores = _scores[p];
        scores.assign(target_original, inf);
        if (! domain_of[p] || assigned[p] != -1)
            continue;
        auto & d = *domain_of[p];
        for (unsigned i = 0; i < candidates[p].size(); ++i) {
            unsigned t = candidates[p][i];
            long long th = theta[p][i];
            bool remove = is_inf(th);
            if (! remove && add(bound, th - alpha[p] + beta[column_of[t]]) >= upper_bound) {
                remove = true;
                any_by_bound = true;
            }
            if (remove) {
                d.values.reset(t);
                --d.count;
                ++_removals;
                changed = true;
                if (_want_certificates)
                    _certificate.removed.emplace_back(p, t);
            }
            else
                scores[t] = th;
        }
        if (0 == d.count) {
            if (any_by_bound)
                record_certificate();
            return false;
        }
    }

    if (any_by_bound)
        record_certificate();

    return true;
}
