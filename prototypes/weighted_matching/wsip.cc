// Prototype, not part of the solver: weighted subgraph isomorphism on the graph3 scene-graph
// data. See README.md in this directory. Minimise the sum of -log(weight) over pattern edges
// (and vertices), injective, label-preserving, where every pattern edge (u, v, label,
// directedness) must map to a target edge with the same label and directedness. Parallel edges
// with different labels between one pair are allowed, which is why this does not use InputGraph.
//
//   wsip <mode> pattern.csv target.csv [time limit in seconds] [dual iterations]
//
// prints: status, objective, cost of the name-preserving mapping, root bound, nodes, seconds,
// and how many pattern vertices were mapped to the target vertex of the same name.
//
// Modes that search, differing only in the bound:
//   none : partial cost only (what the collaborators' modified solver does)
//   gl   : Gilmore-Lawler style: per-vertex cheapest completion, pair costs split in half,
//          no injectivity across vertices; with cost-based domain filtering
//   lap  : as gl, but the per-vertex costs are combined by a min-cost assignment
//          (Hungarian), which adds injectivity; reduced-cost filtering
//   dual : as lap, but the half/half split of each pair's cost is replaced by a dual-ascent
//          reparametrisation, recomputed at every node
//
// Modes that write a pseudo-Boolean model to stdout instead, costs scaled by 10^4 and rounded
// (variable names need renum.py before RoundingSat will read them):
//   opb-weak   : y >= x + x' - 1 linearisation of each pair cost
//   opb-strong : local-polytope encoding, sum_s y(t, s) = x(u, t) and symmetrically

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

constexpr double INF = numeric_limits<double>::infinity();

struct Graph
{
    vector<string> names, labels;
    vector<double> vweight;
    map<string, int> index;
    // edges: from, to, label, directed, weight
    struct E
    {
        int f, t;
        string label;
        bool directed;
        double w;
    };
    vector<E> edges;

    int vertex(const string & n)
    {
        auto it = index.find(n);
        if (it != index.end()) return it->second;
        int i = names.size();
        index.emplace(n, i);
        names.push_back(n);
        labels.push_back("");
        vweight.push_back(1.0);
        return i;
    }
};

auto split(const string & s, char c) -> vector<string>
{
    vector<string> r;
    string cur;
    for (char x : s) {
        if (x == c) {
            r.push_back(cur);
            cur.clear();
        }
        else
            cur += x;
    }
    r.push_back(cur);
    return r;
}

auto read_graph(const string & fn) -> Graph
{
    Graph g;
    ifstream in{fn};
    if (! in) {
        cerr << "can't open " << fn << endl;
        exit(1);
    }
    string line;
    while (getline(in, line)) {
        if (! line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        auto p = split(line, ',');
        if (p.size() >= 4 && p[1].empty()) {
            // vertex: name,,label,weight[,D/B]
            int v = g.vertex(p[0]);
            g.labels[v] = p[2];
            g.vweight[v] = stod(p[3]);
        }
        else if (p.size() == 3 && p[0].find('>') != string::npos) {
            auto q = split(p[0], '>');
            g.edges.push_back({g.vertex(q[0]), g.vertex(q[1]), p[1], true, stod(p[2])});
        }
        else if (p.size() == 4) {
            g.edges.push_back({g.vertex(p[0]), g.vertex(p[1]), p[2], false, stod(p[3])});
        }
        else {
            cerr << "bad line " << line << endl;
            exit(1);
        }
    }
    return g;
}

// A pattern pair {u < v} and the target cost matrix over (x, y) = (f(u), f(v)).
struct Pair
{
    int u, v;
    vector<double> c; // nt * nt, INF where infeasible (includes x == y)
};

struct Problem
{
    int np, nt;
    vector<double> unary; // np * nt
    vector<Pair> pairs;
    vector<vector<pair<int, int>>> incident; // per pattern vertex: (pair index, 0 if u is .u else 1)
    vector<vector<int>> dom; // initial domains
};

struct Stats
{
    long long nodes = 0, sols = 0;
};

struct Solver
{
    const Problem & P;
    string bound;
    double ub = INF;
    vector<int> best;
    Stats st;
    chrono::steady_clock::time_point start;
    double limit;
    bool timed_out = false;
    int dual_iters = 5;

    Solver(const Problem & p, string b, double lim) : P(p), bound(b), limit(lim)
    {
    }

    auto C(const Pair & pr, int x, int y) const -> double
    {
        return pr.c[x * P.nt + y];
    }

    // Given a partial assignment (-1 = unassigned) and domains, compute for each unassigned u and
    // x in D(u) a lower bound contribution h(u, x); also returns g, the cost of the assigned part.
    // Domains may be pruned for infeasibility (h = INF).
    auto compute_h(const vector<int> & f, vector<vector<int>> & D, vector<vector<double>> & h, double & g) -> void
    {
        int np = P.np, nt = P.nt;
        g = 0;
        for (int u = 0; u < np; ++u)
            if (f[u] != -1) g += P.unary[u * nt + f[u]];
        for (auto & pr : P.pairs)
            if (f[pr.u] != -1 && f[pr.v] != -1) g += C(pr, f[pr.u], f[pr.v]);

        h.assign(np, vector<double>(nt, INF));
        for (int u = 0; u < np; ++u) {
            if (f[u] != -1) continue;
            for (int x : D[u])
                h[u][x] = P.unary[u * nt + x];
        }

        // Reparametrised pair costs for pairs with both endpoints unassigned; pairs with one
        // endpoint assigned fold straight into the unary term.
        for (size_t pi = 0; pi < P.pairs.size(); ++pi) {
            auto & pr = P.pairs[pi];
            bool au = f[pr.u] != -1, av = f[pr.v] != -1;
            if (au && av) continue;
            if (au) {
                for (int y : D[pr.v])
                    h[pr.v][y] += C(pr, f[pr.u], y);
            }
            else if (av) {
                for (int x : D[pr.u])
                    h[pr.u][x] += C(pr, x, f[pr.v]);
            }
        }

        if (bound == "none") return;

        if (bound == "gl" || bound == "lap") {
            for (auto & pr : P.pairs) {
                if (f[pr.u] != -1 || f[pr.v] != -1) continue;
                for (int x : D[pr.u]) {
                    double m = INF;
                    for (int y : D[pr.v])
                        m = min(m, C(pr, x, y));
                    h[pr.u][x] += m / 2;
                }
                for (int y : D[pr.v]) {
                    double m = INF;
                    for (int x : D[pr.u])
                        m = min(m, C(pr, x, y));
                    h[pr.v][y] += m / 2;
                }
            }
        }
        else if (bound == "dual") {
            // Local-polytope dual (min-sum), block coordinate ascent in the style of MPLP:
            // theta_u(x) = h(u, x) plus messages; each free pair keeps a residual matrix.
            vector<size_t> free_pairs;
            for (size_t pi = 0; pi < P.pairs.size(); ++pi)
                if (f[P.pairs[pi].u] == -1 && f[P.pairs[pi].v] == -1) free_pairs.push_back(pi);
            // messages m_pu[k][x] (into u from pair k), m_pv[k][y]
            vector<vector<double>> mu(free_pairs.size(), vector<double>(nt, 0.0)), mv(free_pairs.size(), vector<double>(nt, 0.0));
            vector<vector<double>> th = h; // unary with all messages included
            for (int it = 0; it < dual_iters; ++it) {
                for (size_t k = 0; k < free_pairs.size(); ++k) {
                    auto & pr = P.pairs[free_pairs[k]];
                    // remove this pair's messages
                    for (int x : D[pr.u])
                        if (! isinf(th[pr.u][x])) th[pr.u][x] -= mu[k][x];
                    for (int y : D[pr.v])
                        if (! isinf(th[pr.v][y])) th[pr.v][y] -= mv[k][y];
                    // MPLP update: new mu(x) = 1/2 min_y [c(x,y) + th_v(y)] - 1/2 th_u(x), etc
                    vector<double> a(nt, INF), b(nt, INF);
                    for (int x : D[pr.u])
                        for (int y : D[pr.v]) {
                            double c = C(pr, x, y);
                            a[x] = min(a[x], c + th[pr.v][y]);
                            b[y] = min(b[y], c + th[pr.u][x]);
                        }
                    for (int x : D[pr.u]) {
                        mu[k][x] = isinf(th[pr.u][x]) ? 0.0 : (isinf(a[x]) ? INF : 0.5 * a[x] - 0.5 * th[pr.u][x]);
                        th[pr.u][x] += mu[k][x];
                    }
                    for (int y : D[pr.v]) {
                        mv[k][y] = isinf(th[pr.v][y]) ? 0.0 : (isinf(b[y]) ? INF : 0.5 * b[y] - 0.5 * th[pr.v][y]);
                        th[pr.v][y] += mv[k][y];
                    }
                }
            }
            // After MPLP, the bound is sum_u min_x th_u(x) + sum_k min_{x,y} [c - mu - mv]. Each
            // pair residual min is added in, split evenly between endpoints so h stays unary.
            for (size_t k = 0; k < free_pairs.size(); ++k) {
                auto & pr = P.pairs[free_pairs[k]];
                double r = INF;
                for (int x : D[pr.u])
                    for (int y : D[pr.v]) {
                        double c = C(pr, x, y);
                        if (isinf(c) || isinf(mu[k][x]) || isinf(mv[k][y]) || isinf(th[pr.u][x]) || isinf(th[pr.v][y])) continue;
                        r = min(r, c - mu[k][x] - mv[k][y]);
                    }
                if (isinf(r)) r = 0; // domain empty handled elsewhere
                for (int x : D[pr.u])
                    th[pr.u][x] += r / 2;
                for (int y : D[pr.v])
                    th[pr.v][y] += r / 2;
            }
            h = th;
        }
    }

    // Rectangular Hungarian: rows = unassigned pattern vertices, cols = target vertices.
    // Returns cost, and dual potentials for reduced costs.
    auto lap(const vector<int> & rows, const vector<vector<double>> & h, vector<double> & ru, vector<double> & rv) -> double
    {
        int n = rows.size(), m = P.nt;
        // standard O(n^2 m) with potentials; INF handled by big M
        const double BIG = 1e9;
        vector<double> u(n + 1, 0), v(m + 1, 0);
        vector<int> p(m + 1, 0), way(m + 1, 0);
        for (int i = 1; i <= n; ++i) {
            p[0] = i;
            int j0 = 0;
            vector<double> minv(m + 1, INF);
            vector<char> used(m + 1, false);
            do {
                used[j0] = true;
                int i0 = p[j0], j1 = 0;
                double delta = INF;
                for (int j = 1; j <= m; ++j)
                    if (! used[j]) {
                        double c = h[rows[i0 - 1]][j - 1];
                        if (isinf(c)) c = BIG;
                        double cur = c - u[i0] - v[j];
                        if (cur < minv[j]) minv[j] = cur, way[j] = j0;
                        if (minv[j] < delta) delta = minv[j], j1 = j;
                    }
                for (int j = 0; j <= m; ++j)
                    if (used[j])
                        u[p[j]] += delta, v[j] -= delta;
                    else
                        minv[j] -= delta;
                j0 = j1;
            } while (p[j0] != 0);
            do {
                int j1 = way[j0];
                p[j0] = p[j1];
                j0 = j1;
            } while (j0);
        }
        double cost = 0;
        for (int j = 1; j <= m; ++j)
            if (p[j]) {
                double c = h[rows[p[j] - 1]][j - 1];
                if (isinf(c)) return INF;
                cost += c;
            }
        ru.assign(n, 0);
        rv.assign(m, 0);
        for (int i = 0; i < n; ++i)
            ru[i] = u[i + 1];
        for (int j = 0; j < m; ++j)
            rv[j] = v[j + 1];
        return cost;
    }

    // Propagate: prune values by bound; returns false on wipeout / bound failure.
    auto propagate(const vector<int> & f, vector<vector<int>> & D, vector<vector<double>> & h, double & lb) -> bool
    {
        int np = P.np;
        while (true) {
            // injectivity: remove assigned values from other domains
            for (int u = 0; u < np; ++u)
                if (f[u] == -1)
                    for (int w = 0; w < np; ++w)
                        if (f[w] != -1) erase(D[u], f[w]);

            double g;
            compute_h(f, D, h, g);
            bool changed = false;
            vector<int> rows;
            for (int u = 0; u < np; ++u)
                if (f[u] == -1) {
                    vector<int> nd;
                    for (int x : D[u])
                        if (! isinf(h[u][x])) nd.push_back(x);
                    if (nd.size() != D[u].size()) changed = true;
                    D[u] = nd;
                    if (D[u].empty()) return false;
                    rows.push_back(u);
                }
            if (changed) continue;

            if (bound == "none") {
                lb = g;
                return lb < ub;
            }

            if (bound == "gl") {
                vector<double> mins(np, 0);
                double s = 0;
                for (int u : rows) {
                    double m = INF;
                    for (int x : D[u])
                        m = min(m, h[u][x]);
                    mins[u] = m;
                    s += m;
                }
                lb = g + s;
                if (lb >= ub) return false;
                for (int u : rows) {
                    vector<int> nd;
                    for (int x : D[u])
                        if (lb - mins[u] + h[u][x] < ub) nd.push_back(x);
                    if (nd.size() != D[u].size()) changed = true;
                    D[u] = nd;
                    if (D[u].empty()) return false;
                }
                if (! changed) return true;
                continue;
            }

            // lap / dual
            vector<double> ru, rv;
            double c = rows.empty() ? 0 : lap(rows, h, ru, rv);
            lb = g + c;
            if (lb >= ub) return false;
            for (size_t i = 0; i < rows.size(); ++i) {
                int u = rows[i];
                vector<int> nd;
                for (int x : D[u])
                    if (lb + (h[u][x] - ru[i] - rv[x]) < ub - 1e-9) nd.push_back(x);
                if (nd.size() != D[u].size()) changed = true;
                D[u] = nd;
                if (D[u].empty()) return false;
            }
            if (! changed) return true;
        }
    }

    static auto erase(vector<int> & d, int x) -> void
    {
        auto it = find(d.begin(), d.end(), x);
        if (it != d.end()) d.erase(it);
    }

    auto search(vector<int> & f, vector<vector<int>> D) -> void
    {
        if (timed_out) return;
        if ((++st.nodes & 1023) == 0 && chrono::duration<double>(chrono::steady_clock::now() - start).count() > limit) {
            timed_out = true;
            return;
        }
        vector<vector<double>> h;
        double lb;
        if (! propagate(f, D, h, lb)) return;

        int bu = -1;
        for (int u = 0; u < P.np; ++u)
            if (f[u] == -1 && (bu == -1 || D[u].size() < D[bu].size() || (D[u].size() == D[bu].size() && P.incident[u].size() > P.incident[bu].size()))) bu = u;

        if (bu == -1) {
            double g = lb;
            if (g < ub) {
                ub = g;
                best = f;
                ++st.sols;
            }
            return;
        }

        vector<int> vals = D[bu];
        sort(vals.begin(), vals.end(), [&](int a, int b) { return h[bu][a] < h[bu][b]; });
        for (int x : vals) {
            f[bu] = x;
            auto D2 = D;
            D2[bu] = {x};
            search(f, D2);
            f[bu] = -1;
            if (timed_out) return;
        }
    }
};

int main(int argc, char * argv[])
{
    if (argc < 4) {
        cerr << "usage: " << argv[0] << " bound pattern target [timelimit] [dual_iters]" << endl;
        return 1;
    }
    string bound = argv[1];
    auto pg = read_graph(argv[2]);
    auto tg = read_graph(argv[3]);
    double limit = argc > 4 ? stod(argv[4]) : 60;

    Problem P;
    P.np = pg.names.size();
    P.nt = tg.names.size();
    int nt = P.nt;

    // target edge costs keyed by (label, directed) -> dense matrix
    map<pair<string, bool>, vector<double>> tm;
    for (auto & e : tg.edges) {
        auto & m = tm[{e.label, e.directed}];
        if (m.empty()) m.assign(nt * nt, INF);
        double c = -log(e.w);
        m[e.f * nt + e.t] = c;
        if (! e.directed) m[e.t * nt + e.f] = c;
    }

    P.unary.assign(P.np * nt, INF);
    for (int u = 0; u < P.np; ++u)
        for (int x = 0; x < nt; ++x)
            if (pg.labels[u] == tg.labels[x]) P.unary[u * nt + x] = -log(tg.vweight[x]);

    map<pair<int, int>, int> pair_index;
    P.incident.resize(P.np);
    for (auto & e : pg.edges) {
        auto it = tm.find({e.label, e.directed});
        if (e.f == e.t) {
            for (int x = 0; x < nt; ++x)
                P.unary[e.f * nt + x] += (it == tm.end() ? INF : it->second[x * nt + x]);
            continue;
        }
        int a = min(e.f, e.t), b = max(e.f, e.t);
        auto [pit, fresh] = pair_index.emplace(pair{a, b}, P.pairs.size());
        if (fresh) {
            Pair pr{a, b, vector<double>(nt * nt, 0.0)};
            for (int x = 0; x < nt; ++x)
                pr.c[x * nt + x] = INF;
            P.pairs.push_back(std::move(pr));
            P.incident[a].emplace_back(pit->second, 0);
            P.incident[b].emplace_back(pit->second, 1);
        }
        auto & pr = P.pairs[pit->second];
        for (int x = 0; x < nt; ++x)
            for (int y = 0; y < nt; ++y) {
                // x = f(a), y = f(b)
                double c;
                if (it == tm.end())
                    c = INF;
                else if (e.f == a)
                    c = it->second[x * nt + y];
                else
                    c = it->second[y * nt + x];
                pr.c[x * nt + y] += c;
            }
    }

    P.dom.resize(P.np);
    for (int u = 0; u < P.np; ++u)
        for (int x = 0; x < nt; ++x)
            if (! isinf(P.unary[u * nt + x])) P.dom[u].push_back(x);

    // identity cost (pattern vertex name == target vertex name), for reference
    double idcost = 0;
    bool idok = true;
    {
        vector<int> f(P.np);
        for (int u = 0; u < P.np; ++u) {
            auto it = tg.index.find(pg.names[u]);
            if (it == tg.index.end()) {
                idok = false;
                break;
            }
            f[u] = it->second;
        }
        if (idok) {
            for (int u = 0; u < P.np; ++u)
                idcost += P.unary[u * nt + f[u]];
            for (auto & pr : P.pairs)
                idcost += pr.c[f[pr.u] * nt + f[pr.v]];
        }
    }

    if (bound == "opb-weak" || bound == "opb-strong") {
        // x_u_t: pattern u -> target t. y_k_t_s: pair k maps to (t, s). Costs scaled by 10^4.
        bool strong = bound == "opb-strong";
        auto scale = [](double c) { return (long long)llround(c * 10000); };
        stringstream obj, cons;
        int ncons = 0;
        set<string> vars;
        auto X = [&](int u, int t) { auto s = "x" + to_string(u) + "_" + to_string(t); vars.insert(s); return s; };
        for (int u = 0; u < P.np; ++u) {
            for (int t : P.dom[u]) {
                long long c = scale(P.unary[u * nt + t]);
                if (c) obj << " +" << c << " " << X(u, t);
            }
            for (int t : P.dom[u])
                cons << "+1 " << X(u, t) << " ";
            cons << "= 1 ;\n";
            ++ncons;
        }
        for (int t = 0; t < nt; ++t) {
            int cnt = 0;
            stringstream line;
            for (int u = 0; u < P.np; ++u)
                if (find(P.dom[u].begin(), P.dom[u].end(), t) != P.dom[u].end()) line << "-1 " << X(u, t) << " ", ++cnt;
            if (cnt > 1) cons << line.str() << ">= -1 ;\n", ++ncons;
        }
        for (size_t k = 0; k < P.pairs.size(); ++k) {
            auto & pr = P.pairs[k];
            auto Y = [&](int t, int s) { auto n = "y" + to_string(k) + "_" + to_string(t) + "_" + to_string(s); vars.insert(n); return n; };
            for (int t : P.dom[pr.u])
                for (int s : P.dom[pr.v]) {
                    double c = pr.c[t * nt + s];
                    if (isinf(c)) {
                        if (t != s) cons << "+1 ~" << X(pr.u, t) << " +1 ~" << X(pr.v, s) << " >= 1 ;\n", ++ncons;
                        continue;
                    }
                    long long ci = scale(c);
                    if (strong || ci) {
                        auto y = Y(t, s);
                        if (ci) obj << " +" << ci << " " << y;
                        if (! strong) cons << "+1 " << y << " +1 ~" << X(pr.u, t) << " +1 ~" << X(pr.v, s) << " >= 1 ;\n", ++ncons;
                    }
                }
            if (strong) {
                // sum_s y(t, s) = x(u, t) and sum_t y(t, s) = x(v, s)
                for (int t : P.dom[pr.u]) {
                    for (int s : P.dom[pr.v])
                        if (! isinf(pr.c[t * nt + s])) cons << "+1 " << Y(t, s) << " ";
                    cons << "-1 " << X(pr.u, t) << " = 0 ;\n", ++ncons;
                }
                for (int s : P.dom[pr.v]) {
                    for (int t : P.dom[pr.u])
                        if (! isinf(pr.c[t * nt + s])) cons << "+1 " << Y(t, s) << " ";
                    cons << "-1 " << X(pr.v, s) << " = 0 ;\n", ++ncons;
                }
            }
        }
        cout << "* #variable= " << vars.size() << " #constraint= " << ncons << "\n";
        cout << "min:" << obj.str() << " ;\n"
             << cons.str();
        return 0;
    }

    Solver s{P, bound, limit};
    if (argc > 5) s.dual_iters = stoi(argv[5]);
    s.start = chrono::steady_clock::now();

    // root bound, before any search
    double rootlb = 0;
    {
        vector<int> f(P.np, -1);
        auto D = P.dom;
        vector<vector<double>> h;
        s.propagate(f, D, h, rootlb);
    }
    vector<int> f(P.np, -1);
    s.search(f, P.dom);
    double t = chrono::duration<double>(chrono::steady_clock::now() - s.start).count();

    int matches = 0;
    if (! s.best.empty())
        for (int u = 0; u < P.np; ++u)
            if (tg.names[s.best[u]] == pg.names[u]) ++matches;

    cout << (s.timed_out ? "timeout" : "optimal") << "\t" << s.ub << "\t" << idcost << "\t" << rootlb << "\t" << s.st.nodes << "\t" << t << "\t" << matches << "/" << P.np << endl;
}
