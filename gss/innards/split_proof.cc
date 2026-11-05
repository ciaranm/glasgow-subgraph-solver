#include <gss/counting_statement.hh>
#include <gss/innards/homomorphism_proofs.hh>
#include <gss/innards/proof.hh>
#include <gss/innards/split_proof.hh>

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <map>
#include <numeric>
#include <optional>
#include <ranges>
#include <tuple>

using namespace gss;
using namespace gss::innards;

using std::make_unique;
using std::map;
using std::optional;
using std::pair;
using std::shared_ptr;
using std::string;
using std::to_string;
using std::tuple;
using std::unique_ptr;
using std::vector;

namespace
{
    auto g(int i) -> string
    {
        return "g" + to_string(i);
    }

    auto y(int u, int t) -> string
    {
        return counting_variable_name(u, t);
    }

    auto a(int i) -> string
    {
        return counting_activity_name(i);
    }

    auto t(int i, int j) -> string
    {
        return counting_tag_name(i, j);
    }

    // A node's gate says that its assignment set holds: all of a cube's literals, all of an
    // and's children, some child of an or. Its definition, as (the rows that say the gate
    // implies its definition, the row that says the definition implies the gate), and the
    // definition as one constraint, for the preserved set.
    struct GateDefinition
    {
        optional<string> forward;
        string backward, constraint;
    };

    // A distinct node's gate is the and of its sub-gates (see SplitProof::distinct), which
    // are given.
    auto gate_definition(int i, const CountingNode & node, const vector<string> & sub_gates) -> GateDefinition
    {
        vector<string> literals;
        if (node.kind == CountingNodeKind::Cube)
            for (auto & [u, v] : node.assignment)
                literals.push_back(y(u, v));
        else if (node.kind == CountingNodeKind::Distinct)
            literals = sub_gates;
        else
            for (auto & c : node.children)
                literals.push_back(g(c));

        auto k = to_string(literals.size());
        string forward, backward, constraint;
        if (node.kind == CountingNodeKind::Or) {
            forward = "1 ~" + g(i);
            backward = k + " " + g(i);
            for (auto & l : literals) {
                forward += " 1 " + l;
                backward += " 1 ~" + l;
                constraint += "1 " + l + " ";
            }
            return GateDefinition{forward + " >= 1", backward + " >= " + k, constraint + ">= 1"};
        }

        if (literals.empty())
            return GateDefinition{std::nullopt, "1 " + g(i) + " >= 1", ">= 0"};
        forward = k + " ~" + g(i);
        backward = "1 " + g(i);
        for (auto & l : literals) {
            forward += " 1 " + l;
            backward += " 1 ~" + l;
            constraint += "1 " + l + " ";
        }
        return GateDefinition{forward + " >= " + k, backward + " >= 1", constraint + ">= " + k};
    }

    auto ids(const vector<long> & lines) -> string
    {
        string result;
        for (auto & l : lines)
            result += " " + to_string(l);
        return result;
    }
}

struct SplitProof::Imp
{
    shared_ptr<Proof> proof;
    HomomorphismProofs * hom_proofs = nullptr;
    int pattern_size = 0, target_size = 0;
    std::string statement_file;
    std::ofstream statement_stream;
    bool finished = false;
    optional<CountingStatementWriter> writer;
    // the two halves of each copy y_p_t <-> x_p_t, and the at-most-one rows over y
    vector<vector<long>> y_implies_x, x_implies_y;
    vector<long> y_at_most_one;
    int number_of_assumptions = 0;

    // The statement so far, for the end of the proof, with each node's gate rows, and for
    // each or, the vertex that tells its children apart and their values there.
    CountingStatement statement;
    vector<long> forward, backward;
    map<int, pair<int, vector<int>>> separators;

    map<pair<vector<int>, vector<int>>, int> distincts;

    // A distinct node's gate is the and of sub-gates, each defined over y, and shared
    // between nodes: "u is at one of these values", and "at most one of these vertices is
    // at t". Each at-most-one sub-gate holds, by the input's injectivity, and those facts,
    // with the at-most-one rows over y they come from, are deleted before injectivity is.
    struct SubGate
    {
        string name, constraint;
        long forward = 0, backward = 0;
    };
    vector<SubGate> sub_gates;
    map<pair<int, vector<int>>, int> in_set_gates;
    map<pair<vector<int>, int>, int> at_most_one_gates;
    map<pair<vector<int>, int>, long> y_at_most_one_over;
    vector<long> injectivity_lemmas;

    // For each distinct node, its sub-gates, by index: one per vertex, then, with two or
    // more vertices, one per value. The observer takes the next node's from pending.
    map<int, vector<int>> distinct_sub_gates;
    vector<int> pending_sub_gates;

    auto sub_gate_names(int node) const -> vector<string>
    {
        vector<string> result;
        if (auto d = distinct_sub_gates.find(node); d != distinct_sub_gates.end())
            for (auto & i : d->second)
                result.push_back(sub_gates[i].name);
        return result;
    }

    auto define_sub_gate(const string & name, const string & forward, const string & backward, const string & constraint) -> int
    {
        SubGate gate{name, constraint};
        proof->at_level_zero([&] {
            gate.forward = proof->emit_proof_line("red " + forward + " : " + name + " -> 0 ;");
            gate.backward = proof->emit_proof_line("red " + backward + " : " + name + " -> 1 ;");
            proof->emit_proof_directive("core id " + to_string(gate.forward) + " " + to_string(gate.backward) + " ;");
        });
        sub_gates.push_back(gate);
        return sub_gates.size() - 1;
    }

    // u is at one of these values: an or over y
    auto in_set_gate(int u, const vector<int> & values) -> int
    {
        auto [where, inserted] = in_set_gates.emplace(pair{u, values}, -1);
        if (inserted) {
            auto name = "ds" + to_string(sub_gates.size());
            auto d = to_string(values.size());
            string forward = "1 ~" + name, backward = d + " " + name, constraint;
            for (auto & t : values) {
                forward += " 1 " + y(u, t);
                backward += " 1 ~" + y(u, t);
                constraint += "1 " + y(u, t) + " ";
            }
            where->second = define_sub_gate(name, forward + " >= 1", backward + " >= " + d, constraint + ">= 1");
        }
        return where->second;
    }

    // At most one of these vertices is at t, over y, from the input's injectivity row on t
    // and the copies y -> x, weakening away the other vertices.
    auto y_at_most_one_of(const vector<int> & vertices, int t) -> long
    {
        auto [where, inserted] = y_at_most_one_over.emplace(pair{vertices, t}, 0);
        if (inserted) {
            string line = "pol " + proof->injectivity_labels().at(t);
            for (auto & u : vertices)
                line += " " + to_string(y_implies_x[u][t]) + " +";
            for (int v = 0; v < pattern_size; ++v)
                if (! std::ranges::binary_search(vertices, v))
                    line += " x" + proof->variable_name(v, t) + " w";
            proof->at_level_zero([&] { where->second = proof->emit_proof_line(line + " ;"); });
            injectivity_lemmas.push_back(where->second);
        }
        return where->second;
    }

    // At most one of these vertices is at t: a reified cardinality constraint, which holds.
    auto at_most_one_gate(const vector<int> & vertices, int t) -> int
    {
        auto [where, inserted] = at_most_one_gates.emplace(pair{vertices, t}, -1);
        if (inserted) {
            auto name = "ds" + to_string(sub_gates.size());
            auto k = to_string(vertices.size() - 1);
            string forward = k + " ~" + name, backward = "2 " + name, constraint;
            for (auto & u : vertices) {
                forward += " 1 ~" + y(u, t);
                backward += " 1 " + y(u, t);
                constraint += "1 ~" + y(u, t) + " ";
            }
            where->second = define_sub_gate(name, forward + " >= " + k, backward + " >= 2", constraint + ">= " + k);
            auto amo = y_at_most_one_of(vertices, t);
            proof->at_level_zero([&] {
                injectivity_lemmas.push_back(proof->emit_proof_line(
                    "pol " + to_string(sub_gates[where->second].backward) + " " + to_string(amo) + " + s ;"));
            });
        }
        return where->second;
    }
};

SplitProof::SplitProof(const shared_ptr<Proof> & proof, HomomorphismProofs * hom_proofs, int pattern_size, int target_size,
    const string & statement_file) :
    _imp(make_unique<Imp>())
{
    _imp->proof = proof;
    _imp->hom_proofs = hom_proofs;
    _imp->pattern_size = pattern_size;
    _imp->target_size = target_size;
    _imp->statement_file = statement_file;
    _imp->statement_stream.open(statement_file);
    if (! _imp->statement_stream)
        throw ProofError{"Error writing counting statement to '" + statement_file + "'"};
    _imp->statement.domain_sizes.assign(pattern_size, target_size);
    _imp->writer.emplace(_imp->statement_stream, _imp->statement.domain_sizes);
    // Without a proof, the statement is only written, as a compressed form of the solutions.
    if (! _imp->proof)
        return;
    // Each node's gate is defined, both ways, into the core, while it is still fresh: VeriPB
    // checks a definition against the derived constraints too, so the rows that say what a
    // gate implies cannot wait until the search has derived things about it.
    _imp->writer->set_observer([this](int i, const CountingNode & node) {
        auto & proof = *_imp->proof;
        if (node.kind == CountingNodeKind::Distinct) {
            _imp->distinct_sub_gates.emplace(i, std::move(_imp->pending_sub_gates));
            _imp->pending_sub_gates.clear();
        }
        auto definition = gate_definition(i, node, _imp->sub_gate_names(i));
        proof.at_level_zero([&] {
            long forward = 0;
            if (definition.forward)
                forward = proof.emit_proof_line("red " + *definition.forward + " : " + g(i) + " -> 0 ;");
            auto backward = proof.emit_proof_line("red " + definition.backward + " : " + g(i) + " -> 1 ;");
            proof.emit_proof_directive("core id" + (forward ? " " + to_string(forward) : string{}) + " " + to_string(backward) + " ;");
            _imp->forward.push_back(forward);
            _imp->backward.push_back(backward);
        });
        _imp->statement.nodes.push_back(CountingNode{node.kind, 0, node.assignment, node.children, node.variables, node.values});
    });
}

SplitProof::~SplitProof()
{
    // an incomplete statement, or none at all because there were no solutions, is not left
    // lying around to be mistaken for a real one
    if (! _imp->finished) {
        _imp->statement_stream.close();
        std::remove(_imp->statement_file.c_str());
    }
}

auto SplitProof::proving() const -> bool
{
    return bool(_imp->proof);
}

auto SplitProof::begin() -> void
{
    if (! _imp->proof)
        return;
    auto & proof = *_imp->proof;

    // Each y_p_t is a copy of x_p_t, and joins the preserved set. (The x side leaves it in
    // finish, once nothing but the copies mentions x.)
    proof.emit_proof_directive("% the statement's variables are copies of the encoding's");
    _imp->y_implies_x.resize(_imp->pattern_size);
    _imp->x_implies_y.resize(_imp->pattern_size);
    for (int p = 0; p < _imp->pattern_size; ++p)
        for (int v = 0; v < _imp->target_size; ++v) {
            auto x = "x" + proof.variable_name(p, v);
            auto forward = proof.emit_proof_line("red 1 ~" + y(p, v) + " 1 " + x + " >= 1 : " + y(p, v) + " -> 0 ;");
            auto backward = proof.emit_proof_line("red 1 " + y(p, v) + " 1 ~" + x + " >= 1 : " + y(p, v) + " -> 1 ;");
            proof.emit_proof_directive("core id " + to_string(forward) + " " + to_string(backward) + " ;");
            proof.emit_proof_directive("preserved_add " + y(p, v) + " : 1 " + x + " >= 1 ;");
            _imp->y_implies_x[p].push_back(forward);
            _imp->x_implies_y[p].push_back(backward);
        }

    // The exactly-one rows over y, which the tagged formula has: at least one by unit
    // propagation through the copies, and at most one by adding the y -> x halves of the
    // copies to the encoding's at most one over x.
    for (int p = 0; p < _imp->pattern_size; ++p) {
        string at_least = "rup", at_most = "pol " + proof.at_most_one_value_label(p);
        for (int v = 0; v < _imp->target_size; ++v) {
            at_least += " 1 " + y(p, v);
            at_most += " " + to_string(_imp->y_implies_x[p][v]) + " +";
        }
        auto al = proof.emit_proof_line(at_least + " >= 1 ;");
        auto am = proof.emit_proof_line(at_most + " ;");
        _imp->y_at_most_one.push_back(am);
        proof.emit_proof_directive("core id " + to_string(al) + " " + to_string(am) + " ;");
    }
}

auto SplitProof::cube(vector<pair<int, int>> assignment) -> int
{
    // the search's assignments can list a vertex more than once, always with the same value
    std::sort(assignment.begin(), assignment.end());
    assignment.erase(std::unique(assignment.begin(), assignment.end()), assignment.end());
    return _imp->writer->add_cube(assignment);
}

auto SplitProof::conjunction(vector<int> children) -> int
{
    if (1 == children.size())
        return children[0];
    std::sort(children.begin(), children.end());
    return _imp->writer->add_and(children);
}

auto SplitProof::disjunction(const vector<int> & children, int vertex, const vector<int> & values) -> int
{
    if (1 == children.size())
        return children[0];

    // the statement lists children in increasing order, and the tags follow that order
    vector<unsigned> order(children.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](unsigned i, unsigned j) { return children[i] < children[j]; });
    vector<int> sorted_children, sorted_values;
    for (auto & i : order) {
        sorted_children.push_back(children[i]);
        sorted_values.push_back(values.at(i));
    }

    auto id = _imp->writer->add_or(sorted_children);
    if (_imp->proof)
        _imp->separators.emplace(id, pair{vertex, move(sorted_values)});
    return id;
}

auto SplitProof::distinct(const vector<int> & vertices, const vector<int> & values) -> int
{
    auto [where, inserted] = _imp->distincts.emplace(pair{vertices, values}, -1);
    if (inserted) {
        if (_imp->proof) {
            _imp->pending_sub_gates.clear();
            for (auto & u : vertices)
                _imp->pending_sub_gates.push_back(_imp->in_set_gate(u, values));
            if (vertices.size() >= 2)
                for (auto & t : values)
                    _imp->pending_sub_gates.push_back(_imp->at_most_one_gate(vertices, t));
        }
        where->second = _imp->writer->add_distinct(vertices, values);
    }
    return where->second;
}

auto SplitProof::distinct_lemmas(const vector<string> & premises, const vector<int> & vertices, const vector<int> & values) -> vector<long>
{
    vector<long> result;
    if (! _imp->proof)
        return result;
    for (auto & u : vertices) {
        auto h = _imp->in_set_gates.at(pair{u, values});
        auto literals = premises;
        literals.push_back(_imp->sub_gates[h].name);
        result.push_back(rup(literals, true));
    }
    return result;
}

auto SplitProof::distinct_impossible(const vector<string> & premises, const vector<int> & vertices, const vector<int> & values) -> vector<long>
{
    vector<long> result;
    if (! _imp->proof)
        return result;
    // each vertex is at one of the values, and each value has at most one of them, so
    // adding everything up, k vertices fit into fewer than k values, which they cannot
    for (auto & u : vertices) {
        auto literals = premises;
        for (auto & t : values)
            literals.push_back(y(u, t));
        result.push_back(rup(literals, true));
    }
    if (! values.empty()) {
        string line = "pol";
        for (unsigned i = 0; i < result.size(); ++i)
            line += " " + to_string(result[i]) + (0 == i ? "" : " +");
        for (auto & t : values)
            line += " " + to_string(_imp->y_at_most_one_of(vertices, t)) + " +";
        line += " " + to_string(vertices.size()) + " d ;";
        _imp->proof->at_level_zero([&] { result.push_back(_imp->proof->emit_proof_line(line)); });
    }
    return result;
}

auto SplitProof::gate(int node) const -> string
{
    return g(node);
}

auto SplitProof::statement_variable(int v, int t) const -> string
{
    return y(v, t);
}

auto SplitProof::assumption(const vector<pair<int, int>> & excluded) -> string
{
    auto name = "sa" + to_string(++_imp->number_of_assumptions);
    if (! _imp->proof)
        return name;
    auto & proof = *_imp->proof;
    auto k = to_string(excluded.size());
    string all = "red " + k + " ~" + name, not_all = "red 1 " + name;
    for (auto & [v, t] : excluded) {
        all += " 1 ~x" + proof.variable_name(v, t);
        not_all += " 1 x" + proof.variable_name(v, t);
    }
    proof.at_level_zero([&] {
        if (! excluded.empty())
            proof.emit_proof_line(all + " >= " + k + " : " + name + " -> 0 ;");
        proof.emit_proof_line(not_all + " >= 1 : " + name + " -> 1 ;");
    });
    return name;
}

auto SplitProof::rup(const vector<string> & literals, bool at_level_zero) -> long
{
    if (! _imp->proof)
        return 0;
    string line = "rup";
    for (auto & l : literals)
        line += " 1 " + l;
    line += " >= 1 ;";
    long result = 0;
    if (at_level_zero)
        _imp->proof->at_level_zero([&] { result = _imp->proof->emit_proof_line(line); });
    else
        result = _imp->proof->emit_proof_line(line);
    return result;
}

auto SplitProof::derive_signatures(const vector<std::tuple<int, int, int>> & lemmas) -> void
{
    // Each follows by unit propagation through the gates' definitions: every cube below the
    // node, by every path, puts the vertex at the value, or lies below a node with a lemma
    // that says so.
    if (! _imp->proof)
        return;
    auto & proof = *_imp->proof;
    proof.at_level_zero([&] {
        for (auto & [node, v, t] : lemmas)
            proof.emit_proof_line("rup 1 ~" + gate(node) + " 1 " + statement_variable(v, t) + " >= 1 ;");
    });
}

auto SplitProof::forget(const vector<long> & lines) -> void
{
    if (lines.empty() || ! _imp->proof)
        return;
    _imp->proof->at_level_zero([&] { _imp->proof->emit_proof_directive("del id" + ids(lines) + " ;"); });
}

namespace
{
    // Lemmas "node m's gate implies a fact about vertex v", each by one hinted RUP, built
    // per node bottom up: from a cube's forward row, from one child's lemma for an and, or
    // from every child's lemma for an or. A fact is y_v_w ("is"), ~y_v_w ("not"), or "v is
    // at none of these values" ("outside", with the values as a set by number).
    class VertexLemmas
    {
    public:
        enum class Fact
        {
            Is,
            Not,
            Outside
        };

        // the fact, the vertex, and the value or the set of values
        using Key = tuple<Fact, int, int>;

    private:
        Proof & _proof;
        const CountingStatement & _statement;
        const vector<long> & _forward;
        const vector<long> & _y_at_most_one;
        int _target_size;
        vector<vector<int>> _value_sets;
        map<vector<int>, int> _value_set_numbers;

        // how each (node, fact) is refuted: -1 a cube, -2 an or, -3 not at all, or the child
        // an and uses
        map<pair<int, Key>, int> _how;
        map<pair<int, Key>, long> _lemmas;

    public:
        vector<pair<long, vector<long>>> lemmas_and_hints;

        VertexLemmas(Proof & proof, const CountingStatement & statement, const vector<long> & forward,
            const vector<long> & y_at_most_one, int target_size) :
            _proof(proof),
            _statement(statement),
            _forward(forward),
            _y_at_most_one(y_at_most_one),
            _target_size(target_size)
        {
        }

        // v is at none of these values, which must be sorted
        auto outside(int v, const vector<int> & values) -> Key
        {
            auto [it, inserted] = _value_set_numbers.emplace(values, _value_sets.size());
            if (inserted)
                _value_sets.push_back(values);
            return Key{Fact::Outside, v, it->second};
        }

        auto supports(const Key & key, int value) const -> bool
        {
            auto & [fact, _, w] = key;
            switch (fact) {
            case Fact::Is: return value == w;
            case Fact::Not: return value != w;
            case Fact::Outside: return ! std::ranges::binary_search(_value_sets[w], value);
            }
            return false;
        }

        auto row(int node, const Key & key) const -> string
        {
            auto & [fact, v, w] = key;
            string result = "1 ~" + g(node);
            switch (fact) {
            case Fact::Is: result += " 1 " + y(v, w); break;
            case Fact::Not: result += " 1 ~" + y(v, w); break;
            case Fact::Outside:
                for (int value = 0; value < _target_size; ++value)
                    if (! std::ranges::binary_search(_value_sets[w], value))
                        result += " 1 " + y(v, value);
                break;
            }
            return result + " >= 1";
        }

        auto plan(int n, const Key & key) -> void
        {
            auto v = std::get<1>(key);
            vector<pair<int, bool>> stack{{n, false}};
            while (! stack.empty()) {
                auto [m, expanded] = stack.back();
                stack.pop_back();
                if (_how.contains({m, key}))
                    continue;
                auto & node = _statement.nodes[m];
                if (node.kind == CountingNodeKind::Distinct) {
                    // never needed: an or's children are told apart by vertices outside
                    // any distinct node below them
                    _how.emplace(pair{m, key}, -3);
                    continue;
                }
                if (node.kind == CountingNodeKind::Cube) {
                    bool ok = std::ranges::any_of(node.assignment, [&](const auto & uv) { return uv.first == v && supports(key, uv.second); });
                    _how.emplace(pair{m, key}, ok ? -1 : -3);
                    continue;
                }
                if (! expanded) {
                    stack.emplace_back(m, true);
                    for (auto & c : node.children)
                        if (! _how.contains({c, key}))
                            stack.emplace_back(c, false);
                    continue;
                }
                if (node.kind == CountingNodeKind::And) {
                    // any child that refutes it will do, but a cube is quickest
                    int chosen = -3;
                    for (auto & c : node.children)
                        if (_how.at({c, key}) != -3 && (chosen == -3 || _statement.nodes[c].kind == CountingNodeKind::Cube))
                            chosen = c;
                    _how.emplace(pair{m, key}, chosen);
                }
                else
                    _how.emplace(pair{m, key}, std::ranges::all_of(node.children, [&](int c) { return _how.at({c, key}) != -3; }) ? -2 : -3);
            }
        }

        // The lemma for node n, deriving it and those below it as needed, into the core.
        auto get(int n, const Key & key) -> long
        {
            if (auto l = _lemmas.find({n, key}); l != _lemmas.end())
                return l->second;

            plan(n, key);
            vector<long> fresh;
            if (-3 == _how.at({n, key})) {
                // it does not hold, so there is no proof of it: leave VeriPB to say so
                auto line = _proof.emit_proof_line("rup " + row(n, key) + " ;");
                _lemmas.emplace(pair{n, key}, line);
                lemmas_and_hints.emplace_back(line, vector<long>{});
                fresh.push_back(line);
            }
            else {
                vector<pair<int, bool>> stack{{n, false}};
                while (! stack.empty()) {
                    auto [m, expanded] = stack.back();
                    stack.pop_back();
                    if (_lemmas.contains({m, key}))
                        continue;
                    auto how = _how.at({m, key});
                    vector<int> below;
                    if (-2 == how)
                        below = _statement.nodes[m].children;
                    else if (how >= 0)
                        below = {how};
                    if (! expanded) {
                        stack.emplace_back(m, true);
                        for (auto & c : below)
                            if (! _lemmas.contains({c, key}))
                                stack.emplace_back(c, false);
                        continue;
                    }
                    vector<long> hints{_forward[m]};
                    for (auto & c : below)
                        hints.push_back(_lemmas.at({c, key}));
                    if (-1 == how && std::get<0>(key) == Fact::Not)
                        hints.push_back(_y_at_most_one[std::get<1>(key)]);
                    auto line = _proof.emit_proof_line("rup " + row(m, key) + " : ~" + ids(hints) + " ;");
                    _lemmas.emplace(pair{m, key}, line);
                    lemmas_and_hints.emplace_back(line, move(hints));
                    fresh.push_back(line);
                }
            }
            _proof.emit_proof_directive("core id" + ids(fresh) + " ;");
            return _lemmas.at({n, key});
        }
    };
}

auto SplitProof::finish(int root) -> void
{
    _imp->writer->finish(root);
    _imp->statement_stream.close();
    _imp->finished = true;
    if (! _imp->proof)
        return;

    auto & proof = *_imp->proof;
    auto x = [&](int p, int t) { return "x" + proof.variable_name(p, t); };
    auto n = _imp->pattern_size;

    proof.emit_proof_directive("setlvl 0;");
    auto root_line = proof.emit_proof_line("rup 1 " + gate(root) + " >= 1 ;");
    proof.emit_proof_directive("core id " + to_string(root_line) + " ;");

    // what the distinct nodes' sub-gates took from injectivity goes before it does
    if (! _imp->injectivity_lemmas.empty())
        proof.emit_proof_directive("del id" + ids(_imp->injectivity_lemmas) + " ;");

    // Every adjacency row goes by unit propagation from the root: under the negated row,
    // every cube that puts its first vertex where the row says has its second vertex
    // somewhere the row rules out.
    proof.emit_proof_directive("% the root implies adjacency");
    for (auto & label : _imp->hom_proofs->model_adjacency_labels())
        proof.emit_proof_directive("del id " + label + " ;");

    // Each injectivity row: for each pair of pattern vertices, not both at t, by unit
    // propagation from the root. Within a cube that is immediate. Across the components of
    // a split, two vertices that could both be at t are recorded at t, and the split never
    // pairs two signatures that both put a vertex there. Then recover the at most one
    // from the pairs.
    proof.emit_proof_directive("% the root implies injectivity");
    for (auto & [t, label] : proof.injectivity_labels()) {
        if (n < 2) {
            proof.emit_proof_directive("del id " + label + " ;");
            continue;
        }

        proof.begin_deletion_subproof("del id " + label + " : : subproof");
        auto negated = "@sneg" + to_string(t);
        proof.emit_proof_line(negated + " pol -1 ;");

        auto pair_label = [&](int u, int v) { return "@sp" + to_string(t) + "_" + to_string(u) + "_" + to_string(v); };
        for (int u = 0; u < n; ++u)
            for (int v = u + 1; v < n; ++v)
                proof.emit_proof_line(pair_label(u, v) + " rup 1 ~" + x(u, t) + " 1 ~" + x(v, t) + " >= 1 ;");

        // at most one, from the pairs, a vertex at a time
        proof.emit_proof_line("@samo" + to_string(t) + "_2 pol " + pair_label(0, 1) + " ;");
        for (int k = 2; k < n; ++k) {
            string line = "@samo" + to_string(t) + "_" + to_string(k + 1) + " pol @samo" + to_string(t) + "_" + to_string(k) + " " + to_string(k - 1) + " *";
            for (int i = 0; i < k; ++i)
                line += " " + pair_label(i, k) + " +";
            proof.emit_proof_line(line + " " + to_string(k) + " d ;");
        }
        proof.emit_proof_line("pol @samo" + to_string(t) + "_" + to_string(n) + " " + negated + " + ;");
        proof.emit_proof_directive("qed : -1 ;");
    }

    // The exactly-one rows over x follow from the ones over y and the copies.
    proof.emit_proof_directive("% the exactly-one rows over x follow from the ones over y");
    for (int p = 0; p < n; ++p) {
        proof.emit_proof_directive("del id " + proof.at_least_one_value_label(p) + " ;");
        string at_most = "pol -1 " + to_string(_imp->y_at_most_one[p]) + " +";
        for (int t = 0; t < _imp->target_size; ++t)
            at_most += " " + to_string(_imp->x_implies_y[p][t]) + " +";
        proof.begin_deletion_subproof("del id " + proof.at_most_one_value_label(p) + " : : subproof");
        proof.emit_proof_line(at_most + " ;");
        proof.emit_proof_directive("qed : -1 ;");
    }

    // Now x is only in the copies: it leaves the preserved set, and then the copies go,
    // with x -> y as the witness.
    proof.emit_proof_directive("% x leaves the preserved set");
    for (int p = 0; p < n; ++p)
        for (int t = 0; t < _imp->target_size; ++t)
            proof.emit_proof_directive("preserved_rm " + x(p, t) + " : 1 " + y(p, t) + " >= 1 ;");
    for (int p = 0; p < n; ++p)
        for (int t = 0; t < _imp->target_size; ++t) {
            proof.emit_proof_directive("del id " + to_string(_imp->y_implies_x[p][t]) + " : " + x(p, t) + " -> " + y(p, t) + " ;");
            proof.emit_proof_directive("del id " + to_string(_imp->x_implies_y[p][t]) + " : " + x(p, t) + " -> " + y(p, t) + " ;");
        }

    finish_tagged(root, root_line);

    proof.emit_proof_directive("output EQUIENUMERABLE FILE ;");
    proof.emit_proof_directive("conclusion NONE ;");
    proof.emit_proof_directive("end pseudo-Boolean proof ;");
}

auto SplitProof::finish_tagged(int root, long root_line) -> void
{
    auto & proof = *_imp->proof;
    auto & statement = _imp->statement;
    statement.root = root;
    auto & nodes = statement.nodes;
    int size = nodes.size();
    auto & forward = _imp->forward;
    auto & backward = _imp->backward;
    auto & y_at_most_one = _imp->y_at_most_one;

    // Each node's parents, as (parent, position), in order of parent, so in the order of
    // its selectors.
    vector<vector<pair<int, int>>> parents(size);
    for (int p = 0; p < size; ++p)
        for (unsigned j = 0; j < nodes[p].children.size(); ++j)
            parents[nodes[p].children[j]].emplace_back(p, j);
    auto selector = [&](int p, int j) { return CountingNodeKind::And == nodes[p].kind ? a(p) : t(p, j); };
    auto selector_index = [&](int c, int p, int j) {
        auto index = std::ranges::find(parents[c], pair{p, j}) - parents[c].begin();
        if (index == std::ssize(parents[c]))
            throw ProofError{"node " + to_string(c) + " has no parent " + to_string(p) + " at position " + to_string(j)};
        return index;
    };
    auto selectors = [&](int i) {
        string result;
        for (auto & [p, j] : parents[i])
            result += "1 " + selector(p, j) + " ";
        return result + ">= 1";
    };

    auto & sub_gates = _imp->sub_gates;
    proof.emit_proof_directive("% the tagged formula: gates, and distinct nodes' sub-gates, join the preserved set, for now");
    for (auto & sub : sub_gates)
        proof.emit_proof_directive("preserved_add " + sub.name + " : " + sub.constraint + " ;");
    for (int i = 0; i < size; ++i)
        proof.emit_proof_directive("preserved_add " + g(i) + " : " + gate_definition(i, nodes[i], _imp->sub_gate_names(i)).constraint + " ;");

    // Activity and tags, top down, each defined from the gates. The root is active when
    // its gate holds, any other node when one of its selectors is, and a tag when its or
    // is active and its child's gate holds. Deleting each tag's "if" row, later, is where
    // the or's children are shown to be disjoint.
    proof.emit_proof_directive("% activity and tags, top down");
    long root_zero = 0, root_one = 0;
    vector<long> zero(size);
    vector<vector<long>> one(size), tag_zero(size), tag_child(size), tag_if(size);
    for (int i = size - 1; i >= 0; --i) {
        if (i == root) {
            root_zero = proof.emit_proof_line("red 1 ~" + a(i) + " 1 " + g(i) + " >= 1 : " + a(i) + " -> 0 ;");
            root_one = proof.emit_proof_line("red 1 " + a(i) + " 1 ~" + g(i) + " >= 1 : " + a(i) + " -> 1 ;");
            proof.emit_proof_directive("core id " + to_string(root_zero) + " " + to_string(root_one) + " ;");
            proof.emit_proof_directive("preserved_add " + a(i) + " : 1 " + g(i) + " >= 1 ;");
        }
        else {
            string row = "1 ~" + a(i);
            for (auto & [p, j] : parents[i])
                row += " 1 " + selector(p, j);
            zero[i] = proof.emit_proof_line("red " + row + " >= 1 : " + a(i) + " -> 0 ;");
            for (auto & [p, j] : parents[i])
                one[i].push_back(proof.emit_proof_line("red 1 " + a(i) + " 1 ~" + selector(p, j) + " >= 1 : " + a(i) + " -> 1 ;"));
            proof.emit_proof_directive("core id " + to_string(zero[i]) + ids(one[i]) + " ;");
            proof.emit_proof_directive("preserved_add " + a(i) + " : " + selectors(i) + " ;");
        }

        if (CountingNodeKind::Or == nodes[i].kind)
            for (unsigned j = 0; j < nodes[i].children.size(); ++j) {
                auto c = nodes[i].children[j];
                tag_zero[i].push_back(proof.emit_proof_line("red 1 ~" + t(i, j) + " 1 " + a(i) + " >= 1 : " + t(i, j) + " -> 0 ;"));
                tag_child[i].push_back(proof.emit_proof_line("red 1 ~" + t(i, j) + " 1 " + g(c) + " >= 1 : " + t(i, j) + " -> 0 ;"));
                tag_if[i].push_back(proof.emit_proof_line("red 1 " + t(i, j) + " 1 ~" + a(i) + " 1 ~" + g(c) + " >= 1 : " + t(i, j) + " -> 1 ;"));
                proof.emit_proof_directive("core id " + to_string(tag_zero[i][j]) + " " + to_string(tag_child[i][j]) + " " + to_string(tag_if[i][j]) + " ;");
                proof.emit_proof_directive("preserved_add " + t(i, j) + " : 1 " + a(i) + " 1 " + g(c) + " >= 2 ;");
            }
    }

    // An active node's gate holds: through a parent's tag, or an and parent's own activity
    // and gate.
    proof.emit_proof_directive("% activity implies the gate, top down");
    vector<long> active(size);
    active[root] = root_zero;
    for (int i = size - 1; i >= 0; --i) {
        if (i == root)
            continue;
        vector<long> hints;
        for (auto & [p, j] : parents[i]) {
            if (CountingNodeKind::And == nodes[p].kind) {
                hints.push_back(forward[p]);
                hints.push_back(active[p]);
            }
            else
                hints.push_back(tag_child[p][j]);
        }
        hints.push_back(zero[i]);
        active[i] = proof.emit_proof_line("rup 1 ~" + a(i) + " 1 " + g(i) + " >= 1 : ~" + ids(hints) + " ;");
        proof.emit_proof_directive("core id " + to_string(active[i]) + " ;");
    }

    proof.emit_proof_directive("% an active cube's literals hold");
    vector<vector<long>> cube_rows(size);
    for (int i = 0; i < size; ++i)
        if (CountingNodeKind::Cube == nodes[i].kind) {
            for (auto & [u, v] : nodes[i].assignment)
                cube_rows[i].push_back(proof.emit_proof_line("rup 1 ~" + a(i) + " 1 " + y(u, v) + " >= 1 : ~ " + to_string(active[i]) + " " + to_string(forward[i]) + " ;"));
            if (! cube_rows[i].empty())
                proof.emit_proof_directive("core id" + ids(cube_rows[i]) + " ;");
        }

    // An active distinct node's vertices are each at one of its values, from its sub-gate
    // "in the set", and at most one is at each value, from its at-most-one sub-gate, which
    // is a cutting planes step.
    proof.emit_proof_directive("% an active distinct node's rows hold");
    vector<vector<long>> distinct_in(size), distinct_at_most(size);
    vector<long> scratch;
    for (int i = 0; i < size; ++i)
        if (CountingNodeKind::Distinct == nodes[i].kind) {
            auto & sub = _imp->distinct_sub_gates.at(i);
            auto & variables = nodes[i].variables;
            auto k = variables.size();
            for (unsigned j = 0; j < k; ++j) {
                string row = "1 ~" + a(i);
                for (auto & value : nodes[i].values)
                    row += " 1 " + y(variables[j], value);
                distinct_in[i].push_back(proof.emit_proof_line("rup " + row + " >= 1 : ~ " + to_string(active[i]) + " " + to_string(forward[i]) + " " + to_string(sub_gates[sub[j]].forward) + " ;"));
            }
            if (k >= 2)
                for (unsigned j = 0; j < nodes[i].values.size(); ++j) {
                    auto & m = sub_gates[sub[k + j]];
                    auto holds = proof.emit_proof_line("rup 1 ~" + a(i) + " 1 " + m.name + " >= 1 : ~ " + to_string(active[i]) + " " + to_string(forward[i]) + " ;");
                    scratch.push_back(holds);
                    distinct_at_most[i].push_back(proof.emit_proof_line("pol " + to_string(holds) + " " + to_string(k - 1) + " * " + to_string(m.forward) + " + ;"));
                }
            vector<long> rows = distinct_in[i];
            rows.insert(rows.end(), distinct_at_most[i].begin(), distinct_at_most[i].end());
            proof.emit_proof_directive("core id" + ids(rows) + " ;");
        }
    if (! scratch.empty())
        proof.emit_proof_directive("del id" + ids(scratch) + " ;");

    // For each or, each child's lemma about the vertex that tells the children apart, and
    // for a child at none of the others' values, that it is at none of them one by one.
    proof.emit_proof_directive("% what each or's children do to the vertex that tells them apart");
    VertexLemmas lemmas(proof, statement, forward, y_at_most_one, _imp->target_size);
    map<int, vector<long>> vertex_lemma;
    map<int, vector<vector<long>>> not_lemmas;
    map<int, vector<int>> named_values;
    for (auto & [i, separator] : _imp->separators) {
        auto & [v, values] = separator;
        vector<int> named;
        for (auto & value : values)
            if (-1 != value)
                named.push_back(value);
        std::ranges::sort(named);
        named_values[i] = named;
        for (unsigned j = 0; j < values.size(); ++j) {
            auto c = nodes[i].children[j];
            not_lemmas[i].emplace_back();
            if (-1 != values[j])
                vertex_lemma[i].push_back(lemmas.get(c, {VertexLemmas::Fact::Is, v, values[j]}));
            else {
                vertex_lemma[i].push_back(lemmas.get(c, lemmas.outside(v, named)));
                for (auto & value : named)
                    not_lemmas[i].back().push_back(lemmas.get(c, {VertexLemmas::Fact::Not, v, value}));
            }
        }
    }

    // The tagged formula's own rows: the root is active, and each active or has exactly
    // one tag. (The rest are already the definitions' rows.) The at-most-one is one
    // cutting planes step, over each tag's child's lemma and the at-most-one on the vertex.
    proof.emit_proof_directive("% the root is active, and each active or has exactly one tag");
    auto root_active = proof.emit_proof_line("rup 1 " + a(root) + " >= 1 : ~ " + to_string(root_line) + " " + to_string(root_one) + " ;");
    proof.emit_proof_directive("core id " + to_string(root_active) + " ;");
    vector<long> at_least_one_tag(size);
    for (int i = 0; i < size; ++i) {
        if (CountingNodeKind::Or != nodes[i].kind)
            continue;
        auto & children = nodes[i].children;
        string row = "1 ~" + a(i);
        for (unsigned j = 0; j < children.size(); ++j)
            row += " 1 " + t(i, j);
        at_least_one_tag[i] = proof.emit_proof_line("rup " + row + " >= 1 : ~ " + to_string(active[i]) + " " + to_string(forward[i]) + ids(tag_if[i]) + " ;");
        proof.emit_proof_directive("core id " + to_string(at_least_one_tag[i]) + " ;");

        auto & [v, values] = _imp->separators.at(i);
        auto & named = named_values.at(i);
        string sum = "pol";
        bool rest = false;
        for (unsigned j = 0; j < children.size(); ++j) {
            string tag_row = "1 ~" + t(i, j);
            if (-1 != values[j])
                tag_row += " 1 " + y(v, values[j]);
            else {
                rest = true;
                for (int value = 0; value < _imp->target_size; ++value)
                    if (! std::ranges::binary_search(named, value))
                        tag_row += " 1 " + y(v, value);
            }
            auto line = proof.emit_proof_line("rup " + tag_row + " >= 1 : ~ " + to_string(tag_child[i][j]) + " " + to_string(vertex_lemma[i][j]) + " ;");
            sum += " " + to_string(line) + (0 == j ? "" : " +");
        }
        sum += " " + to_string(y_at_most_one[v]) + " +";
        if (! rest)
            for (int value = 0; value < _imp->target_size; ++value)
                if (! std::ranges::binary_search(named, value))
                    sum += " " + y(v, value) + " w";
        proof.emit_proof_directive("core id " + to_string(proof.emit_proof_line(sum + " ;")) + " ;");
    }

    // Everything else goes. Each deletion is checked, with a one-line subproof.
    auto del = [&](long line, const vector<long> & hints) {
        proof.begin_deletion_subproof("del id " + to_string(line) + " : : subproof");
        proof.emit_proof_line("rup >= 1 : -1" + ids(hints) + " ;");
        proof.emit_proof_directive("qed : -1 ;");
    };

    // An active or whose child holds has tagged that child, since the tag it has chosen is
    // on a child that puts the vertex somewhere else.
    proof.emit_proof_directive("% each or's children are disjoint");
    for (auto & [i, separator] : _imp->separators) {
        auto & [v, values] = separator;
        for (unsigned j = 0; j < values.size(); ++j) {
            vector<long> hints;
            if (-1 != values[j])
                hints = {vertex_lemma[i][j], y_at_most_one[v]};
            else
                hints = not_lemmas[i][j];
            for (unsigned k = 0; k < values.size(); ++k)
                if (k != j) {
                    hints.push_back(vertex_lemma[i][k]);
                    hints.push_back(tag_child[i][k]);
                }
            hints.push_back(at_least_one_tag[i]);
            del(tag_if[i][j], hints);
        }
    }

    proof.emit_proof_directive("% the root's gate, and its activity's definition");
    del(root_line, {root_active, root_zero});
    del(root_one, {root_active});

    proof.emit_proof_directive("% the lemmas about vertices, top down");
    for (auto & [line, hints] : lemmas.lemmas_and_hints | std::views::reverse)
        del(line, hints);

    proof.emit_proof_directive("% a tag implies its child's gate");
    for (int i = 0; i < size; ++i)
        for (unsigned j = 0; j < tag_child[i].size(); ++j) {
            auto c = nodes[i].children[j];
            del(tag_child[i][j], {one[c][selector_index(c, i, j)], active[c]});
        }

    // Top down, so that each node's children still have theirs: an active cube's literals
    // hold, so its gate does; an active and's children are active; an active or has a tag,
    // whose child is active.
    proof.emit_proof_directive("% activity implies the gate, top down");
    for (int i = size - 1; i >= 0; --i) {
        vector<long> hints;
        switch (nodes[i].kind) {
        case CountingNodeKind::Cube:
            hints = cube_rows[i];
            hints.push_back(backward[i]);
            break;
        case CountingNodeKind::And:
            for (unsigned j = 0; j < nodes[i].children.size(); ++j) {
                auto c = nodes[i].children[j];
                hints.push_back(one[c][selector_index(c, i, j)]);
                hints.push_back(active[c]);
            }
            hints.push_back(backward[i]);
            break;
        case CountingNodeKind::Or:
            hints.push_back(backward[i]);
            for (unsigned j = 0; j < nodes[i].children.size(); ++j) {
                auto c = nodes[i].children[j];
                hints.push_back(active[c]);
                hints.push_back(one[c][selector_index(c, i, j)]);
            }
            hints.push_back(at_least_one_tag[i]);
            break;
        case CountingNodeKind::Distinct: {
            // each sub-gate holds, from the node's rows, and then so does its gate
            auto & sub = _imp->distinct_sub_gates.at(i);
            auto k = nodes[i].variables.size();
            auto negated = proof.begin_deletion_subproof("del id " + to_string(active[i]) + " : : subproof");
            vector<long> lemmas;
            for (unsigned j = 0; j < k; ++j)
                lemmas.push_back(proof.emit_proof_line("rup 1 ~" + a(i) + " 1 " + sub_gates[sub[j]].name + " >= 1 : ~ " +
                    to_string(distinct_in[i][j]) + " " + to_string(sub_gates[sub[j]].backward) + " ;"));
            for (unsigned j = 0; j < distinct_at_most[i].size(); ++j)
                lemmas.push_back(proof.emit_proof_line("pol " + to_string(sub_gates[sub[k + j]].backward) + " " + to_string(distinct_at_most[i][j]) + " + s ;"));
            proof.emit_proof_line("rup >= 1 : " + to_string(negated) + ids(lemmas) + " " + to_string(backward[i]) + " ;");
            proof.emit_proof_directive("qed : -1 ;");
            continue;
        }
        }
        del(active[i], hints);
    }

    // The activity variables are determined by the tags, and the gates by y, so neither
    // needs to be preserved. Activity goes bottom up, since each is defined by its parents'
    // selectors, and the root's is always true; then the gates, top down; then their
    // definitions, with witnesses.
    proof.emit_proof_directive("% activity and gates leave the preserved set");
    for (int i = 0; i < size; ++i)
        if (i != root)
            proof.emit_proof_directive("preserved_rm " + a(i) + " : " + selectors(i) + " ;");
    proof.emit_proof_directive("preserved_rm " + a(root) + " : >= 0 ;");
    for (int i = size - 1; i >= 0; --i)
        proof.emit_proof_directive("preserved_rm " + g(i) + " : " + gate_definition(i, nodes[i], _imp->sub_gate_names(i)).constraint + " ;");
    for (auto & sub : sub_gates | std::views::reverse)
        proof.emit_proof_directive("preserved_rm " + sub.name + " : " + sub.constraint + " ;");
    for (int i = size - 1; i >= 0; --i) {
        proof.emit_proof_directive("del id " + to_string(backward[i]) + " : " + g(i) + " -> 1 ;");
        if (forward[i])
            proof.emit_proof_directive("del id " + to_string(forward[i]) + " : " + g(i) + " -> 0 ;");
    }
    for (auto & sub : sub_gates | std::views::reverse) {
        proof.emit_proof_directive("del id " + to_string(sub.backward) + " : " + sub.name + " -> 1 ;");
        proof.emit_proof_directive("del id " + to_string(sub.forward) + " : " + sub.name + " -> 0 ;");
    }
}
