#include <gss/counting_statement.hh>
#include <gss/innards/homomorphism_proofs.hh>
#include <gss/innards/proof.hh>
#include <gss/innards/split_proof.hh>

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <optional>

using namespace gss;
using namespace gss::innards;

using std::make_unique;
using std::optional;
using std::pair;
using std::shared_ptr;
using std::string;
using std::to_string;
using std::unique_ptr;
using std::vector;

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
};

SplitProof::SplitProof(const shared_ptr<Proof> & proof, HomomorphismProofs & hom_proofs, int pattern_size, int target_size,
    const string & statement_file) :
    _imp(make_unique<Imp>())
{
    _imp->proof = proof;
    _imp->hom_proofs = &hom_proofs;
    _imp->pattern_size = pattern_size;
    _imp->target_size = target_size;
    _imp->statement_file = statement_file;
    _imp->statement_stream.open(statement_file);
    if (! _imp->statement_stream)
        throw ProofError{"Error writing counting statement to '" + statement_file + "'"};
    _imp->writer.emplace(_imp->statement_stream, vector<int>(pattern_size, target_size));
    // Each node's gate is defined, both ways, into the core, while it is still fresh: VeriPB
    // checks a definition against the derived constraints too, so the rows that say what a
    // gate implies cannot wait until the search has derived things about it.
    _imp->writer->set_observer([this](int i, const CountingNode & node) {
        auto & proof = *_imp->proof;
        proof.at_level_zero([&] {
            auto rows = counting_node_definition(i, node);
            auto g = counting_gate_name(i);
            string ids;
            for (unsigned r = 0; r < rows.size(); ++r) {
                // every row but the last says the gate implies its definition, so holds
                // with the gate false; the last says the definition implies the gate
                bool last = r + 1 == rows.size();
                auto row = rows[r].substr(0, rows[r].size() - 2);
                ids += " " + to_string(proof.emit_proof_line("red " + row + " : " + g + " -> " + (last ? "1" : "0") + " ;"));
            }
            proof.emit_proof_directive("core id" + ids + " ;");
        });
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

auto SplitProof::begin() -> void
{
    auto & proof = *_imp->proof;

    // Each y_p_t is a copy of x_p_t, and joins the preserved set. (The x side leaves it in
    // finish, once nothing but the copies mentions x.)
    proof.emit_proof_directive("% the statement's variables are copies of the encoding's");
    _imp->y_implies_x.resize(_imp->pattern_size);
    _imp->x_implies_y.resize(_imp->pattern_size);
    for (int p = 0; p < _imp->pattern_size; ++p)
        for (int t = 0; t < _imp->target_size; ++t) {
            auto y = counting_variable_name(p, t), x = "x" + proof.variable_name(p, t);
            auto forward = proof.emit_proof_line("red 1 ~" + y + " 1 " + x + " >= 1 : " + y + " -> 0 ;");
            auto backward = proof.emit_proof_line("red 1 " + y + " 1 ~" + x + " >= 1 : " + y + " -> 1 ;");
            proof.emit_proof_directive("core id " + to_string(forward) + " " + to_string(backward) + " ;");
            proof.emit_proof_directive("preserved_add " + y + " : 1 " + x + " >= 1 ;");
            _imp->y_implies_x[p].push_back(forward);
            _imp->x_implies_y[p].push_back(backward);
        }

    // The exactly-one rows over y, which the counting formula has: at least one by unit
    // propagation through the copies, and at most one by adding the y -> x halves of the
    // copies to the encoding's at most one over x.
    for (int p = 0; p < _imp->pattern_size; ++p) {
        string at_least = "rup", at_most = "pol " + proof.at_most_one_value_label(p);
        for (int t = 0; t < _imp->target_size; ++t) {
            at_least += " 1 " + counting_variable_name(p, t);
            at_most += " " + to_string(_imp->y_implies_x[p][t]) + " +";
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

auto SplitProof::conjunction(const vector<int> & children) -> int
{
    return 1 == children.size() ? children[0] : _imp->writer->add_and(children);
}

auto SplitProof::disjunction(const vector<int> & children) -> int
{
    return 1 == children.size() ? children[0] : _imp->writer->add_or(children);
}

auto SplitProof::gate(int node) const -> string
{
    return counting_gate_name(node);
}

auto SplitProof::statement_variable(int v, int t) const -> string
{
    return counting_variable_name(v, t);
}

auto SplitProof::assumption(const vector<pair<int, int>> & excluded) -> string
{
    auto & proof = *_imp->proof;
    auto name = "sa" + to_string(++_imp->number_of_assumptions);
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
    auto & proof = *_imp->proof;
    proof.at_level_zero([&] {
        for (auto & [node, v, t] : lemmas)
            proof.emit_proof_line("rup 1 ~" + gate(node) + " 1 " + statement_variable(v, t) + " >= 1 ;");
    });
}

auto SplitProof::forget(const vector<long> & lines) -> void
{
    if (lines.empty())
        return;
    string ids;
    for (auto & id : lines)
        ids += " " + to_string(id);
    _imp->proof->at_level_zero([&] { _imp->proof->emit_proof_directive("del id" + ids + " ;"); });
}

auto SplitProof::finish(int root) -> void
{
    auto & proof = *_imp->proof;
    auto x = [&](int p, int t) { return "x" + proof.variable_name(p, t); };
    auto y = [](int p, int t) { return counting_variable_name(p, t); };
    auto n = _imp->pattern_size;

    _imp->writer->finish(root);
    _imp->statement_stream.close();
    _imp->finished = true;

    proof.emit_proof_directive("setlvl 0;");
    auto root_line = proof.emit_proof_line("rup 1 " + gate(root) + " >= 1 ;");
    proof.emit_proof_directive("core id " + to_string(root_line) + " ;");

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

        proof.emit_proof_directive("del id " + label + " : : subproof");
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
        proof.emit_proof_directive("del id " + proof.at_most_one_value_label(p) + " : : subproof");
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

    proof.emit_proof_directive("output EQUIENUMERABLE FILE ;");
    proof.emit_proof_directive("conclusion NONE ;");
    proof.emit_proof_directive("end pseudo-Boolean proof ;");
}
