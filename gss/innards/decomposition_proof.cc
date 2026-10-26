#include <gss/innards/decomposition_proof.hh>
#include <gss/innards/homomorphism_proofs.hh>
#include <gss/innards/proof.hh>

#include <algorithm>
#include <set>
#include <string>
#include <utility>

using namespace gss;
using namespace gss::innards;

using std::make_shared;
using std::make_unique;
using std::optional;
using std::pair;
using std::set;
using std::shared_ptr;
using std::string;
using std::to_string;
using std::unique_ptr;
using std::vector;

struct DecompositionProof::Imp
{
    const InputGraph * target = nullptr;
    shared_ptr<Proof> proof;
    unique_ptr<HomomorphismProofs> hom_proofs;
    int number_of_folds = 0;
    set<int> refuted_leaves;
};

namespace
{
    auto unit(const string & literal) -> string
    {
        return "rup 1 " + literal + " >= 1 ;";
    }
}

DecompositionProof::DecompositionProof(const ProofOptions & options, const InputGraph & pattern, const InputGraph & target) :
    _imp(make_unique<Imp>())
{
    _imp->target = &target;
    _imp->proof = make_shared<Proof>(options);
    _imp->hom_proofs = make_unique<HomomorphismProofs>(_imp->proof, pattern, target);
    auto & proof = *_imp->proof;

    // the solver's own encoding of the whole instance, with its preserved set
    HomomorphismParams model_params;
    model_params.count_solutions = true;
    _imp->hom_proofs->emit_model(pattern, target, model_params);

    // Each y_p_t is a copy of x_p_t, and joins the preserved set. (The x side leaves it in
    // phase 5, once nothing but the copies mentions x.)
    proof.emit_proof_directive("% the statement's variables are copies of the encoding's");
    vector<vector<long>> y_implies_x(pattern.size());
    for (int p = 0; p < pattern.size(); ++p)
        for (int t = 0; t < target.size(); ++t) {
            auto y = counting_variable_name(p, t), x = "x" + proof.variable_name(p, t);
            auto forward = proof.emit_proof_line("red 1 ~" + y + " 1 " + x + " >= 1 : " + y + " -> 0 ;");
            auto backward = proof.emit_proof_line("red 1 " + y + " 1 ~" + x + " >= 1 : " + y + " -> 1 ;");
            proof.emit_proof_directive("core id " + to_string(forward) + " " + to_string(backward) + " ;");
            proof.emit_proof_directive("preserved_add " + y + " : 1 " + x + " >= 1 ;");
            y_implies_x[p].push_back(forward);
        }

    // The exactly-one rows over y, which the counting formula has: at least one by unit
    // propagation through the copies, and at most one by adding the y -> x halves of the
    // copies to the encoding's at most one over x.
    for (int p = 0; p < pattern.size(); ++p) {
        string at_least = "rup", at_most = "pol " + proof.at_most_one_value_label(p);
        for (int t = 0; t < target.size(); ++t) {
            at_least += " 1 " + counting_variable_name(p, t);
            at_most += " " + to_string(y_implies_x[p][t]) + " +";
        }
        auto al = proof.emit_proof_line(at_least + " >= 1 ;");
        auto am = proof.emit_proof_line(at_most + " ;");
        proof.emit_proof_directive("core id " + to_string(al) + " " + to_string(am) + " ;");
    }
}

DecompositionProof::~DecompositionProof() = default;

auto DecompositionProof::define_gate(int i, const CountingNode & node) -> void
{
    auto & proof = *_imp->proof;
    proof.at_level_zero([&] {
        auto rows = counting_node_definition(i, node);
        auto g = counting_gate_name(i);
        string ids;
        for (unsigned r = 0; r < rows.size(); ++r) {
            // every row but the last says the gate implies its definition, so holds with
            // the gate false; the last says the definition implies the gate
            bool last = r + 1 == rows.size();
            auto row = rows[r].substr(0, rows[r].size() - 2);
            ids += " " + to_string(proof.emit_proof_line("red " + row + " : " + g + " -> " + (last ? "1" : "0") + " ;"));
        }
        proof.emit_proof_directive("core id" + ids + " ;");
    });
}

auto DecompositionProof::search_component(const InputGraph & component, const vector<int> & vertices,
    HomomorphismParams params, CountingStatementWriter & writer) -> SearchedComponent
{
    auto & proof = *_imp->proof;
    SearchedComponent result;

    // The fold: guard_k <-> guard_{k-1} or cube_k, over the cubes found so far. Its
    // definitions stay derived, and go once the leaf is derived.
    vector<long> fold_lines;
    optional<string> guard;

    NestedSolve nested;
    nested.on_solution = [&](const vector<pair<int, int>> & solution) -> string {
        vector<pair<int, int>> assignment;
        vector<int> embedding(component.size());
        for (auto & [p, t] : solution) {
            assignment.emplace_back(vertices.at(p), t);
            embedding.at(p) = t;
        }
        std::sort(assignment.begin(), assignment.end());
        result.embeddings.push_back(std::move(embedding));

        // the writer's observer defines the cube's gate
        auto id = writer.add_cube(assignment);
        result.cubes.push_back(id);
        auto cube = counting_gate_name(id);

        auto fold = "fold" + to_string(++_imp->number_of_folds);
        proof.at_level_zero([&] {
            if (guard) {
                fold_lines.push_back(proof.emit_proof_line("red 1 ~" + fold + " 1 " + *guard + " 1 " + cube + " >= 1 : " + fold + " -> 0 ;"));
                fold_lines.push_back(proof.emit_proof_line("red 2 " + fold + " 1 ~" + *guard + " 1 ~" + cube + " >= 2 : " + fold + " -> 1 ;"));
            }
            else {
                fold_lines.push_back(proof.emit_proof_line("red 1 ~" + fold + " 1 " + cube + " >= 1 : " + fold + " -> 0 ;"));
                fold_lines.push_back(proof.emit_proof_line("red 1 " + fold + " 1 ~" + cube + " >= 1 : " + fold + " -> 1 ;"));
            }
        });
        guard = fold;
        return fold;
    };

    params.proof_options = ProofOptions{};
    params.nested_proof = _imp->proof;
    params.enumerate_callback = nullptr;

    proof.emit_proof_directive("% searching a component of the pattern");
    proof.begin_nested_solve(std::move(nested));
    auto solved = solve_homomorphism_problem(component, *_imp->target, params);
    auto conclusion = proof.end_nested_solve();
    result.complete = solved.complete && conclusion.has_value();
    if (! result.complete || result.cubes.empty())
        return result;

    // the last guard is the or of every cube, and so is the leaf
    result.leaf = writer.add_or(result.cubes);
    proof.emit_proof_line(unit(counting_gate_name(*result.leaf)));
    fold_lines.push_back(*conclusion);
    string ids;
    for (auto & id : fold_lines)
        ids += " " + to_string(id);
    proof.emit_proof_directive("del id" + ids + " ;");
    return result;
}

auto DecompositionProof::derive(const string & literal) -> void
{
    _imp->proof->at_level_zero([&] { _imp->proof->emit_proof_line(unit(literal)); });
}

auto DecompositionProof::refute_glued_leaf(int leaf, const vector<int> & cubes) -> void
{
    if (! _imp->refuted_leaves.insert(leaf).second)
        return;
    for (auto & c : cubes)
        derive("~" + counting_gate_name(c));
    derive("~" + counting_gate_name(leaf));
}

auto DecompositionProof::finish_pattern_bigger_than_target() -> void
{
    _imp->hom_proofs->failure_due_to_pattern_bigger_than_target();
    _imp->proof->finish_enumeration_proof(0, true);
}

auto DecompositionProof::finish_no_solutions() -> void
{
    _imp->proof->finish_enumeration_proof(0, true);
}

auto DecompositionProof::finish() -> void
{
    _imp->proof->finish_unknown_proof();
}
