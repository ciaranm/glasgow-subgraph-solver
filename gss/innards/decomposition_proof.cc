#include <gss/innards/decomposition_proof.hh>
#include <gss/innards/homomorphism_proofs.hh>
#include <gss/innards/proof.hh>

#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <utility>

using namespace gss;
using namespace gss::innards;

using std::make_shared;
using std::make_unique;
using std::map;
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
    const InputGraph * pattern = nullptr;
    const InputGraph * target = nullptr;
    // the two halves of each copy y_p_t <-> x_p_t, and the exactly-one rows over y
    vector<vector<long>> y_implies_x, x_implies_y;
    vector<long> y_at_most_one;
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
    _imp->pattern = &pattern;
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
    auto & y_implies_x = _imp->y_implies_x;
    y_implies_x.resize(pattern.size());
    _imp->x_implies_y.resize(pattern.size());
    for (int p = 0; p < pattern.size(); ++p)
        for (int t = 0; t < target.size(); ++t) {
            auto y = counting_variable_name(p, t), x = "x" + proof.variable_name(p, t);
            auto forward = proof.emit_proof_line("red 1 ~" + y + " 1 " + x + " >= 1 : " + y + " -> 0 ;");
            auto backward = proof.emit_proof_line("red 1 " + y + " 1 ~" + x + " >= 1 : " + y + " -> 1 ;");
            proof.emit_proof_directive("core id " + to_string(forward) + " " + to_string(backward) + " ;");
            proof.emit_proof_directive("preserved_add " + y + " : 1 " + x + " >= 1 ;");
            y_implies_x[p].push_back(forward);
            _imp->x_implies_y[p].push_back(backward);
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
        _imp->y_at_most_one.push_back(am);
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

auto DecompositionProof::derive(const string & literal) -> long
{
    long result = 0;
    _imp->proof->at_level_zero([&] { result = _imp->proof->emit_proof_line(unit(literal)); });
    return result;
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

auto DecompositionProof::finish(const Structure & structure) -> void
{
    auto & proof = *_imp->proof;
    auto & pattern = *_imp->pattern;
    auto & target = *_imp->target;
    auto g = [](int i) { return counting_gate_name(i); };
    auto x = [&](int p, int t) { return "x" + proof.variable_name(p, t); };
    auto y = [](int p, int t) { return counting_variable_name(p, t); };
    auto join = [](const vector<string> & literals) {
        string result = "rup";
        for (auto & l : literals)
            result += " 1 " + l;
        return result + " >= 1 ;";
    };

    proof.emit_proof_directive("core id " + to_string(structure.root_line) + " ;");

    // Every adjacency row goes by unit propagation from the root: it reaches the leaf of
    // the component, every cube of which the negated row falsifies.
    proof.emit_proof_directive("% the root implies adjacency");
    for (auto & label : _imp->hom_proofs->model_adjacency_labels())
        proof.emit_proof_directive("del id " + label + " ;");

    // Each injectivity row: for each pair of pattern vertices, not both at t. Within a
    // component, every cube separates them. Across components, the collision lemma: if
    // they collide, then whichever pair of cubes holds completes a cube of some overlap ρ,
    // so ρ's product holds, and so either ρ does, or one of ρ's own overlaps does. Every one
    // of those is an overlap of the root, which the root says do not hold, so unit
    // propagation finds the contradiction once the pair of cubes is fixed. Then recover
    // the at most one from the pairs.
    vector<int> component_of(pattern.size());
    vector<int> position_of(pattern.size());
    for (unsigned c = 0; c < structure.components.size(); ++c)
        for (unsigned i = 0; i < structure.components[c].vertices.size(); ++i) {
            component_of[structure.components[c].vertices[i]] = c;
            position_of[structure.components[c].vertices[i]] = i;
        }

    proof.emit_proof_directive("% the root implies injectivity");
    for (int t = 0; t < target.size(); ++t) {
        auto label = proof.injectivity_label(t);
        if (pattern.size() < 2) {
            proof.emit_proof_directive("del id " + label + " ;");
            continue;
        }

        proof.emit_proof_directive("del id " + label + " : : subproof");
        auto negated = "@neg" + to_string(t);
        proof.emit_proof_line(negated + " pol -1 ;");
        auto b = structure.overlaps ? vector<string>{g(*structure.overlaps)} : vector<string>{};
        auto with_b = [&](vector<string> literals) {
            literals.insert(literals.end(), b.begin(), b.end());
            return literals;
        };

        auto pair_label = [&](int u, int v) { return "@p" + to_string(t) + "_" + to_string(u) + "_" + to_string(v); };
        for (int u = 0; u < pattern.size(); ++u)
            for (int v = u + 1; v < pattern.size(); ++v) {
                if (component_of[u] != component_of[v]) {
                    auto & k1 = structure.components[component_of[u]];
                    auto & k2 = structure.components[component_of[v]];
                    for (unsigned i = 0; i < k1.cubes.size(); ++i) {
                        if (k1.embeddings[i][position_of[u]] != t)
                            continue;
                        for (unsigned j = 0; j < k2.cubes.size(); ++j) {
                            if (k2.embeddings[j][position_of[v]] != t)
                                continue;
                            // the overlap these two cubes complete, gluing every pair
                            // of their vertices that share an image
                            map<int, int> other_end;
                            for (unsigned a = 0; a < k1.vertices.size(); ++a)
                                for (unsigned c = 0; c < k2.vertices.size(); ++c)
                                    if (k1.embeddings[i][a] == k2.embeddings[j][c])
                                        other_end[k1.vertices[a]] = k2.vertices[c];
                            vector<vector<int>> partition;
                            set<int> glued;
                            for (auto & [w1, w2] : other_end) {
                                partition.push_back({std::min(w1, w2), std::max(w1, w2)});
                                glued.insert(w1);
                                glued.insert(w2);
                            }
                            for (int w = 0; w < pattern.size(); ++w)
                                if (! glued.contains(w))
                                    partition.push_back({w});
                            std::sort(partition.begin(), partition.end());
                            if (! structure.overlaps_by_partition.contains(partition))
                                throw ProofError{"two cubes collide in an overlap that has no node"};
                            proof.emit_proof_line(join(with_b({"~" + g(structure.product), "~" + g(k1.cubes[i]), "~" + g(k2.cubes[j])})));
                        }
                        proof.emit_proof_line(join(with_b({"~" + g(structure.product), "~" + g(k1.cubes[i]), "~" + y(v, t)})));
                    }
                    proof.emit_proof_line(join(with_b({"~" + g(structure.product), "~" + y(u, t), "~" + y(v, t)})));
                }
                proof.emit_proof_line(pair_label(u, v) + " " + join({"~" + x(u, t), "~" + x(v, t)}));
            }

        // at most one, from the pairs, a vertex at a time
        proof.emit_proof_line("@amo" + to_string(t) + "_2 pol " + pair_label(0, 1) + " ;");
        for (int k = 2; k < pattern.size(); ++k) {
            string line = "@amo" + to_string(t) + "_" + to_string(k + 1) + " pol @amo" + to_string(t) + "_" + to_string(k) + " " + to_string(k - 1) + " *";
            for (int i = 0; i < k; ++i)
                line += " " + pair_label(i, k) + " +";
            proof.emit_proof_line(line + " " + to_string(k) + " d ;");
        }
        proof.emit_proof_line("pol @amo" + to_string(t) + "_" + to_string(pattern.size()) + " " + negated + " + ;");
        proof.emit_proof_directive("qed : -1 ;");
    }

    // The exactly-one rows over x follow from the ones over y and the copies.
    proof.emit_proof_directive("% the exactly-one rows over x follow from the ones over y");
    for (int p = 0; p < pattern.size(); ++p) {
        proof.emit_proof_directive("del id " + proof.at_least_one_value_label(p) + " ;");
        string at_most = "pol -1 " + to_string(_imp->y_at_most_one[p]) + " +";
        for (int t = 0; t < target.size(); ++t)
            at_most += " " + to_string(_imp->x_implies_y[p][t]) + " +";
        proof.emit_proof_directive("del id " + proof.at_most_one_value_label(p) + " : : subproof");
        proof.emit_proof_line(at_most + " ;");
        proof.emit_proof_directive("qed : -1 ;");
    }

    // Now x is only in the copies: it leaves the preserved set, and then the copies go,
    // with x -> y as the witness.
    proof.emit_proof_directive("% x leaves the preserved set");
    for (int p = 0; p < pattern.size(); ++p)
        for (int t = 0; t < target.size(); ++t)
            proof.emit_proof_directive("preserved_rm " + x(p, t) + " : 1 " + y(p, t) + " >= 1 ;");
    for (int p = 0; p < pattern.size(); ++p)
        for (int t = 0; t < target.size(); ++t) {
            proof.emit_proof_directive("del id " + to_string(_imp->y_implies_x[p][t]) + " : " + x(p, t) + " -> " + y(p, t) + " ;");
            proof.emit_proof_directive("del id " + to_string(_imp->x_implies_y[p][t]) + " : " + x(p, t) + " -> " + y(p, t) + " ;");
        }

    proof.emit_proof_directive("output EQUIENUMERABLE FILE ;");
    proof.emit_proof_directive("conclusion NONE ;");
    proof.emit_proof_directive("end pseudo-Boolean proof ;");
}
