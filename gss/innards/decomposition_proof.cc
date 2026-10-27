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
using std::tuple;
using std::unique_ptr;
using std::vector;

struct DecompositionProof::Imp
{
    const InputGraph * pattern = nullptr;
    const InputGraph * target = nullptr;
    vector<int> hubs;
    vector<bool> is_hub;
    // the two halves of each copy y_p_t <-> x_p_t, and the exactly-one rows over y
    vector<vector<long>> y_implies_x, x_implies_y;
    vector<long> y_at_most_one;
    shared_ptr<Proof> proof;
    unique_ptr<HomomorphismProofs> hom_proofs;
    int number_of_folds = 0;
    set<int> refuted_leaves;
    set<string> derived_units;

    // each cube written so far, by its assignment
    map<vector<pair<int, int>>, int> cubes_by_assignment;

    // the scaffolding for hubs: off_hubs[v] is v's literal and the line of its definition
    // that says it puts v on no hub, and the rest is by what it says
    map<int, pair<string, long>> off_hubs;
    map<vector<string>, string> conjunctions;
    set<tuple<int, int, vector<int>>> degree_restrictions;
    set<pair<int, int>> hubs_too_small, ruled_out;
    int number_of_guess_folds = 0;
    optional<string> guess_fold;
    set<long> scaffolding;
};

namespace
{
    auto unit(const string & literal) -> string
    {
        return "rup 1 " + literal + " >= 1 ;";
    }
}

DecompositionProof::DecompositionProof(const ProofOptions & options, const InputGraph & pattern, const InputGraph & target,
    const vector<int> & hubs) :
    _imp(make_unique<Imp>())
{
    _imp->pattern = &pattern;
    _imp->target = &target;
    _imp->hubs = hubs;
    _imp->is_hub.resize(target.size(), false);
    for (auto h : hubs)
        _imp->is_hub[h] = true;
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
    HomomorphismParams params, CountingStatementWriter & writer, const optional<string> & assumption) -> SearchedComponent
{
    auto & proof = *_imp->proof;
    SearchedComponent result;

    // The fold: guard_k <-> guard_{k-1} or cube_k, over the cubes found so far. Its
    // definitions stay derived, and go once the leaf is derived.
    optional<string> guard;

    NestedSolve nested;
    nested.assumption = assumption;
    nested.on_solution = [&](const vector<pair<int, int>> & solution) -> string {
        vector<pair<int, int>> assignment;
        vector<int> embedding(component.size());
        for (auto & [p, t] : solution) {
            assignment.emplace_back(vertices.at(p), t);
            embedding.at(p) = t;
        }
        std::sort(assignment.begin(), assignment.end());
        result.embeddings.push_back(std::move(embedding));

        // the writer's observer defines the cube's gate, if it is new
        auto [where, inserted] = _imp->cubes_by_assignment.emplace(assignment, 0);
        if (inserted)
            where->second = writer.add_cube(assignment);
        result.cubes.push_back(where->second);
        auto cube = counting_gate_name(where->second);

        auto fold = "fold" + to_string(++_imp->number_of_folds);
        proof.at_level_zero([&] {
            if (guard) {
                result.fold_lines.push_back(proof.emit_proof_line("red 1 ~" + fold + " 1 " + *guard + " 1 " + cube + " >= 1 : " + fold + " -> 0 ;"));
                result.fold_lines.push_back(proof.emit_proof_line("red 2 " + fold + " 1 ~" + *guard + " 1 ~" + cube + " >= 2 : " + fold + " -> 1 ;"));
            }
            else {
                result.fold_lines.push_back(proof.emit_proof_line("red 1 ~" + fold + " 1 " + cube + " >= 1 : " + fold + " -> 0 ;"));
                result.fold_lines.push_back(proof.emit_proof_line("red 1 " + fold + " 1 ~" + cube + " >= 1 : " + fold + " -> 1 ;"));
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
    result.conclusion = proof.end_nested_solve();
    result.complete = solved.complete && result.conclusion.has_value();
    if (result.cubes.empty() && result.conclusion)
        _imp->scaffolding.insert(*result.conclusion);
    return result;
}

auto DecompositionProof::conclude_component(const SearchedComponent & searched, int leaf, const optional<string> & assumption) -> void
{
    // the last guard is the or of every cube found, all of which are in the leaf
    auto & proof = *_imp->proof;
    proof.emit_proof_line("rup" + (assumption ? " 1 ~" + *assumption : string{}) + " 1 " + counting_gate_name(leaf) + " >= 1 ;");
    auto lines = searched.fold_lines;
    lines.push_back(*searched.conclusion);
    forget(lines);
}

auto DecompositionProof::derive(const string & literal) -> long
{
    long result = 0;
    _imp->proof->at_level_zero([&] { result = _imp->proof->emit_proof_line(unit(literal)); });
    return result;
}

auto DecompositionProof::derive_clause(const vector<string> & literals) -> long
{
    string line = "rup";
    for (auto & l : literals)
        line += " 1 " + l;
    long result = 0;
    _imp->proof->at_level_zero([&] { result = _imp->proof->emit_proof_line(line + " >= 1 ;"); });
    _imp->scaffolding.insert(result);
    return result;
}

auto DecompositionProof::forget(const vector<long> & lines) -> void
{
    if (lines.empty())
        return;
    string ids;
    for (auto & id : lines) {
        ids += " " + to_string(id);
        _imp->scaffolding.erase(id);
    }
    _imp->proof->at_level_zero([&] { _imp->proof->emit_proof_directive("del id" + ids + " ;"); });
}

auto DecompositionProof::forget_scaffolding() -> void
{
    forget(vector<long>(_imp->scaffolding.begin(), _imp->scaffolding.end()));
}

auto DecompositionProof::assignment(int v, int t) const -> string
{
    return "x" + _imp->proof->variable_name(v, t);
}

auto DecompositionProof::off_hubs(int v) -> string
{
    if (auto o = _imp->off_hubs.find(v); o != _imp->off_hubs.end())
        return o->second.first;

    auto & proof = *_imp->proof;
    auto name = "offhub" + to_string(v);
    auto c = to_string(_imp->hubs.size());
    string none = "red " + c + " ~" + name, some = "red 1 " + name;
    for (auto h : _imp->hubs) {
        none += " 1 ~x" + proof.variable_name(v, h);
        some += " 1 x" + proof.variable_name(v, h);
    }
    long none_line = 0;
    proof.at_level_zero([&] {
        none_line = proof.emit_proof_line(none + " >= " + c + " : " + name + " -> 0 ;");
        _imp->scaffolding.insert(none_line);
        _imp->scaffolding.insert(proof.emit_proof_line(some + " >= 1 : " + name + " -> 1 ;"));
    });
    _imp->off_hubs.emplace(v, pair{name, none_line});
    return name;
}

auto DecompositionProof::conjunction(vector<string> literals) -> string
{
    std::sort(literals.begin(), literals.end());
    if (auto c = _imp->conjunctions.find(literals); c != _imp->conjunctions.end())
        return c->second;

    auto & proof = *_imp->proof;
    auto name = "hguess" + to_string(_imp->conjunctions.size() + 1);
    auto k = to_string(literals.size());
    string all = "red " + k + " ~" + name, not_all = "red 1 " + name;
    for (auto & l : literals) {
        all += " 1 " + l;
        not_all += " 1 ~" + l;
    }
    proof.at_level_zero([&] {
        _imp->scaffolding.insert(proof.emit_proof_line(all + " >= " + k + " : " + name + " -> 0 ;"));
        _imp->scaffolding.insert(proof.emit_proof_line(not_all + " >= 1 : " + name + " -> 1 ;"));
    });
    _imp->conjunctions.emplace(literals, name);
    return name;
}

auto DecompositionProof::restrict_by_degree(int v, int t, const vector<int> & w) -> void
{
    if (! _imp->degree_restrictions.emplace(v, t, w).second)
        return;

    // The adjacency rows for v at t, one per neighbour in w, say that each of them goes to
    // a neighbour of t. Injectivity on each non-hub neighbour of t, and each neighbour's
    // being off the hubs, cancel all of those, and what is left says that v is not at t,
    // or some neighbour is on a hub. (A target loop at t lets t count as a neighbour of
    // itself here, but injectivity on t also cancels v being there, so it comes out even.)
    auto & proof = *_imp->proof;
    auto x = [&](int p, int u) { return "x" + proof.variable_name(p, u); };
    string pol = "pol";
    set<long> permitted;
    bool first = true;
    for (auto n : w) {
        auto [label, row] = _imp->hom_proofs->model_adjacency(v, n, t);
        pol += " " + label + (first ? "" : " +");
        first = false;
        permitted.insert(row.begin(), row.end());
    }
    bool some_hub = false;
    for (auto u : permitted) {
        if (_imp->is_hub[u])
            some_hub = true;
        else
            pol += " " + proof.injectivity_label(u) + " +";
    }
    string conclusion = "ia 1 ~" + x(v, t);
    for (auto n : w) {
        if (some_hub)
            pol += " " + to_string(_imp->off_hubs.at(n).second) + " +";
        conclusion += " 1 ~" + _imp->off_hubs.at(n).first;
    }
    proof.at_level_zero([&] {
        auto sum = proof.emit_proof_line(pol + " s ;");
        _imp->scaffolding.insert(proof.emit_proof_line(conclusion + " >= 1 : " + to_string(sum) + " ;"));
        _imp->scaffolding.insert(sum);
    });
}

auto DecompositionProof::rule_out(int v, int t) -> void
{
    if (_imp->ruled_out.emplace(v, t).second)
        _imp->proof->at_level_zero([&] { _imp->scaffolding.insert(_imp->proof->emit_proof_line(unit("~" + assignment(v, t)))); });
}

auto DecompositionProof::hub_too_small(int v, int h) -> void
{
    if (! _imp->hubs_too_small.emplace(v, h).second)
        return;

    // as for restrict_by_degree, but with every neighbour, and every neighbour of h
    auto & proof = *_imp->proof;
    auto & pattern = *_imp->pattern;
    string pol = "pol";
    set<long> permitted;
    bool first = true;
    for (int n = 0; n < pattern.size(); ++n)
        if (n != v && pattern.adjacent(v, n)) {
            auto [label, row] = _imp->hom_proofs->model_adjacency(v, n, h);
            pol += " " + label + (first ? "" : " +");
            first = false;
            permitted.insert(row.begin(), row.end());
        }
    for (auto u : permitted)
        pol += " " + proof.injectivity_label(u) + " +";
    proof.at_level_zero([&] {
        auto sum = proof.emit_proof_line(pol + " s ;");
        _imp->scaffolding.insert(proof.emit_proof_line("ia 1 ~x" + proof.variable_name(v, h) + " >= 1 : " + to_string(sum) + " ;"));
        _imp->scaffolding.insert(sum);
    });
}

auto DecompositionProof::too_many_off_hubs(const vector<int> & vertices) -> void
{
    // each vertex goes somewhere, on a hub only if it is not off the hubs, and each
    // non-hub takes at most one of them
    auto & proof = *_imp->proof;
    string pol = "pol";
    bool first = true;
    for (auto v : vertices) {
        pol += " " + proof.at_least_one_value_label(v) + (first ? "" : " +");
        first = false;
        pol += " " + to_string(_imp->off_hubs.at(v).second) + " +";
    }
    for (int t = 0; t < _imp->target->size(); ++t)
        if (! _imp->is_hub[t])
            pol += " " + proof.injectivity_label(t) + " +";
    proof.at_level_zero([&] { _imp->scaffolding.insert(proof.emit_proof_line(pol + " ;")); });
}

auto DecompositionProof::fold_guess(int gate) -> string
{
    auto & proof = *_imp->proof;
    auto g = counting_gate_name(gate);
    auto fold = "hfold" + to_string(++_imp->number_of_guess_folds);
    proof.at_level_zero([&] {
        if (_imp->guess_fold) {
            _imp->scaffolding.insert(proof.emit_proof_line("red 1 ~" + fold + " 1 " + *_imp->guess_fold + " 1 " + g + " >= 1 : " + fold + " -> 0 ;"));
            _imp->scaffolding.insert(proof.emit_proof_line("red 2 " + fold + " 1 ~" + *_imp->guess_fold + " 1 ~" + g + " >= 2 : " + fold + " -> 1 ;"));
        }
        else {
            _imp->scaffolding.insert(proof.emit_proof_line("red 1 ~" + fold + " 1 " + g + " >= 1 : " + fold + " -> 0 ;"));
            _imp->scaffolding.insert(proof.emit_proof_line("red 1 " + fold + " 1 ~" + g + " >= 1 : " + fold + " -> 1 ;"));
        }
    });
    _imp->guess_fold = fold;
    return fold;
}

auto DecompositionProof::refute_glued_leaf(int leaf, const vector<int> & cubes) -> void
{
    if (! _imp->refuted_leaves.insert(leaf).second)
        return;
    for (auto & c : cubes)
        if (_imp->derived_units.insert("~" + counting_gate_name(c)).second)
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

    // Every adjacency row goes by unit propagation from the root: under the negated row,
    // every guess fails, either because it puts a vertex somewhere the row rules out, or
    // because the row falsifies every cube of some component's leaf.
    proof.emit_proof_directive("% the root implies adjacency");
    for (auto & label : _imp->hom_proofs->model_adjacency_labels())
        proof.emit_proof_directive("del id " + label + " ;");

    // Each injectivity row: for each pair of pattern vertices, not both at t. Every guess
    // puts at most one vertex on each hub and keeps the rest off them, and within a
    // component, every cube separates them. That leaves two vertices in different
    // components of what some guess leaves off the hubs, which is the collision lemma: if
    // they collide, then whichever pair of cubes holds completes a cube of some overlap ρ,
    // so ρ's product holds, and so either ρ does, or one of ρ's own overlaps does. Every
    // one of those is an overlap of the guess, which the guess says do not hold, so unit
    // propagation finds the contradiction once the pair of cubes is fixed. Then recover
    // the at most one from the pairs.
    struct Where
    {
        vector<int> component_of, position_of;
    };
    vector<Where> where;
    for (auto & guess : structure.guesses) {
        Where w{vector<int>(pattern.size(), -1), vector<int>(pattern.size(), -1)};
        for (unsigned c = 0; c < guess.components.size(); ++c)
            for (unsigned i = 0; i < guess.components[c].vertices.size(); ++i) {
                w.component_of[guess.components[c].vertices[i]] = c;
                w.position_of[guess.components[c].vertices[i]] = i;
            }
        where.push_back(std::move(w));
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

        auto pair_label = [&](int u, int v) { return "@p" + to_string(t) + "_" + to_string(u) + "_" + to_string(v); };
        for (int u = 0; u < pattern.size(); ++u)
            for (int v = u + 1; v < pattern.size(); ++v) {
                for (unsigned gi = 0; gi < structure.guesses.size(); ++gi) {
                    auto & guess = structure.guesses[gi];
                    auto & w = where[gi];
                    if (w.component_of[u] == -1 || w.component_of[v] == -1 || w.component_of[u] == w.component_of[v])
                        continue;
                    auto not_guess = "~" + g(guess.gate);
                    auto & k1 = guess.components[w.component_of[u]];
                    auto & k2 = guess.components[w.component_of[v]];
                    for (unsigned i = 0; i < k1.cubes.size(); ++i) {
                        if (k1.embeddings[i][w.position_of[u]] != t)
                            continue;
                        for (unsigned j = 0; j < k2.cubes.size(); ++j) {
                            if (k2.embeddings[j][w.position_of[v]] != t)
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
                            for (auto & k : guess.components)
                                for (auto z : k.vertices)
                                    if (! glued.contains(z))
                                        partition.push_back({z});
                            std::sort(partition.begin(), partition.end());
                            if (! guess.overlaps_by_partition.contains(partition))
                                throw ProofError{"two cubes collide in an overlap that has no node"};
                            proof.emit_proof_line(join({not_guess, "~" + g(k1.cubes[i]), "~" + g(k2.cubes[j])}));
                        }
                        proof.emit_proof_line(join({not_guess, "~" + g(k1.cubes[i]), "~" + y(v, t)}));
                    }
                    proof.emit_proof_line(join({not_guess, "~" + y(u, t), "~" + y(v, t)}));
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
