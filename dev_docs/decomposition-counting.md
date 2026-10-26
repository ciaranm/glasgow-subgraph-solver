# Certified counting by decomposition

A design plan for counting subgraph embeddings by decomposition, following Section 3.6 of
Jessica Ryan's PhD thesis, and for certifying those counts with VeriPB without enumerating
the solutions. This is forward-looking: it describes where the code is going, not where it
is. The first piece, the c = 0 case without proofs, exists as a prototype in
`gss/sip_decomposer.cc`; everything else here is a plan. For the current architecture see
[architecture.md](architecture.md); for proof logging see
[proof-logging.md](proof-logging.md).

Throughout, "count" means the number of **labelled** embeddings, which is what the
certified SIP encoding and the solver both count. The thesis goes on to divide by
#Aut(H) to count unlabelled copies; that step is deliberately out of scope.

## References

- J. L. Ryan, *Parameterised algorithms for counting subgraphs, matchings, and
  monochromatic partitions*, PhD thesis, University of Glasgow, 2023,
  <https://theses.gla.ac.uk/83568/>. Section 3.6 is the algorithm; Lemma 3.8 is the
  inclusion–exclusion recursion.
- C. McCreesh, J. Nordström, A. Oertel and Y. K. Tan, *Proof Logging for Projected
  Enumeration (and Counting?) Problems in VeriPB*, CP 2026,
  <https://ciaranm.github.io/papers/cp2026-enumeration.pdf>. Section 4 is the
  `EQUIENUMERABLE` output and the preserved-set rules.

## Why

### A single search is not FPT here

Ryan's algorithm counts embeddings of a pattern H (m vertices) into a host G (n vertices)
that has at most c vertices of degree above Δ, in time O(c^c n² m^(2m+c) Δ^m). In outline:

1. Guess which pattern vertices go to the high-degree vertices, and where: at most
   m^c · c^c guesses.
2. Split what is left of the pattern into connected components. Each component, with
   lists of allowed targets imposed by the guess, embeds into the bounded-degree part,
   and a connected-ordering search counts those embeddings in O(n Δ^(m−1)) nodes
   (Lemmas 3.1 and 3.2).
3. Multiply the per-component counts, then correct for embeddings in which components
   overlap, by inclusion–exclusion over the ways of gluing components together (Lemmas
   3.3 to 3.10).

The FPT bound comes from step 3, not from the vertex ordering. A smallest-domain-first
heuristic dominates the ordering property of Lemma 3.1 at every node where that property
applies, so for c = 0 and a connected pattern a single search already has the folklore
bound. But a single backtracking search *multiplies* across components where the thesis
*adds*, so with ℓ components it costs about n^ℓ Δ^m · f(m), which is XP. Two examples:

- **Counting**: G = K_{1,n−1} (so c = 1 and Δ = 1) and H = K_{1,k}. There are Θ(n^k)
  embeddings, and an enumerating solver needs one leaf per solution.
- **Deciding, unsatisfiable**: G is a hub with n/2 pendant leaves plus a disjoint
  bipartite grid on n vertices, and H is K_{1,k} plus a disjoint C_5. After the hub is
  placed, the leaves' domains (about n/2) are smaller than the C_5's (about n), so all k
  leaves are branched on first and the C_5 fails again under each of (n/2)^k leaf
  assignments. Restart nogoods are decision prefixes, so they do not carry across.

### The recursion we use

Let the components (after the guess) be C_1, …, C_ℓ. A partition σ of their vertices is
*valid* if no block contains two vertices of the same component. Valid partitions are
exactly the ways that embeddings of the components can overlap, since each embedding is
injective on its own component. H/σ glues each block into one vertex; it has no new
loops, because adjacent vertices share a component.

Write none(σ) for the number of injective embeddings of H/σ, and product(σ) for the
product, over the connected components of H/σ, of each one's injective embedding count.
A tuple counted by product(σ) is an injective embedding of H/τ for exactly one τ that
coarsens σ by gluing blocks from *distinct* components of H/σ. Hence

    none(σ) = product(σ) − Σ_{τ ∈ Up(σ), τ ≠ σ} none(τ)

which is Lemma 3.8, and the answer is none of the discrete partition. Every term is
non-negative, so a zero product prunes everything above it. Rearranged, this is a
*positive* identity, which is what certification relies on:

    product(σ) = Σ_{τ ∈ Up(σ)} none(τ)

For example, two disjoint edges into the path P4: 36 = 6 · 6 = 8 (disjoint) + 4 · 4 (four
ways to glue one vertex, each giving a P3 with 4 embeddings) + 2 · 6 (two ways to glue both
vertices, each giving one edge with 6 embeddings).

### Is it worth it?

Only on instances where the n^ℓ term bites: patterns that are disconnected, or that fall
apart once a few vertices are placed on hubs, in hosts with a heavy degree drop-off. The
benchmark sets mostly contain instances we can already solve or nearly solve, which is a
selection bias against this. Deciding when to decompose should eventually be a runtime
choice.

## Where the prototype is

`solve_sip_by_decomposition` (`--decomposition`) does the c = 0 case when counting
injective, non-induced, undirected embeddings. Vertex labels, edge labels and loops are
supported. In brief:

- Isolated pattern vertices use a closed form: a falling factorial for each vertex label.
- Everything else goes through the recursion above. Gluing vertices with different labels,
  or gluing two edges with different labels onto each other, gives zero.
- Each glued component is counted by one call to `solve_homomorphism_problem`, cached up
  to isomorphism: colour refinement, then brute force within colour classes, falling back
  to the exact graph above 5040 orderings.
- none() is memoised on the quotient's isomorphism class, which is the sorted multiset of
  its components' keys.
- It reports `decomposition_components`, `decomposition_solver_calls`,
  `decomposition_quotients` and `decomposition_coarsenings_considered`.
- It refuses proof logging, and produces no mapping when counting.

`gss/sip_decomposer_test.cc` compares it with the plain solver on 200 random disconnected
patterns across loops and labels.

Counting on a 40×40 grid, where the plain solver runs out of time on all of these:

| Pattern | Decomposition |
|---|---|
| 3 disjoint edges | 2.4×10¹¹ in 64 ms, 5 solver calls |
| 3 × P3 | 6.0×10¹² in 3.3 s, 51 calls |
| 6 disjoint edges | 5.7×10²² in 12.7 s, 52 calls, 1.8M coarsenings |
| 4 × P4 | times out at 120 s: 812 distinct glued trees, each counted by enumeration |

Known limits: the coarsening enumeration does not use the symmetry between identical
components, which dominates the six-edge case; large glued components are still counted
by enumeration; and there is no c > 0 case yet.

## The counting language

The aim is to move as much of the argument as possible into VeriPB rather than into an
external combiner. The solver states its argument in a small language whose statements
are efficiently countable, and VeriPB checks that the problem is equienumerable with a PB
translation of that statement.

### Nodes

Every node is a **reified gate**: an extension variable r_n defined in both directions
over its children. Gate variables are therefore functionally determined by the original
variables, which is what makes negation, and hence subtraction, cheap.

Variables are **multi-valued**: each pattern vertex u is the exactly-one group x_{u,·}. A
literal fixes one value and implicitly sets the rest of the group false. Without this,
smoothing over Boolean variables would cost a factor of n.

Each node has a *scope*, a set of multi-valued variables, and a count: the number of
assignments to its scope, respecting the exactly-one groups, under which r_n is true.

| Node | Definition | Side condition | Count |
|---|---|---|---|
| cube | r ↔ ∧ of literals, one value for each variable in scope | — | 1 |
| ⊗ | r ↔ a ∧ b | scopes disjoint (syntactic) | #a · #b |
| ⊕ (n-ary) | r ↔ a ∨ b ∨ … | equal scopes; children pairwise disjoint | #a + #b + … |
| ⊖ | r ↔ a ∧ ¬b | equal scopes; b → a | #a − #b |

⊗ and ⊕ are the product and sum nodes of the CPOG framework; ⊖ is the only addition. Any
inclusion–exclusion identity fits these nodes, and constant coefficients, if ever needed,
are an exactly-one over that many fresh variables.

### Translation

A statement translates to a linear-size PB formula: the gate definitions, the exactly-one
groups, and the unit constraint r_root. Because every gate is determined by x, the
formula's count projected onto x is the root's count. For example, ⊖ is

```
1 ~r 1 a >= 1 ;        % r → a
1 ~r 1 ~b >= 1 ;       % r → ¬b
1 r 1 ~a 1 b >= 1 ;    % a ∧ ¬b → r
```

### Why there is no selector-style union

An alternative ⊕ uses selector variables: f_i → a_i with an exactly-one over the f_i. It
needs no disjointness condition, because the selector distinguishes the copies. But
selectors only ever *add* solutions, and the selector is not determined by the scope, so
such a node cannot be negated, and so cannot appear under ⊖. Subtraction has to be
A ∧ ¬B, which is only efficiently encodable when B's auxiliary variables are
functionally determined. So the language uses gate ⊕ throughout and pays for it with a
disjointness side condition.

### Side conditions affect the count

The side conditions say that the gate structure is valid over the definitions *alone*,
not just on the instance's solutions. Equienumerability does not check them. A statement
whose ⊕ children overlap still has a gate formula with the right solution set, so the
main proof still verifies, but the evaluator double-counts the overlap. They must be
checked before any count is reported.

- Most are syntactic. A leaf is a decision structure over cubes (it is the solver's own
  enumeration tree), so sibling cubes conflict on a variable. Leaves are almost all of the
  nodes.
- The rest are about f(m) · m^c · c^c instances: disjointness between hub guesses,
  containment for each Exact_σ, and disjointness of the exact-pattern classes.
- Each of those is emitted by the translator as a small PB instance over the statement's
  own definitions, such as defs ∧ r_b ∧ ¬r_a, which needs an ordinary UNSAT proof. Because
  the instance contains only definitions, it cannot accidentally depend on the problem's
  constraints, and no VeriPB extension is needed.
- Later, these could move inside the main proof with a narrowly scoped VeriPB feature: a
  lemma that may use only the output formula's definitions.

## Compiling Ryan's algorithm

- **Root** = ⊕ over hub guesses (S, f), with |S| ≤ c and f an injective, edge-preserving
  map from S into V_c(G). Guesses that fail those conditions are left out, and the main
  proof shows they are impossible.
- **Guess** = cube(S ↦ f) ⊗ Rest_{S,f}, where Rest is over V(H) \ S.
- **Rest_{S,f}** = Exact_0̂ for the components of H \ S.
- **Exact_σ** = Prod_σ ⊖ (⊕ of Exact_τ over τ ∈ Up(σ) \ {σ}), shared across parents.
  This is the memo table.
- **Prod_σ** = ⊗ of Leaf_K over the components K of H/σ.
- **Leaf_K** = ⊕ of cubes, one per embedding of K into G \ V_c(G) that is consistent with
  f, with glued vertices taking the same value.

Lists need no special treatment: they determine which cubes exist, and a glued vertex's
list is the intersection automatically. Scopes line up (a cube over S ⊗ a node over
V(H) \ S covers V(H)). The size is about m^c · c^c guesses × f(m) partitions × Σ #Emb(K)
cubes, the same FPT budget as the algorithm. Exact nodes cannot be shared between guesses,
because their leaves differ, and isomorphic glued components in different places are
separate nodes. A renaming node, count(ρ(n)) = count(n) for a variable bijection ρ, would
recover the isomorphism sharing; it is deferred.

## The verification model

1. The solver writes a VeriPB proof from the certified SIP encoding into a counting-PB OPB
   file, plus a counting-language statement describing its argument.
2. A separate translator and counter read the statement, write their own copy of the
   counting-PB OPB file, and print the count.
3. The solver's OPB file is discarded and replaced by the counter's.
4. VeriPB checks that the SIP encoding and the counter's OPB file are equienumerable, using
   the solver's proof.
5. The translator also writes the side-condition instances, and VeriPB checks an UNSAT
   proof for each.

Initially, the trust base is the translator and counter (not yet formally verified) plus
the glue between the steps. Both are purely structural, and neither needs to understand
inclusion–exclusion.

### Lining up the files

`EQUIENUMERABLE FILE` needs the output file to match the final database exactly, so:

- **Use one implementation.** The solver emits its gate definitions in the proof by
  calling the translator's code, so the two copies cannot diverge. The independent copy is
  the one checked.
- **Order.** If both walk the DAG in the same order, constraint order matches. If VeriPB
  turns out to compare order and that proves fragile, the options are a VeriPB feature
  request or the planned PERMUTATION output.
- **Names.** x variables follow the certified encoder's naming exactly, as for the CakePB
  encoding. Gate variables are derived from node ids under a prefix the encoder never
  uses.
- **The encoding's own constraints.** Everything from the SIP encoding is deleted except
  what the translator also emits. The overlap is the exactly-one group for each pattern
  vertex: either the translator copies the encoding's exact syntax, or the proof re-derives
  the translator's form and deletes the encoding's.

### The main proof

- **Gates** are introduced by `red`, in both directions, and are not in the preserved set,
  so introducing them preserves the count. No `preserved_add` or `preserved_rm` is needed,
  because everything stays over x.
- **Input ⇒ root.** Each component's embedding lies in its leaf, by that leaf's search
  proof, logged into the one proof file under the guess literals as assumptions. The
  overlaps contradict the input by RUP: every cube under σ ≠ 0̂ sets x_{u,t} and x_{v,t}
  for some glued pair, against the at-most-one constraint for t.
- **Root ⇒ input**, needed to delete the input's constraints as in Section 4.4 of the
  CP 2026 paper. A colliding element of Prod_0̂ has some exact pattern ρ, and to show it is
  in the overlap needs *glued-leaf completeness*. Glued cubes deliberately violate
  injectivity, so those searches must be logged against adjacency plus per-component
  at-most-one constraints, derived from the full ones before they are deleted. The
  deletion order matters, and the proof contains searches over assignments that the input
  rules out. This is the hardest part of the proof engineering.

## Plan

Each phase ends with something that can be checked on its own. Phases 1 and 2 can run in
parallel.

### Phase 0: land the c = 0 prototype

- Open a PR for `decomposition/component-counting`. It also fixes #104 (the
  pattern-bigger-than-target result now reports a complete search); the issue's other
  requests, an `extra_stats` line naming that step and tightening the allowance in
  `random_homomorphism_test.cc`, should go in the same PR.
- Done when: CI is green, with #106 (an unrelated proof-sweep failure) noted.

### Phase 1: emit statements, count them independently

- Define a file format for counting-language statements: nodes, scopes, the variable and
  value names, and the gate-name convention.
- Have the decomposer emit the statement for its argument (c = 0 for now), with leaves
  coming from the per-component enumerations.
- Write a standalone translator and counter: statement → PB file and count.
- Check syntactic side conditions in the counter.
- Done when: the counter reproduces `solution_count` on the random tests and the grid
  instances, and the translator's PB formula has the same count as the input under brute
  force on tiny instances.

### Phase 2: VeriPB feasibility, by hand

Before building proof plumbing, settle the questions VeriPB has to answer:

- Hand-write the P4 example's full certificate: the main proof from a small SIP encoding
  into the gate formula, plus the side-condition instances and their UNSAT proofs.
- Find out exactly what `EQUIENUMERABLE FILE` compares, including constraint order.
- See how large the exact-pattern disjointness refutation is. It is the side condition
  most likely to be expensive, and the mathematics needs it whatever the encoding.
- Done when: the hand-written certificate verifies, or there is a precise list of what
  VeriPB would need.

### Phase 3: c > 0, without proofs

- Add an API for restricting initial domains (useful elsewhere too), so a component can be
  searched with its hub neighbours pinned and its lists applied.
- Implement the hub-guess layer and lists, choosing c (equivalently Δ) at runtime from the
  degree drop-off. Any choice is correct; it only affects speed.
- Extend the statement emission and the tests, including the star-into-hub instance that
  defeats a single search.
- Done when: counts match the plain solver on random instances, and the counter reproduces
  them.

### Phase 4: proof plumbing, input ⇒ root

- Let subcalls log into one shared proof over the original variables, under assumptions,
  instead of writing their own models.
- Emit the gate definitions through the translator's code, and derive the root.
- Done when: VeriPB accepts the proof up to the point of deriving r_root, on small
  instances.

### Phase 5: root ⇒ input, and deletion

- Derive the per-component at-most-one constraints, run and log the glued-leaf
  completeness searches, then delete the input constraints in a valid order.
- Done when: `EQUIENUMERABLE` verifies end to end against the counter's file, for c = 0
  and then c > 0.

### Phase 6: side conditions, and the pipeline

- Have the solver produce UNSAT proofs for the translator's side-condition instances.
- Write the pipeline script: run, replace the OPB file, verify everything, report the
  count only if every step passes.
- Add it to the randomised proof sweep (`test-instances/random_proof_sweep.bash`).

### Later

- Move the side conditions into the main proof, via a definitions-only lemma in VeriPB.
- Formally verify the translator and counter, alongside the CakePB encoding work.
- Renaming nodes, to recover isomorphism sharing in certificates.
- Use the symmetry between identical components when enumerating coarsenings.
- A decision version, using the observation that a component with more than about
  m² Δ^(m−1) list embeddings can always be placed last.

## Open questions

- How expensive are the exact-pattern disjointness proofs in practice?
- Does `EQUIENUMERABLE FILE` care about constraint order?
- How big do certificates get on instances where decomposition actually wins, compared
  with enumeration proofs of the same counts?
- When should the solver choose to decompose at all?
