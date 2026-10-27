# Certified counting by decomposition

A design plan for counting subgraph embeddings by decomposition, following Section 3.6 of
Jessica Ryan's PhD thesis, and for certifying those counts with VeriPB without enumerating
the solutions. This is forward-looking: it describes where the code is going, not where it
is. The algorithm exists in `gss/sip_decomposer.cc`, including the hub guesses (phase 3).
It can state its argument as a counting statement (phase 1), and prove that the solver's
encoding of the instance is equienumerable with that statement, with or without hubs (phases 4
and 5). The side-condition proofs and the pipeline (phase 6) are still a plan. For the current architecture see
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

`solve_sip_by_decomposition` (`--decomposition`) does this when counting injective,
non-induced, undirected embeddings. Vertex labels, edge labels and loops are supported. In
brief:

- First, the plain solver's checks before search run on the whole instance
  (`HomomorphismModel::prepare`). They can rule it out where every component would fit on
  its own, by comparing degree sequences for instance, and decomposing never sees the whole
  instance. On the LV benchmark instances with disconnected patterns, the plain solver
  rejects 749 of 1159 this way, where decomposition used to grind through the inclusion–
  exclusion until it timed out. `DecompositionOptions::check_whole_instance` turns this off,
  which the tests do, so that they still exercise the decomposition on those instances. It
  does not yet happen under proof logging, which would need the plain solver's proof of it.
- Isolated pattern vertices use a closed form: a falling factorial for each vertex label.
- Everything else is guessed onto the hubs in every possible way, and what is left of each
  guess goes through the recursion above, with lists. Gluing vertices with different labels,
  or gluing two edges with different labels onto each other, gives zero, and so does gluing
  vertices whose lists have nothing in common.
- The hubs are the c highest-degree target vertices, ties going to the lower-numbered
  vertex. `DecompositionOptions::hubs` (the hidden `--decomposition-hubs N`) sets c.
  Otherwise c is chosen to minimise a rough estimate of the work, (number of guesses) × P
  + E, worked out in logarithms. Here:
  - P is the number of valid partitions of the pattern's vertices, which bounds each
    guess's inclusion–exclusion;
  - E estimates the embeddings that enumeration visits: of each component, and of the
    biggest glued component, which glues one vertex of greatest degree from each component
    together;
  - a graph with e edges and degrees g(v) gets Π_v F(g(v)) / F(1)^e, where F(j) =
    Σ_t d(t)_j over the non-hubs t, d(t) is the degree once the hubs are gone, and d_j is a
    falling factorial. This is the expected number of embeddings into a random graph with
    those degrees (the configuration model).

  The estimate is exact for a star, which is what makes a heavy-tailed degree distribution
  worth guessing on: the top few vertices dominate F(j) for large j. It charges little for
  a dense pattern in a sparse target. An earlier version used Σ d(t)^(k−1) for the biggest
  glued component of k vertices, which is the same for a star but about 10⁶⁰ for a
  45-vertex component. That made any number of guesses look cheap, so it picked c in the
  tens and then never finished guessing. A grid gets c = 0 on all of the patterns in the
  table below. The number of guesses only grows with c, so the search for c stops once the
  guesses alone cost more than the best so far. Hubs need an undirected target.
- A guess puts some of the pattern vertices onto distinct hubs, respecting labels, loops
  and the edges between them, and only onto hubs with at least their degree. Every other
  vertex gets a list. It contains the non-hubs that have:
  - an edge, with the right label, to every hub that its neighbours were guessed onto;
  - at least as many non-hub neighbours as it has neighbours not on hubs.

  Both hold of any embedding, because it is injective and every vertex not on a hub goes to
  a non-hub. A guess that leaves any vertex an empty list is skipped. The degree condition
  is what prunes guesses that leave a high-degree pattern vertex off the hubs when only the
  hubs could take it.
- The guesses and the coarsenings check the timeout themselves, because there can be
  enough of them to outlast it without any solver call noticing. The lists reach the solver through
  `HomomorphismParams::restrict_initial_domains`, which any caller can use.
- Each glued component is counted by one call to `solve_homomorphism_problem`, cached up
  to isomorphism: colour refinement, then brute force within colour classes, falling back
  to the exact graph above 5040 orderings. A vertex's list is part of its colour, so the
  cache is shared between guesses wherever their lists agree. Lists are numbered as they
  are first seen, so that a key stays short.
- none() is memoised on the quotient's isomorphism class, which is the sorted multiset of
  its components' keys. This is shared between guesses too.
- A quotient with more blocks of some label than the target has vertices with that label
  has no injective embedding, and is skipped. Without this, a pattern bigger than its
  target explores every way of gluing itself together.
- It reports `decomposition_components`, `decomposition_hubs`, `decomposition_guesses`,
  `decomposition_solver_calls`, `decomposition_quotients` and
  `decomposition_coarsenings_considered`, and `decomposition_rejected_before_search = true`
  when the checks on the whole instance ruled it out.
- It produces no mapping when counting. It allows proof logging only when writing a
  statement, and then makes the same choice of hubs as without it.
- With `--decomposition-statement FILE` it also writes its argument as a counting
  statement; see phase 1 below.

`gss/sip_decomposer_test.cc` compares it with the plain solver on 200 random disconnected
patterns across loops and labels, with the decomposer's own choice of c and with c from 0 to
3.

Counting on a 40×40 grid, where the plain solver runs out of time on all of these:

| Pattern | Decomposition |
|---|---|
| 3 disjoint edges | 2.4×10¹¹ in 64 ms, 5 solver calls |
| 3 × P3 | 6.0×10¹² in 3.3 s, 51 calls |
| 6 disjoint edges | 5.7×10²² in 12.7 s, 52 calls, 1.8M coarsenings |
| 4 × P4 | times out at 120 s: 812 distinct glued trees, each counted by enumeration |

Known limits:
- The coarsening enumeration does not use the symmetry between identical components, which
  dominates the six-edge case.
- Large glued components are still counted by enumeration.
- The choice of c is a rough heuristic.
- Inclusion–exclusion is at its worst when the count is zero or small but the pieces have
  many embeddings. Then it works through every quotient to reach a total that a plain search
  of the whole pattern would refute quickly. For instance, one hub splits a 24-vertex
  pattern into four small dense graphs that each fit in LV g13 many times. That case takes
  38,635 solver calls and times out at 100 s, where the plain solver refutes the whole
  pattern in 1.5 s.

## The counting language

The aim is to move as much of the argument as possible into VeriPB rather than into an
external combiner. The solver states its argument in a small language whose statements
are efficiently countable, and VeriPB checks that the problem is equienumerable with a PB
translation of that statement. The file format and the exact translation are specified in
[counting-statements.md](counting-statements.md); this section is about why the language
is the way it is.

### Nodes

Every node is a **reified gate**: an extension variable r_n defined in both directions
over its children. Gate variables are therefore functionally determined by the original
variables, which is what makes negation, and hence subtraction, cheap.

Variables are **multi-valued**: each is an exactly-one group, such as y_{u,·} for pattern
vertex u. A literal fixes one value and implicitly sets the rest of the group false. Without
this, smoothing over Boolean variables would cost a factor of n.

The statement's variables are its own, not the input's. The counting formula mentions
nothing from the input encoding, and the proof moves the preserved set across to it, as in
Section 4 of the CP 2026 paper. That keeps the verified counter and translator independent
of the problem: they would serve a constraint solver whose preserved set is the bits of its
integer variables, or an AND/OR search, as well as this. It also leaves the solver free to
choose its variables. One per pattern vertex is a renaming; one per component, with one
value per embedding, would make each component's leaf just its exactly-one group. Phase 2
measured the renaming as a small, constant cost per input variable.

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
groups, and the unit constraint r_root. Because every gate is determined by the variables,
the formula's count projected onto them is the root's count. For example, ⊖ is

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
- The translator emits all of those as one PB instance over the statement's own
  definitions, with a selector variable for each condition saying that it is violated
  (for a ⊖ node, r_b ∧ ¬r_a), and a clause saying that one of them is. That needs an
  ordinary UNSAT proof. Because the instance contains only definitions, it cannot
  accidentally depend on the problem's constraints, and no VeriPB extension is needed.
- Phase 2 confirmed that this is load-bearing. A statement that lists one of its glued
  leaves twice under the same ⊕ translates to a formula whose main proof still verifies
  as `EQUIENUMERABLE`, because the duplicate does not change which assignments the root
  accepts. But the counter adds the duplicate's count twice, and prints 4 for an instance
  with 8 solutions. Only the side-condition instance, which then has a solution, catches
  it.
- Later, these could move inside the main proof with a narrowly scoped VeriPB feature: a
  lemma that may use only the output formula's definitions.

## Compiling Ryan's algorithm

- **Root** = ⊕ over hub guesses (S, f), with |S| ≤ c and f an injective, edge-preserving
  map from S into V_c(G). Guesses that fail those conditions are left out, and the main
  proof shows they are impossible. So are guesses whose Rest has no node, because some
  component has no embeddings, or by pigeonhole. With only one guess left, the root is
  that guess.
- **Guess** = cube(S ↦ f) ⊗ Rest_{S,f}, where Rest is over V(H) \ S. When S is empty,
  this is just Rest, and when S is everything, just the cube. With c = 0, the only guess
  has S empty, so the statement is the same as before there were hubs.
- **Rest_{S,f}** = Exact_0̂ for the components of H \ S.
- **Exact_σ** = Prod_σ ⊖ (⊕ of Exact_τ over τ ∈ Up(σ) \ {σ}), shared across parents.
  This is the memo table.
- **Prod_σ** = ⊗ of Leaf_K over the components K of H/σ.
- **Leaf_K** = ⊕ of cubes, one per embedding of K into G \ V_c(G) that is consistent with
  f, with glued vertices taking the same value.

Lists need no special treatment: they determine which cubes exist, and a glued vertex's
list is the intersection automatically. Scopes line up (a cube over S ⊗ a node over
V(H) \ S covers V(H)). The size is about m^c · c^c guesses × f(m) partitions × Σ #Emb(K)
cubes, the same FPT budget as the algorithm. Exact nodes are shared between guesses only
where their blocks have the same lists, which is what the memo is keyed on. Isomorphic glued
components in different places are separate nodes. A renaming node, count(ρ(n)) = count(n) for a variable bijection ρ, would
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

What `EQUIENUMERABLE FILE` compares, from `veripb-checker/src/rules/output.rs` and checked
in phase 2: every constraint in the output file must be in the database and in core, and
every core constraint must be in the output file, looked up after normalisation. So it is a
set comparison. Order does not matter, a duplicated row is harmless, and `-1 x … >= -1`
matches `1 x … <= 1`. Derived constraints are ignored. The preserved sets must be equal as
sets. Checked deletion must never have been downgraded, no solution may have been logged,
and no witness may touch a preserved variable.

- **Use one implementation.** The solver emits its gate definitions in the proof by
  calling the translator's code, so the two copies cannot diverge. The independent copy is
  the one checked.
- **Names.** The counting formula uses only its own names (`y`, `g`), which the encoder
  never uses, so nothing has to match the encoder's naming. That would not have been
  true of the formula sharing the input's x variables, which would have tied the
  translator to one encoder's names and broken for inputs with named vertices.
- **The encoding's own constraints.** Everything from the SIP encoding is deleted,
  including its exactly-one rows, once the y versions have been derived from them.

### The main proof

- **Variables.** Each y variable is defined from the input's variables and added to the
  preserved set, and at the end each x variable is removed from it, and the links are
  deleted with an x ↦ y witness. The steps are listed in
  [counting-statements.md](counting-statements.md), "From an input to the counting formula".
- **Gates** are introduced by `red`, in both directions, and are not in the preserved set,
  so introducing them preserves the count.
- **Input ⇒ root.** Each component's embedding lies in its leaf, by a search of that
  component logged into the one proof. The search runs under a guard that it does not
  know in advance. Until it finds a solution there is no guard. Each solution's cube is
  then written to the statement and its gate defined, and a fold literal is defined:
  f_k ↔ f_{k−1} ∨ cube_k. From then on every backtracking nogood carries "∨ f_k". So the
  search concludes f_N, which implies the leaf, defined only afterwards over all the cubes.
  Every definition is fresh when it is made, so nothing has to be known before the search
  starts. Propagation failures need no guard: no solution extends them. The overlaps
  contradict the input by RUP: every cube under σ ≠ 0̂ sets x_{u,t} and x_{v,t} for some
  glued pair, against the at-most-one constraint for t.
- **Root ⇒ input**, needed to delete the input's constraints as in Section 4.4 of the
  CP 2026 paper. A colliding element of Prod_0̂ has some exact pattern ρ, and to show it is
  in the overlap needs *glued-leaf completeness*. This was expected to need searches logged
  against adjacency plus per-component at-most-one constraints. Phase 2 found otherwise:
  the adjacency rows go by unit propagation, and the injectivity rows by the collision
  lemma, which needs only the statement's own definitions (see phase 2's findings). This
  is still the hardest part of the proof engineering.

## Plan

Each phase ends with something that can be checked on its own. Phases 1 and 2 can run in
parallel.

### Phase 0: land the c = 0 prototype

- PR #107 for `decomposition/component-counting`. The complete-flag and proof-sweep
  failures it ran into (#104, #106) are fixed on `main`.
- Done when: CI is green and the PR is merged.

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
- Status: done, except for the larger grid instances, which need the renaming node.
  - The format and the translation are specified in
    [counting-statements.md](counting-statements.md). It is an s-expression in the same
    small dialect as the constraint solver's `.scp` files, chosen so that a verified
    parser stays simple: numbers and a few keywords only, no names, strings or comments,
    fixed sections, and node ids that must count up from zero with children defined first.
  - `gss/counting_statement.{hh,cc}` reads, checks, counts and translates statements, and
    writes them as they are built. `src/count_statement.cc` is the standalone counter and
    translator. `--decomposition-statement FILE` makes the solver write one.
  - Statements build the c = 0 recursion over every pattern vertex: isolated vertices are
    components like any other, and there is one node per quotient, not per isomorphism
    class. Isomorphic components still share one enumeration, renamed into place.
  - There is no node for "nothing". If some component has no embeddings, or the pattern
    cannot fit by pigeonhole, no statement is written, and an ordinary unsatisfiability
    proof is the right certificate. A count of zero for any other reason is expressible,
    as a ⊖ whose two sides are equal.
  - A pigeonholed quotient gets no node, so it is left out of its parents' ⊕ nodes. That
    is sound for the count, but the certificate will have to refute those quotients, by a
    counting argument from the at-most-one constraints.
  - `gss/counting_statement_test.cc` tests the language on its own: counts, every way of
    being malformed, and the exact translation.
  - `gss/sip_decomposer_test.cc` checks statement counts against the plain solver on
    random disconnected patterns. On tiny instances it evaluates every assignment against
    the meaning of every node and the side conditions. On the solutions, and a sample of the
    other assignments, it also checks the translator's rows: every definition holds, the
    root row holds exactly when the root does, flipping any one gate breaks one of its own
    rows, and no violation selector in the side-condition instance can be switched on.
    Planting any of six bugs, three in the translator and three in the decomposer, makes it
    fail.
  - Phase 2's certificate is a regression test (see below), and the translator reproduces
    an independently written Python translator's output byte for byte on it.

  Measured sizes:

  | Instance | Nodes | File | Build | Count |
  |---|---|---|---|---|
  | 3 edges, 6×6 grid | 34k | 1.9 MB | 0.1 s | 0.1 s |
  | 3 × P3, 6×6 grid | 2.1M | 175 MB | 3.8 s | 4.0 s, 640 MB |
  | 3 edges, 40×40 grid | 2.3M | 161 MB | 3.1 s | 3.5 s, 470 MB |
  | 6 edges, 6×6 grid | stopped at 300 s, after 15.7 GB of statement | | | |
  | 4 × P4, 6×6 grid | stopped after 3 GB of statement | | | |

  The solver writes statements as it goes, so it builds them in under 25 MB, but they get
  big quickly. Two things blow up: there is one node per partition rather than per
  isomorphism class (6 edges have 113 classes, and vastly more partitions), and every glued
  component's embeddings are written out as cubes. The renaming node (see "Later") fixes
  the first, and is a prerequisite for anything that is not small. The second is inherent
  to leaves being plain enumerations. The counter holds the whole statement in memory, as
  a verified one probably would.
- Instances are in `test-instances/decomposition/`, as CSV:
  - `grid6` and `grid40` are 6×6 and 40×40 grid graphs;
  - `edges2`, `edges3` and `edges6` are disjoint unions of edges;
  - `p3x3` is three disjoint P3s, and `p4x4` is four disjoint P4s;
  - `path4` is the path on four vertices.

  Run experiments under `ulimit -v` and `ulimit -f`, and use the solver's own `--timeout`
  rather than `timeout`, so that it can remove an incomplete statement.

### Phase 2: VeriPB feasibility, by hand

Before building proof plumbing, settle the questions VeriPB has to answer:

- Hand-write the P4 example's full certificate: the main proof from a small SIP encoding
  into the gate formula, plus the side-condition instances and their UNSAT proofs.
- Find out exactly what `EQUIENUMERABLE FILE` compares, including constraint order.
- See how large the exact-pattern disjointness refutation is. It is the side condition
  most likely to be expensive, and the mathematics needs it whatever the encoding.
- Done when: the hand-written certificate verifies, or there is a precise list of what
  VeriPB would need.
- Status: done for two disjoint edges into P4, which is in
  `test-instances/counting/two_edges_into_path4/` and checked by the
  `counting_certificate_two_edges_into_path4` test. It has the statement, the solver's
  encoding of the instance, the translator's two files, and two proofs:
  - `main.pbp`, 758 lines: `veripb --force-checked-deletion input.opb main.pbp output.opb`
    says `s VERIFIED OUTPUT EQUIENUMERABLE`;
  - `side.pbp`, 126 lines: `veripb side.opb side.pbp` says `s VERIFIED UNSATISFIABLE`.

  How the main proof goes:
  1. Each y variable is defined as a copy of its x variable and added to the preserved set,
     and the exactly-one rows over y are derived. Then every gate is introduced by `red`
     in both directions, autoproved, and moved to core.
  2. Input ⇒ root. Each component's leaf holds, by a search over that component that
     closes off each solution with the clause "not this assignment, or the leaf gate"
     instead of `solx`. Every cube of a glued leaf is refuted by RUP against one
     injectivity row. Then the ⊕ and the ⊖ follow by RUP.
  3. Root ⇒ input, by deleting the encoding's rows.
     - Each adjacency row is deleted with no subproof at all. Unit propagation from the
       root reaches the component's leaf, and every cube in it is falsified, either by
       the at-most-one or by the missing neighbour.
     - Each injectivity row is deleted with a subproof. It derives the pairwise clauses
       ¬x_{u,t} ∨ ¬x_{v,t}, recovers the at-most-one from them with `pol`, and
       contradicts the negated row. For u and v in the same component, the pairwise
       clause is one RUP step, because every cube is injective. Across components, it
       is a case split on the other vertices' values: each case completes a cube of some
       glued leaf, and the root says that no glued leaf holds.
  4. The x exactly-one rows are deleted, each x variable leaves the preserved set, and the
     links are deleted with x ↦ y witnesses.
  5. `output EQUIENUMERABLE FILE`, against a file that mentions no x at all.

  Findings:
  - Moving the preserved set to the statement's own variables costs 144 of the 758 lines,
    a constant number per input variable, and VeriPB proves every one of those steps by
    itself apart from four one-line `pol` subproofs. An earlier version kept the x
    variables and their exactly-one rows in the output, which needed none of this (VeriPB
    normalises before comparing, so the rows matched), but tied the translator to the
    encoder.
  - Root ⇒ input does not need the adjacency rows. The proof still verifies with every
    adjacency row deleted *before* the injectivity rows. The completeness it needs comes
    from the glued leaves' own cubes, not from searches against adjacency. So the plan's
    "glued-leaf completeness searches logged against adjacency plus per-component
    at-most-one" are unnecessary, at least here.
  - The side conditions are load-bearing, as the plan said: see "Side conditions affect
    the count".
  - Here every Exact under the root's ⊕ is a leaf, because every quotient of two edges is
    connected. Disjointness was one RUP per cube per pair of leaves, then `pol`.

  With three or more components, an Exact node under a ⊕ is a ⊖ in its own right, and
  disjointness and root ⇒ input turn out to need the same lemma. Call it the *collision
  lemma* at τ: Prod_τ ∧ x_{u,t} ∧ x_{v,t} → B_τ, for u and v in different components of
  H/τ, where B_τ is the ⊕ under Exact_τ.
  - Two different Exacts under one ⊕ have partitions that differ on some pair u, v,
    glued in one (so every cube there has u and v equal) but not in the other.
  - If u and v are in the same component of the other quotient, its cubes all separate
    them. That gives one RUP step per value t, with no case split.
  - Otherwise the other Exact implies, by its collision lemma, that u and v differ.
  - Deleting injectivity row t is the collision lemma at the root.

  The lemma follows by cases on pairs of cubes, one from each of the two components that
  collide. A colliding pair completes a cube of the glued leaf of some ρ above τ, and so
  Prod_ρ by unit propagation, and Prod_ρ implies Exact_ρ ∨ B_ρ. Every node under B_ρ is
  also under B_τ, so B_ρ → B_τ is one RUP step. The number of cases is the number of
  colliding cube pairs, which is about the size of the glued leaves, so this is polynomial
  in what the solver already enumerated.

  Checked by generating both proofs for statements the solver writes, with a Python
  generator that reads the partition structure back off the statement's leaves. It handles
  unlabelled, loopless LAD instances with at least as many target vertices as pattern
  vertices, so nothing is pigeonholed. VeriPB accepts the main proof as
  `OUTPUT EQUIENUMERABLE` and the side-condition proof as `UNSATISFIABLE` on:
  - all 98 random instances with a statement among 100 small ones (3 or 4 components, at
    most 6 pattern vertices, the other two having count 0);
  - 15 larger ones (3 or 4 components, up to 7 pattern vertices, up to 39k statement
    nodes), which are all of those that VeriPB finished; of the rest, 10 had their proofs
    generated but not checked within 30 minutes, 3 hit a 2 GB cap on the statement, and
    1 had count 0;
  - three disjoint edges into P6.

  Disjointness has two possible shapes, and the first one tried does not scale:
  - Pairwise, as in the P4 certificate, is quadratic in the number of Exacts under a ⊕, and
    under the root's ⊕ that is one per way of gluing. It was two thirds of the side proof.
  - A decision tree over vertex pairs is linear in that number. Split the Exacts on a pair
    u, v that some glue and some do not. Introduce a gate a ↔ ∨ (each half), and exclude
    the two halves from each other with one RUP step per value. Then combine with the same
    `pol`-and-divide trick that recovers an at-most-one from pairwise clauses. The
    generator uses this, and it halved the side proof on three edges into P6.

  Sizes, per statement node, over the 98 small instances: the main proof is about 5 lines;
  the side-condition proof is 9 with three components and 18 with four. The side proof
  is now spread over restriction lemmas for containment, equality and separation lemmas
  for single leaves, the collision lemmas' cases, and the decision trees, with none of
  them dominating. On larger instances the side proof grows to about 16 to 21 lines per
  node. Checking it is the slow part. On one 7-vertex instance with 39k nodes and 1008
  solutions:
  - VeriPB took 17 s for the 160k-line main proof;
  - it took 34 minutes for the 640k-line side proof;
  - the solver's own enumeration proof took 0.05 s.

  The side proof never deletes anything, which is the first thing to look at; this hasn't
  been investigated yet. These are not instances where decomposition should win, since
  their counts are small. The statements themselves remain the hard limit: random
  patterns with 10 or 11 vertices hit the 2 GB cap on the statement alone.

  The generator is in `tmp/gss-decomp/phase2b/` on the shared cluster, not in this
  repository. The structure it recovers from leaves is what the solver will know directly
  when it writes these proofs itself (phases 4 to 6).

### Phase 3: c > 0, without proofs

- Add an API for restricting initial domains (useful elsewhere too), so a component can be
  searched with its hub neighbours pinned and its lists applied.
- Implement the hub-guess layer and lists, choosing c (equivalently Δ) at runtime from the
  degree drop-off. Any choice is correct; it only affects speed.
- Extend the statement emission and the tests, including the star-into-hub instance that
  defeats a single search.
- Done when: counts match the plain solver on random instances, and the counter reproduces
  them.
- Status: done, as described in "Where the prototype is" above.
  - `HomomorphismParams::restrict_initial_domains` is the API: a predicate on (pattern
    vertex, target vertex), applied when the domains are first built. The target-loop and
    clique shortcuts step aside for it, and proof logging refuses it for now.
    `random_homomorphism_test` checks it against the brute-force oracle in every family
    and mode.
  - The statement for a guess is cube(S ↦ f) ⊗ Rest, and the root is the or of the
    guesses that have embeddings. With c = 0, statements are as they were.
  - `sip_decomposer_test` checks the counts and the statements with the decomposer's own
    choice of c and with c from 0 to 3, including the brute-force semantics and side
    conditions with c up to 2. A star with four leaves into a star with 199 leaves (1.5 ×
    10⁹ solutions) takes c = 1 and counts at once. A wider sweep of 400 random instances
    into small heavy-tailed targets, with c up to 6, also agreed with the plain solver
    and `count_statement` every time.
  - Proofs with hubs came later, in phase 5. The side conditions still have to cover the
    guesses' disjointness, which is phase 6.
  - Counting, on a 300-vertex preferential-attachment graph (maximum degree 43), and on a
    120-vertex G(n, 0.06):

    | Target | Pattern | Plain | c = 0 | Chosen c | Time |
    |---|---|---|---|---|---|
    | preferential | K₁,₅ | 44.6 s | 42.8 s | 10 | 0.3 s |
    | preferential | K₁,₄ + P₃ | > 120 s | > 300 s | 5 | 13.6 s |
    | preferential | 3 × P₃ | > 120 s | > 300 s | 4 | 26.9 s |
    | preferential | 6 disjoint edges | | > 300 s | 1 | 124 s |
    | G(n, p) | 3 × P₃ | | 37.8 s | 3 | 23.0 s |
    | G(n, p) | K₁,₄ + P₃ | | 12.8 s | 5 | 5.4 s |

    Guesses are cheap: 3 × P₃ took 21.5 s with 163k guesses at c = 8, its best. On the
    40 × 40 grid a needless hub costs 1.2 to 1.5 times as long, which the heuristic
    avoids for the patterns in the table above. 4 × P₄ on the preferential graph times
    out at every c tried.

### Phase 4: proof plumbing, input ⇒ root

- Let subcalls log into one shared proof over the original variables, under assumptions,
  instead of writing their own models.
- Emit the gate definitions through the translator's code, and derive the root.
- Done when: VeriPB accepts the proof up to the point of deriving r_root, on small
  instances.
- Status: done.
  - `--prove NAME` now works with `--decomposition` when `--decomposition-statement` is
    given (and only then), for unlabelled instances.
  - The proof:
    1. writes the solver's usual encoding of the whole instance;
    2. moves the preserved set to the statement's variables;
    3. defines every node's gate, into the core, as the statement writer writes it,
       through `counting_node_definition`;
    4. searches each component of the pattern as a nested solve under the fold guard
       above;
    5. derives the root as a unit, and stops with `conclusion NONE`.

    When a component has no embeddings, its search derives a contradiction, and the proof
    concludes that there are no solutions instead.
  - `gss/innards/decomposition_proof.{hh,cc}` has the proof's side of this. The nested
    solve is `Proof::begin_nested_solve`, which `HomomorphismParams::nested_proof` points a
    solve at. A nested solve's pattern must name its vertices as the outer pattern does.
    Then its variables and its at-most-one, at-least-one and adjacency rows are the outer
    model's and are not written again. Its injectivity rows are derived per target, by
    adding the other pattern vertices' literals to the outer row.
  - Component searches do not share enumerations by isomorphism, because each proof is
    about its own variables. This costs at most a factor of the number of components.
  - For now, component searches run with `--no-supplementals --no-nds`. Glued components,
    which need no search proof, run as before.
  - `test-instances/decomposition_proof_sweep.bash` (ctest `decomposition_proof_sweep`)
    checks 30 random instances, with and without loops. For each it checks the count,
    `count_statement`, VeriPB, and that the root is derived. Before it was cut down for
    ctest, the same checks passed on 100 small instances and 40 loopy ones of up to 7
    pattern vertices.
  - Planting bugs shows which steps matter:
    - dropping the guard from backtracking nogoods breaks every instance;
    - so does dropping it from solution clauses, or leaving the previous guard out of a
      fold;
    - dropping it from propagation failures changes nothing, which is why they do not
      carry one;
    - skipping the refutation of every glued cube breaks almost every instance, but
      skipping the last of a leaf's cubes does not: unit propagation finds that one by
      itself.

### Phase 5: root ⇒ input, and deletion

- Delete the adjacency rows, which unit propagation from the root justifies, and the
  injectivity rows, by the collision lemma at the root (see phase 2). Pigeonholed
  quotients, which have no nodes, need a counting argument instead.
- Done when: `EQUIENUMERABLE` verifies end to end against the counter's file, for c = 0
  and then c > 0.
- Status for c = 0: done.
  - The proof now ends with `output EQUIENUMERABLE FILE` against `count_statement --opb`'s
    file. VeriPB accepts it with `--force-checked-deletion`.
  - After the root (phase 4):
    1. The root goes into the core.
    2. The adjacency rows are deleted with no subproof.
    3. Each injectivity row is deleted with a subproof:
       - for each pair of pattern vertices, "not both at t";
       - within a component, one RUP step;
       - across components, the collision lemma at the root, by cases on pairs of cubes;
       - then the at-most-one is recovered with `pol`.
    4. The x exactly-one rows are deleted, x leaves the preserved set, and the copies are
       deleted with x ↦ y witnesses.
  - The collision lemma needs no `B_ρ → B_τ` lemmas, unlike the Python generator's side
    proofs. At the root, ¬B_0 falsifies every overlap of the root by unit propagation, and
    ρ's own overlaps are all among them, so B_ρ fails too. Planting their removal changed
    nothing, which is how this came to light. Without the per-pair-of-cubes cases, 4 of the
    30 sweep instances fail: the rest have at most one candidate cube, which unit
    propagation finds by itself.
  - The overlap that a pair of colliding cubes completes is found by its partition, which
    the decomposer hands over with the root's structure.
  - Pigeonholed quotients never come up here. Every overlap of the root has fewer blocks
    than the root, and the root has passed the pigeonhole test, or there would be no
    statement. They will come up in the side-condition proofs, at inner quotients.
  - `decomposition_proof_sweep` checks `OUTPUT EQUIENUMERABLE` (or zero solutions) on its 30
    instances. The same checks pass on:
    - 100 small instances;
    - 40 loopy ones;
    - 22 larger ones, of up to 8 pattern vertices and 371k statement nodes.

    Checking time is mostly VeriPB's:
    - a 39k-node statement (1008 solutions) checks in 19 s, where the Python generator's
      side proof alone took 34 minutes;
    - the 371k-node one takes 18.5 minutes for 1.5M proof lines;
    - three more, of 3M to 14M proof lines, did not finish checking within 30 minutes.

    The side conditions are not part of this proof; they are phase 6.
  - Labels are now defined once each: `Proof` throws if a line reuses one. Before this, the
    loop-cancelled adjacency rows took over the model rows' `@adj` labels, which VeriPB
    allows. A nested solve's derivations then leaked into the outer proof: deleting
    `@adj…` deleted the derived row, not the model's, and the output check failed on
    loopy instances. The loop-cancelled rows now have labels of their own (`…_loopfree`),
    which the adjacency cache points at.

- Status for c > 0: done, but certificates on the instances where hubs win are too big to
  check (see the sizes below).
  - The proof walks the same tree of guesses as `for_each_guess`. The decisions at a node
    are literals, each saying that a pattern vertex is on some hub (x) or on none of them
    (`offhubv`, defined as the and of v's negated hub literals). Each node derives "not
    all of these, or the fold so far", where the fold is over the guesses' gates in the
    order they are written, f_k ↔ f_{k−1} ∨ G_k, as in a component search.
    - At a leaf with embeddings, this is unit propagation through the guess's gate: its
      cube from the decisions, its rest's leaves from their searches, and its rest's
      overlaps from their glued leaves, which are refuted as for c = 0.
    - At a leaf without, it is unit propagation from why there are none: a component's
      search failed, or some vertex has no list, or more vertices are off the hubs than
      there are non-hubs. The last needs a `pol`: each such vertex's at-least-one row and
      `offhub` definition, plus injectivity on every non-hub.
    - At any other node, it is unit propagation from the node's children. The vertex it
      decides goes off the hubs, or onto a hub that the node either tried or can rule out.
      A hub already in use, or not adjacent to the hub of a neighbour, goes by unit
      propagation from the decisions. A hub without a loop, for a vertex with one, needs
      its own one-line `rup`, because the adjacency row only conflicts once the
      assignment is assumed. A hub with fewer neighbours than the vertex needs a degree
      `pol`.

    The root of the tree is the last fold as a unit, and the statement's root follows from
    it.
  - Each component of a guess's rest is searched under an *assumption*: a literal defined
    as the and of its vertices' `offhub` literals and the decisions that put its
    neighbours on hubs. The search depends on nothing else, so guesses that agree on it
    share the search and its leaf.
    - The nested solve carries the assumption's negation on every line it writes: its
      nogoods, its solutions, its propagation failures, any value it rules out before
      searching, and its conclusion (`NestedSolve::assumption`). It may then restrict its
      initial domains to the lists, which is otherwise refused under proof.
    - Unit propagation from the assumption gives the hub and hub-adjacency parts of the
      lists. The degree part needs a `pol` for each vertex and target vertex it removes.
      It sums the adjacency rows for the vertex's off-hub neighbours, injectivity on the
      target vertex's non-hub neighbours, and the neighbours' `offhub` definitions. That
      gives "not there, or some neighbour is on a hub", which guesses also share.
  - All of this is scaffolding: derived, never in the core, and deleted once the root is
    derived.
  - Root ⇒ input is as for c = 0, except that the collision lemma is per guess, with the
    guess's gate in place of the root. Under a negated adjacency row, unit propagation
    still falsifies every guess, so those rows still need no subproof.
  - A loop-free adjacency row that an earlier nested solve derived is reused rather than
    derived again under the same label, now that one component can be searched several
    times.
  - `decomposition_proof_sweep` now runs each instance twice: with the decomposer's own
    choice of hubs, and with two hubs. A wider sweep agreed with the plain solver and
    `count_statement`, and VeriPB accepted every proof. It had 480 instances: 60 seeds,
    with and without loops, with c from 1 to 3 and the decomposer's own choice. Every kind
    of new step came up in it.
  - Each of ten planted bugs makes some of those 480 fail:
    - no degree restrictions: 70;
    - no pigeonhole `pol`: 194;
    - no loop units: 7;
    - no hub-degree `pol`: 2;
    - propagation failures without the assumption: 20;
    - nogoods without it: 341;
    - no collision cases: 95;
    - an `offhub` definition missing a hub: 360;
    - leaving the hubs' literals out of a search's assumption: 271.

    So does dropping a guess from the statement. On every instance where that changes the
    count, VeriPB rejects the proof, not just `count_statement`'s disagreement with the
    solver.
  - Sizes, on the instances where decomposition beat the plain solver in the generated
    sweep (2 blocks of 3 to 5 vertices, joined through 2 to 4 vertices):
    - s17 into LV g25 (3 hubs, 748 guesses, 31,424,939 solutions) takes 242 s with proof
      logging, against 7.4 s without. Its statement has 24.1M nodes, 1.7 GB, all but
      8.6k of them cubes, and its proof is 15 GB. Of that:
      - 6.1 GB defines the statement's gates;
      - about 8 GB is the component searches, which found 12.7M embeddings between them
        under the guesses' various assumptions;
      - 1.5 GB is the collision cases.
    - s3, s12 and s17 into LV g23 (2 × 10⁸ to 1.6 × 10⁹ solutions) each exceed a 20 GB
      cap on the proof.

    The cubes are the glued leaves' embeddings, written out one by one, so the
    certificate is about as large as the enumeration the decomposition does, which here
    is close to the number of solutions. That is inherent to leaves being plain
    enumerations. The renaming node would share isomorphic leaves, but with lists those
    are few: counting without a statement needed 2481 isomorphism classes where the
    statement has 4990 quotients. A certificate that is small where decomposition wins
    needs leaves whose counts are themselves certified some other way.

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

- How expensive are the exact-pattern disjointness proofs in practice? Phase 2 found them
  to be the same problem as deleting the injectivity constraints (the collision lemma).
  With a decision tree instead of pairwise exclusion, the side proof is about 10 to 20
  lines per statement node on small random instances.
- How big do certificates get on instances where decomposition actually wins, compared
  with enumeration proofs of the same counts? On the generated instances where hubs win,
  about as big as an enumeration (phase 5), because the glued leaves are enumerations.
- When should the solver choose to decompose at all?
