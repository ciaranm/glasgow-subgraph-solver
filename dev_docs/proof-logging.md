# Proof logging

The solver can emit a machine-checkable certificate of its answer, verifiable with
[VeriPB](https://gitlab.com/MIAOresearch/software/VeriPB). This document covers which option
combinations are supported, how to run it, and the current known limitations. The encoding follows
the constraint-programming-to-pseudo-Boolean approach from the project's papers (see the
[README references](../README.md#references)); for the cleanup status of the open bugs see
the [GitHub issues](https://github.com/ciaranm/glasgow-subgraph-solver/issues).

## What it produces

Passing `--prove NAME` writes two files:

- `NAME.opb` — the pseudo-Boolean **model**: an OPB encoding of the problem instance (one set of
  Boolean variables per CP variable, injectivity constraints, adjacency constraints, …).
- `NAME.pbp` — the **proof log**: the sequence of VeriPB rules justifying every inference,
  backtrack, and the final conclusion.

You then check the pair with VeriPB:

```shell session
$ veripb NAME.opb NAME.pbp
```

For an enumeration or counting run (`--count-solutions`, `--enumerate`, or
`--print-all-solutions`), the model also declares a `preserved:` set — the assignment
variables — so the proof's solution count is in terms of the high-level mapping rather
than any auxiliary encoding variables.

Useful companion flags:

- `--verbose-proofs` writes extra `*` comment lines into the log, for tracing.

## Supported option combinations (homomorphism / subgraph isomorphism)

Proof logging for `glasgow_subgraph_solver` is currently incompatible with a number of the solver's
"extra" features. Requesting `--prove` together with any of these throws an
`UnsupportedConfiguration` (see `gss/homomorphism.cc`, the guard block near the top of
`solve_homomorphism_problem`):

| Requires | Because |
| --- | --- |
| a single thread (no `--parallel`, `--threads 1`) | proof logging is not thread-safe yet |
| `--no-clique-detection` | the clique-detection shortcut is not yet logged |
| no less-than / occurs-less symmetry constraints | not yet logged |
| unlabelled graphs (no vertex or edge labels) | labels are not yet encoded |
| not `--staged` together with `--count-solutions` | the stage transition is a restart, and an enumeration proof cannot yet survive one (same restriction as counting with restarts) |

`--staged` on its own *is* supported: the supplemental graphs are derived mid-proof at the level-0
restart boundary, and an instance that concludes in the cheap first round emits no supplemental
derivations at all.

Injective and non-injective proofs support supplemental graphs, distance-3 (`--distance3`),
neighbourhood degree sequences, and clique-size constraints (`--cliques`), on both loopless and
loopy instances.

Local injectivity (`--locally-injective`) is encoded by *neighbourhood*-injectivity constraints
(`@linj`): for each pattern vertex `v` and target `t`, at most one of `v`'s neighbours maps to `t`
(its closed neighbourhood if `v` has a self-loop). The degree, NDS and exact-path-graph derivations
use these in place of the global injectivity constraints — wherever those derivations would sum "at
most one pattern vertex maps to this target", local injectivity instead sums "at most one of the
neighbours of some `v` maps to this target", which is exactly what the argument needs (the relevant
pattern vertices are all neighbours of `v`). For the exact-path (distance-2 supplemental) graphs:
dropping the "`q` maps to `t`" term uses the neighbourhood-injectivity of a common neighbour of `p`
and `q` (both are in its neighbourhood), and the insufficient-paths pigeonhole uses the
neighbourhood-injectivity of `p` (the intermediate vertices are all neighbours of `p`). Loopy
instances fall back to plain adjacency + local-injectivity propagation, with the degree/NDS/exact-path
filters disabled ([issue #58]), so those derivations only run loopless. Distance-3 and `--k4` shape
graphs are only built under full injectivity, so they never arise under local injectivity.

[issue #58]: https://github.com/ciaranm/glasgow-subgraph-solver/issues/58

Loops used to be incompatible with supplemental graphs ([issue #56], now fixed). The adjacency
constraint keeps the target's self-loop term (so a loop→loop mapping satisfies the model, see
[issue #49]), but that term is a stray when the constraint is summed into a pseudo-Boolean
derivation. So before any such derivation, `HomomorphismProofs::derive_loop_fixed_adjacencies` derives the
loop-cancelled form of each loop-bearing adjacency constraint — `~x_p_t` together with the neighbours of `t` other
than `t` itself, which follows from the constraint plus injectivity on `t` — and the degree,
supplemental-graph and distance-3 pols sum *that* in its place. The induced encoding additionally
forbids a non-loopy pattern vertex from mapping to a loopy target (the `q == p` case of induced
non-edge preservation), which the model previously left out.

The induced non-edge constraint (for non-adjacent `p`, `q`: if `p` maps to `t` then `q` maps to a
non-neighbour of `t`) takes the *full* set of non-neighbours of `t`, including `t` itself when `t`
has no self-loop. Under full injectivity `q` cannot share `t` with `p`, so leaving `t` out was
harmless; under local injectivity `p` and `q` may both map to a loopless `t`, and since `t` is not
adjacent to itself that is a legitimate induced non-edge, so the model has to keep `t` in the set.

[issue #56]: https://github.com/ciaranm/glasgow-subgraph-solver/issues/56

A minimal worked example:

```shell session
$ ./build/glasgow_subgraph_solver --induced --no-supplementals --no-clique-detection --no-nds \
    --prove myproof --format lad pattern target
$ veripb myproof.opb myproof.pbp
```

Clique and maximum-common-subgraph proof logging also exist (`CliqueParams::proof_options`,
`CommonSubgraphParams::proof_options`), and the common-subgraph reduction extends a clique proof
internally.

## Conclusions

The proof ends with one of the following conclusions, depending on what the solver was asked and
whether the search completed:

| Run | Outcome | Conclusion |
| --- | --- | --- |
| decision | a mapping found | `SAT` (the mapping is logged with one `solx`) |
| decision | no mapping, search exhausted | `UNSAT` |
| counting / enumeration | search exhausted | `ENUMERATION_COMPLETE <n>` |
| counting / enumeration | stopped early (timeout or `--solution-limit`) | `ENUMERATION_PARTIAL <n>` |
| `--minimise-cost` | search exhausted, a mapping found | `BOUNDS c c` (each improvement logged with `soli`) |
| `--minimise-cost` | no mapping at all | `BOUNDS INF INF` (VeriPB does not accept `UNSAT` with an objective) |
| `--minimise-cost` | stopped early | `NONE` |

Each solution is logged with the `solx` rule at the *top* proof level, so the blocking constraint it
introduces survives the `wiplvl` cleanup of the search subtree on backtrack — this is what keeps the
solution count sound. VeriPB checks the claimed count `<n>` against the number of `solx` rules.

Each backtrack nogood is moved into the core (`core id`); on backtracking out of a level, the blocking
constraints (and the now-subsumed deeper core nogoods) recorded at that level are checked-deleted
(`del id`), re-deriving by RUP from the subsuming nogood, and the nogood itself is deleted when we
backtrack past its own level. This matters most for counting, where the per-solution blocking
constraints would otherwise accumulate and make the proof linear in the number of solutions; deleting
them keeps it linear in the search *depth* instead.

## Minimising cost and multigraphs

A multigraph, or a target with costs when minimising, is reified before search (see
`gss/innards/reification.hh`): every edge becomes a vertex. Under `--prove`, minimising always
reifies, even with vertex costs only. The OPB model is **not** written from the reified graphs.
`HomomorphismProofs::emit_reified_model` writes it from the *original* graphs, independently of the
reification code, so that a reification bug shows up as a proof that does not check rather than
being faithfully encoded. The model has:

- a variable `x` for each original pattern vertex and each target vertex with a compatible label,
  with the usual `@al1`, `@am1` and `@inj` constraints. These are exactly the search's values for
  original vertices;
- a variable `z` for each pair of adjacent pattern vertices and each ordered pair of target vertices
  that carries every edge the pattern pair needs (each label, in the right direction), so at most
  one `z` for each pair of a pattern edge and a target edge. A missing `z` is what forbids a pair of
  images: there are no adjacency constraints;
- linking equalities in both directions, `sum_y z(a, b, x, y) = x(a, x)` and likewise for `b`,
  labelled `@lnk...ge` and `@lnk...le`, which make each `z` the conjunction of its two `x`. This is
  the local-polytope encoding, whose linear relaxation is the one the cost bound works in;
- a pattern loop as a requirement on its vertex's image;
- when minimising, an objective of the target costs of the `x` and `z` variables.

Variables are named by index (`xp3_t17`, `zp0_p1_t4_t9`), since vertex names may contain anything.
The search's edge-vertices have no variables: their values follow from the `z`, so no proof line may
mention one. Three things keep it that way, and each is sound on its own terms:

- Branching is on original vertices only, and a solution (`solx` or `soli`) lists only their values.
  The case where search would branch on an edge-vertex, a pattern edge without a label with
  several parallel target edges to choose from, is refused under proof.
- All-different, and the initial Hall check, look only at original vertices; injectivity on
  edge-vertices follows from it.
- Under proof, the simple injectivity propagation skips edge-vertex assignments too. Two pattern
  edges wanting the one target edge is a pigeonhole argument over the original vertices in the
  model, which unit propagation cannot see, so search is made to reach it through a Hall
  violator over the original vertices instead.

The degree and NDS filters, and the whole-instance degree check, are off for a reified instance
under proof, since their derivations cite adjacency constraints this model does not have. So is the
pattern-bigger-than-target refutation, unless the *original* pattern is bigger: a reified pattern can
outnumber a reified target when the originals do not, and the model only has injectivity on
originals. Search refutes those instead. (Reaching search with more pattern vertices than target
vertices exposed a sizing bug in `cheap_all_different`, whose bucket arrays were sized by the
target; they are now sized by the number of domains.)

Each new best mapping is logged with `soli` at the top level, and the cost bound
(`gss/innards/cost_bound.hh`) prunes against the objective-improving constraint it adds. **Each
call of the bound that fails or removes values by the bound emits one `pol`**, which adds up:

- the latest objective-improving constraint, `obj <= U - 1`;
- each original vertex's `@al1` (or `@am1`, for a negative multiplier) times its row potential
  from the assignment step, plus the residual of every pair where it is the smaller vertex;
- each target vertex's `@inj` times its column potential, which is never negative;
- each linking inequality, `ge` or `le` according to sign, times the dual-ascent message on that
  value, plus the pair's residual on the smaller vertex's side.

The result is `sum(-reduced cost * variable) >= bound - (U - 1)`. Every live candidate's
coefficient is minus its reduced cost, never positive, and every other variable with a positive
coefficient is false at the node by propagation. So at the node the constraint conflicts, when
`bound >= U`, or propagates exactly the values with `bound + reduced cost >= U`, which are the
ones the bound removed. Nothing else is written: the search's own backtrack nogoods then follow by
RUP. Values the bound removes because no pair of images supports them need no derivation, since
propagation over the linking equalities finds them. An assignment step with no finite solution at
all writes the Hall violator the Hungarian algorithm's alternating tree gives.

Two things about the bound are there to make that exact. Every original vertex is a row, with a
single candidate once it has a value, and every pattern pair is one pairwise term, so each
quantity in the bound is a multiplier on one model-B constraint. And the assignment step solves
the square problem padded with zero-cost rows, whose dual, shifted, is an exact dual of the
rectangular one; the unpadded Hungarian algorithm's potentials fall short of the assignment's cost
when there are more columns than rows, and cost plus a reduced cost is then not even a valid bound.
The bound is the dual's value, so it is exactly what the `pol` proves. The messages are rounded down
to integers, so every multiplier is an integer.

This was developed by first writing each conclusion as an `a` (assumption) rule, to check that the
rest of the proof holds together, then replacing each with its derivation followed by a RUP of the
same conclusion, and finally dropping those RUPs once they had all checked.

With the bound certified, proofs have the same search as without them. On the graph3 benchmark, the
largest supplied pattern (10 people, 76 edges) proves in 11 nodes with a 3.4 MB proof that
verifies in under three seconds; a random 10-person query (781 nodes) writes 238 MB, almost all of
it these `pol` lines, averaging about 700 terms, and takes three minutes to verify. Their size is
the next thing to work on.

`test-instances/weighted` has fixed instances for each case above, registered as `proof_weighted_*`
and `proof_multigraph_*`. `test-instances/weighted_proof_sweep.py`, registered as
`proof_weighted_random_sweep` when Python is available, checks 150 random instances across
directedness, loops, labels, parallel edges and vertex and edge costs, negative costs included.
Planting a bug in the reifier (a wrong cost on one label, or a reversed direction on one label in
the target only) made the proofs of the affected instances fail to verify, even with the driver's
own cost and solution checks switched off.

## Current status and known limitations

- **Refutation (UNSAT) proofs verify.** Proving that *no* mapping exists works end to end with
  VeriPB 3.0.2 (`s VERIFIED UNSATISFIABLE`).
- **Solution, counting and enumeration proofs verify**, including loop-preserving mappings. Decision
  proofs conclude `s VERIFIED SATISFIABLE`; counting/enumeration proofs conclude
  `s VERIFIED {COMPLETE,PARTIAL} ENUMERATION OF n SOLUTIONS`. The adjacency constraint keeps the
  target self-loop term in its neighbour sum, so a loop→loop solution satisfies the model (this was
  [issue #49], now fixed).
- **Enumeration proof size is linear in the search depth.** Each solution's `solx` blocking
  constraint is checked-deleted once we backtrack out of the level it was found at, by moving the
  subsuming backtrack nogoods into the core (see above). This depends on an upstream VeriPB fix to a
  `core id` monotonicity bug ([issue #59]); an earlier attempt was reverted while that bug was open.

- **Supplemental derivations are emitted lazily, and only the strongest per head.** The adjacency
  constraints for the supplemental graphs nest (`distance3` ⊇ `exact-path-1` ⊇ `exact-path-2` …,
  same head), so only the set-minimal one per head is derived; anything that needs an elided one
  re-derives it as a one-step weakening and deletes it again. On top of that, the kept derivations
  are not emitted during model build but *materialised on first use* — a root degree/NDS read, an
  assignment of `p → t`, or a forward-check removal of `t` from `dom(p)` — so a head search never
  touches is never derived. Together these cut the proof by roughly 2–3x on supplemental-heavy
  configurations, with the OPB and the search tree unchanged. The reasoning for why deferral is
  admissible where omission is not, and the maintenance obligation it creates, is in
  [preprocessor-refactor.md](preprocessor-refactor.md#lazy-supplemental-emission). Pass
  `--no-proof-supplemental-subsumption` to emit everything, for studying the effect.

[issue #49]: https://github.com/ciaranm/glasgow-subgraph-solver/issues/49
[issue #59]: https://github.com/ciaranm/glasgow-subgraph-solver/issues/59

## Verifying with the CakePB checker

The proof can also be checked end to end by the formally verified CakePB checker, `cake_pb_iso`.
VeriPB *elaborates* the user-friendly proof down to a kernel subset, checking it against cake's own
OPB encoding (which it derives directly from the LAD files); `cake_pb_iso` then checks that
elaborated proof:

```shell session
$ ./build/glasgow_subgraph_solver --format lad --no-clique-detection --prove myproof pattern target
$ cake_pb_iso pattern target > myproof.cakeopb
$ veripb myproof.cakeopb myproof.pbp --elaborate myproof.core.pbp
$ cake_pb_iso pattern target myproof.core.pbp
```

This checks the solver's proof against cake's own (independently, formally verified) OPB encoding
rather than the solver's, giving an end-to-end formally verified result.

## Tests

The `ctest` suite includes proof-verification tests (registered only when `veripb` is on the
`PATH`): they run a solver with `--prove` on small instances and check the proof with VeriPB. Each
test states the exact feature combination it exercises (the harness adds no flags of its own), so
between them they trigger every conditional proof-writing path: refutation, decision, complete and
partial enumeration; the supplemental-graph, distance-3, neighbourhood-degree-sequence and
clique-size derivations; the clique and common-subgraph solvers; and loop-preserving instances. See
`src/CMakeLists.txt` and `test-instances/verify_proof.bash`.

If `cake_pb_iso` is found (point CMake at it with `-DCAKE_PB_ISO_EXECUTABLE=/path/to/cake_pb_iso`),
the suite additionally registers `cake_*` tests that run the whole verified pipeline above, including
loop cases. See `test-instances/verify_cake_pipeline.bash`.

```shell session
$ ctest --preset release -R proof
```

### Coverage

To check that the tests actually trigger each conditional path in the proof-writing code, build with
the `coverage` preset and run the `coverage` target, which runs the suite under `gcov` instrumentation
and reports branch/decision coverage of `gss/innards/proof.cc` (requires `gcovr`):

```shell session
$ cmake --preset coverage && cmake --build --preset coverage
$ cmake --build build-coverage --target coverage   # writes build-coverage/proof-coverage.{txt,html}
```

CI runs this and fails if the coverage regresses, so a new conditional proof feature has to arrive
with a test that triggers it. (The one path no test currently reaches is the connected-clique
connectivity backtrack, `backtrack_from_binary_variables`.)
