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
derivation. So before any such derivation, `Proof::loop_fix_adjacencies` derives the loop-cancelled
form of each loop-bearing adjacency constraint — `~x_p_t` together with the neighbours of `t` other
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

Each solution is logged with the `solx` rule at the *top* proof level, so the blocking constraint it
introduces survives the `wiplvl` cleanup of the search subtree on backtrack — this is what keeps the
solution count sound. VeriPB checks the claimed count `<n>` against the number of `solx` rules.

Each backtrack nogood is moved into the core (`core id`); on backtracking out of a level, the blocking
constraints (and the now-subsumed deeper core nogoods) recorded at that level are checked-deleted
(`del id`), re-deriving by RUP from the subsuming nogood, and the nogood itself is deleted when we
backtrack past its own level. This matters most for counting, where the per-solution blocking
constraints would otherwise accumulate and make the proof linear in the number of solutions; deleting
them keeps it linear in the search *depth* instead.

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
