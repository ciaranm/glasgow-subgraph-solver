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
| no lackey (`--send-to-lackey` / `--receive-from-lackey`) | external propagation is not logged |
| no less-than / occurs-less symmetry constraints | not yet logged |
| injective or non-injective only (not `--locally-injective`) | only these two are encoded |
| unlabelled graphs (no vertex or edge labels) | labels are not yet encoded |

The README's worked example also passes `--no-supplementals` and `--no-nds`; that is the
known-good combination this project tests against:

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
introduces survives the deletions that clean up the search subtree on backtrack — this is what keeps
the solution count sound. VeriPB checks the claimed count `<n>` against the number of `solx` rules.

## Current status and known limitations

- **Refutation (UNSAT) proofs verify.** Proving that *no* mapping exists works end to end with
  VeriPB 3.0.2 (`s VERIFIED UNSATISFIABLE`).
- **Solution, counting and enumeration proofs verify**, including loop-preserving mappings. Decision
  proofs conclude `s VERIFIED SATISFIABLE`; counting/enumeration proofs conclude
  `s VERIFIED {COMPLETE,PARTIAL} ENUMERATION OF n SOLUTIONS`. The adjacency constraint keeps the
  target self-loop term in its neighbour sum, so a loop→loop solution satisfies the model (this was
  [issue #49], now fixed).

[issue #49]: https://github.com/ciaranm/glasgow-subgraph-solver/issues/49

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
`PATH`): they run the solver with `--prove` on small instances and check the proof with VeriPB,
covering refutation, decision, complete enumeration and partial enumeration, for both loopless and
loop-preserving instances. See `src/CMakeLists.txt` and `test-instances/verify_proof.bash`.

If `cake_pb_iso` is found (point CMake at it with `-DCAKE_PB_ISO_EXECUTABLE=/path/to/cake_pb_iso`),
the suite additionally registers `cake_*` tests that run the whole verified pipeline above, including
loop cases. See `test-instances/verify_cake_pipeline.bash`.

```shell session
$ ctest --preset release -R proof
```
