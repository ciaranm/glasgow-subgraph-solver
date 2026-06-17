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

Useful companion flags:

- `--verbose-proofs` writes extra `*` comment lines into the log, for tracing.
- `--recover-proof-encoding` emits the encoding in the form expected by verified (CakeML) encoders.

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

## Current status and known limitations

- **Refutation (UNSAT) proofs verify.** Proving that *no* mapping exists works end to end with
  VeriPB 3.0.2. The example above (an unsatisfiable induced instance) verifies
  `s VERIFIED UNSATISFIABLE`.
- **Solution (SAT) proofs do not currently verify.** There are two independent open bugs, both
  tracked in the [GitHub issues](https://github.com/ciaranm/glasgow-subgraph-solver/issues):
  1. The solver logs solutions with the VeriPB `solx` rule, which is the *enumeration* rule and
     requires the OPB to declare a `preserved:` set of variables — which is not emitted. A plain
     decision proof should use `sol` instead.
  2. The adjacency-constraint encoding drops the loop→loop term, so a valid
     loop-preserving solution falsifies the model.

  These are "not yet implemented / to be fixed" issues, not fundamental restrictions, and are being
  reworked in coordination with the verified CakeML encoder so the two encodings agree on loops and
  solution logging.

## Tests

The `ctest` suite includes proof-verification tests (registered only when `veripb` is on the `PATH`):
they run the solver with `--prove` on small instances and check the proof with VeriPB. The currently
broken solution proofs are registered as `WILL_FAIL` tripwires, so the suite stays green while the
bugs exist and will flag the moment they are fixed. See `src/CMakeLists.txt` and
`test-instances/verify_proof.bash`.

```shell session
$ ctest --preset release -R proof
```
