# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code
in this repository.

It is deliberately short. Almost everything a developer needs here is project
documentation rather than agent guidance, and it lives in `dev_docs/`, `README.md`
or `CONTRIBUTING.md`. A fact copied into two files goes stale in one of them: an
earlier version of this file described a `lackey.cc` and a `symmetries.cc` that do
not exist, a GAP integration the README says was removed, a CMake version that was
wrong, and Boost.Program_options in place of cxxopts — all of it plausible, none of
it true. So this file routes rather than restates: what to read, the mistakes that
have actually been made in this tree, and what to run before committing.

## Read the relevant document first

A C++20 solver for subgraph isomorphism, maximum clique and maximum common subgraph:
backtracking search with constraint propagation, which can also emit a
VeriPB-checkable proof of its answer. That proof-logging ability shapes a lot of the
design and very little of it is guessable from the code alone. Read the document
covering what you are about to change *before* you change it.

| If you are... | Start with |
|---|---|
| getting oriented anywhere in the tree | [`dev_docs/architecture.md`](dev_docs/architecture.md) — the layer map, and what each component is for |
| touching a format reader, a writer, or `InputGraph` | [`dev_docs/file-formats.md`](dev_docs/file-formats.md) |
| writing or debugging proof logging | [`dev_docs/proof-logging.md`](dev_docs/proof-logging.md) — including which option combinations are supported, which is not all of them |
| touching an option, a filter, or the traits layer | [`dev_docs/option-compatibility.md`](dev_docs/option-compatibility.md) — which combinations are legal, which are silently disabled rather than refused, and what the sweeps cover |
| changing anything the solver does before search | [`dev_docs/preprocessor-refactor.md`](dev_docs/preprocessor-refactor.md) — where that code is going, which is not where it is |
| looking for what an option does, or a file format | [`README.md`](README.md) |
| about to commit | [`CONTRIBUTING.md`](CONTRIBUTING.md) |

`gss/*.hh` is the public API — `HomomorphismParams` and `solve_homomorphism_problem`
are the main entry point, with `clique.hh` and `common_subgraph.hh` alongside.
`gss/innards/` is everything that is not part of it. `src/` is the command-line
drivers, which use cxxopts.

## Mistakes that have been made here before

Each of these is a real bug or a real wasted afternoon, not a hypothetical. The
reasoning is in the linked document; the instruction is here so it is visible
without following the link.

- **Do not infer a graph-level property from how a graph happened to be built.**
  Adding a label to an undirected CSV file used to make it directed, because the
  labelled path went through `add_directed_edge`, which set the flag (#86, #87).
  Directedness is now declared to the `InputGraph` constructor and
  `add_directed_edge` requires it. See
  [`file-formats.md`](dev_docs/file-formats.md).
- **After changing what a reader reports, grep for everything that branches on it.**
  Fixing the above broke the unwritten invariant that every edge-labelled graph was
  also directed, and two callers had been relying on it: the model allocated the
  forward/reverse target rows only for directed patterns while the searcher always
  reads them once there are edge labels (a segfault), and the SIP decomposer
  discarded every edge label when rebuilding its reduced pattern (silently no
  solutions). Both were unreachable before and neither was caught by a type error.
- **Do not compare user-supplied label text against magic strings.** A filter
  skipping target edges labelled `"unlabelled"` outlived the representation it was
  guarding by five years, and turned a satisfiable instance into `status = false`
  for anyone who used that word as a label (#88).
- **Loops and the induced / locally-injective modes are where the bugs live.** Four
  separate ones: the self-loop adjacency term dropped from solution proofs (#49),
  supplemental-graph derivations not verifying on loopy graphs (#56),
  locally-injective enumeration over-pruning (#58), and the loop/clique shortcuts
  being wrong for induced mappings into loopy targets. When you touch propagation or
  a proof derivation, test a loopy instance and the induced and locally-injective
  combinations, not just the default path. `test-instances/small` and `large` have a
  self-loop for this reason.
- **Do not "tidy" the `FORCE` off the per-configuration flags in `CMakeLists.txt`.**
  CMake pre-creates `CMAKE_CXX_FLAGS_<CONFIG>` as empty cache entries during
  compiler detection in `project()`, so a plain `set(... CACHE ...)` there is a
  silent no-op and the build type gets none of its flags — which for `sanitize`
  means a build containing no sanitizers that still passes.
- **The empty `//` markers in the cxxopts option blocks in `src/` are
  load-bearing.** Each option block is one chained expression, and with
  `ColumnLimit: 0` clang-format packs it onto a single line hundreds of characters
  long; a trailing comment cannot be joined with what follows it, so the marker pins
  the break. Keep them, and see the comment in `.clang-format`.
- **A filter that needs k *distinct* images needs injectivity, or a loopless target.**
  Nothing collapses two adjacent pattern vertices onto one image unless that image
  carries a self-loop — so the clique-size filter and the clique reduction are sound
  without injectivity right up until a target loop turns up, which is why both went
  years without the premise being stated (#91, #94). The supplemental graphs need a
  different premise again. Before adding a filter, write its condition down in
  [`option-compatibility.md`](dev_docs/option-compatibility.md) next to the others.
- **A builder that reads adjacency rows must ask whether *either* graph is directed.**
  A row is a vertex's out-neighbours, so anything written for undirected rows is wrong
  on an asymmetric one, whichever graph it came from: the k4 and extra-shape builders
  marked pairs in an orientation they never checked (#97, #99), the clique solver held
  "cliques" together with one-way arcs (#93), and local injectivity asked about common
  successors when it meant common predecessors (#96). `ProcessedGraphsData::directed` follows the *pattern*,
  because that is what picks the searcher's propagation path; `either_graph_directed` is
  the question a builder wants.
- **The randomised oracle now covers directed and labelled graphs, and the sweep covers
  the filtering options.** `random_homomorphism_test.cc` enumerates all sixteen
  combinations of loops × directed × vertex labels × edge labels against a brute-force
  verifier, and `option_sweep_test.cc` multiplies a pairwise array over the options by
  nineteen instance families. Green there is real evidence — for everything except
  threads, proof logging (covered by `test-instances/random_proof_sweep.bash`),
  `--shape`, and instances whose pattern and target have different shapes. For those,
  add a fixed instance.
- **A stale build can lie to you, particularly when checking that a test catches a
  bug.** Reverting a fix with `git stash push <file>` and rebuilding does not
  reliably recompile the restored file, and `touch` is not always enough either; the
  reliable move is to delete the object file
  (`build/gss/CMakeFiles/glasgow_subgraphs.dir/<path>.cc.o`) or to reconfigure.
  This cost two rounds of confusion in one session, once presenting a working fix as
  broken and once the reverse. Verifying a fix by reverting it is a good habit — just
  confirm the compiler actually ran.

## Before committing

[`CONTRIBUTING.md`](CONTRIBUTING.md) is the contract: the policy on AI-assisted
contributions (declare it, with a `Co-Authored-By:` trailer naming the tool), the
formatting requirement, and what a contribution should pass. Work on a branch and
open a pull request; do not commit to `main`.

```shell
# Format. CI checks this, and pins clang-format 21.1.8 because major versions
# disagree; `pip install clang-format==21.1.8` gets exactly it.
git ls-files '*.cc' '*.hh' | xargs clang-format -i

j=$(nproc 2>/dev/null || sysctl -n hw.logicalcpu)

cmake --preset release  && cmake --build --preset release  && ctest --preset release  -j $j
cmake --preset sanitize && cmake --build --preset sanitize && ctest --preset sanitize -j $j

# Not registered with ctest, so neither the presets above nor CI run it. Needs the
# release preset's build/ directory, which is the one it looks in.
./run-tests.bash
```

Run a single test binary directly (`./build/homomorphism_test`, and Catch2 filtering
works: `./build/homomorphism_test "*loop*"`, `-l` to list). The proof-verification
tests are only registered when `veripb` is on the `PATH`, so a green `ctest` on a
machine without it has checked considerably less than it looks.

For a change that is meant to be a pure refactor, `test-instances/proof_metrics.bash`
emits a TSV of supplemental-graph counts, proof sizes, solution counts and node
counts over a fixed matrix of option combinations; diffing it against a `main`-built
binary is how this tree checks that nothing moved.
[`preprocessor-refactor.md`](dev_docs/preprocessor-refactor.md) ("Reproducing the
guardrails") has the two-worktree recipe.
