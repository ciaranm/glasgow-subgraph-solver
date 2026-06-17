# Architecture

A developer's-eye map of the Glasgow Subgraph Solver: how the code is laid out, what the main
components are, and how they fit together. For *using* the solvers see the [README](../README.md);
for the proof-logging machinery see [proof-logging.md](proof-logging.md).

## What's in the box

One static library (`gss/`, target `glasgow_subgraphs`) and a handful of command-line drivers
(`src/`). The library solves three related problems:

| Problem | Public entry point | Header |
| --- | --- | --- |
| Subgraph isomorphism / homomorphism | `solve_homomorphism_problem` | `gss/homomorphism.hh` |
| Clique (decision / maximum) | `solve_clique_problem` | `gss/clique.hh` |
| Maximum common (connected) subgraph | `solve_common_subgraph_problem` | `gss/common_subgraph.hh` |

The drivers are thin: they parse options with `cxxopts`, read graphs, fill in a `*Params` struct,
call the matching `solve_*` function, and print the `*Result`.

## Layers

```
            src/glasgow_subgraph_solver.cc   glasgow_clique_solver.cc   glasgow_common_subgraph_solver.cc
                          │                          │                          │   (also create_random_graph, convert_to_lad)
                          ▼                          ▼                          ▼
  ┌───────────────────────────────────────────────────────────────────────────────────────┐
  │ Public API  (gss/*.hh)                                                                   │
  │   homomorphism · clique · common_subgraph · sip_decomposer                               │
  │   HomomorphismParams / CliqueParams / CommonSubgraphParams   +   *Result structs         │
  │   restarts · timeout · value_ordering · loooong · vertex_to_vertex_mapping               │
  └───────────────────────────────────────────────────────────────────────────────────────┘
                          │ uses
                          ▼
  ┌───────────────────────────────────────────────────────────────────────────────────────┐
  │ Implementation detail  (gss/innards/*)                                                    │
  │   homomorphism_model · homomorphism_searcher · homomorphism_domain · homomorphism_traits  │
  │   cheap_all_different · graph_traits · watches (nogoods) · svo_bitset                      │
  │   proof (VeriPB logging) · verify · lackey (external solver) · symmetries (GAP) · threads  │
  └───────────────────────────────────────────────────────────────────────────────────────┘
                          │ uses
                          ▼
  ┌───────────────────────────────────────────────────────────────────────────────────────┐
  │ Graph input  (gss/formats/*)            Utilities  (gss/utils/*)                          │
  │   InputGraph  +  readers: csv, lad, dimacs, vfmcs, read_file_format                       │
  │   graph_file_error                       hashing_utils · vertex_name_map                   │
  └───────────────────────────────────────────────────────────────────────────────────────┘
```

The convention is that everything under `gss/innards/` is private to the library — callers should
only need the headers directly under `gss/`. A few innards types do currently leak through the public
API (the `Lackey` handle and `ProofOptions` on the `*Params` structs, and `SVOBitset` in a
`CliqueParams` callback); these are deliberate extension points rather than accidents.

## The homomorphism solver in more detail

This is the largest and most general engine; the subgraph isomorphism and (non-)injective
homomorphism variants are all configurations of it.

- **`HomomorphismModel`** (`innards/homomorphism_model.{hh,cc}`) turns the two `InputGraph`s plus the
  `HomomorphismParams` into the internal constraint model: it re-encodes the graphs as bitset
  adjacency, builds the *supplemental graphs* (exact-path / distance-3 / k4 / extra shapes) used for
  stronger filtering, precomputes degrees and neighbourhood-degree sequences, and seeds the initial
  variable domains. The heavy lifting happens in `prepare()`.
- **`HomomorphismDomain`** (`innards/homomorphism_domain.hh`) is one CP variable's domain: an
  `SVOBitset` of still-possible target vertices plus bookkeeping.
- **`HomomorphismSearcher`** (`innards/homomorphism_searcher.{hh,cc}`) is the backtracking engine:
  variable/value ordering, constraint propagation (adjacency, all-different via
  `cheap_all_different`, less-than / occurs-less symmetry constraints, optional lackey propagation),
  restarts, and nogood recording via `Watches`. The propagation hot path is templated on
  `<directed, has_edge_labels, induced, verbose_proofs>` so the per-node inner loop has no runtime
  branches on those flags — a deliberate performance choice.
- **`solve_homomorphism_problem`** (`homomorphism.cc`) wires these together: it sets up optional
  proof logging, applies cheap early-exit heuristics (pattern bigger than target, the loop shortcut,
  clique detection), builds the model, and dispatches to a sequential or threaded solver.

`sip_decomposer` offers an alternative top level that solves subgraph isomorphism by decomposing the
pattern into biconnected components.

## The clique and common-subgraph solvers

- **Clique** (`clique.{hh,cc}`) is a branch-and-bound maximum-clique solver using a greedy colouring
  as both the bound and the branching heuristic, over `SVOBitset` adjacency, with the same
  `Watches`-based nogoods and restarts. `CliqueParams::decide` switches it to a decision problem.
- **Common subgraph** (`common_subgraph.{hh,cc}`) finds the maximum common *induced* subgraph by a
  partition-refinement branch and bound. It can optionally be solved by reducing to a clique on the
  association graph (`CommonSubgraphParams::clique`), and supports a connected variant. The clique
  solver is reused for the reduction, which is why `clique.hh` exposes the `connected` callback and
  the proof-extension hooks.

## Cross-cutting pieces

- **`InputGraph`** (`formats/input_graph.{hh,cc}`) is the format-agnostic graph the readers produce
  and the solvers consume: a pimpl over an edge map, with vertex names/labels and edge labels. It is
  deliberately not performance-tuned — the solvers re-encode it into bitsets.
- **Format readers** (`formats/`) all take a `std::istream` and a filename and return an
  `InputGraph`; `read_file_format` dispatches by name and can auto-detect. Supported: CSV, LAD,
  directed/vertex-labelled/labelled LAD, DIMACS, and VFMCS.
- **`SVOBitset`** (`innards/svo_bitset.hh`) is the small-vector-optimised bitset used everywhere for
  domains and adjacency: inline storage for up to 16 64-bit words, heap beyond that.
- **`loooong`** (`gss/loooong.hh`) is a thin GMP `mpz_t` wrapper for solution counts, which overflow
  64 bits readily.
- **`Watches`** (`innards/watches.hh`) is a generic two-watched-literal nogood store shared by the
  homomorphism and clique searchers.
- **`RestartsSchedule`** (`gss/restarts.hh`) and **`Timeout`** (`gss/timeout.hh`) are the search
  control knobs. Restart policies: none, Luby, geometric, timed, and a thread-synchronised variant.
- **`Proof`** (`innards/proof.{hh,cc}`) emits the VeriPB model (`.opb`) and proof log (`.pbp`).
- **`Lackey`** (`innards/lackey.{hh,cc}`) talks to an external constraint solver over named pipes.
- **`find_symmetries`** (`innards/symmetries.{hh,cc}`) shells out to the GAP computer algebra system
  to find pattern/target symmetries.

## Build and tests

CMake with three presets — `release`, `debug`, and `sanitize` (ASan + UBSan) — see
`CMakePresets.json`. The library and drivers need GMP; Catch2 and cxxopts are fetched via
`FetchContent`. Unit tests live next to the code they cover (`gss/**/<name>_test.cc`, registered in
`gss/CMakeLists.txt`) and run under `ctest`. Proof-verification tests (`src/CMakeLists.txt`,
`test-instances/verify_proof.bash`) run the solver under VeriPB and are only registered when `veripb`
is found. `run-tests.bash` is a small end-to-end smoke test over the binaries.
