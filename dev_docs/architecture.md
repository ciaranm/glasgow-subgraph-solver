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
                          │                          │                          │   (also create_random_graph, convert_to_lad, convert_to_json)
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
  │   solve_state · supplemental_graphs · processed_graphs_data · clique_size_constraints      │
  │   cheap_all_different · graph_traits · watches (nogoods) · svo_bitset · filter_activations │
  │   homomorphism_proofs · proof (VeriPB logging) · verify · threads                          │
  └───────────────────────────────────────────────────────────────────────────────────────┘
                          │ uses
                          ▼
  ┌───────────────────────────────────────────────────────────────────────────────────────┐
  │ Graph input  (gss/formats/*)            Utilities  (gss/utils/*)                          │
  │   InputGraph  +  readers: csv, json, lad, dimacs, vfmcs, read_file_format                 │
  │   graph_file_error                       hashing_utils · vertex_name_map                   │
  └───────────────────────────────────────────────────────────────────────────────────────┘
```

The convention is that everything under `gss/innards/` is private to the library — callers should
only need the headers directly under `gss/`. A few innards types do currently leak through the public
API (`ProofOptions` on the `*Params` structs, and `SVOBitset` in a `CliqueParams` callback); these
are deliberate extension points rather than accidents.

## The homomorphism solver in more detail

This is the largest and most general engine; the subgraph isomorphism and (non-)injective
homomorphism variants are all configurations of it.

- **`HomomorphismModel`** (`innards/homomorphism_model.{hh,cc}`) turns the two `InputGraph`s plus the
  `HomomorphismParams` into the internal constraint model: it re-encodes the graphs as bitset
  adjacency, owns the *supplemental graphs* (exact-path / distance-3 / k4 / extra shapes) used for
  stronger filtering, precomputes degrees and neighbourhood-degree sequences, and seeds the initial
  variable domains. `prepare()` drives this; the graph construction itself is free functions in
  `innards/supplemental_graphs.{hh,cc}` writing into the `ProcessedGraphsData`
  (`innards/processed_graphs_data.hh`) the model hands them, from one `ShapeGraphSpec` plan that is
  the single source of truth for which slot holds what. `build_supplemental_graphs()` and
  `tighten_domains_with_supplementals()` are split out so staged solving can defer them past a first
  search round. Clique-size constraints live in `innards/clique_size_constraints.{hh,cc}`.
- **`HomomorphismDomain`** (`innards/homomorphism_domain.hh`) is one CP variable's domain: an
  `SVOBitset` of still-possible target vertices plus bookkeeping.
- **`HomomorphismSearcher`** (`innards/homomorphism_searcher.{hh,cc}`) is the backtracking engine:
  variable/value ordering, constraint propagation (adjacency, all-different via
  `cheap_all_different`, less-than / occurs-less symmetry constraints),
  restarts, and nogood recording via `Watches`. The propagation hot path is templated on
  `<directed, has_edge_labels, induced, track_removals>` so the per-node inner loop has no runtime
  branches on those flags — a deliberate performance choice, and a measurable one: attributing a
  removal to the graph pair that made it costs a popcount per graph pair, so the two things that
  want that attribution (verbose proof comments and filter-activation recording) share the one
  instantiation rather than adding a branch.
- **`solve_homomorphism_problem`** (`homomorphism.cc`) wires these together as a **pipeline of
  `SolveStep`s** over a shared `SolveContext`, run in registration order and stopping at the first
  step that concludes: emit the OPB model, the pattern-bigger-than-target refutation, the
  target-loop shortcut, the clique-pattern reduction, and then `MainSolveStep`, which builds the
  model and searches. Search is simply the terminal step — there is no preprocess/search
  distinction in the control flow. `SolveState` (`innards/solve_state.hh`) is what the pipeline
  carries: the model, the root domains, and the nogood store, so that steps can grow the model and
  accumulate nogoods between rounds. `--staged` uses exactly that: a cheap first round (original
  graph only, degree + Hall, no NDS or supplementals) under a bounded restart schedule, and only if
  it does not conclude are the supplemental graphs built, the domains re-filtered, and the search
  resumed unbounded, with the nogoods carried across. Sequential only; see
  [preprocessor-refactor.md](preprocessor-refactor.md).

- **Multigraphs and costs.** `solve_homomorphism_problem` first **reifies** a multigraph, or a
  target with edge costs when minimising (`innards/reification.{hh,cc}`): each edge becomes a
  vertex of its own, so the pipeline above only ever sees simple graphs. `--minimise-cost` adds
  **`CostBound`** (`innards/cost_bound.{hh,cc}`), which the searcher runs in `propagate()`: it
  rebuilds pairwise costs from the edge-vertices' domains, tightens them by a local-polytope
  dual ascent, and combines them with a minimum-cost assignment whose reduced costs remove
  values. Every mapping search reaches is a new incumbent. See
  [option-compatibility.md](option-compatibility.md) for what this can be combined with.

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
  deliberately not performance-tuned — the solvers re-encode it into bitsets. Directedness is
  *declared* to the constructor, not inferred from which edges get added: `add_directed_edge`
  requires it, and `add_edge` keeps an undirected graph undirected however its edges are labelled.
  See [file-formats.md](file-formats.md) for why that matters.
- **Format readers** (`formats/`) all take a `std::istream` and a filename and return an
  `InputGraph`; `read_file_format` dispatches by name and can auto-detect. Supported: CSV, LAD,
  directed/vertex-labelled/labelled LAD, DIMACS, VFMCS, and the `gss-graph` JSON format. Only the
  last of those declares its graph-level properties rather than inferring them, and only it has a
  writer; see [file-formats.md](file-formats.md).
- **`SVOBitset`** (`innards/svo_bitset.hh`) is the small-vector-optimised bitset used everywhere for
  domains and adjacency: inline storage for up to 16 64-bit words, heap beyond that.
- **`loooong`** (`gss/loooong.hh`) is a thin GMP `mpz_t` wrapper for solution counts, which overflow
  64 bits readily.
- **`Watches`** (`innards/watches.hh`) is a generic two-watched-literal nogood store shared by the
  homomorphism and clique searchers.
- **`FilterActivations`** (`innards/filter_activations.hh`) counts what each filter actually
  removed, when `HomomorphismParams::record_filter_activations` asks for it, and reports it in
  the result's extra stats. Nothing in a normal solve touches it; it is what lets the option
  sweep tell a filter that is correct here from one that did nothing here. See
  [option-compatibility.md](option-compatibility.md).
- **`RestartsSchedule`** (`gss/restarts.hh`) and **`Timeout`** (`gss/timeout.hh`) are the search
  control knobs. Restart policies: none, Luby, geometric, timed, and a thread-synchronised variant.
- **`Proof`** (`innards/proof.{hh,cc}`) emits the VeriPB model (`.opb`) and proof log (`.pbp`). It
  holds the generic pseudo-Boolean machinery — variable naming, model and proof line emission,
  levels, the dedup caches — plus the derivations the clique and common-subgraph solvers share.
  **`HomomorphismProofs`** (`innards/homomorphism_proofs.{hh,cc}`) sits between the homomorphism
  solver and `Proof` and owns everything homomorphism-exclusive: the OPB model emission, the
  adjacency / exact-path / distance-3 / extra-shape derivations, the degree / NDS / Hall filter
  proofs, and the proof-size economies (subsumption elision, deferral, and the lazy
  materialisation the searcher drives). Keeping it out of `Proof` is what stops the bottom layer
  from having to know what an exact-path graph is.

## Build and tests

CMake with four presets — `release`, `debug`, `sanitize` (ASan + UBSan) and `coverage` (gcov) — see
`CMakePresets.json`. The library and drivers need GMP; Catch2, cxxopts and nlohmann/json are fetched
via `FetchContent` when not already installed. Unit tests live next to the code they cover (`gss/**/<name>_test.cc`, registered in
`gss/CMakeLists.txt`) and run under `ctest`. Proof-verification tests (`src/CMakeLists.txt`,
`test-instances/verify_proof.bash`) run the solver under VeriPB and are only registered when `veripb`
is found. `run-tests.bash` is a small end-to-end smoke test over the binaries.

Two of the tests are sweeps rather than fixed cases, and between them they are what most of the
correctness confidence rests on. `random_homomorphism_test` checks the solver against a
brute-force oracle, exhaustively over the problem axes, on all sixteen combinations of loops ×
directed × vertex labels × edge labels. `option_sweep_test` needs no oracle — its reference is
the same instance solved with the filtering off — so it can afford instances big enough for the
filters to fire, and runs a pairwise covering array over the option space against nineteen
instance families, checking its per-cell outcomes against `gss/option_sweep_golden.tsv`. Which
option combinations are legal, which are silently disabled, and what neither sweep covers is in
[option-compatibility.md](option-compatibility.md).
