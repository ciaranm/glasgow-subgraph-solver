# GSS bitset PR benchmark kit (temporary)

Not for merging. This directory exists so that other machines can measure #135, #136 and #137,
plus a rejected fused-popcount experiment, on CPUs other than the Zen 3 they were developed on.
The tracking issue, #138, says what to run and what to post.

Compares seven builds of the solver, each compiled with GCC and with clang (Release, so
`-O3 -march=native`), at pinned commits (see the top of `run-kit.bash`):

| arm | what |
|---|---|
| `main` | `main` at `01ef6ca` |
| `pr1` | #135, `svo-bitset-iteration` (iteration API) |
| `pr2` | #136, `cheap-all-different-hall-early-stop` |
| `pr3` | #137, `recount-only-filtered-domains`: recount a domain only when adjacency filtered it |
| `v1` | `main` + `fused-variants.patch`, `-DGSS_COUNT_VARIANT=1`: apply all row intersections, then `count()` |
| `v2` | same, variant 2: the last intersection is fused with the popcount |
| `v3` | same, variant 3: one pass over the words does every intersection and the popcount |

The `v` arms are the fused-popcount experiment. It lost on Zen 3 (`v3` was 10–12% slower than
`v1`) and is not proposed as a PR. It is here because the `popcnt` false output dependency on
some Intel cores could change the answer. `v1` already contains a form of the recount change
(close to `pr3`), so `v2` and `v3` against `v1` measure the fusion itself.

## Files

- `run-kit.bash [REPS]`: clones the repository into `work/`, builds the fourteen binaries,
  then runs the 59 instances listed in `instances.txt`, one solver at a time. Within each
  instance it rotates the order of the arms. It is resumable. `SMOKE=1` does a quick
  three-instance, one-rep check into a separate results file. Other knobs are documented at
  the top of the script.
- `analyse-campaign.py RESULTS META [metrics]`: checks that every run of an instance did
  identical work (same nodes and status), then reports geometric-mean ratios against `main`
  (and of `v2`/`v3` against `v1`), per instance family.
- `make-comment.py RESULTS MACHINE META`: the markdown comment to post on the issue.
- `instances/`: 93 graph files (6.5 MB) from the standard subgraph isomorphism benchmark
  families: LV, images, meshes, scale-free, si, and phase-transition.
  `instances-meta.txt` gives each one's family and its runtime on `main` on a 2.3 GHz Zen 3,
  in ms.
- `fused-variants.patch`: the experimental variants, applied to `main` for the `v` arms.

Timing on the development machine (EPYC 7643, boost off, 2.3 GHz): about 3.5 hours for
`./run-kit.bash 4`, plus about ten minutes of building.
