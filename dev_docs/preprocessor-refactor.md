# Refactoring the solve pipeline

A design plan for restructuring everything the homomorphism solver does "before search"
into a uniform, reorderable pipeline. This is forward-looking — it describes where the
code is going, not where it is. For the current architecture see
[architecture.md](architecture.md); for proof logging see
[proof-logging.md](proof-logging.md). Scope is the **homomorphism** solver only; the
clique and common-subgraph solvers are deliberately left alone.

## Why

Four things are awkward today, and they share a root cause — preprocessing is a fixed,
implicit sequence split across four sites with proof emission inlined into it:

- **S1 — adding a feature is a 5–7 site change.** A new supplemental graph or filter
  touches `calculate_n_shape_graphs`, a `_build_*` method, a `supports_*` trait, an
  `if`-block in `prepare()` with manual `next_pattern_supplemental` / `next_target_supplemental`
  index threading, a parallel inline proof loop, the slot-index sanity check, and the
  `_check_*` chain plus its NDS-style proof second pass. There is no single registration
  point.
- **S2 — proofs are bigger than they need to be.**
  - S2a: supplemental-graph derivations are emitted to the proof for *every* `(g, p, q, t)`
    triple regardless of whether anything ever filters on that edge, plus the eager
    loop-incompatibility block and `loop_fix_adjacencies` over every loopy adjacency.
  - S2b: the full OPB model is emitted before any refutation check runs — even
    `pattern > target` writes the entire model first.
- **S3 — no reordering or staging.** `prepare()` builds every applicable supplemental
  unconditionally; there is no way to reorder steps, nor to "run fast preprocessing,
  search for a bit, and only then run full preprocessing".
- **S4 — the sequence is hard to follow.** Control flow is spread across the
  `HomomorphismModel` constructor, `solve_homomorphism_problem`, `prepare()`, and
  `initialise_domains()`, with proof scaffolding inlined into the build loops
  (`prepare()` is ~320 lines, about half of it proof emission).

A telling observation: the current code *already* concludes the problem during
"preprocessing" — `global_degree_is_preserved` Hall checking returns UNSAT, and the
loop shortcut and clique reduction return SAT — they are just special-cased as early
returns in `solve_homomorphism_problem`. The preprocess/search distinction is artificial.

## Two economies that are NOT available

Before the design, two tempting ideas that do not work, so we do not attempt them:

- **Minimal / lazy OPB — no.** The OPB is the trusted encoding. A checker has no
  independent way to confirm that a partial OPB is still a faithful, complete model of
  the instance, so omitting any model constraint is indistinguishable from cheating
  (silently making the instance over-permissive). The OPB always contains the complete
  encoding — adjacency, (local-)injectivity, domains — no more and no less. (Counterpart
  rule: do not *add* derived consequences to the OPB either; derive them in the proof.)
- **Lazy-on-cite proof emission — no.** RUP steps record nothing about which constraints
  they consume; VeriPB unit-propagates over the whole database. Once a derived constraint
  exists, any later RUP — every search nogood included — can silently rely on it. We can
  never prove a derived constraint is unused, so we cannot defer or omit it on demand.

So **S2a is achieved by static subsumption / inertness elimination**: do not *generate*
the useless supplementals in the first place. This is decidable, result-preserving, and
a search-speed win (fewer model bitset rows), not only a smaller proof:

- A whole supplemental-graph pair `(Pˢ, Gˢ)` is inert — skip it entirely, no slot, no
  build, no derivations — when, by direct bitset checks: `Gˢ` is complete (the adjacency
  constraint is vacuous and degree/NDS on it are maximal so cannot fire), or `Pˢ` is
  empty (no constraints generated), or `(Pˢ, Gˢ)` duplicates an already-registered slot.
  Dense instances are where this bites: distance-3 / exact-path graphs often saturate to
  complete (maximal proof, zero filtering), and duplicates between e.g. distance-2 and
  exact-path-1 are common.
- Within a graph we do keep, do not derive a constraint subsumed by one already derived
  (cite the stronger one where a `pol` step would otherwise reference the weaker).

No justification step is needed: the OPB is untouched (supplementals were never in it),
we simply decline to derive and to allocate inert ones.

## The abstraction: a unified `SolveStep`

One step type, no preprocess/search distinction. A step is **a partial solver: a
transformer of shared state that may also conclude the problem.** It may (a) tighten
domains, (b) produce data later steps consume (a supplemental graph, degrees, NDS,
clique sizes, nogoods), (c) conclude (SAT / UNSAT / enumeration-complete), or (d) make
progress and hand off.

```cpp
struct SolveStep {
    virtual auto name() const -> std::string;
    virtual auto cost() const -> Cost;                                  // Fast / Medium / Expensive
    virtual auto applicable(const HomomorphismParams &, const SolveState &) const -> bool;
    virtual auto run(SolveState &) -> StepOutcome;   // Continue | Concluded(result) | Aborted
};
```

- **`SolveState`** is the shared working state, spanning what is today split between
  `HomomorphismModel` and the searcher: domains; the model graph rows plus slot
  allocation (`max_graphs` *grows* as builders commit) plus degree / NDS / clique caches;
  accumulated nogoods / watches; the best mapping and solution count; the `Proof` object;
  and the timeout / budget.
- **`applicable()`** subsumes the `supports_*` traits *and* the mode gates — e.g. the
  loop shortcut is "not applicable when counting" (today's `! params.count_solutions`) —
  so a step's possible findings are always mode-appropriate by construction.
- **`StepOutcome`** is `Continue` (drove state forward, not done), `Concluded(result)`
  (the driver emits the matching proof conclusion and stops), or `Aborted` (budget /
  timeout → unknown). This outcome shape is provisional; revisit if it proves awkward.

Steps are then: setup (recode, emit OPB model), filters (degree / label / loop /
global-degree / NDS / clique / Hall), graph builders (exact-path / distance-2 / distance-3
/ k4 / extra-shape, each of which *commits or discards* a candidate slot), shortcuts
(pattern > target, loop, clique reduction), and **search** — a bounded or unbounded
search engine wrapped as a step. The registered step list *is* the documented,
reorderable solve strategy.

**Staging (S3) is then native, not bolted on.** "Run fast preprocessing, search a bit,
then full preprocessing" is just a step list: Fast filters → cheap builders → a *bounded*
search step (which returns `Continue` after recording nogoods if it does not conclude) →
Expensive builders (which re-filter against the new graphs and accumulated nogoods) → an
unbounded search step. Nogoods persist in `SolveState` across rounds and feed the later
search.

## Staged execution

Each phase is an independently mergeable PR.

- **Phase 0 — Instrument, no behaviour change.** Add PBP line-count, OPB size,
  `max_graphs`, and solution-count columns to `test-instances/random_proof_sweep.bash`
  plus a ctest, so every later win is measurable and Phases 1–3 can be proven
  byte-identical.
- **Phase 1 — Extract the solve pipeline (S4; S1 foundation).** Introduce `SolveStep` and
  a pipeline driver; move the existing phases into steps **in the same order, with no
  logic change** — the cheap checks and shortcuts become steps, and the existing
  sequential search becomes the **terminal step**. Keep the `HomomorphismModel` / searcher
  boundary and have the search step wrap the existing engine; do **not** unify `SolveState`
  yet. Guardrail: byte-identical proofs (sweep + cake).
- **Phase 1b — Unify `SolveState`.** *Done (commit 504c96c).* The constraint model, the
  (root) domains and the nogood store now live in one `SolveState` (`gss/innards/solve_state.hh`)
  carried through the pipeline, instead of as locals scattered across the model build and the two
  solvers. The model is built into the carried state — still lazily by the main solve step, after
  the cheap shortcut steps, so the `max_graphs` guard cannot fire before a shortcut concludes — and
  is "frozen" during a search step. The searcher no longer owns its nogood store: it takes a
  `Watches &` from the caller (the sequential path supplies `state.watches`; the threaded terminal
  search keeps a per-thread `vector<Watches>`, as threaded staging is deferred). This is the
  prerequisite for Phase 6's growable-model-across-rounds and persistent nogoods. Guardrail held:
  274/274 proof files byte-identical, 49/49 ctest, clean under ASan/UBSan, threaded runs match the
  sequential results.
- **Phase 2 — Self-allocating slots + commit/discard (S1).** Replace
  `calculate_n_shape_graphs`, the manual `next_*_supplemental` threading, and the
  slot-index sanity check with context-driven slot allocation; `max_graphs` becomes a
  result of registration. The discard path exists but always commits (no behaviour change
  yet). Proofs unchanged.
- **Phase 3 — Co-locate proof emission (S1, S4).** Move each inline proof loop out of
  `prepare()` into its step, behind a narrow proof interface. Byte-identical. After this,
  a new graph or filter is genuinely one new file.
- **Phase 4 — Supplemental constraint-level subsumption (S2a).** *Done (commit cc72e5f).* The
  originally-planned graph-level inert/duplicate elimination was **shelved**: measurement showed the
  supplemental adjacency derivations are ~92–100% never *explicitly* cited but remain RUP-load-bearing in
  search, so dropping whole graphs (or eliding by apparent deadness) is unsound — only a post-hoc proof
  trimmer can safely act on dynamic deadness. What *is* sound and instance-independent is **constraint-level
  subsumption**: the adjacency constraints nest (`distance3 ⊇ exact-path-1 ⊇ exact-path-2 …`, same head), so
  emit only the set-minimal (strongest) constraint per head; a degree/NDS pigeonhole that needs an elided
  per-graph constraint re-derives it on demand as a one-step weakening of the kept one (`ia <wider> : @kept`)
  and deletes it again. Hidden `--no-proof-supplemental-subsumption` (default on) restores the full emission
  for studying the effect (e.g. on a proof trimmer). Guardrail: identical solution counts (correctness sweep);
  smaller proofs that still verify (proof sweep, 128/128); ~25–39% PBP shrink on preprocessing-heavy configs.
  (Graph rows / `max_graphs` / search-side filtering are deliberately untouched.)
- **Phase 5 — Cost ordering: defer the adjacency derivations (S2b).** *Done.* The cheap
  concluding steps (pattern-too-big, target-loop, clique) already run before the supplemental
  *builders* (those live in `MainSolveStep` via `prepare()`), so an early conclusion never paid
  for a supplemental derivation. What it *did* pay for was the **loop-cancelled adjacency
  derivations** (`loop_fix_adjacencies`): these are PBP derivations the degree / supplemental /
  distance-3 pols later cite, but they were emitted at the tail of `emit_model` — the very first
  step — so even a trivial refutation carried them. Under `--induced` *every* non-edge constraint
  has the target in its permitted set (a vertex is its own non-neighbour), so all of them get
  rewritten: the induced `c3 -> trident` refutation (pattern bigger than target, a one-line
  injectivity pigeonhole that cites no adjacency constraint) emitted **2430** dead `@adj` rup
  lines before concluding. Phase 5 moves the derivation out of `emit_model` into a new
  `HomomorphismProofs::derive_loop_fixed_adjacencies()`, called at the start of `MainSolveStep`
  (after the cheap steps, before `prepare()`'s first `@adj`-citing pol). Nothing is emitted to the
  proof between the two points for an instance that reaches search, so its proof is **byte-identical**;
  an instance that concludes early simply omits the derivation. The OPB is untouched (it stays
  complete and up front). Guardrail held: 273/274 proof-metrics + sweep files byte-identical, the
  one change being `induced_unsat` PBP **2439 -> 9** (still verifies); all 49 ctests pass (release +
  ASan/UBSan), including the random proof sweep and the cake pipeline. (A *further* S2b lever — when
  search does start but the cheap original-graph degree/Hall filter already refutes at the root,
  before the supplementals it built were needed — is left as a follow-up; it is the non-staged
  cousin of staged solving and would change the filtering order, so it overlaps with "make `--staged`
  the default" rather than this byte-identical reordering.)
- **Phase 6 — Staged solving (S3, first-class).** *v1 done (commit 52185f3; flag-gated, sequential, decision).*
  `--staged` runs a cheap first round (original graph only: degree + Hall, no NDS, no
  supplemental graphs), and only if a small backtrack budget passes without concluding does it
  build the supplemental graphs and NDS, re-filter the domains, and continue unbounded. The
  budget is the transition trigger: `prepare()` is split into the Stage-1 work and a new
  `build_supplemental_graphs()`, and the driver in `SequentialSolver` runs Stage 1 under its own
  fixed-budget restart schedule, then at the first restart calls `build_supplemental_graphs()` +
  `tighten_domains_with_supplementals()` and switches to the user's schedule. Nogoods persist in
  the carried `state.watches` (Phase 1b) across the transition. **Full proof support:** the
  supplemental graphs are derived mid-proof at the level-0 restart boundary (`back_up_to_top`),
  and the re-filter emits only the new degree/NDS prunings; staged proofs verify under VeriPB,
  and easy instances that conclude in Stage 1 emit *no* supplemental derivations (e.g. trident
  PBP 2168 → 136). Guardrails: unstaged proofs byte-identical; staged proofs verify (sweep
  160/160); staging is satisfiability-neutral vs unstaged (random correctness sweep). **Deferred
  to follow-ups:** staged *counting* (a mid-enumeration restart would recount explored subtrees
  — needs solution-blocking nogoods, so `--staged --count-solutions` is rejected for now);
  threaded staging (the threaded engine's barriers / nogood sharing make a *bounded threaded*
  round harder — the unbounded threaded search stays terminal); finer multi-tier schedules;
  expressing the bounded round as a first-class composable `SolveStep`; making `--staged` the
  default. (Aside, found while testing: the decision-mode clique and target-loop shortcut steps
  ignore `--induced` on loopy targets — a *pre-existing* bug, orthogonal to staging.)

### Strand → phase

| Strand | Delivered by |
| --- | --- |
| S1 (extensibility) | Phases 1–3 |
| S2a (don't emit unused) | Phase 4 (constraint subsumption; graph-level elision shelved as unsound) |
| S2b (shortcut obvious unsat) | Phase 5 |
| S3 (reorder + stage) | Phase 5 (reorder) + Phase 6 (stage) |
| S4 (sequencing clarity) | Phases 1, 3 |

## Risks and invariants

- The `<directed, has_edge_labels, induced, verbose_proofs>` templated propagation hot
  loop stays out of scope — only the engine *wrapping* it becomes a step.
- Proof byte-stability is the contract for Phases 1–3; Phases 4 onward change proofs
  deliberately, verified by re-running the sweeps, not by diff.
- The OPB stays complete. `finalise_model` is unavoidably up front — even a cheap
  refutation's proof cites OPB constraints — so the S2b / S3 deferral is of *derivations*,
  not of the model.
- Model-growth vs threading: under staging the model can only grow at a barrier between
  threaded-search rounds (it is shared `const` across threads), hence sequential staging
  first.
- Proof under staging: deriving supplementals mid-proof (after a search round) is legal
  but interacts with the existing "no counting with restarts under proof" restriction;
  staged + counting + proof will likely need the same guard initially.
- Preserve the locally-injective / loops gating (issues #56, #58) exactly when the
  `supports_*` traits become `applicable()`. `extra_shapes` re-enters
  `solve_homomorphism_problem` and is simply another step.
