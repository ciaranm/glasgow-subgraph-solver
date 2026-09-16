# Option compatibility

Which combinations of solver options are legal, which are quietly turned off, and which are
simply untested. For what the options *do* see the [README](../README.md); for proof logging's
own restrictions see [proof-logging.md](proof-logging.md), which this file does not duplicate.

This exists because the answer used to be spread across four places — the proof guard block in
`homomorphism.cc`, the throws in `solve_homomorphism_problem`, the predicates in
`innards/homomorphism_traits.cc`, and the conditions on the individual pipeline steps — and
none of them said what the others did. Six wrong answers were found in one sitting by writing
the combinations down and trying them (#91 to #97).

## Three kinds of "no"

A combination that does not work can fail in three quite different ways, and telling them
apart matters:

1. **Refused.** `UnsupportedConfiguration` is thrown, the caller sees a message. Used where
   carrying on would give a wrong answer and there is no sensible fallback.
2. **Silently disabled.** The option is accepted and then ignored, because the filtering it
   asks for is not sound here. This is what the traits layer does, deliberately: `--cliques`
   is a hint, and honouring it where it is valid while dropping it where it is not is more
   useful than refusing the whole run. The cost is that nothing tells you it happened —
   `--record-filter-activations` is the way to see it, and the sweep's golden table records
   it as a `vacuous` cell.
3. **Untested.** Nothing stops you and nothing checks you. The list is at the end.

## What is refused

| Combination | Message |
|---|---|
| `--staged` with threads | Staged solving requires sequential search, use `--threads 1` |
| threads without restarts | Threaded search requires restarts |
| proof logging with threads | Proof logging cannot yet be used with threads |
| proof logging with `--staged --count-solutions` | …with staged counting |
| proof logging with clique detection | …with clique detection, use `--no-clique-detection` |
| proof logging with less-constraints | …with less-constraints |
| proof logging on a labelled pattern | …on labelled graphs |
| proof logging, counting, with restarts | …when counting with restarts, use `--restarts none` |
| more than 8 graph pairs (`--n-exact-path-graphs` too large, with `--distance3` and `--k4`) | Supplemental graphs won't fit in the chosen bitset size |
| a cycle in the `--pattern-less-than` constraints | Pattern less than constraints form a loop |

## What is silently disabled

Every row here is a predicate in `innards/homomorphism_traits.{hh,cc}`, which is the single
place to look. `has_loops` means *either* graph has a self-loop.

| Filter | Enabled only when | Because |
|---|---|---|
| exact-path graphs | supplementals on, not non-injective, not (locally injective and loops) | #58 |
| distance-2 graph | supplementals on, no exact-path graphs, not (locally injective and loops) | as above |
| distance-3 graph (`--distance3`) | supplementals on, fully injective | |
| k4 graph (`--k4`) | supplementals on, not non-injective, not (locally injective and loops), **neither graph directed** | #97 |
| degree and NDS | not non-injective, not (locally injective and loops) | #58 |
| degree and NDS *exactly* | induced, equal sizes, fully injective | a bijection is forced only then |
| whole-instance degree sequence | fully injective | |
| clique-size constraints (`--cliques`) | fully injective, **or** no loops | #91 |
| …on supplementals | as above, and not non-injective | #91 |
| clique detection (`--clique-detection`) | not counting, no proof, and fully injective **or** a loopless target | #94 |
| nogood recording | `--staged`, or the restart schedule might restart | nothing consults a nogood without a restart |

Two conditions recur, and it is worth seeing why they are the same argument twice. Both the
clique-size filter and the clique reduction need k pattern vertices to reach k *distinct*
target vertices. Injectivity gives that outright. So does a loopless target, because a
homomorphism cannot collapse two adjacent vertices onto one image without that image carrying
a self-loop — which is why these survived so long without the premise being stated, and why
adding a loop to a test instance is worth doing (`test-instances/small` and `large` have one).

The supplemental graphs need a *different* premise, and `--cliques-on-supplementals` needs
both. Adjacency in the distance-2 graph means "within distance two", and two such vertices may
share an image in a perfectly loopless target — `build_exact_path_graphs` sets that graph's
diagonal precisely so that propagation allows it.

## The pipeline steps' own guards

The concluding steps in `homomorphism.cc` each have conditions that are not in the traits
layer, because they are about the shape of the instance rather than about an option:

- **`PatternBiggerThanTargetStep`** needs a non-shrinking (fully injective) mapping.
- **`TargetLoopShortcutStep`** needs non-injective, non-induced, an unlabelled pattern, a
  loopy target, and no counting or enumeration callback.
- **`CliqueShortcutStep`** needs `can_use_clique()`, a pattern that `is_simple_clique()`
  accepts — no labels, no loops, not directed — and a target that is neither directed (#93)
  nor loopy-with-induced.

`is_simple_clique()`'s list is exactly what the clique solver's view of a graph leaves out: it
has no notion of a label, a loop, or an edge direction. Anything else handed to that solver
needs the same check.

One ordering used to be load-bearing and is no longer: `TargetLoopShortcutStep` runs before
`CliqueShortcutStep` and used to be what kept the clique reduction away from the non-injective
loopy case. An API caller that set an enumerate callback without `count_solutions` skipped the
first step and reached the wrong answer, so `can_use_clique()` now states the premise itself
(#94). Do not reintroduce reasoning of the form "the earlier step catches it".

## What an instance means when the two graphs disagree

The model reads its directedness and its label-ness off the **pattern**:
`HomomorphismModel`'s constructor branches on `pattern.directed()` and
`pattern.has_edge_labels()`, and builds the target's label tables only when the pattern has
labels at all. The resulting semantics are consistent, but they are nowhere else written
down:

- An **undirected graph is a digraph with both arcs**. A directed pattern arc may therefore be
  mapped onto an undirected target edge, and an undirected pattern edge needs *both* arcs
  present in a directed target.
- **Labels are pattern-driven.** An unlabelled pattern ignores the target's labels entirely,
  which is the sensible reading — an unlabelled pattern vertex matches anything — but it does
  mean a labelled target is silently unlabelled when the pattern says nothing.

A supplemental-graph builder cannot rely on the first of these, because it reads rows rather
than going through the searcher: a builder written for undirected rows is wrong on an
asymmetric row whichever graph it came from, which is why the k4 guard asks whether *either*
graph is directed (#97). `ProcessedGraphsData::directed` follows the pattern, because that is
what selects the searcher's propagation path; `either_graph_directed` is the question a
builder should ask.

Mixed instances are not swept: every family in the sweep gives the pattern and the target the
same shape.

## What the sweep covers

`gss/option_sweep_test.cc` is a pairwise covering array over fifteen option columns, times
nineteen instance families (loops × directed × vertex labels × edge labels, plus a clique
pattern, a pattern bigger than its target, and an equal-sized pair). Its reference is the same
instance solved with the filtering off, and it checks that satisfiability and — where
observable — the solution count are unchanged, and that any mapping it gets back verifies.

Structure is a multiplier rather than another column on purpose. Pairwise over options alone
would not catch #58, which needs local injectivity *and* loops *and* supplemental graphs
together; multiplying gives strength 3 for every (option, option, graph property) triple,
which is the shape of every condition in the table above. The known-dangerous triples
— each filter alone, against each injectivity mode — are added explicitly on top, because
whether the array happens to try a filter in isolation is otherwise an accident of how the
greedy construction came out.

`gss/random_homomorphism_test.cc` is the other half: a brute-force oracle over the problem
axes, exhaustively, on instances small enough to enumerate every function between them.
Metamorphic testing cannot catch a bug the baseline shares, which is exactly why that test
stays separate rather than being folded in.

Not covered, and worth deciding on separately:

- **Threads.** Nondeterministic, so a cell's outcome would not be reproducible.
- **Proof logging.** Much slower, needs VeriPB, and covered separately by
  `test-instances/random_proof_sweep.bash`.
- **`--shape`.** `build_extra_shape()` builds its master graph with `add_edge` and sets vertex
  labels of its own, so it ignores edge direction and labels the same way `build_k4_graphs()`
  did before #97. Sweeping the directed and labelled families would be testing a path already
  known to be incomplete; it probably wants the same guard k4 got.
- **Timed restarts.** Wall-clock dependent.
- **`--decomposition`** (`sip_decomposer`), the clique solver and the common-subgraph solver:
  different top levels, each with its own option space.
- **`--solution-limit`**, which deliberately returns an incomplete search.
- **Mixed instances**, as above.

## The golden table

`gss/option_sweep_golden.tsv` records, for every (configuration, family) cell, whether it came
out `ok`, `unsupported` or `diverges`, and whether any optional filter actually removed
anything (`active`) or the cell was vacuous. A divergence fails the build whatever the table
says; what the table adds is that a *newly* unsupported or *newly* vacuous cell shows up as a
reviewable diff rather than as nothing at all. About two cells in five are vacuous, which
is not a defect: a filter with nothing to do on a family, or one the traits layer switched off,
is exactly what the activation column is there to make visible.

After reviewing a diff:

```shell
GSS_SWEEP_REGENERATE=1 ./build/option_sweep_test
```

`GSS_SWEEP_INSTANCES=n` raises the instances per cell from the default four, for hunting; in
that mode the activation column is read as a floor, since extra instances can wake a filter up
but never silence one. A thousand per cell takes about half a minute, and was clean once #97 was
fixed.

A table that only reproduces on the machine that made it would be worse than no table, so two
things are deliberate and worth keeping if you add a family or a column. The instance generator
uses raw `rng()` arithmetic and integer percentages: `std::mt19937` is specified exactly, but
`std::uniform_real_distribution` is not, and libstdc++ and libc++ consume the generator
differently. And the activation reading comes from a probe solve with the value ordering forced
to `None`, because `Biased` — the default — and `Random` drive the search from
`std::uniform_int_distribution` and `std::shuffle`. The *answers* do not depend on search
order, which is why no outcome ever differed between platforms; how much each filter removes
along the way does.
