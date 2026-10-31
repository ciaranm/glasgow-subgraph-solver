# Splitting a decision search into components

`HomomorphismParams::split_decision` (the hidden `--split-decision`) is the decision-search
counterpart of [split-counting.md](split-counting.md). Wherever the pattern vertices still
unassigned at a search node fall into two or more components, it can search each component
on its own. If one component has no solutions, neither has the node. It is off by default.

## What it does

The components are the same as for counting: connected components of the pattern among
the unassigned vertices, with edges in either direction joining a component
(`HomomorphismSearcher::find_components`, shared by both). Every constraint between an
assigned and an unassigned vertex has already been propagated into the domains. Induced
mappings, local injectivity and less-constraints are refused. So the only constraint
between two components is injectivity, and a component's solutions depend only on its own
domains.

Searching a component on its own drops the injectivity between it and the others, so it
is a relaxation, and the split works like this:

- **Refutation.** If one component has no solutions on its own, the node has none. This
  is what the split is for. A plain search that branches on the vertices of a satisfiable
  component first refutes the unsatisfiable one again under every solution of the first,
  so each refutation is repeated a multiplicative number of times. That is where the
  big wins below come from.
- **Witnesses.** A component's search that finds a solution keeps it as a *witness*. For
  each vertex, `_witness_image` records its image and `_witness_id` records which
  solution the image came from. At a later node, a component still has its witness if:
  - every one of its vertices has a recorded image still in its domain;
  - all those images came from the same solution.

  After a backtrack, the vertices of one component can come from two earlier components
  whose solutions were found separately. Those solutions need not agree on the edges
  between them, so such a mixed component must be searched again (`decision_split_mixed`
  counts these).
- **Combining.** If every component has a solution and no two components' images
  collide (non-injective mappings cannot collide), the witnesses are a solution together.
  Otherwise the split gives no answer, and the node is searched as usual. Components are
  searched smallest first, by the sum of the logs of their domain sizes, so an empty one
  is found as cheaply as possible.
- **Nesting.** A component's search is the same searcher, so a node inside it can split
  again. While `_in_component` is set, a full assignment of the component is not a
  solution of the instance, so the leaf returns it without logging it.

## When it looks

The split is not tried at every node. A node that splits only after its first child has
failed, by a propagation failure or by an unsatisfiable subtree, and at most once:

- **The cost of looking everywhere.** Looking at every node where the pattern falls apart
  costs a component search per node. On a satisfiable instance that the plain search
  would have walked straight into, this is mostly wasted. Over the 14,621 standard pairs,
  looking at every split node gained 28 finishes and lost 70. Almost all of the losses
  were satisfiable instances that default settings solve quickly.
- **What waiting costs.** A search heading for a solution never pays. A node that keeps
  failing finds an empty component after one failure, rather than proving it again
  under every value of the others.

## Restarts and staged search

A component's search uses the node's restart schedule. A restart inside it is propagated
out as a restart of the whole search. The node posts restart nogoods for the values it has
tried, which includes the value whose failure set off the split, just as when a child
returns a restart. Staged search works unchanged, since its stage transition is a restart.

## Proofs

Refuting a node by searching a component first is a choice of which vertex to branch on.
The proof is the ordinary search proof of that branching order, so no new kind of
reasoning is needed. Only one piece of level bookkeeping is new:

- **Where it goes wrong.** A component search that finds a solution returns without
  backtracking out of its levels. It leaves core nogoods registered for deletion at levels
  deeper than the node's. A later, unrelated `forget_level` on those levels would then
  delete them before the node's own backtrack nogood subsumes them, and the checked
  deletion fails.
- **The fix.** `Proof::keep_deeper_lines_at_active_level` copies each such line to the
  node's level (`pol id ;`, `core id new ;`, `del id old ;`) and registers the copy there.
- **How often it matters.** Without it, 8 of the 298 instances that split in a
  2,000-instance random proof sweep fail to verify. The three fixed instances in
  `test-instances/split-decision` are from those 8.

A solution made by combining witnesses is posted as an ordinary solution, and VeriPB
checks it.

## Restrictions

| Combination | Why |
|---|---|
| counting or enumerating | `--split-components` is the counting version; this one stops at the first solution |
| `--locally-injective` | combining witnesses would need local injectivity's constraints between components checked as well as collisions, and it is untested |
| `--induced` | non-edges relate the components |
| threads | the witnesses are per searcher, and it is untested |
| less-constraints | they relate the components |
| `--decomposition` | a different top level |

## Testing

- **`random_homomorphism_test.cc`.**
  - The oracle runs `split_decision` on every instance, with and without restarts
    (restarting after every backtrack), plain and staged. It checks that induced
    mappings and local injectivity are refused.
  - "Splitting a decision search into components matches plain search" compares against
    the plain solver on sparse patterns and on blocks joined through hubs, in every
    family, injective and non-injective, without restarts, restarting, and staged. It
    checks satisfiability and that the mapping verifies, and has floors on the numbers of
    split nodes, refutations, combined solutions, reused witnesses and mixed components.
- **`test-instances/decision_split_proof_sweep.bash`** (ctest `decision_split_proof_sweep`).
  - It checks the answers and the proofs on 360 random hub and sparse instances, under six
    option sets, and on the three fixed instances.
  - It has floors on the numbers of instances that split, that a component refuted, and
    that were solved from components.
  - It takes about 25 seconds.

Each of these mutations is caught by one or both:

| Mutation | Random test | Proof sweep |
|---|---|---|
| no collision check when combining | caught | caught |
| reuse a witness whose image has left the domain | caught | caught |
| reuse a mixed witness | caught | |
| no `keep_deeper_lines_at_active_level` | | caught (the fixed instances) |
| treat a component's full assignment as a solution | | caught |
| restart nogoods on every value, not just those tried | caught | |

## Results

These come from the standard 14,621 pairs, run in decision mode at 100 s. Each split run
is compared with default settings, restarts on, from the same commit, with 80 runs in
parallel. Comparing two sweeps of the family that never splits shows run-to-run noise of
about 15%. No run disagreed with another on satisfiability.

**Overall.** The split gains 18 finishes and loses 22. 1,398 pairs split at least once.
Among the pairs where either run takes at least 100 ms, 45 are more than twice as fast and
90 more than twice as slow.

**Gains.** The gains are the thrashing that the split was meant for:

| Pairs | Default | With the split |
|---|---|---|
| g9-g44, g9-g45, g9-g29, g24-g29/30/31 (LV) | over 100 s | 4–6 ms |
| g38-g60, g38-g73 | over 100 s | 41 and 97 ms |
| g28 LV pairs, several `meshes` pattern3 and two `phase` instances | over 100 s | under 100 s |

**Losses.** Seven of the 22 lost pairs never searched a component, so they are noise at
the timeout. In the rest, components that do have solutions are expensive to solve on
their own, so the work is duplicated. This happens on unsatisfiable instances as well as
satisfiable ones: g35-g93 and g43-g93 search hundreds of components without a single
refutation. g35-g68 goes from 6.7 s to 91 s with 225 refutations in 8,485 searches.

Searching only the components other than the branching vertex's was tried and was no
better. It loses g9-g44, which goes back to 2.8 s.

## Not done

- **Satisfiable components that are expensive to solve.** This is where the remaining
  losses come from. A budget on component searches is the obvious fix and an unappealing
  one; there is no better idea yet.
- **Components recomputed at every failing node.** They are found from scratch each time,
  which costs time linear in the pattern's edges. Counting needs this much more often, and
  an incremental structure shared by both would help counting more.
- **Default settings.** It is off by default.
