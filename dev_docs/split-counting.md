# Counting by splitting into components

`HomomorphismParams::split_components` (the hidden `--split-components S`) counts solutions
without multiplying out independent parts of the pattern, inside the ordinary search. With
`--split-statement FILE`, it writes a counting statement as it goes, in the format of
[counting-statements.md](counting-statements.md). With `--prove` too, it proves that the
instance has exactly as many solutions as the statement's tagged formula, so a count is
certified by one VeriPB check (see [Proofs](#proofs)). Either way, the statement is also a
compressed form of the solutions themselves (see
[The statement as the solutions](#the-statement-as-the-solutions)).

## A worked example

The pattern is a triangle A B C and a pentagon V W X Y Z, joined by the edge C–D, the square
D G U H, and the edge U–V. The target has 17 vertices:

```
pattern: A-B A-C B-C C-D D-G D-H G-U H-U U-V V-W W-X X-Y Y-Z Z-V
target:  5-6 5-7 6-7 2-6 2-7 2-17 1-2 1-3 1-4 2-4 3-4 3-8 3-12 3-14 8-9 9-10 8-10
         10-11 11-12 8-12 4-13 13-17 13-14 14-15 15-16 16-17
```

Counting non-induced injective maps, the whole instance has 56 solutions. Take the part of
the search where G is at 1 and H is at 4, which has 28 of them.

**The node.** D is adjacent to both G and H, so it must go to a common neighbour of 1 and 4,
which is 2 or 3, and so must U. Every path from the triangle's side of the pattern to the
pentagon's went through G or H, so the unassigned vertices now fall into two components,
A B C D and U V W X Y Z. Only injectivity ties them together.

**Each component on its own.** Each is searched from the node's domains:

| Component | Choice | Solutions |
|---|---|---:|
| A B C D | D at 2: C at 6, with A and B at 5 and 7, or C at 7, with A and B at 5 and 6 | 4 |
| A B C D | D at 3: C at 8, with A and B at 9 and 10 | 2 |
| U V W X Y Z | U at 3: V at 8 or 12, round the cycle 8–9–10–11–12 either way, or V at 14, round 14–13–17–16–15 either way | 6 |
| U V W X Y Z | U at 2: V at 17, round 17–13–14–15–16 either way | 2 |

Multiplying 6 by 8 would give 48, but some of those combinations use a target vertex twice.

**Signatures.** The triangle's component can only reach 2, 3, 5–10, and the pentagon's only
2, 3 and 8–17, so only 2, 3, 8, 9 and 10 can ever clash. Each solution is tallied by the ones
of those it uses. The four triangle solutions with D at 2 all have signature {2}, because 5,
6 and 7 are out of the pentagon's reach, so they are one entry with a count of 4. Combining
is a sum of products over the pairs that do not clash:

| Triangle | Count | Pentagon | Count | Clash | Product |
|---|---:|---|---:|---|---:|
| {2} | 4 | {2} | 2 | 2 | 0 |
| {2} | 4 | {3, 8, 9, 10} | 4 | | 16 |
| {2} | 4 | {3} | 2 | | 8 |
| {3, 8, 9, 10} | 2 | {2} | 2 | | 4 |
| {3, 8, 9, 10} | 2 | {3, 8, 9, 10} | 4 | 3, 8, 9, 10 | 0 |
| {3, 8, 9, 10} | 2 | {3} | 2 | 3 | 0 |
| | | | | | 28 |

Nothing is subtracted, and every number is a sum of products of counts of solutions.

**The statement.** The solver splits a little later than this: at this node it branches on
C first, which fixes D and U by propagation, and then splits. Its counting statement, as
`count_statement --render` writes it, is:

```
root [28] = (G->1 H->4) * one of
  | (C->6 D->2 U->3) * n2 * n18  [12]
  | (C->7 D->2 U->3) * n18 * n24  [12]
  | (C->8 D->3 U->2) * {A B} = {9 10} * n37  [4]
n2 [2] = {A B} = {5 7}  [2]
n18 [6] = (V W X Y Z) in {(8 9 10 11 12) (8 12 11 10 9) (12 8 9 10 11) (12 11 10 9 8)
    (14 13 17 16 15) (14 15 16 17 13)}  [6]
n24 [2] = {A B} = {5 6}  [2]
n37 [2] = (V->17) * (W X Y Z) in {(13 14 15 16) (16 15 14 13)}  [2]
```

`C->8` is one assignment, and `{A B} = {9 10}` is A and B at 9 and 10 in either order, so two
ways. The last alternative is the table's row worth 4. The first two, 12 each, are its rows
worth 16 and 8, regrouped by where C goes. `n18` is written once and used twice, since with
C at 6 or at 7 the pentagon's component faces the same subproblem, so the statement is a DAG
rather than a tree.

The files are `test-instances/split/worked_example_*.csv`, the `pinned` pair fixing G and H
by vertex labels. `test-instances/split_worked_example.bash` (ctest `split_worked_example`)
checks the count and the rendering above, and `split_proof_sweep` checks the whole
instance's proof. To run it by hand, with the value ordering fixed, since the default one is
random and differs between standard libraries:

```shell session
$ ./build/glasgow_subgraph_solver --count-solutions --value-ordering none --split-components 1000 \
    --split-statement pinned.stmt --format csv test-instances/split/worked_example_pinned_{pattern,target}.csv
$ ./build/count_statement --render --format csv --pattern test-instances/split/worked_example_pinned_pattern.csv \
    --target test-instances/split/worked_example_pinned_target.csv pinned.stmt
```

Labelled instances cannot be proof-logged, so for a proof run the whole instance, as in
[Proofs](#proofs).

## What it does

At each search node, the unassigned pattern vertices fall into connected components of the
pattern (edges in either direction join a component). Every constraint between an assigned
and an unassigned vertex has already been propagated into the domains. So the solutions
below the node are the combinations of one solution per component that satisfy the
constraints between components.

Without induced mappings or less-constraints, those constraints are only injectivity: a
pair of vertices in different components may not share an image. Under non-injectivity
there are none at all, and the components' counts just multiply. Call a pair that could
share an image, because its domains intersect, a *collision pair*.

A split counts each component on its own, and combines the counts:

- **Signatures.** Each component is searched on its own, by the same searcher, from the
  node's domains. Rather than counting its solutions, the search tallies them by their
  *signature*: the images of the component's vertices that are in collision pairs. An image
  is recorded only if a partner could also take it, and otherwise as −1, which collides
  with nothing.
- **Combining.** The tallies are combined one component at a time: a sum of products over
  pairs of signatures that do not collide. The key keeps only the images that a later
  component could still collide with.
- **Order.** Components are searched smallest first, by the product of their domain sizes,
  and a component with no solutions ends the split at once. So a split whose answer is zero
  costs at most a search of the components up to the empty one, which the plain search
  would also have had to refute.
- **Nesting.** A node inside a component's search can split again, and then its result goes
  into the enclosing component's tally rather than the count. So the tally also records the
  images the enclosing split needs, which are part of the signature at every level.

There is no inclusion–exclusion and no subtraction, so every number in the computation is a
sum of products of counts of solutions. Nor are there any new models or solver calls. A
component is searched with the model, the supplemental graphs and the propagation the plain
search already has. The supplemental graphs stay sound inside a component: they are implied
by the original constraints, and any combination of component solutions that is a solution
of the whole satisfies them.

## When it splits

A node with at least two components is split if the signatures have at most S possible
values: the product, over the recorded vertices, of one more than the number of images
recorded for each. Every table the split builds is keyed by some of those images, so S
bounds each tally, each intermediate table, and the pairs each combining step looks at.

The number of collision pairs is the wrong thing to bound. A bound on it allows four nearly
isolated vertices sharing about 4,400 target vertices, whose second combining step looks at
10¹¹ pairs, and on the benchmarks below it ran out of memory on 179 instances and overran
the timeout on 43 more.

Most nodes are not split, but every node has to ask, and on an instance whose unassigned
vertices stay in two components throughout, every node is refused for the same reason. So
the question is answered apart from the split itself, by `worth_splitting`, which finds the
components and counts the colliding vertices and the signatures in one pass over the
unassigned vertices. It works on bare words at the graphs' real widths, since a `SVOBitset`
works on all of its inline words whatever its size, in a workspace kept from one node to
the next, so that a node it refuses allocates nothing. It finds the images a partner of v
could take without looking at pairs: they are v's domain intersected with the images some
other component could take, which are those outside v's own component's union that some
component could take, together with those inside it that two components could take. Only a
split it lets through builds its components again and lists its collision pairs.

On the 311 instances of the sweep in [Results](#results) that never split, search the same
number of nodes either way, and take at least 1 s plainly, splitting costs 1.03× (p90 1.07×,
at worst 1.18×).

## Caching

The result of a split depends only on the domains of its unassigned vertices and on which of
their images the enclosing tally records. So it is cached on exactly that, and so is each
component's tally. This turned out to matter more than anything else. On si4_r001_m200.09,
10,373 splits had only 28 distinct subproblems, and 53,002 component searches had only 51.
In a sparse pattern, what is left far from the assigned vertices sees the same domains
whichever way the search got there. This is the component caching of #SAT solvers.

The caches are emptied when they reach about `HomomorphismParams::split_cache_limit` bytes
(1 GiB by default, the hidden `--split-cache-limit`), which can only cost time. Without
this, some instances filled 12 GB in well under the timeout. The searcher reports these
counts:

- `splits` and `nested_splits`;
- `split_cache_hits` and `split_cache_entries`;
- `component_cache_hits` and `component_cache_entries`;
- `split_cache_clears`.

## Interchangeable vertices

Under injectivity, a split's components often include single vertices whose neighbours are
all assigned, several with the same domain: the leaves of a hub, say, once the hub is placed.
Such vertices are interchangeable. Searched and combined one at a time, five of them sharing
nine values make 15,120 entries, one for each way of placing them, and on si2_b06m_m400.07
that was every pair the combining looked at.

So a split first collects its **free sets**: the vertices that are components on their own,
whose images the enclosing tally does not need, with the same domain. It keeps some free
sets, largest first, whose domains are pairwise disjoint, and does not search those
vertices. Free sets whose domains overlap a kept one's are searched as ordinary components.
It searches and combines the rest as before, keeping the images that could collide with a
free set's vertices, until the last combining step. There, for each compatible pair, each
kept free set's k vertices take distinct values from its domain less the pair's images. That
is d (d − 1) ⋯ (d − k + 1) ways, for d values left, and none if d < k. The new table is then
keyed only by what the enclosing tally needs. Not under non-injective mappings, where nothing
collides anyway.

A vertex alone in its component has exactly its domain as solutions. Every constraint
between it and the assigned vertices has already been propagated into its domain, and there
is nothing else in the component.

The free sets are placed pair by pair, inside the last combining step, rather than in a step
of their own after it. A separate step would need the last combining step to keep every pair
as its own entry, keyed by the images the free sets need: on gen-s17-LV-g45 that makes
counting 2.6× as long as without free sets, against 1.4× pair by pair. A free set of one vertex is kept too: it gains little for counting, but a
statement gets one node for it rather than one per value.

On the 269 instances of the benchmark families in [Results](#results) that split, against
no free sets:

- si2_b06m_m400.07: counting 0.42 s → 0.06 s, writing the statement 8.5 s → 0.13 s and
  175 MB → 1.4 MB;
- statements, in total: 5.7× smaller in si2_b06, 25× in si2_r01, 14× in si4_r00, 7.5× in
  si6_r00, 1.8× in si2_m4D and 1.2× in LV, where most combining is of components of two or
  three vertices; median per instance 1.06×;
- time: unchanged at the median with and without a statement. Counting is at worst 1.41×
  slower (gen-s17-LV-g45; three instances above 1.1×), and writing the statement is 10%
  faster in total.

## The statement as the solutions

A statement is a DAG of cubes, `and`s whose children have disjoint scopes, and `or`s whose
children are told apart by one pattern vertex: a decision-DNNF. The search writes it so that
its root's objects are the instance's solutions, one each, so it serves as a compressed listing
of them, and since every node knows its count, it answers a lot without being multiplied out.
The proof certifies the count, not which assignments the objects are. That is as for any
VeriPB or CakePB proof: even for a satisfiable decision problem with a witness, the checker
certifies that there is a solution, not what it is, and anything more would have to come
from inside the proof.

`count_statement --render` (`gss/counting_statement_render.hh`) writes a statement for people
to read, as definitions. A node gets a name if it has more than one object and is either
shared or multiplied by another such node. Everything else is multiplied out in place. What all of a definition's alternatives assign alike is written once
in front of them. Alternatives that differ only in their values become a table, and a table
is factored into the product of its projections where it is one, with a block that is every
ordering of one set of values written as `{u v} = {3 4}`, and one that is every way of
giving its vertices distinct values from a larger set as `{u v} distinct in {3 4 5}`. A
`distinct` node is written the same way in place, or as `u in {3 4 5}` for one vertex.
`--pattern` and `--target` give vertex names.

None of this needs a proof, and none of it is part of what checks a count: it reads a
statement that has already been checked for well-formedness. The statement is written the
same, byte for byte, with or without `--prove`. Given the search's order, it is also the same
on every platform: the tallies are hash tables, so wherever nodes or proof lines are written
while going over one, it is gone over in order of key. The default value ordering is random,
though, and standard libraries turn the same random numbers into different choices, so for
the same statement everywhere, fix it, as with `--value-ordering none`.

For example, si2_b09m_m400.08 has 6,451,200 solutions. Its statement is 136 KB, and it renders
in 20 lines, starting:

```
root [6451200] = (0->0 1->4 2->15 3->18 ...) * one of
  | (7->329) * {12 16} = {137 182} * {23 25 28} = {86 109 155} * {32 34} = {8 78} *
        {45 49} = {134 186} * {67 69} = {135 188} * {68 71} = {144 279} * {70 73} = {211 307} *
        n2  [1935360]
  | ...
n2 [5040] = {47 51 52 54 55} = {156 200 203 235 306} * {78 79} distinct in {1 41 90 98 104 216
    284}  [5040]
```

Without free sets (see [Interchangeable vertices](#interchangeable-vertices)) it renders in 64
lines: in place of `n2`, seven alternatives, one per value of vertex 78, each with 79 at one of
the other six and a further node for the permutation.

## Restrictions and testing

It is refused, rather than silently disabled, in these cases:

- with induced mappings, less-constraints, threads, restarts, staging or an enumerate
  callback;
- with local injectivity, because there are no proofs for it. Splitting would work: two
  vertices in different components could collide only if they have a common in-neighbour.
  Everything splitting does, though, is meant to be provable;
- when not counting;
- with proof logging, unless there is a counting statement to prove.

See [option-compatibility.md](option-compatibility.md).

`random_homomorphism_test.cc` checks it against the brute-force oracle on the oracle's
instances. It also checks it against plain counting on 1,600 sparser instances:

- with up to 8 pattern vertices;
- in all sixteen families, both injective and non-injective;
- with S from 1 (split only where nothing is recorded) to unlimited;
- with caches that are emptied at every insertion.

Sixteen hand-made mutations of the splitting code are all caught by the second of these
checks. Each gives a wrong count, except one, which crashes. The mutations are:

- ignoring collisions in the combining, or the −1 exception to them;
- recording every image as −1, or leaving out the images a partner could share;
- recording only one end of each collision pair;
- dropping an image before the last component that could collide with it, or ignoring the
  enclosing tally's need for it;
- leaving the assigned vertices or the enclosing relevance out of a nested result;
- following only out-edges when finding components;
- ordering by component number rather than by search position;
- swapping where a combined key comes from;
- leaving the domains or the relevance out of either cache key.

The weakest of these is leaving the enclosing relevance out of the split cache's key, which
only two of the test's checks catch.

## Proofs

With a statement, the proof does not log solutions. Instead, the search builds the statement
as it goes (`gss/innards/split_proof.cc`), and derives that every solution lies in its root.
Then a tail turns the core into exactly the statement's tagged formula, which
`count_statement --opb` writes, and the proof ends `OUTPUT EQUIENUMERABLE`:

```shell session
$ ./build/glasgow_subgraph_solver --count-solutions --no-clique-detection --split-components 1000000 \
    --split-statement s.stmt --prove p --format lad pattern target
$ ./build/count_statement s.stmt --opb tagged.opb
$ veripb --force-checked-deletion p.opb p.pbp tagged.opb
```

The count is `count_statement`'s, with no further conditions.

### The statement

The statement mirrors the search.

- **Search nodes.** A search node's node is the `or`, over its children, of the `and` of a
  cube for what the child's decision assigned (the decision and its propagation) with the
  child's node. A solution is just its cube. With a tally, a node has one such `or` per
  signature.
- **Components under assumptions.** A split's component searches are proved under an
  assumption of their own: a fresh literal `sa<k>`, defined to mean that none of the values
  outside the component's domains are taken, in place of the decisions above the split.
  So a component's result is a level-0 fact, `~sa<k>` or one of its entries' nodes, and a
  cached result can be used again wherever the domains are the same. The split itself has
  such an assumption for all its vertices, which implies each component's by unit
  propagation.
- **Combining.** Combining two tables makes, for each new entry, the `or` of the `and`s of
  its pairs. The proof then derives, for each left entry, that under the split's assumption
  it goes into one of its new entries: each compatible right entry by the pair's definition,
  and each colliding one by injectivity.

  That last step needs to know where both entries put the vertices they collide on. When
  either entry is a cube or an `and` of cubes, unit propagation sees it. When both are
  `or`s, it cannot, so each entry's node gets a *signature lemma* `~g ∨ y_v_t` for each
  recorded image, derived when the entry is made. Only the two fixed instances in
  `test-instances/split` need these. They were found by searching 2,000 random instances for
  a collision between two `or`s.
- **Free sets.** In the last combining step, each pair's node is the `and` of the two
  entries' nodes and a `distinct` node per kept free set, for its vertices and the values
  the pair leaves them. Its gate is the
  `and` of sub-gates of the proof's own: each vertex is at one of the values, and at most
  one of the vertices is at each value. The second kind are facts, from injectivity by one
  cutting planes step each. The first follow, for each pair, by unit propagation from the
  split's assumption, the entries' signature lemmas and injectivity, and then that the pair
  leads to its new entry. A pair whose free set has fewer values than vertices is refuted by
  pigeonhole, as one cutting planes step adding up each vertex's "at one of the values" and
  each value's at-most-one. With these, each left entry goes into its new entries by unit
  propagation, as before.
- **Regrouping.** A split's result, regrouped by the enclosing tally's signature, goes back
  to the search node as that node's entries.

Every `or` the search writes has its children told apart by one pattern vertex, which the
proof relies on (`SplitProof::disjunction`):

- A search node's `or` uses the branching vertex, which each child puts at its own value.
- Combining and regrouping `or`s are built as tries over the signatures, splitting at the
  first position where they differ.

  A −1 there means "at none of the recorded values that another child has", and at most one
  child of a trie node has it.

### The proof's ending

Every node's gate `g<i>` is defined into the core, both ways, when the node is written.
VeriPB checks a definition's `red` against the derived constraints too, so the downward rows
cannot be added later, once the search has derived things about the gate. The gates belong
to the proof: the tagged formula does not mention them.

At the end, the proof does the following:

1. **Back to the input.** It derives the root, and deletes the adjacency and injectivity
   constraints as implied by it: injectivity a pair at a time, then summed into the
   at-most-ones. It then swaps the `x` variables for the statement's `y` in the preserved
   set.
2. **Activity and tags.** It defines activity and tags from the gates, top down, and adds
   them to the preserved set:
   - the root is active when its gate holds;
   - any other node is active when one of its selectors is;
   - a tag holds when its `or` is active and its child's gate holds.
3. **Lemmas.** It derives these:
   - an active node's gate holds;
   - an active cube's literals hold;
   - for each `or`, what each child does to the vertex that tells the children apart.

   The vertex lemmas are built per node, bottom up, one small RUP each: from a cube's own
   row, from one child of an `and`, or from every child of an `or`.
4. **The tagged formula's rows.** It derives the tagged formula's own rows, among them each
   `or`'s at-most-one tag. That is one cutting-planes step over the children's vertex
   lemmas and the at-most-one on the vertex.
5. **Deletions.** It deletes everything else, each deletion checked with a one-line subproof
   that is a hinted RUP:
   - Deleting a tag's "active and child holds implies tag" row is where VeriPB checks that
     the `or`'s children are disjoint.
   - The activity variables and the gates then leave the preserved set, and the gates'
     definitions go with witnesses.

Every RUP in the ending is hinted, so no step propagates over the whole database. With the
deletions left to unhinted RUP, this part of the check takes 9.3 hours on si4_r001_m200.09,
against 81 s with hints.

With no solutions there is no statement, and the proof ends as an ordinary enumeration of
zero solutions. The statement file is then removed, and so it is on a timeout.

### Testing

`test-instances/split_proof_sweep.bash` (ctest `split_proof_sweep`) checks the count, the
statement's count and the proof on 150 small random instances, half of which split, and on
the two fixed ones.

Mutations that leave out a child's cube or the signature lemmas are caught. So is an `or`
given an extra child that overlaps another, while the search stays valid, so that the count
goes up. For 36 such mutants, VeriPB rejects every proof:

- at the disjointness check, when the copy claims its original's value;
- at the false vertex lemma, when it claims none of the others'.

### Cost

These are six si instances, with targets of 200 vertices, solved with proofs and checked
with `veripb --force-checked-deletion`. Every one verifies `OUTPUT EQUIENUMERABLE` with the
plain count:

| Instance | Count | Solve | Check | Proof |
|---|---:|---:|---:|---:|
| si4_r001_m200.06 | 6,434,208 | 3.7 s | 27 s | 115 MB |
| si6_r001_m200.06 | 21,504,000 | 6.6 s | 27 s | 190 MB |
| si4_r001_m200.04 | 283,557,120 | 3.1 s | 36 s | 103 MB |
| si4_r001_m200.02 | 15,454,800 | 3.9 s | 44 s | 123 MB |
| si6_r001_m200.01 | 70,447,104 | 6.2 s | 57 s | 197 MB |
| si4_r001_m200.09 | 153,135,360 | 4.0 s | 301 s | 131 MB |

Checking time follows the number of search nodes rather than the proof's size. Plain proofs
of instances like these can be as slow, where the supplemental graphs make each RUP step
expensive. On top of that, here a search step's unit propagation runs through every gate
definition written so far, so checking is quadratic in the size of the search. Deleting
injectivity needs a step per pair of pattern vertices per target vertex. Free sets account
for much of what is left: a build without them writes a 1.0 GB proof for si4_r001_m200.09,
which takes 25 s to solve and two hours to check.

## Complexity

Each component is searched once per distinct subproblem, so the components' counts add
rather than multiply. The combining is bounded by S per split, so for a fixed S it is
polynomial. It is not FPT in the number of collision pairs, as Ryan's inclusion–exclusion
recursion is, but that recursion could be applied to the signatures instead.

## Results

These are all 15,356 counting instances of the standard benchmark families plus generated
sparse patterns into LV targets (gen-s, gen-m), without proofs, at S = 10⁶, with a 100 s
timeout, 80 at a time on one machine. Plain counting is the solver without
`--split-components`, from a sweep on the same machine with a 1,000 s timeout, and counted
here as finished only within 100 s. The geometric mean is of split time over plain time, over the
instances both finish where either takes at least 100 ms.

| Family | Instances | Plain finishes | Split finishes | Gained | Lost | Geomean |
|---|---:|---:|---:|---:|---:|---:|
| si | 1,170 | 950 | 1,018 | 69 | 1 | 0.44 |
| gen-s | 427 | 144 | 190 | 46 | 0 | 0.22 |
| gen-m | 308 | 30 | 39 | 9 | 0 | 0.28 |
| LV | 3,831 | 3,189 | 3,191 | 4 | 2 | 0.97 |
| images | 6,302 | 6,302 | 6,302 | 0 | 0 | 1.01 |
| meshes | 3,018 | 2,902 | 2,902 | 0 | 0 | 1.01 |
| scalefree | 100 | 100 | 100 | 0 | 0 | 1.01 |
| phase | 200 | 0 | 0 | 0 | 0 | |
| all | 15,356 | 13,617 | 13,742 | 128 | 3 | 0.91 |

Every count agrees with plain counting, on every instance both finish. The
three lost instances all take over 98 s plainly. Of the 128 gained, 69 are not finished
plainly within 1,000 s either. The largest gains are on the very sparse si patterns:

- si4_r001_m200.04 goes from 74 s to 3 ms;
- si4_r001_m200.03 goes from 81 s to 13 ms;
- si4_r001_m200.09, the caching example above, goes from 41 s to 17 ms;
- g6-g59-normed goes from 83 s to 75 ms.

On the 245 instances that split at least once and take at least 1 s either way, the
geometric mean is 0.13. The worst of them is gen-s3-LV-g47, at 1.7×. On instances that
never split, the cost is that of asking, as in [When it splits](#when-it-splits).
