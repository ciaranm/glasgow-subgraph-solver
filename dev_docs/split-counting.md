# Counting by splitting into components

`HomomorphismParams::split_components` (the hidden `--split-components S`) counts solutions
without multiplying out independent parts of the pattern, inside the ordinary search. With
`--prove` and `--split-statement FILE`, it writes a counting statement as it goes, in the
format of [counting-statements.md](counting-statements.md). It proves that the instance has
exactly as many solutions as the statement's tagged formula, so a count is certified by one
VeriPB check (see [Proofs](#proofs)).

## What it does

At each search node, the unassigned pattern vertices fall into connected components of the
pattern (edges in either direction join a component). Every constraint between an assigned
and an unassigned vertex has already been propagated into the domains. So the solutions
below the node are the combinations of one solution per component that satisfy the
constraints between components.

Without induced mappings or less-constraints, those constraints are only injectivity: a
pair of vertices in different components may not share an image. Under local injectivity
that holds only for a pair with a common in-neighbour, and under non-injectivity for no pair
at all. Call a pair that could share an image, because its domains intersect, a *collision
pair*.

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

A first version bounded the number of collision pairs instead, which is the wrong measure.
It allowed four nearly isolated vertices sharing about 4,400 target vertices, whose second
combining step would have looked at 10¹¹ pairs. It ran out of memory on 179 benchmark
instances and overran the timeout on 43 more.

Finding the components costs a pass over the unassigned vertices at every node, and the
collision pairs a pass over pairs of them wherever there are two or more components. On an
instance that never splits, that costs about 20%.

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

## Restrictions and testing

It is refused, rather than silently disabled, in these cases:

- with induced mappings, less-constraints, threads, restarts, staging or an enumerate
  callback;
- when not counting;
- with proof logging, unless there is a counting statement to prove;
- with proof logging under local injectivity.

See [option-compatibility.md](option-compatibility.md).

`random_homomorphism_test.cc` checks it against the brute-force oracle on the oracle's
instances. It also checks it against plain counting on 1,600 sparser instances:

- with up to 8 pattern vertices;
- in all sixteen families and all three injectivity modes;
- with S from 1 (split only where nothing is recorded) to unlimited;
- with caches that are emptied at every insertion.

Seventeen hand-made mutations of the splitting code are all caught by the second of these
checks. Each gives a wrong count, except one, which crashes. The mutations are:

- ignoring collisions in the combining, or the −1 exception to them;
- recording every image as −1, or leaving out the images a partner could share;
- recording only one end of each collision pair;
- dropping an image before the last component that could collide with it, or ignoring the
  enclosing tally's need for it;
- leaving the assigned vertices or the enclosing relevance out of a nested result;
- ignoring local injectivity's common in-neighbours;
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

Every RUP in the ending is hinted, so no step propagates over the whole database. An
earlier version, which left the deletions to unhinted RUP, took 9.3 hours for this part of
the check on si4_r001_m200.09, against 81 s with hints.

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

These are six si instances, with targets of 200 or 280 vertices, solved with proofs and
checked with `veripb --force-checked-deletion`. Every one verifies `OUTPUT EQUIENUMERABLE`
with the plain count:

| Instance | Count | Solve | Check | Proof |
|---|---:|---:|---:|---:|
| si4_r001_m200.06 | 6,434,208 | 3.5 s | 25 s | 116 MB |
| si6_r001_m200.06 | 21,504,000 | 6.5 s | 27 s | 192 MB |
| si4_r001_m200.02 | 15,454,800 | 4.3 s | 102 s | 164 MB |
| si6_r001_m200.01 | 70,447,104 | 6.7 s | 101 s | 224 MB |
| si4_r001_m200.04 | 283,557,120 | 3.9 s | 106 s | 142 MB |
| si4_r001_m200.09 | 153,135,360 | 25 s | 7,220 s | 1.0 GB |

Checking time follows the number of search nodes rather than the proof's size. Plain proofs
of instances like these can be as slow, where the supplemental graphs make each RUP step
expensive. On top of that, here a search step's unit propagation runs through every gate
definition written so far, so checking is quadratic in the size of the search. Deleting
injectivity needs a step per pair of pattern vertices per target vertex.

Not done yet: proofs under local injectivity.

## Complexity

Each component is searched once per distinct subproblem, so the components' counts add
rather than multiply. The combining is bounded by S per split, so for a fixed S it is
polynomial. It is not FPT in the number of collision pairs, as Ryan's inclusion–exclusion
recursion is, but that recursion could be applied to the signatures instead.

## Results

These are 2,153 counting instances (si, LV, and generated sparse patterns into LV targets),
without proofs, with a 100 s timeout, at S = 10⁶. The geometric mean is of split time over
plain time. It is over the instances both finish where either takes at least 100 ms.

| Family | Instances | Plain finishes | Split finishes | Geomean |
|---|---|---|---|---|
| si | 380 | 159 | 226 | 0.08 |
| gen-s | 427 | 142 | 188 | 0.24 |
| gen-m | 308 | 30 | 39 | 0.34 |
| LV | 722 | 79 | 83 | 0.80 |

Neither finishes any meshes or phase-transition instance. In all, 410 instances finish
plainly and 536 with splitting; at S = 10⁴, 506 do. The largest gains are on the very sparse
si patterns:

- si4_r001_m200.04 goes from 75 s to 4 ms;
- si4_r001_m200.03 goes from 81 s to 36 ms;
- si4_r001_m200.09, the caching example above, goes from 41 s to 35 ms;
- g6-g59-normed goes from 87 s to 0.12 s.

Every count agrees with plain counting on the 410 instances both finish, and with S = 10⁴ on
the 506 both of those finish.

On the 63 finished instances that never split, it is 1.20× slower (at worst 1.36×), which is
the cost of looking for components. That cost loses the one instance only plain counting
finishes, si2_r005_m600.08, which takes 91 s plainly.

Where it does split, the worst case of note is gen-m15-LV-g23 at 1.7×. Nearly all of its 2.5
million splits are cache hits, so the time goes on building cache keys. 24 finished
instances emptied their caches, at most five times.

Over all 14,621 standard pattern–target pairs, with a 100 s timeout, there were no
disagreements. The number of finished instances went from 13,441 to 13,509: 73 more, all but
a few of them si, and 5 fewer.
