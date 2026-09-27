The Glasgow Subgraph Solver
===========================

This is a solver for subgraph isomorphism (induced and non-induced) problems, based upon a series of
papers by subsets of Blair Archibald, Ciaran McCreesh, Patrick Prosser and James Trimble at the
University of Glasgow, and Fraser Dunlop and Ruth Hoffmann at the University of St Andrews. A clique
decision / maximum clique solver is also included.

If you use this software for research, please cite [icgt/McCreeshP020]. If you use this solver in a
non-research setting, please get in touch if you can. This software is an output of taxpayer funded
research, and it is very helpful for us if we can demonstrate real-world impact when we write grant
applications.

Please contact [Ciaran McCreesh](mailto:ciaran.mccreesh@glasgow.ac.uk) with any queries.

Compiling
---------

To build, you will need:
- C++20 compiler, such as GCC 10.3.
- gmp, installation instructions can be found [here](https://gmplib.org/)

```shell
cmake -S . -B build
cmake --build build
```

Running
-------

To run:

```shell session
$ ./build/glasgow_subgraph_solver pattern-file target-file
```

If you would like induced subgraph isomorphisms rather than non-induced (that is, if non-adjacent
vertices must be mapped to non-adjacent vertices), you must request it:

```shell session
$ ./build/glasgow_subgraph_solver --induced pattern-file target-file
```

The default mode is to display the first found solution, or to prove unsatisfiability if no solution
exists. To count or print all solutions, use one of:

```shell session
$ ./build/glasgow_subgraph_solver [ --count-solutions | --print-all-solutions ] pattern-file target-file
```

Note that printing all solutions can be exponentially slower than counting solutions.

The solver supports parallel search. Usually you should enable this, as follows:

```shell session
$ ./build/glasgow_subgraph_solver --parallel ...
```

Note that parallel search, in its default configuration, is non-deterministic.

Preprocessing is normally done in full before search starts. As an experimental alternative,
staged solving does only the cheap filtering first, and builds the more expensive supplemental
graphs only if a short first round of search does not solve the instance:

```shell session
$ ./build/glasgow_subgraph_solver --staged pattern-file target-file
```

This can help a lot on instances that are easy, and costs a little on instances that are not. It
is currently sequential only (so it cannot be combined with `--parallel`), and with proof logging
it cannot yet be combined with `--count-solutions`.

If the target's vertices or edges carry costs, the solver can find a mapping of least total cost,
where a mapping costs the sum of the costs of the target vertices and target edges it uses:

```shell session
$ ./build/glasgow_subgraph_solver --minimise-cost --format json pattern.json target.json
```

It reports `cost = ...`, and `optimal = true` if the search finished. Costs, and pairs of vertices
joined by several edges with different labels, can only be given in the [JSON
format](#the-gss-graph-json-format). A multigraph works without `--minimise-cost` too, as a
decision problem. Both need an injective, non-induced mapping; minimising also needs sequential
search without restarts (the default when minimising), and cannot be combined with counting. Both
can be proof logged (`--prove`, with `--no-clique-detection`), including the cost bound; see
[the proof logging notes](dev_docs/proof-logging.md#minimising-cost-and-multigraphs). The supplemental graphs are not used on these instances, since every edge
becomes a vertex in the graphs the solver searches.

`tools/scene_graph_csv_to_json.py` converts scene graphs in the CSV dialect of the graph3
benchmark (a fourth column of confidences, and parallel edges with different labels) into this
JSON format, with costs of `round(-log(confidence) × 10^6)`.

File Formats
------------

We try to auto-detect the input format, but it's best to specify it using, for example:

```shell session
$ ./build/glasgow_subgraph_solver --format lad pattern-file target-file
```

In particular, note that auto-detection can easily fail if, for example, the first vertex in the
graph has no neighbours.  We can read LAD, Labelled LAD (labels on vertices, and optionally also on
edges), CSV, DIMACS 2, and [the gss-graph JSON format](#the-gss-graph-json-format) formatted graphs. [The LAD
format](https://perso.liris.cnrs.fr/christine.solnon/SIP.html) is a nice simple choice. If you need
to support named vertices, labels on vertices and / or edges, or directed edges, consider using the
CSV format. To specify a directed edge, use a greater-than sign rather than a comma as the delimiter
between the first two columns.  To specify an edge label, include a third column in the file. To
specify a vertex label, leave the second column empty and use the third column for the label. For
example, the following describes a graph with four vertices, with colours for edge labels and shapes
for vertex labels.

```
first>second,red
second>first,blue
first,third,purple
first>first,green
first,,circle
second,,circle
third,,square
fourth,,square
```

Labels are all or nothing, separately for vertices and for edges: if any vertex in a file is
labelled then every vertex must be, and likewise for edges, so a file labelling only some of its
vertices or only some of its edges is rejected. Labelling the vertices does not oblige you to label
the edges, or the other way around. There is no way to write "this element has no label constraint"
in a file that uses labels of that kind, and treating an undeclared label as the empty string would
quietly restrict that element to unlabelled target elements instead.

Labelling the edges of an undirected graph does not make it directed: use the greater-than delimiter
for that.

### The gss-graph JSON format

Every format above works out what kind of graph it is holding from whichever syntax happens to turn
up in the file. That is unavoidably ambiguous in places: whether a graph is directed depends on
whether any line used `>`, an isolated vertex needs the trailing-comma idiom, and auto-detection
cannot always tell LAD from Labelled LAD. The `json` format instead *declares* those properties, and
validates strictly, so that anything ambiguous is an error rather than a guess:

```json
{"format": "gss-graph", "version": 1, "directed": false,
 "vertices": 4, "edges": [[0, 1], [1, 2], [2, 3]]}
```

Named and labelled, with a self-loop:

```json
{
  "format": "gss-graph", "version": 1,
  "directed": false,
  "vertices": [
    {"name": "c1", "label": "C"},
    {"name": "n1", "label": "N"},
    {"name": "o1", "label": "O"}
  ],
  "edges": [
    ["c1", "n1", "single"],
    ["n1", "o1", "double"],
    ["c1", "c1", "single"]
  ]
}
```

`vertices` is either a count, giving that many anonymous vertices, or an array whose entries are
names (`["a", "b"]`, exact sugar for `[{"name": "a"}, {"name": "b"}]`) or objects with an optional
`name`, `label` and `cost`. `edges` entries are either `[from, to]` / `[from, to, label]` or objects
with `from`, `to` and an optional `label` and `cost`. The array form is fixed at two or three elements permanently:
every key added in future will live only in the object form, so the compact spelling can never come
to mean something new.

The rules, each of which is a hard failure naming what was wrong:

- `format` and `version` are required. A reader refuses a version above the one it knows, naming both.
- `directed` is required and never inferred. This is the rule that stops a labelled undirected graph
  turning into a directed one.
- In an undirected graph `[u, v]` and `[v, u]` are the same edge, so giving both is a duplicate
  rather than being quietly idempotent as it is in CSV. A self-loop `[v, v]` is exactly one edge in
  directed and undirected graphs alike, never implicitly doubled.
- Labels are all or nothing per element type, as for CSV above, and `""` is a real label rather than
  an absent one.
- An endpoint addresses a vertex by index when it is an integer and by name when it is a string, and
  a file has to pick one and keep to it. This is where JSON beats every text format here: the vertex
  named `"1"` and the vertex at index `1` are different tokens, so numeric vertex names stop being a
  trap.
- Any unrecognised key is an error, so a misspelled `"wieght"` fails loudly instead of being
  silently dropped. Keys beginning `x-` are reserved for third-party annotation and always ignored.
- Vertex names are unique, and an isolated vertex is simply one listed in `vertices` and absent from
  `edges`, with no special syntax.

`"multigraph": true` allows parallel edges with different labels, by switching the uniqueness rule
from the `(from, to)` pair to the `(from, to, label)` triple. Several edges with the *same* label
between one pair would need an edge object's `multiplicity`, which this build refuses rather than
quietly merging them. Only the homomorphism solver accepts a multigraph.

A `cost` is a 64-bit signed integer. Like labels, costs are all or nothing per element type, since an
absent cost is not the same as 0, and a cost can only appear in the object form. Costs are what
`--minimise-cost` minimises; without it they are read and ignored. A scene graph whose weights are
probabilities, for instance, would give each element the cost `round(-log(p) × 10^6)`.

`convert_to_json` writes any graph the other readers accept into this format, and unlike
`convert_to_lad` it has nothing to refuse, since the format carries directedness, names and both
kinds of label:

```shell session
$ ./build/convert_to_json --format lad my-graph.lad > my-graph.json
```

Its output is canonical, so converting it again gives the same bytes. Reading a file and writing it
back preserves every property the format carries, which is what makes the claim that the format is
unambiguous something that can be tested rather than just asserted.

Symmetries
----------

Symmetry elimination support is currently very experimental, and is probably only useful for solution
counting. Symmetry-breaking constraints can be supplied manually: a "less than" constraint forces one
pattern vertex to map below another, and an "occurs less than" constraint orders how often target
vertices are used. Vertices are named as in the input files.

```shell session
$ ./build/glasgow_subgraph_solver --pattern-less-than 'a<b' --target-occurs-less-than '0<1' \
    --count-solutions pattern-file target-file
```

Automatic detection of these constraints (it previously shelled out to the GAP computer algebra
system) has been removed pending a more robust replacement.

Proof Logging
-------------

As a highly experimental feature, the solver can output a proof log. First, install the following
program:

* VeriPB from https://gitlab.com/MIAOresearch/software/VeriPB

And then you can produce and verify a log like this:

```shell session
$ ./build/glasgow_subgraph_solver --induced --no-supplementals --no-clique-detection --no-nds \
    --prove myproof --format lad pattern-file target-file
$ veripb myproof.opb myproof.pbp
```

This writes the pseudo-Boolean model to `myproof.opb` and the proof to `myproof.pbp`. Refutation
(unsatisfiable), decision (satisfiable), and counting/enumeration proofs (`--count-solutions`,
`--enumerate`, `--print-all-solutions`) all verify, including loop-preserving mappings. Most other
features are not yet supported with proof logging — this is a "not yet implemented" problem, not a
fundamental restriction. See [dev_docs/proof-logging.md](dev_docs/proof-logging.md) for the supported
option combinations, the conclusions produced, and how to check proofs with the formally verified
CakePB checker.

Clique Solving
--------------

To run the clique solver, use:

```shell session
$ ./build/glasgow_clique_solver graph-file
```

Details on the Algorithms
-------------------------

The subgraph solver is a constraint programming style backtracker, which recursively builds up a
mapping from pattern vertices to target vertices. It includes inference based upon paths (not just
adjacency) and neighbourhood degree sequences, has a fast all-different propagator, and uses
sophisticated variable- and value-ordering heuristics to direct a slightly-random restarting search.

Chronologically, our first subgraph isomorphism solver is [cp/McCreeshP15]. We introduced new
variants of this solver in [lion/KotthoffMS16], and described a refactored version (which can solve
an optimisation variant of the problem) in [aaai/HoffmannMR17]. We also investigated search ordering
heuristics in more detail in [jair/McCreeshPST18], and [cpaior/ArchibaldDHMPT19] describes its new
restarting search algorithm. There is currently no paper describing the entire algorithm, but
[icgt/McCreeshP020] summarises the main aspects of it.

The clique solver (with its default configuration) is a branch and bound solver that uses a greedy
colouring both as the bound function, and as a branching heuristic. It is based upon the "domains of
size two first" variant described in [cp/McCreeshP14], which is in turn derived from the "MCSa1"
algorithm described by [algorithms/Prosser12] combined with the bit-parallelism techniques discussed
by [ol/SegundoMRH13]; this in turn is a simplification of "MCS" described by [walcom/TomitaSHTW10].
The solver also incorporates the fast clique detection technique described by [jco/BatsynGMP14].

Funding Acknowledgements
------------------------

This work was supported by the Engineering and Physical Sciences Research Council (grant numbers
EP/P026842/1, EP/M508056/1, and EP/N007565). This work used the Cirrus UK National Tier-2 HPC
Service at EPCC (http://www.cirrus.ac.uk) funded by the University of Edinburgh and EPSRC
(EP/P020267/1).

References
----------

* [walcom/TomitaSHTW10]: https://dblp.org/rec/html/conf/walcom/TomitaSHTW10
  **walcom/TomitaSHTW10**:
  Etsuji Tomita, Yoichi Sutani, Takanori Higashi, Shinya Takahashi, Mitsuo Wakatsuki:
  A Simple and Faster Branch-and-Bound Algorithm for Finding a Maximum Clique. WALCOM 2010: 191-203.
  DBLP: [walcom/TomitaSHTW10]

* [algorithms/Prosser12]: https://dblp.org/rec/html/journals/algorithms/Prosser12
  **algorithms/Prosser12**:
  Patrick Prosser: Exact Algorithms for Maximum Clique: A Computational Study. Algorithms 5(4):
  545-587 (2012). DBLP: [algorithms/Prosser12].

* [ol/SegundoMRH13]: https://dblp.org/rec/html/journals/ol/SegundoMRH13
  **ol/SegundoMRH13**:
  Pablo San Segundo, Fernando Matía, Diego Rodríguez-Losada, Miguel Hernando: An improved bit
  parallel exact maximum clique algorithm. Optimization Letters 7(3): 467-479 (2013). DBLP:
  [ol/SegundoMRH13].

* [cp/McCreeshP14]: https://dblp.org/rec/html/conf/cp/McCreeshP14
  **cp/McCreeshP14**:
  Ciaran McCreesh, Patrick Prosser: Reducing the Branching in a Branch and Bound Algorithm for the
  Maximum Clique Problem. CP 2014: 549-563. DBLP: [cp/McCreeshP14].

* [jco/BatsynGMP14]: https://dblp.org/rec/html/journals/jco/BatsynGMP14
  **jco/BatsynGMP14**:
  Improvements to MCS algorithm for the maximum clique problem. J. Comb. Optim. 27(2): 397-416
  (2014). DBLP: [jco/BatsynGMP14]

* [cp/McCreeshP15]: https://dblp.org/rec/html/conf/cp/McCreeshP15
  **cp/McCreeshP15**:
  Ciaran McCreesh, Patrick Prosser: A Parallel, Backjumping Subgraph Isomorphism Algorithm Using
  Supplemental Graphs. CP 2015: 295-312. DBLP: [cp/McCreeshP15].

* [lion/KotthoffMS16]: https://dblp.org/rec/html/conf/lion/KotthoffMS16
  **lion/KotthoffMS16**:
  Lars Kotthoff, Ciaran McCreesh, Christine Solnon: Portfolios of Subgraph Isomorphism Algorithms.
  LION 2016: 107-122. DBLP: [lion/KotthoffMS16].

* [aaai/HoffmannMR17]: https://dblp.org/rec/html/conf/aaai/HoffmannMR17
  **aaai/HoffmannMR17**:
  Ruth Hoffmann, Ciaran McCreesh, Craig Reilly: Between Subgraph Isomorphism and Maximum Common
  Subgraph. AAAI 2017: 3907-3914. DBLP: [aaai/HoffmannMR17].

* [jair/McCreeshPST18]: https://dblp.org/rec/html/journals/jair/McCreeshPST18
  **jair/McCreeshPST18**:
  Ciaran McCreesh, Patrick Prosser, Christine Solnon, James Trimble: When Subgraph Isomorphism is
  Really Hard, and Why This Matters for Graph Databases. J. Artif. Intell. Res. 61: 723-759 (2018).
  DBLP: [jair/McCreeshPST18].

* [cpaior/ArchibaldDHMPT19]: http://dblp.org/rec/html/conf/cpaior/ArchibaldDHMP019
  **cpaior/ArchibaldDHMPT19**:
  Blair Archibald, Fraser Dunlop, Ruth Hoffmann, Ciaran McCreesh, Patrick Prosser and James Trimble:
  Sequential and Parallel Solution-Biased Search for Subgraph Algorithms. CPAIOR 2019: 20-38.
  DBLP: [cpaior/ArchibaldDHMPT19].

* [icgt/McCreeshP020]: http://dblp.org/rec/html/conf/gg/McCreeshP020
  **icgt/McCreeshP020**:
  Ciaran McCreesh, Patrick Prosser, James Trimble:
  The Glasgow Subgraph Solver: Using Constraint Programming to Tackle Hard Subgraph Isomorphism
  Problem Variants. ICGT 2020: 316-324.
  DBLP: [icgt/McCreeshP020].

<!-- vim: set tw=100 spell spelllang=en : -->
