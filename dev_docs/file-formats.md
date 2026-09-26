# File formats

How graphs get into the solver, what each format can and cannot say, and the design rule the newer
parts of this code follow. For *using* the formats see the [README](../README.md#file-formats); for
where the readers sit in the tree see [architecture.md](architecture.md).

## The rule: declare, do not infer

Every format here except JSON works out what kind of graph it is holding from whichever syntax
happened to turn up in the file. Whether a CSV graph is directed depends on whether any line used
`>`; whether it has edge labels depends on whether any label was non-empty; an isolated vertex needs
the trailing-comma idiom; and auto-detection sometimes cannot tell LAD from Labelled LAD at all.

That is not merely untidy. It has produced real bugs, twice in the same week:

- **An undirected graph with edge labels loaded as directed** (#86, fixed in #87). The labelled
  branch of `read_csv` built each edge from two `add_directed_edge` calls, and that function set the
  directed flag unconditionally, so adding a label to a file flipped `directed()`.
- **Partial vertex labelling became "labelled with the empty string"** (#86, fixed in #87).
  `vertex_labels` was indexed with `operator[]`, so an undeclared label default-constructed to `""`,
  and since label matching is exact that silently confined the vertex to unlabelled targets.

The rule that came out of it: **a graph-level property is declared, never reconstructed from how the
graph happened to be built.** Concretely, in `InputGraph`:

- Directedness is a constructor argument. It defaults to `false`, and `add_directed_edge` *requires*
  it rather than setting it, throwing `std::logic_error` otherwise. So an asymmetric edge set can
  never sit behind a `directed()` of `false`, and a directed graph with no edges — which nothing
  about an edge list can tell you about — is still directed.
- `add_edge(a, b, label)` adds a labelled *undirected* edge. Use it rather than a pair of
  `add_directed_edge` calls, which is what caused #86.
- `vertex_has_name(v)` distinguishes a name that was set from the index fallback `vertex_name(v)`
  returns. A writer needs this so that it does not invent names for anonymous vertices.

When adding a property to `InputGraph`, apply the same test the JSON format's spec uses: a property
may be optional only if its default is the restrictive reading, the one that causes an error rather
than a silent reinterpretation. `directed` fails that test — both readings parse cleanly and give
different answers — which is why it is mandatory in the file format and why its C++ default is
backed by a throw.

## What each format can express

| Format | Names | Vertex labels | Edge labels | Directed | Writer |
|---|---|---|---|---|---|
| `csv` | yes | yes | yes | per-edge `>` | no |
| `lad` | via `0=name` | no | no | no | `convert_to_lad`, unlabelled only |
| `vertexlabelledlad` | via `0=name` | yes | no | no | no |
| `labelledlad` | via `0=name` | yes | yes | always | no |
| `directedlad` | via `0=name` | no | no | always | no |
| `dimacs` | 1-indexed | no | no | no | no |
| `vfmcs` | no | some variants | no | per-variant | no |
| `json` | optional, per vertex | yes | yes | **declared** | `convert_to_json` |

Only JSON can express a multigraph (parallel edges with different labels) or integer costs on
vertices and edges, and both are declared: `"multigraph"` at the top level, and costs all or nothing
per element type, as labels are.

Two things worth knowing about the older readers. LAD and DIMACS name *every* vertex, LAD with its
own index as a string, so a partly named graph can only come from JSON or from building an
`InputGraph` by hand. And `labelledlad` is directed because the format is, not because it has
labels — `read_any_lad` branches on directedness explicitly to keep those two facts separate.

## Labels are all or nothing

In both CSV and JSON, if any vertex carries a label then every vertex must, and likewise for edges;
otherwise the file is a `GraphFileError`. The two kinds are independent — labelling the vertices does
not oblige you to label the edges.

The reason is that an absent label and `""` are different things, `""` being a real label, and label
matching is exact. Reading an undeclared label as `""` therefore silently constrains that element to
match only unlabelled ones, which is a wrong answer rather than a loud failure. This was the second
half of #86.

## The `gss-graph` JSON format

`gss/formats/json_graph.{hh,cc}`, `--format json`, auto-detected on a leading `{`. The
[README](../README.md#the-gss-graph-json-format) documents the format for users and lists the eight
validation rules; the design rationale is in #85. Points that matter when changing the code:

- **The positional edge form `[from, to]` / `[from, to, label]` is frozen at two or three elements
  permanently.** Every key added in future lives only in the object form. Do not extend the array.
- **The writer's output is canonical**, and that is load-bearing rather than cosmetic: it is what
  lets `json_graph_test` assert that reading and writing is a fixed point in both the graph and the
  bytes, over graphs from the CSV and LAD readers as well as JSON's own. That round trip is the only
  real evidence that the format is unambiguous, so if you add a property to the format, add it to
  the writer in the same change or the test silently stops covering it.
- **`"multigraph": true` switches the uniqueness key from `(from, to)` to `(from, to, label)`.**
  `InputGraph` holds such a graph by keeping a short list of edges per endpoint pair, so
  `adjacent()` and `degree()` are unchanged (degree counts neighbours), and `edge_label()` throws
  on a multigraph rather than choosing one of several answers. Anything that reads edge labels
  pairwise must therefore either refuse a multigraph or not see one: the homomorphism solver
  rewrites a multigraph into a simple graph before anything else looks at it (see
  `gss/innards/reification.hh`), and the clique and common-subgraph solvers refuse one.
- **`multiplicity` is still a known key that this build refuses** for any value but 1. It is how
  the format would say "several edges with the same label", which `InputGraph` cannot hold.
- **Costs are integers and are all or nothing**, for the same reason labels are: reading an absent
  cost as 0 would make that element free, which is a wrong answer rather than a loud failure.
  `InputGraph` enforces this too: a costed graph's edges must be added with a cost, and reading a
  vertex cost that was never set throws. The array form cannot carry a cost, so the writer uses
  the object form for every edge of a graph with edge costs, and sorts edges by
  `(from, to, label)` so that the output stays canonical whatever order parallel edges were added
  in.
- **Unknown keys are an error, except `x-`.** A misspelled key must fail loudly rather than being
  dropped.
- `nlohmann/json` is included only by `json_graph.cc` and linked `PRIVATE`. Keep it out of headers:
  it is slow to compile and should cost one translation unit.

## Who depends on these properties

Changing what a reader reports is not a local change. `directed()` and `has_edge_labels()` both
steer the solver, and they used to be entangled: **every edge-labelled graph was also directed**,
because the only way to get edge labels was through `add_directed_edge`. Fixing #86 broke that
invariant and immediately exposed two callers that had been relying on it:

- `HomomorphismModel` allocated `forward_target_graph_rows` / `reverse_target_graph_rows` only for a
  directed pattern, while `HomomorphismSearcher` always instantiates the directed specialisation of
  `propagate_adjacency_constraints` once there are edge labels, since it has to compare each edge's
  forward and reverse labels separately. An edge-labelled undirected graph then indexed an empty
  vector — a segfault, caught only because the new tests solved such an instance.
- `solve_sip_by_decomposition` declines to decompose a *directed* pattern but not a labelled one, and
  rebuilt the reduced pattern with `add_directed_edge(..., "")`, discarding every edge label.
  Previously unreachable; afterwards `--decomposition` on a labelled pattern silently returned no
  solutions.

So: after changing a graph-level property, grep for every use of it and check whether the code is
branching on it for a *representational* reason or a *correctness* one.
`grep -rn 'directed()\|has_edge_labels()' gss src` is the whole list, and at around twenty call
sites outside the tests it can be read in one sitting.

A related trap, since it was in this code for five years: `homomorphism_model.cc` used to skip target
edges whose label was the string `"unlabelled"`, a sentinel from a representation deleted in the very
commit that introduced the filter. A graph using it as a real label lost those edges from directed
propagation and came back unsatisfiable (#88). Do not compare user-supplied label text against magic
strings.

## Adding a format

1. `gss/formats/<name>.{hh,cc}`, taking `std::istream &&` and a filename, returning `InputGraph`,
   throwing `GraphFileError` with the offending line or key named.
2. Register the source in `gss/CMakeLists.txt` and the name in `read_file_format`, plus auto-detection
   in `detect_format` if the format has an unambiguous signature. Note that the CSV regex is very
   permissive, so anything new has to be tested before it.
3. Declare graph-level properties; do not infer them. If the format cannot express one, say so in the
   table above rather than guessing on the file's behalf.
4. Add the format to the `--format` help in all three solvers and both converters.
5. Tests in `gss/formats/formats_test.cc`, or a dedicated binary if the format needs enough rejection
   tests to swamp that file (as JSON did). A writer, if you add one, earns a round-trip fixed-point
   test — see `json_graph_test.cc`.
