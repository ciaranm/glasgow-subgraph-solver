# Counting statements

The file format in which `--decomposition` states its counting argument, and the PB
translation of it that VeriPB checks against the input. This is the reference for anything
that reads or writes the format. For why the format exists and where it fits in the
certification pipeline, see [decomposition-counting.md](decomposition-counting.md).

`glasgow_subgraph_solver --count-solutions --decomposition --decomposition-statement FILE`
writes one, and `count_statement FILE [--opb OUT] [--side-conditions OUT]` reads it, checks it,
prints its count, and writes its translations. A statement describes a pseudo-Boolean formula
whose number of solutions is cheap to compute. Three things are done with it, and one day the
programs doing them will be formally verified:

- the **counter** computes the count;
- the **translator** writes the formula as an OPB file, which VeriPB checks is
  equienumerable with the solver's encoding of the input;
- the translator also writes the **side-condition instance**, an OPB file that must be
  unsatisfiable for the count to be right.

The format is therefore designed to keep a verified parser and checker as small as
possible: no names, no strings, no comments, no optional parts, no forward references, and
nothing the reader has to infer that it could instead check.

## Lexical structure

A statement is a single s-expression, in the same deliberately small dialect as the `.scp`
files of the Glasgow Constraint Solver. Whitespace (space, tab, newline, carriage return,
form feed, vertical tab) and the parentheses `(` and `)` are the only delimiters; every
other maximal run of characters is an atom. There are no comments and no quoted strings.

Two kinds of atom are meaningful:

- **keywords**, compared as exact strings: `version`, `domain`, `nodes`, `root`, `cube`,
  `and`, `or`, `minus`;
- **numbers**: a non-empty run of the digits `0`–`9`, with no sign and no leading zero
  (so `0` is the only number starting with `0`), read as a natural number in decimal.

Anything else is an error.

## Grammar

```
statement ::= "(" version domain nodes root ")"
version   ::= "(" "version" "1" ")"
domain    ::= "(" "domain" NAT NAT ")"
nodes     ::= "(" "nodes" node* ")"
node      ::= "(" NAT "cube" pair* ")"
           |  "(" NAT "and" NAT* ")"
           |  "(" NAT "or" NAT NAT* ")"
           |  "(" NAT "minus" NAT NAT ")"
pair      ::= "(" NAT NAT ")"
root      ::= "(" "root" NAT ")"
```

The four sections are all required and come in this order. The only version is `1`, and a
reader rejects any other, so that a future change to the format fails cleanly instead of
being misread.

`(domain m n)` declares m **variables**, numbered `0` to m − 1, each taking one of n
**values**, numbered `0` to n − 1. For subgraph isomorphism the variables are the pattern
vertices and the values are the target vertices, both numbered as the solver numbers them
internally, which for a LAD file is the file's own numbering.

Each node starts with its **id**. Ids must be `0, 1, 2, …` in the order the nodes appear,
and every child named by a node must have a smaller id. So the nodes form a DAG listed in
topological order, which the reader checks in a single pass. A child may be shared by
several parents.

## Meaning

An **assignment** is a function f from the variables to the values. Every node has a
**scope**, a set of variables, and is true or false of each assignment, depending only on
what the assignment does to its scope. Its **count** is the number of functions from its
scope to the values for which it is true.

| Node | Scope | True of f when | Count |
|---|---|---|---|
| `(i cube (u₁ t₁) … (uₖ tₖ))` | {u₁, …, uₖ} | f(uⱼ) = tⱼ for every j | 1 |
| `(i and c₁ … cₖ)` | the union of the cᵢ's scopes | every cⱼ is true | the product of the cⱼ's counts |
| `(i or c₁ … cₖ)` | the cⱼ's common scope | some cⱼ is true | the sum of the cⱼ's counts |
| `(i minus a b)` | a's and b's common scope | a is true and b is false | a's count minus b's count |

The **count of the statement** is the count of the node named by `(root r)`.

The counts in the last column are only right under conditions on the nodes, which come in
two kinds: well-formedness, which the reader checks as it parses, and side conditions,
which are about every assignment and so need a proof.

### Well-formedness

A reader must reject a statement unless all of these hold.

- Every variable is less than m and every value is less than n.
- In a cube, the variables are strictly increasing. (So a cube names each of its variables
  once, and a cube is determined by its list of pairs.)
- The children of an `and` have pairwise disjoint scopes.
- An `or` has at least one child, and all its children have the same scope.
- The two children of a `minus` have the same scope.
- The root's scope is every variable, `0` to m − 1.

An `and` may have no children, and a cube no pairs: either is then true of everything, has
an empty scope, and counts 1.

There is no way to write a node that is never true. That is deliberate: an empty `or` would
have no children to take its scope from, and giving it one would need either a scope
annotation or a wildcard scope, both of which the checker would then have to handle
everywhere. A count of zero is still expressible whenever something can be subtracted from
itself, and when the solver finds no way to place some component of the pattern at all,
the instance is better certified by an ordinary unsatisfiability proof anyway.

### Side conditions

- The children of an `or` must be **pairwise disjoint**: no assignment makes two of them
  true. Otherwise the sum over-counts.
- For `(i minus a b)`, **b must imply a**: every assignment that makes b true makes a true.
  Otherwise the difference under-counts, or goes negative.

These are about all mn assignments, not just the ones that are solutions to the instance,
and the main VeriPB proof does not check them. A statement that violates one can still
translate to a formula equienumerable with the input. The output check passes, and the
count printed is wrong.

Some are checked syntactically by the counter instead of being proved:

- An `or` all of whose children are cubes is disjoint if no two of its children are the
  same cube. Its children have equal scopes, so two different cubes assign some variable
  different values. This covers every leaf of a decomposition, which is almost all of its
  nodes.

Every other `or`, and every `minus`, goes into the side-condition instance below. A counter
that finds a `minus` whose second child's count exceeds its first's has found a violated
side condition, and must say so rather than print a negative count.

## Translation to pseudo-Boolean form

The translator turns a well-formed statement into an OPB file whose solutions, projected
onto the x variables, correspond one-to-one with the assignments that make the root true.
Every name and every row is fixed by the statement, so two correct translators produce the
same file, byte for byte.

### Names

- `x<u>_<t>`: variable u takes value t. This is the name the solver's proof encoding uses
  for "pattern vertex u maps to target vertex t" when vertices are named by their numbers,
  as they are in LAD files, and it is the name that `cake_pb_iso`, which reads only LAD,
  uses too. For a format with vertex names, such as CSV, the solver's encoding names its
  variables after the names instead, so the two do not line up; that is a problem for the
  proof plumbing (phase 4) to solve, not for this format.
- `g<i>`: the gate of node i, true exactly when node i is true.
- `v<i>`: in the side-condition instance only, "node i's side condition is violated".

Every variable in the solver's encoding and proofs is written with a leading `x`, so these
cannot clash with anything it writes.

### The counting formula

In this order:

1. `preserved: x0_0 x0_1 … x0_<n−1> x1_0 … x<m−1>_<n−1> ;`, every x variable, ordered by
   variable and then by value.
2. For each variable u in increasing order, the exactly-one group, as two rows:
   `1 xu_0 1 xu_1 … 1 xu_<n−1> >= 1 ;` and then the same terms `<= 1 ;`.
3. For each node i in order, its gate definition. Writing k for the number of literals or
   children, ℓⱼ for the literals of a cube (`xu_t`) or the gates of the children (`gc`) of
   an `and`:

   | Node | Rows |
   |---|---|
   | cube or `and`, k ≥ 1 | `k ~gi 1 ℓ₁ … 1 ℓₖ >= k ;` then `1 gi 1 ~ℓ₁ … 1 ~ℓₖ >= 1 ;` |
   | cube or `and`, k = 0 | `1 gi >= 1 ;` |
   | `or` over c₁ … cₖ | `1 ~gi 1 gc₁ … 1 gcₖ >= 1 ;` then `k gi 1 ~gc₁ … 1 ~gcₖ >= k ;` |
   | `minus a b` | `1 ~gi 1 ga >= 1 ;` then `1 ~gi 1 ~gb >= 1 ;` then `1 gi 1 ~ga 1 gb >= 1 ;` |

   In each case the rows before the last say "the gate implies its definition", and the
   last row says "the definition implies the gate".
4. `1 g<root> >= 1 ;`

Numbers are written in decimal, terms are separated by single spaces, and every row ends
`;` and a newline. There is no `*` header line.

VeriPB compares an output file with its database as a set of normalised constraints, so the
order of the rows, and the choice between `>= -1` with negated coefficients and `<= 1`,
does not matter to the check. They are fixed anyway, so that the file is reproducible and
diffable.

### The side-condition instance

The rows of item 2 and item 3 above, with no `preserved:` line and no root row, then for
each node i whose side condition is not discharged syntactically, in node order:

| Node | Rows |
|---|---|
| `or` over c₁ … cₖ | `2 ~vi 1 gc₁ … 1 gcₖ >= 2 ;` |
| `minus a b` | `1 ~vi 1 gb >= 1 ;` then `1 ~vi 1 ~ga >= 1 ;` |

and finally `1 vi₁ 1 vi₂ … >= 1 ;` over all of those i. (If there are none, the translator
writes no instance, since there is nothing to prove.)

The instance must be unsatisfiable. Every assignment extends to exactly one assignment of
the gates that satisfies items 2 and 3, so if some assignment violated node i's condition,
setting vᵢ true and the other v's false would satisfy the instance.

## Why the count is right

Suppose the statement is well-formed, every syntactically discharged condition holds, and
the side-condition instance is unsatisfiable. Write F for the counting formula and X for
its x variables.

1. The rows of item 2 make an assignment to X that satisfies them the same thing as a
   function f from variables to values.
2. Each gate's rows say gi ↔ (its definition), and the definitions only mention x variables
   and gates of earlier nodes. So by induction on i, each gᵢ is determined by f, and is true
   exactly when node i is true of f. Hence the solutions of F, projected onto X, are exactly
   the functions f that make the root true.
3. By induction on i, node i's count is the number of functions from its scope to the
   values that make it true: a cube fixes its whole scope; an `and`'s children depend on
   disjoint parts of its scope, so the functions that satisfy all of them are tuples of
   functions satisfying each; an `or`'s children are disjoint and depend on the same scope,
   so the sets add; and for a `minus`, b's set is a subset of a's, so the sets subtract.
4. The root's scope is all of the variables, so its count is the number of functions f
   that make it true, which by step 2 is the number of solutions of F projected onto X.

VeriPB then checks, with `EQUIENUMERABLE`, that F has as many solutions projected onto X as
the solver's encoding of the instance does.

## An example

Two disjoint edges, 0–1 and 2–3, into the path 0–1–2–3. Six embeddings of each edge, so
36 pairs, of which 28 overlap. Abbreviated:

```
(
  (version 1)
  (domain 4 4)
  (nodes
    (0 cube (0 0) (1 1))
    (1 cube (0 1) (1 0))
    …
    (6 or 0 1 2 3 4 5)
    …
    (13 or 7 8 9 10 11 12)
    (14 and 6 13)
    (15 cube (0 1) (1 0) (2 1) (3 2))
    …
    (49 or 19 24 29 34 41 48)
    (50 minus 14 49)
  )
  (root 50)
)
```

Nodes 6 and 13 are the two edges' embeddings, node 14 is all 36 pairs, nodes 19 to 48 are the
six ways the edges can overlap (gluing one pair of endpoints gives a path of three, which
embeds 4 ways; gluing both gives one edge, 6 ways), node 49 is the 28 overlapping pairs, and
node 50 is the 8 that do not overlap. The side conditions left to prove are node 49's
disjointness and node 50's containment.
