# Counting statements

The file format in which a solver states a counting argument, and the PB translation of it
that VeriPB checks against the input. This is the reference for anything that reads or
writes the format. For why the format exists and where it fits in the certification
pipeline, see [split-counting.md](split-counting.md).

`glasgow_subgraph_solver --count-solutions --split-components N --split-statement FILE --prove
P` writes one, and `count_statement FILE [--opb OUT]` reads it, checks it, prints its count,
and writes its translation. A statement describes a pseudo-Boolean formula whose number of
solutions is cheap to compute. Two things are done with it, and one day the programs doing
them will be formally verified:

- the **counter** checks the statement and computes its count;
- the **translator** writes the formula as an OPB file, the **tagged formula**, over the
  statement's own fresh variables.

VeriPB then checks one proof, from the solver's encoding of the input to the tagged formula,
ending `OUTPUT EQUIENUMERABLE`. That proof is the only thing the solver has to get right.
Whatever the statement says, the tagged formula has exactly as many solutions as the
counter's count. So a statement that is wrong about the input fails the check, rather than
passing it with a wrong count.

The format is designed to keep a verified parser and checker as small as possible: no
names, no strings, no comments, no optional parts, no forward references, and nothing the
reader has to infer that it could instead check. Nor does it know anything about the problem
it came from. Its variables are its own, and connecting them to the input's encoding is the
job of the proof (see "From an input to the tagged formula"). So the same counter and
translator can serve subgraph counting, constraint programming, or anything else that can
state its argument this way.

## Lexical structure

A statement is a single s-expression, in the same deliberately small dialect as the `.scp`
files of the Glasgow Constraint Solver. Whitespace (space, tab, newline, carriage return,
form feed, vertical tab) and the parentheses `(` and `)` are the only delimiters, and every
other maximal run of characters is an atom. There are no comments and no quoted strings.

Two kinds of atom are meaningful:

- **keywords**, compared as exact strings: `version`, `domain`, `nodes`, `root`, `cube`,
  `and`, `or`;
- **numbers**: a non-empty run of the digits `0`–`9`, with no sign and no leading zero, so
  `0` is the only number starting with `0`. A number is read as a natural number in
  decimal, with no bound on its size, since counts get large.

Anything else is an error.

## Grammar

```
statement ::= "(" version domain nodes root ")"
version   ::= "(" "version" "2" ")"
domain    ::= "(" "domain" NAT* ")"
nodes     ::= "(" "nodes" node* ")"
node      ::= "(" NAT NAT "cube" pair* ")"
           |  "(" NAT NAT "and" NAT* ")"
           |  "(" NAT NAT "or" NAT NAT* ")"
pair      ::= "(" NAT NAT ")"
root      ::= "(" "root" NAT ")"
```

The four sections are all required and come in this order. The only version is `2`, and a
reader rejects any other, so that a change to the format fails cleanly instead of being
misread.

`(domain n₀ n₁ … nₘ₋₁)` declares m **variables**, numbered `0` to m − 1, where variable u
takes one of the **values** `0` to nᵤ − 1. What they stand for is up to whoever writes the
statement. The subgraph solver uses one variable per pattern vertex, numbered as the solver
numbers them, each with one value per target vertex.

A node is `(id count kind …)`.

- The **id** must be `0, 1, 2, …` in the order the nodes appear.
- Every child must have a smaller id than its parent, so the nodes form a DAG listed in
  topological order, which the reader checks in a single pass.
- A child may be shared by several parents.
- The **count** is the node's count, as defined below. It is redundant: the reader
  computes it and rejects the statement if the two differ. It is there so that a person,
  or a debugger, reading a statement can see what each node contributes without
  recomputing anything.

## Meaning

An **assignment** is a function f taking each variable u to one of its values. Every node
has a **scope**, a set of variables, and a finite set of **objects**. Each object
determines an assignment to the node's scope, but different objects may determine the same
assignment. The node's **count** is the number of its objects.

| Node | Scope | Objects | Count |
|---|---|---|---|
| `(i n cube (u₁ t₁) … (uₖ tₖ))` | {u₁, …, uₖ} | one, assigning tⱼ to each uⱼ | 1 |
| `(i n and c₁ … cₖ)` | the union of the cⱼ's scopes | a tuple of one object of each cⱼ, assigning what each of them assigns | the product of the cⱼ's counts |
| `(i n or c₁ … cₖ)` | the cⱼ's common scope | a pair (j, o) of a position j and an object o of cⱼ, assigning what o assigns | the sum of the cⱼ's counts |

So an `or` is a **tagged** sum: an assignment that two of its children both allow gives two
objects, one through each child. Nothing needs the children to be disjoint for the count to
be right. If they overlap where the input has only one solution, the tagged formula has
more solutions than the input, and the proof cannot be completed.

The **count of the statement** is the count of the node named by `(root r)`, which is the
number of objects of the root. The statement claims that the input has this many
solutions.

### Well-formedness

A reader must reject a statement unless all of these hold.

- Every variable u is less than m, and every value given to it is less than nᵤ.
- In a cube, the variables are strictly increasing, so a cube names each of its variables
  once.
- In an `and` or an `or`, the children are strictly increasing, so no node is a child of the
  same parent twice.
- The children of an `and` have pairwise disjoint scopes.
- An `or` has at least one child, all its children have the same scope, and that scope is
  not empty.
- Every node's count is as defined above.
- The root's scope is every variable, `0` to m − 1.

A node need not be reachable from the root. Nothing counts it, and its translation keeps it
inactive. A solver that writes its statement as it searches cannot always know in advance
which nodes it will end up using.

An `and` may have no children, and a cube no pairs. Either then has an empty scope and one
object.

There is no way to write a node with no objects. An empty `or` would have no children to
take its scope from. When a solver finds no way to place some part of the pattern at all,
the instance is better certified by an ordinary unsatisfiability proof anyway.

## Translation to pseudo-Boolean form

The translator turns a well-formed statement into the tagged formula. Its solutions,
projected onto its preserved variables, correspond one-to-one with the root's objects. Every
name and every row is fixed by the statement, so two correct translators produce the same
file, byte for byte.

### Names

- `y<u>_<t>`: variable u takes value t.
- `a<i>`: node i is **active**, that is, part of the object.
- `t<i>_<j>`: `or` node i takes its child in position j, counting from 0.

Every variable in the subgraph solver's encoding is written with a leading `x`, so none of
these can clash with the input. A proof may introduce names of its own, such as a gate
`g<i>` per node, and is free to mention both the input's variables and these, which is what
moving the preserved set needs.

The **selectors** of a node i other than the root, one per parent p:

- `a<p>`, if p is an `and`;
- `t<p>_<j>`, if p is an `or` and i is its child in position j.

They are listed in order of p. A node with no parents has no selectors.

### The tagged formula

In this order:

1. `preserved: y0_0 y0_1 … y0_<n₀−1> y1_0 … t… ;`, listing every y variable, ordered by
   variable and then by value, then every t variable, ordered by node and then by
   position.
2. For each variable u in increasing order, the exactly-one group, as two rows:
   `1 yu_0 1 yu_1 … 1 yu_<nᵤ−1> >= 1 ;` and then the same terms `<= 1 ;`.
3. `1 a<r> >= 1 ;` for the root r.
4. For each node i in order:
   1. if i is not the root, with selectors s₁ … sₚ: `1 ~ai 1 s₁ … 1 sₚ >= 1 ;`, then
      `1 ai 1 ~sⱼ >= 1 ;` for each j. With no selectors, this is just `1 ~ai >= 1 ;`, and
      the node is never active.
   2. if i is a cube, for each of its pairs (u, t): `1 ~ai 1 yu_t >= 1 ;`.
   3. if i is an `or` with k children, the at-least-one `1 ~ai 1 ti_0 … 1 ti_<k−1> >= 1 ;`.
      Then `1 ~ti_j 1 ai >= 1 ;` for each j. Then, if k ≥ 2, the at-most-one
      `1 ~ti_0 … 1 ~ti_<k−1> >= k−1 ;`.

Numbers are written in decimal, terms are separated by single spaces, and every row ends
`;` and a newline. There is no `*` header line.

VeriPB compares an output file with its database as a set of normalised constraints, so the
order of the rows, and the choice between `>= -1` with negated coefficients and `<= 1`,
does not matter to the check. They are fixed anyway, so that the file is reproducible and
diffable.

The a variables are not preserved. Each one is determined by the t variables, so leaving
them out changes no count, and the proof does not have to account for them.

## Why the count is right

Suppose the statement is well-formed. Write T for the tagged formula, and fix an assignment
to its y variables that satisfies the rows of item 2, which is the same thing as an
assignment f.

1. **Objects to solutions.** Take an object of the root that determines f. Make active
   every node the object passes through: the root; every child of an active `and`; and the
   chosen child of an active `or`. Set each `t<i>_<j>` true exactly when i is active and the
   object takes position j there.

   This satisfies T:
   - each node is active exactly when one of its selectors is;
   - each active `or` has exactly one tag, and an inactive one has none;
   - each active cube is consistent with f, because the object determines f.
2. **No node is passed through twice.** Passing through a node twice would need two
   children of an active `and` to lead to it. Its scope would then lie in both of their
   scopes, which are disjoint, so it would be empty. A node with an empty scope is not an
   `or`, and has no `or` below it, so there is no choice that could differ between the two
   visits. So the activity of one object is a set of nodes, and the tags record every
   choice it makes. Different objects therefore give different tag assignments.
3. **Solutions to objects.** Take a solution of T that extends f. The root is active. By the
   rows of item 4.1, a node is active exactly when one of its selectors is. By item 4.3,
   each active `or` has exactly one tag true, and each tag implies its `or` is active.
   Following the tags from the root therefore gives one object of the root, and the active
   nodes are exactly the ones it passes through. By item 4.2, its cubes agree with f, so it
   determines f. This inverts step 1, and the a variables are determined by the t
   variables.
4. **The count.** So the solutions of T, projected onto y and t, are in one-to-one
   correspondence with the root's objects. The root's scope is every variable, so each
   object fixes the whole of f. By induction on the nodes, the number of objects is the
   count in the table:
   - an `and`'s objects are tuples;
   - an `or`'s objects are tagged pairs;
   - a cube has one object.

VeriPB checks, with `EQUIENUMERABLE`, that T has as many solutions projected onto its
preserved set as the solver's encoding of the instance has projected onto its own. So the
input's count is the root's count. No property of the statement beyond well-formedness is
assumed anywhere in this argument.

## From an input to the tagged formula

The tagged formula shares no variables with the input, so a proof that the two are
equienumerable has to move the preserved set across, as in Section 4 of the CP 2026 paper.
How a solver does this is its own business: nothing above depends on it. The split
searcher's proof does it in three parts.

1. **The renaming.** When the statement's variables correspond to the input's multi-valued
   variables, as they do for the subgraph solver (pattern vertex u, target vertex t, input
   variable `xu_t`), this is a renaming. It costs a constant number of steps per variable,
   and VeriPB proves every one of them by itself.
   1. For each `xu_t`, define `yu_t` by `red 1 ~yu_t 1 xu_t >= 1 : yu_t -> 0 ;` and
      `red 1 yu_t 1 ~xu_t >= 1 : yu_t -> 1 ;`. Move both to the core, and then
      `preserved_add yu_t : 1 xu_t >= 1 ;`.
   2. Derive the exactly-one rows over y. The at-least-one comes by RUP. The at-most-one
      comes by `pol`, adding the links `yu_t → xu_t` to the input's at-most-one over x.
   3. Do the work below, in terms of either x or y: the links let unit propagation move
      freely between them.
   4. Delete the input's remaining exactly-one rows over x. The at-most-one needs a
      one-line `pol` subproof: the y at-most-one plus the links `xu_t → yu_t`, added to the
      negated row. Then `preserved_rm xu_t : 1 yu_t >= 1 ;` for every x variable, and
      delete the links with the witness `xu_t -> yu_t`.
2. **The argument.** Give each node a gate `g<i>`, defined both ways, saying that its set of
   assignments holds:
   - a cube's literals all hold;
   - all of an `and`'s children hold;
   - some child of an `or` holds.

   From the input, derive the root's gate, and delete the input's other constraints. This
   is the solver's search, and it is the same proof whatever the count turns out to be.
3. **From gates to tags.** Top down:
   1. Define the root's `a` as its gate. Define each other `a` as the or of its selectors.
      Define each tag as `t<i>_<j> ↔ a<i> ∧ g<cⱼ>`. Add each to the preserved set
      temporarily.
   2. Derive the tagged formula's rows. They need the lemmas "active implies its gate" and,
      for each `or`, that its children are disjoint.
   3. Delete the rows `t<i>_<j> ∨ ~a<i> ∨ ~g<cⱼ>`. This is where VeriPB checks
      disjointness: each deletion has to show that an active `or` whose child holds must
      have tagged that child.
   4. Take the gates and the `a`s out of the preserved set, and delete their definitions
      with witnesses.

   The split searcher can show its `or`s are disjoint cheaply. Each `or` it writes is told
   apart by one pattern vertex: each child puts it at a different value, or at most one
   child puts it at none of the others' values. So every step can be a hinted RUP, and
   the cost is a small part of checking the main proof.

## An example

An edge, 0–1, and an isolated vertex, 2, into a triangle, counting injective maps: three
choices for vertex 2, then two ways to place the edge on the other two target vertices.

```
(
  (version 2)
  (domain 3 3 3)
  (nodes
    (0 1 cube (2 2))
    (1 1 cube (2 0))
    (2 1 cube (2 1))
    (3 1 cube (1 0))
    (4 1 cube (1 1))
    (5 1 cube (0 2))
    (6 1 and 3 5)
    (7 1 and 4 5)
    …
    (18 1 and 0 16)
    …
    (23 1 and 2 6)
    (24 6 or 18 19 20 21 22 23)
  )
  (root 24)
)
```

Nodes 18 to 23 each place all three vertices, and node 24 takes any one of them. Its tagged
formula has one tag per position of node 24, `t24_0` to `t24_5`. For each solution of the
input, exactly one of them is true.
