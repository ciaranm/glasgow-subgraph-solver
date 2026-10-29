**Proposal: signed counting expressions with one equienumerability proof**

This is an alternative to the counting statements in PRs #111--#117. It is a design proposal, not an implemented or formally verified frontend. Its central change is to count explicit tagged objects, rather than interpret addition as a union of assignment sets. Subtraction is arithmetic cancellation, not set difference.

The recommended workflow needs a verified wrapper around the source OPB. It produces one augmented source and one target, checked by one ordinary `EQUIENUMERABLE` proof. There are no semantic side instances. An alternative retaining the exact source OPB is described at the end.

**1. The certificate equation**

Let `F` be the original formula, counted over its declared preserved set. The frontend converts a signed expression `E` into two ordinary finite object descriptions, `Pos(E)` and `Neg(E)`, and computes their cardinalities `p` and `n`.

The one proof checks:

```
F disjoint-sum Neg(E)  EQUIENUMERABLE  Pos(E)
```

Here disjoint sum includes a counted branch tag. The verified frontend establishes that the source has exactly `#F + n` projected solutions and the target exactly `p`. Acceptance therefore establishes `#F = p - n`.

The frontend must construct the augmented source from the actual original formula. It must not accept a solver-supplied replacement source. Nor can the extra solutions be introduced by ordinary equienumerability-preserving steps inside a proof starting from `F`: the wrapper's counting theorem is essential.

If `n = 0`, the wrapper can be omitted. Thus ordinary positive AND/OR counting can use the original source directly.

**2. Object semantics**

The nonnegative constructors describe finite objects with explicit identity:

| Constructor | Objects | Cardinality |
| --- | --- | --- |
| `empty` | No objects | 0 |
| `unit` | The empty tuple | 1 |
| `range k` | Integers `0 <= i < k` | k |
| `table domains rows` | A row index, carrying that row as data | Number of rows |
| `sum A ...` | A branch index and an object from that branch | Sum of cardinalities |
| `product A ...` | An independently chosen object from each child | Product of cardinalities |
| `repeat k A` | An index in `range k` and an object from A | k times the cardinality |

Zero-child sum is empty; zero-child product is unit. `range 0` is empty. `repeat` is convenience syntax for a product with a range.

A table's row index is part of object identity. Repeated rows are permitted and count separately. Similarly, `sum A A` contains two copies, and `product A A` contains two independent choices. Reusing a node means reusing a description, not identifying its choices.

Table data is useful for connecting objects to problem assignments inside the proof. It does not determine object identity. No assertion that table rows are distinct embeddings, that sum branches describe disjoint assignments, or that product factors describe independent CP variables is trusted by the frontend. The equienumerability proof must establish the required relationship with the input's projected solutions.

Products therefore do not require disjoint assignment-variable scopes. Their semantics is an actual Cartesian product, even if two payloads happen to mention the same problem variable. A producer that incorrectly multiplies overlapping subproblems ordinarily gets the wrong target cardinality, which the main proof cannot certify.

**3. Signed expressions and syntax**

The expression language adds `difference A B`, interpreted as integer subtraction. All constructors can take signed children where their arithmetic meaning permits it. Intermediate integers may be negative. There is no claim that B is a subset of A.

Use one s-expression, fixed section order, natural-number identifiers, and references only to earlier nodes. Every node carries a mandatory claimed evaluation, giving the shape `(id count expression)`:

```text
(counting
  (version 1)
  (nodes
    (0 6 (table (4 4)
         ((0 1) (1 0) (1 2) (2 1) (2 3) (3 2))))
    (1 4 (table (4 4 4)
         ((1 0 2) (1 2 0) (2 1 3) (2 3 1))))
    (2 36 (product 0 0))
    (3 16 (repeat 4 1))
    (4 12 (repeat 2 0))
    (5 28 (sum 3 4))
    (6 8 (difference 2 5)))
  (root 6))
```

This computes `6*6 - (4*4 + 2*6) = 8`. Node 0 lists oriented edges of a four-vertex path. Node 1 lists embeddings of a three-vertex path, with its centre in the first column. The repetition indices distinguish four single-endpoint overlap patterns and two double-endpoint overlap patterns. The frontend does not need to know those interpretations.

The grammar is:

```text
statement ::= (counting (version 1) (nodes node*) (root NAT))
node      ::= (NAT count expression)
count     ::= NAT | (neg POS)
expression ::= (empty)
             | (unit)
             | (range NAT)
             | (table (NAT*) ((NAT*)*))
             | (sum NAT*)
             | (product NAT*)
             | (repeat NAT NAT)
             | (difference NAT NAT)
```

Node identifiers are consecutive from zero. Child references precede their parent. In `repeat`, the first number is a multiplicity and the second a child identifier. Table domains are positive; each row has the declared width, and entries lie in their column domains. Empty tables and zero-column tables are permitted. Numeric atoms are arbitrary-precision naturals; `POS` is a strictly positive natural. A negative evaluation is written `(neg magnitude)`, avoiding a second numeric token syntax and excluding negative zero. Version 1 needs only parentheses, whitespace, the fixed atoms above, and decimal numbers.

The frontend must check every claimed evaluation against the constructor's arithmetic rule. Leaf counts follow from the leaf data; internal counts follow from already checked child evaluations. Any mismatch is an error, even at a node that the root does not use. This is a single bottom-up pass and does not require counting models of an OPB formula. For example, `(6 9 (difference 2 5))` is rejected with a diagnostic such as `node 6: claimed 9, expected 8 (= 36 - 28)`.

For nonnegative object expressions, the annotation is the number of objects. For signed expressions, it is their integer evaluation: it may be negative and is not itself the cardinality of a collection. For example, `(7 (neg 2) (difference 1 0))` correctly records `4 - 6`. Annotations are checked assertions, not additional assumptions for the proof checker. The frontend may use them as cached values after checking them, but must still compute the positive and negative collection sizes used by the balancing translation; their difference alone does not determine those two sizes.

The trusted reader operates on an existing s-expression representation and checks this shape. Names, comments, producer metadata, and proof hints can remain outside this core format.

**4. Eliminating subtraction without semantic conditions**

Write `A+` and `A-` for the two nonnegative descriptions associated with A. All unions below are tagged sums, and all products are independent Cartesian products.

```text
leaf K:       K+ = K                  K- = empty

A + B:        (A+B)+ = A+ + B+        (A+B)- = A- + B-

A - B:        (A-B)+ = A+ + B-        (A-B)- = A- + B+

A * B:        (A*B)+ = A+*B+ + A-*B-
              (A*B)- = A+*B- + A-*B+
```

N-ary operations can be folded. A multiplicity is multiplication by `range k`. These are identities over integers; they require neither disjointness of payloads nor containment of sets.

Structural simplifications such as removing empty summands and unit factors are useful. A frontend may also cancel an identical object description appearing on both sides, provided that rewrite is included in its verified implementation. This is optional.

The fundamental induction is:

```
evaluate(E) = cardinality(Pos(E)) - cardinality(Neg(E)).
```

A negative final result can be rejected immediately. No nonnegativity restriction is needed on intermediate signed expressions.

**5. Translating nonnegative objects to OPB**

Each object description has a Boolean representation, a validity predicate, and an explicit preserved set of identity bits. The compiler establishes a bijection between objects and satisfying assignments projected onto those bits. All auxiliary gate variables are outside the preserved set.

`range k` uses `ceil(log2(k))` identity bits for `k > 1` and a binary upper bound. Its sole nontrivial constraint is that the encoded unsigned integer is below k. Handle k=0 by false and k=1 by true over no bits. This keeps multiplicities binary-sized.

Tables use the same representation for their row indices. If payload bits are emitted, they are determined by the selected row; they introduce no choices. Their defining circuits should be total, using a fixed default for invalid index bit patterns. Validity excludes those patterns.

A product uses disjoint identity-bit blocks and conjoins the children's validity predicates. Even repeated references get independent blocks.

A sum uses a branch tag. In branch i, its child representation must be valid, and identity bits reserved for other branches must take a canonical value. A convenient layout uses a shared payload buffer of the maximum child width: branch i uses its child's prefix and forces the remaining buffer bits to zero. Alternatively, use separate blocks and zero every inactive block. The compiler must implement and verify one particular convention.

The branch tag is preserved. Auxiliary variables may be existential; inactive *preserved* bits must not float. This prevents the usual accidental multiplicative factors in a selector encoding.

Compile the resulting Boolean circuit with full definitions into OPB and assert its validity gate. The scalar result is not trusted from the producer: the frontend computes both cardinalities itself.

Node sharing is sharing of a template. Memoising a compiled template on a particular identity-bit slice can share gates. Different slices generally require different instantiations. The design does not claim that an arbitrary expression DAG always has a linear-size OPB expansion. In particular, repeated products can describe numbers whose bit lengths already grow exponentially in the DAG size. Resource bounds must account for identity width, numeric bit lengths, and instantiated circuits.

**6. The verified source wrapper**

Let the original formula be `F(X,A)`, with preserved variables X and auxiliary variables A. Let `N(Z,B)` encode the negative objects, with preserved identity bits Z and auxiliaries B. Rename namespaces apart and introduce one preserved tag s.

The source wrapper represents:

```text
(not s and F(X,A) and Z=0)
or
(s and N(Z,B) and X=0).
```

Its preserved set is `{s} union X union Z`. Every original formula constraint is conditional on the original branch. Original auxiliary variables outside that branch are unpreserved. The analogous conditions apply to the negative branch. Canonicalising only the inactive identity variables suffices for projected counting.

The target is the positive-object encoding. The verified wrapper theorem is:

```text
#projected(wrapper(F,N)) = #projected(F) + cardinality(N).
```

This holds even when F is unsatisfiable, N is empty, the original preserved set is empty, or the input has many auxiliary extensions of the same projected assignment.

The complete driver is conceptually:

```text
verified_frontend original.opb count.sexpr
    -> balanced-input.opb, positive-output.opb, claimed_count

veripb --force-checked-deletion
    balanced-input.opb proof.pbp positive-output.opb

report claimed_count only after EQUIENUMERABLE acceptance
```

The trusted driver binds both generated files and the reported integer to the same expression and original input. The source wrapper can reuse the existing verified OPB parser or consume the formula already parsed by the checker; it does not need a new counting-language parser beyond s-expressions.

**7. How the running example is proved**

For two disjoint edges into the four-vertex path, the source consists of either an actual injective embedding, or one of the 28 tagged overlap objects. The target consists of an independent pair of the six oriented edge embeddings, giving 36 objects.

Map a source embedding to its pair of edge embeddings. Map an overlap object to the pair obtained by the selected gluing. Conversely, classify each target pair by its endpoint equalities: either it is injective, or it belongs to exactly one of the six tagged overlap classes.

This gives the relation the single proof must establish. Missing or duplicated overlap objects change the source cardinality and cannot be hidden by idempotence of Boolean OR. There is no separate certificate saying the overlap classes are contained in the product or disjoint from one another.

The proof generator will still reason about endpoint equalities and uniqueness. Those arguments now justify a bijection between the two counted encodings, rather than certify premises of a separate arithmetic evaluator.

**8. Ryan's recursion and CP AND/OR search**

Jessica Laurette Ryan's *Parameterised algorithms for counting subgraphs, matchings, and monochromatic partitions*, PhD thesis, University of Glasgow, 2023, Section 3.6 and Lemma 3.8, gives a product of component embedding counts minus a sum of recursively counted overlap classes. This maps directly to `product`, `sum`, and `difference`. Hub guesses give an outer sum; lists determine the component tables. See <https://theses.gla.ac.uk/83568/3/2023RyanPhD.pdf>.

Each quotient/intersection state gets one signed expression node. The pair transformation handles nested subtraction without requiring the signed subexpressions to denote literal subsets of a common assignment space. Isomorphic component tables can be reused as object templates. The proof generator supplies the vertex correspondences when connecting their payloads to embeddings.

This covers the embedding-count recurrence. The thesis's final conversion from embeddings to unlabelled copies by division by the automorphism count is a separate semantic step; equienumerability with an input encoding of labelled embeddings cannot itself certify a smaller unlabelled count.

PR #117 itself combines compatible signature entries using sums of products; the subtraction-based algorithm is in the earlier decomposition PRs. Its positive expressions have no negative objects, so its source can remain unchanged. Decision branches become tagged sums, splits become products and sums over compatible pairs, and cached results become shared templates.

The same construction applies to CP AND/OR search. An OR node records a tagged choice of child; an AND node records a tuple of child objects. A successful leaf contributes unit, potentially carrying its assignment as a singleton table; a failed leaf contributes empty. The main proof establishes that these objects are in bijection with the intended *projected* CP solutions. Arbitrary filtering is deliberately absent from the countable object constructors.

**9. Evidence and remaining proof engineering**

A small proof was checked with VeriPB 3.0.2 and forced checked deletion. The original problem is `x != y` over Boolean variables and has two solutions. The negative objects are the two diagonal pairs; the positive objects are all four pairs. The augmented source has four solutions and verifies as equienumerable with an unconstrained two-bit target.

The bijection uses source tag s and negative identity bit z. With x=y=0 on the negative branch and z=0 on the original branch, the target bits are `u = x or z` and `v = y or z`. The proof introduces these variables, transfers the preserved set, and deletes the source constraints. It uses no solution-exclusion rules and no side proof.

This verifies that the workflow works for a nontrivial subtraction example using the existing checker. It does not establish a polynomial proof-size bound for Ryan's general recursion, nor implement the proposed frontend. Large recursive cancellations may need substantial proof engineering. The positive and negative collections can also be much larger than their difference. Those are explicit costs of this proposal.

The relevant formalisation milestones are the object-encoding bijection, signed-pair arithmetic, the source-wrapper counting theorem, and their composition with the existing equienumerability theorem. A useful next prototype is the 36=8+28 edge example, followed by a recursive quotient example and a cached AND/OR example.

**10. If the original source must remain unchanged**

There is a different construction that retains assignment-set circuits and folds their semantic obligations into the target. It needs a trusted upper bound on the input's projected count. The number p of variables in the actual input preserved set supplies `#F <= 2^p`.

Let `R(Y)` be the ordinary root predicate, and let `Bad(W)` express that some counting rule has a counterexample under an independent assignment W: overlapping union children, or a subtracted child outside its parent. Include all required domain constraints and total definitions in the witness circuit, but neither the root assertion nor the original input constraints.

Introduce a preserved mode bit b and p+1 preserved padding bits T. Compile:

```text
(not b and R(Y) and T=0)
or
(b and Bad(W) and Y=canonical).
```

Preserve b, the identity bits Y, and T; do not preserve W or witness selectors. The normal branch has exactly `#R` projected models. If any counterexample exists, the other branch has exactly `2^(p+1)` projected models. Therefore its count is

```text
#R + 2^(p+1) * indicator(exists W. Bad(W)).
```

A proof of equienumerability with the unchanged input forces the indicator to be zero, since the input has at most `2^p` projected models. The usual counting induction then establishes that the evaluated count is `#R = #F`. The huge number is represented by p+1 free bits, not that many rows.

The trusted driver must obtain p from the actual source preserved set, not an untrusted statement field. Canonical values must exist, and inactive witness formulas must not constrain the normal branch. A careful compiler can enforce these requirements.

This alternative gives one target file and one proof with an unchanged source. It integrates the former side conditions into the target's semantics; it does not eliminate the need to prove them. Its target has a cheaply computed count only after acceptance excludes the counterexample branch. The signed-object proposal above has genuinely countable endpoints before verification and is the recommended design if a verified source wrapper is acceptable.

A third, simpler unchanged-source option is to evaluate the arithmetic expression to N and output a binary interval `0 <= r < N`. That has exactly N models unconditionally. It is sound, but erases the component-object structure and leaves the proof generator to establish an arbitrary ranking into that interval. The balanced construction preserves more of the combinatorial argument.
