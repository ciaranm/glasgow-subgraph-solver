# Weighted matching prototype

A standalone experiment, not part of the solver and not built by CMake. It answers one
question: for the scene-graph matching problem our collaborators are solving with a fork of
the subgraph solver, which lower bound is worth building? Nothing here is meant to be merged
as it is.

## The problem

A target scene graph gives, for every pair of entities, a softmax over the labels in each
relation category (`spatial:`, `distance:`, `interaction:`, `hoi:`), plus sparse `has_*` edges
to attribute vertices. A pattern is a query over entities and relations. Find an injective
mapping that preserves vertex labels, where each pattern edge `(u, v, label)` maps to a target
edge with the same label and directedness, minimising the sum of `-log(weight)` over the
matched edges and vertices. There can be several edges between one pair, with different
labels, and that is why `InputGraph` cannot hold these graphs.

Because the target is complete in every category except the attributes, feasibility prunes
almost nothing, and the problem is a quadratic assignment problem in all but name.

## Files

- `wsip.cc`: the CSV reader, a branch and bound with four bounds (`none`, `gl`, `lap`,
  `dual`), and two PB encoders (`opb-weak`, `opb-strong`). The header comment describes the
  modes and the output columns.
- `gen.py`: generates the harder instances: larger ground-truth patterns (`gt<k>_<rep>`) and
  random queries (`q<k>_<rep>`). Seed 1 gives the instances behind `results/generated*.tsv`.
- `renum.py`: renames variables so that RoundingSat will read the OPB.
- `rs_one.sh`: encodes one instance and runs RoundingSat on it with a 60 second limit.
- `results/`: the raw numbers behind the write-up.

The graph3 data is not in the repository. It belongs to the collaborators.

## Running

```shell
c++ -std=c++20 -O2 -o wsip wsip.cc
./wsip dual ~/graph3/patterns/pat_1.csv ~/graph3/graph3.csv 60

./gen.py ~/graph3 generated 1
./wsip dual generated/q10_1.csv ~/graph3/graph3.csv 120

GRAPH3=~/graph3 ROUNDINGSAT=/path/to/roundingsat ./rs_one.sh ~/graph3/patterns/pat_1.csv pat_1
```

RoundingSat is <https://gitlab.com/MIAOresearch/software/roundingsat>, built in release mode with
its default options.

## Results files

`graph3_patterns.tsv` and `generated.tsv` have one row per (bound, instance): bound, instance,
then the columns `wsip` prints (status, objective, cost of the name-preserving mapping, root
bound, nodes, seconds, name matches). The time limits were 60 s for the supplied patterns and
120 s for the generated ones. The two `*_roundingsat.tsv` files have columns tag, status,
objective scaled by 10^4, and seconds, with a 60 s limit. Runs were six to eight at a time on
one laptop, so timings are indicative, not benchmark quality.

## Findings in brief

- On the 100 supplied patterns, any bound at all is enough. With no bound, 6 time out; with
  `gl`, all 100 solve in 0.2 s in total. `dual` proves optimality at the root on 94.
- The matching (assignment) bound adds almost nothing over `gl`. What matters is how each pair
  cost is split between its endpoints, and the local-polytope dual (`dual`) optimises that.
- RoundingSat is very sensitive to the encoding: 13 s with `opb-weak` against 0.09 s with
  `opb-strong` on `pat_1`, and on `pat_96` over 9 minutes with `opb-weak`, still unproved when
  stopped, against 13 s with `opb-strong`.
- Random queries are the hard case. `dual` solves 10-person queries in about 2 s and
  12-person ones in about 13 s, and RoundingSat solves none above 6 people in 60 s. The root
  gap is 20-35% on these, and more dual iterations do not close it, so it is a limit of the
  local polytope.

## Soundness caveats

Costs are doubles and bound pruning compares them directly. Floating-point error could
therefore prune an optimal solution. The PB encodings round costs to integers, so their optima
are for a slightly different problem. On all 107 instances where both finished, the two
optima agree to within 0.0003, which is about the precision `wsip` prints. A real implementation should use fixed-point integer costs
throughout.
