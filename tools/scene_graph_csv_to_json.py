#!/usr/bin/env python3
# Convert a scene graph in the CSV dialect of the graph3 benchmark into the gss-graph JSON
# format, for glasgow_subgraph_solver --minimise-cost. This is a stopgap for one dataset,
# kept here for now rather than being a supported reader.
#
#   scene_graph_csv_to_json.py --role target|pattern [--scale 1000000] in.csv > out.json
#
# CSV lines: "name,,label,weight" is a vertex, "a>b,label,weight" a directed edge, and
# "a,b,label,weight" an undirected edge. A pair may have several edges, one per label.
#
# Both graphs become directed multigraphs, since one graph mixes directed and
# undirected relations. In the target, an undirected edge becomes two arcs with the same
# cost; in the pattern, it becomes one arc, in the direction written, so that its cost is
# counted once. A target element's cost is round(-log(weight) * scale). The pattern
# carries no costs: its weights are not used.
#
# The "_reified" variants of those files, whose extra fifth column marks vertices that
# stand for edges, are refused: the solver does its own reification.

import argparse
import json
import math
import sys

parser = argparse.ArgumentParser()
parser.add_argument('--role', choices=['target', 'pattern'], required=True)
parser.add_argument('--scale', type=int, default=1000000)
parser.add_argument('csv')
args = parser.parse_args()


def cost(w):
    w = float(w)
    if not (0 < w <= 1):
        sys.exit('weight %r is not a probability' % w)
    return round(-math.log(w) * args.scale)


names, labels, vcosts, edges = [], {}, {}, {}


def vertex(n):
    if n not in labels:
        names.append(n)
        labels[n] = None
    return n


for number, line in enumerate(open(args.csv), 1):
    line = line.rstrip('\r\n')
    if not line:
        continue
    p = line.split(',')
    if len(p) >= 5:
        sys.exit('%s:%d: this looks like a reified file (a fifth column), which this does not read' % (args.csv, number))
    if len(p) == 4 and p[1] == '':
        vertex(p[0])
        labels[p[0]] = p[2]
        vcosts[p[0]] = cost(p[3])
        continue
    if len(p) == 3 and '>' in p[0]:
        a, b = p[0].split('>')
        arcs = [(a, b)]
        label, w = p[1], p[2]
    elif len(p) == 4:
        a, b = p[0], p[1]
        arcs = [(a, b), (b, a)] if args.role == 'target' else [(a, b)]
        label, w = p[2], p[3]
    else:
        sys.exit('%s:%d: cannot parse %r' % (args.csv, number, line))
    vertex(a)
    vertex(b)
    for f, t in arcs:
        key = (f, t, label)
        if key in edges:
            sys.exit('%s:%d: repeats the edge %r' % (args.csv, number, key))
        edges[key] = cost(w)

missing = [n for n in names if labels[n] is None]
if missing:
    sys.exit('vertices with no label line: %r' % missing[:5])

doc = {
    'format': 'gss-graph', 'version': 1, 'directed': True, 'multigraph': True,
    'vertices': [dict({'name': n, 'label': labels[n]}, **({'cost': vcosts[n]} if args.role == 'target' else {})) for n in names],
    'edges': [dict({'from': f, 'to': t, 'label': l}, **({'cost': c} if args.role == 'target' else {})) for (f, t, l), c in edges.items()],
}
json.dump(doc, sys.stdout)
