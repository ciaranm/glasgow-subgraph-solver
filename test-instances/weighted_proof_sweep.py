#!/usr/bin/env python3
# Generate random weighted multigraph instances and check that the solver's proof
# verifies on each, and that a verified optimum is the cost the solver reported. The
# axes are the ones the proof's model has to get right for itself, independently of
# reification: directed or undirected graphs (and the two mixed), loops, vertex labels,
# edge labels, parallel edges, and vertex and edge costs, negative ones included. Most
# instances minimise cost; the rest are multigraph decision problems.
#
# The in-process weighted_homomorphism_test checks the costs against a brute-force
# oracle; this checks the proofs.
#
#   weighted_proof_sweep.py <solver> <veripb> <workdir> [seed] [instances]
#
# Exits 0 iff every proof verifies. Instances the solver refuses (a pattern without edge
# labels on a multigraph target, or a target whose costs are on no edges) are skipped.

import json
import os
import random
import re
import subprocess
import sys

solver, veripb, workdir = sys.argv[1:4]
seed = int(sys.argv[4]) if len(sys.argv) > 4 else 1
instances = int(sys.argv[5]) if len(sys.argv) > 5 else 60
rng = random.Random(seed)


def graph(n, directed, loops, vertex_labels, edge_labels, multigraph, vertex_costs, edge_costs, density):
    names = ['v%d' % i for i in range(n)]
    vertices = []
    for v in names:
        d = {'name': v}
        if vertex_labels:
            d['label'] = rng.choice('AB')
        if vertex_costs:
            d['cost'] = rng.randint(-3, 9)
        vertices.append(d)

    def labels():
        if not edge_labels:
            return [None]
        if not multigraph:
            return [rng.choice('xy')]
        chosen = [l for l in 'xy' if rng.random() < 0.6]
        return chosen or [rng.choice('xy')]

    edges = []
    for i in range(n):
        for j in range(n):
            if (i == j and not loops) or ((not directed) and j < i):
                continue
            if rng.random() < (0.4 if i == j else density):
                for l in labels():
                    e = {'from': names[i], 'to': names[j]}
                    if l is not None:
                        e['label'] = l
                    if edge_costs:
                        e['cost'] = rng.randint(-3, 9)
                    edges.append(e)

    doc = {'format': 'gss-graph', 'version': 1, 'directed': directed, 'vertices': vertices, 'edges': edges}
    if multigraph:
        doc['multigraph'] = True
    return doc


failures = checked = skipped = 0
pattern_file = os.path.join(workdir, 'wps_pattern.json')
target_file = os.path.join(workdir, 'wps_target.json')
proof = os.path.join(workdir, 'wps_proof')

for k in range(instances):
    loops, vertex_labels, edge_labels = rng.random() < 0.4, rng.random() < 0.5, rng.random() < 0.7
    multigraph = edge_labels and rng.random() < 0.6
    vertex_costs, edge_costs = rng.random() < 0.5, rng.random() < 0.7
    if not (vertex_costs or edge_costs):
        vertex_costs = True
    minimise = rng.random() < 0.8
    if not (minimise or multigraph):
        multigraph = edge_labels = True

    pattern = graph(rng.randint(1, 4), rng.random() < 0.5, loops, vertex_labels, edge_labels, multigraph, False, False, rng.random())
    target = graph(rng.randint(2, 6), rng.random() < 0.5, loops, vertex_labels, edge_labels, multigraph, vertex_costs, edge_costs, 0.3 + 0.6 * rng.random())
    with open(pattern_file, 'w') as f:
        json.dump(pattern, f)
    with open(target_file, 'w') as f:
        json.dump(target, f)

    args = [solver, '--no-clique-detection', '--prove', proof, '--format', 'json', pattern_file, target_file]
    if minimise:
        args.insert(1, '--minimise-cost')
    run = subprocess.run(args, capture_output=True, text=True)
    if run.returncode != 0:
        if 'cannot yet be used' in run.stderr or 'needs a target with' in run.stderr:
            skipped += 1
            continue
        print('instance %d: solver failed:\n%s%s' % (k, run.stdout[-500:], run.stderr[-500:]))
        failures += 1
        continue

    cost = re.search(r'^cost = (-?\d+)$', run.stdout, re.M)
    check = subprocess.run([veripb, proof + '.opb', proof + '.pbp'], capture_output=True, text=True)
    status = [l for l in check.stdout.splitlines() if l.startswith('s ')]
    status = status[0] if status else ''

    ok = 'VERIFIED' in status
    if ok and minimise:
        bounds = re.search(r'BOUNDS (\S+) <= obj <= (\S+)', status)
        expected = cost.group(1) if cost else 'INF'
        ok = bool(bounds) and bounds.group(1) == bounds.group(2) == expected

    checked += 1
    if not ok:
        failures += 1
        print('instance %d (%s): %s\n%s\n  pattern: %s\n  target: %s' % (
            k, 'minimise' if minimise else 'decide', status or 'no verdict', check.stdout[-400:] + check.stderr[-400:],
            json.dumps(pattern), json.dumps(target)))

print('weighted proof sweep: %d checked, %d skipped, %d failed' % (checked, skipped, failures))
sys.exit(1 if failures else 0)
