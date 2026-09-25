#!/usr/bin/env python3
# Generate harder instances from the graph3 target, in the same CSV format as its patterns:
#   gt<k>_<rep>.csv : k people chosen at random, keeping each ground-truth edge among them (and
#                     their attribute edges) with probability 0.4-0.7, as the supplied patterns do
#   q<k>_<rep>.csv  : random queries over k people, with relations not taken from the ground
#                     truth: for each pair and each of spatial / distance / interaction, a random
#                     label from that category with probability 1/3
#
#   gen.py <graph3 directory> <output directory> [seed]

import collections
import os
import random
import sys

data, out = sys.argv[1], sys.argv[2]
random.seed(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
os.makedirs(out, exist_ok=True)

V = {}
GT = []
for line in open(os.path.join(data, 'graph3.gt.csv')):
    p = line.strip().split(',')
    if len(p) == 4 and p[1] == '':
        V[p[0]] = p[2]
    else:
        GT.append(p)


def ends(p):
    if len(p) == 3:
        a, b = p[0].split('>')
        return a, b
    return p[0], p[1]


persons = [v for v in V if V[v] == 'person']

# labels by category, and whether each category is directed, taken from the full target
labs = collections.defaultdict(set)
dirs = {}
for line in open(os.path.join(data, 'graph3.csv')):
    p = line.strip().split(',')
    if len(p) == 4 and p[1] == '':
        continue
    lab = p[1] if len(p) == 3 else p[2]
    cat = lab.split(':')[0]
    labs[cat].add(lab)
    dirs[cat] = (len(p) == 3)


def write(fn, vs, es):
    with open(os.path.join(out, fn), 'w') as f:
        for v in sorted(vs):
            f.write('%s,,%s,1.0\n' % (v, V[v]))
        for e in es:
            f.write(','.join(e) + ',1.0\n')


for k in [12, 14, 16, 18]:
    for rep in range(5):
        vs = set(random.sample(persons, k))
        r = random.uniform(0.4, 0.7)
        es = []
        for p in GT:
            a, b = ends(p)
            if a in vs and (b in vs or not V[b] == 'person' and b.startswith('attribute')) and random.random() < r:
                if b not in vs:
                    vs.add(b)
                es.append(p[:-1])
        write('gt%d_%d.csv' % (k, rep), vs, es)

for k in [4, 6, 8, 10, 12, 14, 16]:
    for rep in range(5):
        vs = ['q%d' % i for i in range(k)]
        for v in vs:
            V[v] = 'person'
        es = []
        d = 0.5
        for i in range(k):
            for j in range(i + 1, k):
                for cat in ['spatial', 'distance', 'interaction']:
                    if random.random() < d / 3 * 2:
                        lab = random.choice(sorted(labs[cat]))
                        a, b = (vs[i], vs[j]) if random.random() < 0.5 else (vs[j], vs[i])
                        es.append([a + '>' + b, lab] if dirs[cat] else [a, b, lab])
        write('q%d_%d.csv' % (k, rep), vs, es)
