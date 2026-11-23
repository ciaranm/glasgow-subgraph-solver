#!/usr/bin/env python3
# usage: analyse-campaign.py RESULTS META [metric...]
# RESULTS lines: arm cfg name rep pad instructions cycles nodes status rc   (bench-one.bash)
#           or:  arm cfg name rep pad instructions cycles runtime nodes status rc (kit)
import sys, math, statistics, collections
res, meta = sys.argv[1], sys.argv[2]
fam = {}
for l in open(meta):
    f = l.split(); fam[(f[0], f[1])] = f[2]
runs = collections.defaultdict(list)   # (arm,cfg,name) -> list of dict
bad = []
for l in open(res):
    f = l.split()
    if len(f) == 10: arm, cfg, name, rep, pad, ins, cyc, nodes, st, rc = f; rt = 'NA'
    elif len(f) == 11: arm, cfg, name, rep, pad, ins, cyc, rt, nodes, st, rc = f
    else: continue
    if rc != '0' or st not in ('true', 'false'): bad.append(l.strip()); continue
    num = lambda x: float(x) if x not in ('NA', '<not') else None
    runs[(arm, cfg, name)].append(dict(rep=int(rep), ins=num(ins), cyc=num(cyc), rt=num(rt), nodes=nodes, st=st))
print(f"bad runs: {len(bad)}"); [print("  ", b) for b in bad[:10]]
# Fixed-work check, per compiler. The default value ordering draws from
# std::uniform_int_distribution, whose output differs between standard libraries, so a GCC
# build (libstdc++) and an Apple clang build (libc++) can legitimately search differently.
# Within one compiler every arm and rep must do identical search; across compilers a
# difference only drops that instance from the clang-vs-gcc lines.
comp_of = lambda arm: arm.split('-')[-1]
work = collections.defaultdict(set)
for (arm, cfg, name), rs in runs.items():
    for r in rs: work[(cfg, name, comp_of(arm))].add((r['nodes'], r['st']))
diff = {k: v for k, v in work.items() if len(v) > 1}
names = sorted({(c, n) for c, n, _ in work})
print(f"instances: {len(names)}; with differing nodes/status across arms or reps: {len(diff)}")
for k, v in list(diff.items())[:10]: print("  ", k, v)
ok = {k for k in work if k not in diff}
cross = [(c, n) for c, n in names if (c, n, 'g') in ok and (c, n, 'c') in ok and work[(c, n, 'g')] != work[(c, n, 'c')]]
if cross:
    print(f"searching differently between the gcc and clang builds (expected with different standard libraries; only left out of the clang-vs-gcc lines): {len(cross)}")
    for c, n in cross[:10]: print("  ", (c, n), "gcc", work[(c, n, 'g')], "clang", work[(c, n, 'c')])
order = ['main', 'pr1', 'pr2', 'pr3', 'v0', 'v1', 'v2', 'v3', 'v4']
arms = sorted({k[0] for k in runs}, key=lambda a: (a.split('-')[-1], order.index(a.split('-')[0]) if a.split('-')[0] in order else len(order), a))
metrics = [m for m in sys.argv[3:]] or ['cyc', 'ins']
def med(arm, cfg, name, m):
    xs = [r[m] for r in runs.get((arm, cfg, name), []) if r[m] is not None]
    return statistics.median(xs) if xs else None
def spread(arm, cfg, name, m):
    xs = [r[m] for r in runs.get((arm, cfg, name), []) if r[m] is not None]
    return max(xs) / min(xs) if len(xs) > 1 else None
insts_for = lambda comp: [(c, n) for c, n in names if (c, n, comp) in ok]
same = [(c, n) for c, n in names if (c, n, 'g') in ok and (c, n, 'c') in ok and (c, n) not in cross]
def gm(xs): return math.exp(sum(math.log(x) for x in xs) / len(xs)) if xs else float('nan')
for m in metrics:
    print(f"\n=== metric {m}: geometric mean ratio of per-instance medians vs baseline, [min, max] ===")
    for comp in ('g', 'c'):
        for base in ('main', 'v1'):
            b = f"{base}-{comp}"
            if not any(k[0] == b for k in runs): continue
            print(f"-- baseline {b}")
            for a in arms:
                if not a.endswith('-' + comp) or a == b: continue
                if base == 'v1' and a.split('-')[0] not in ('v0', 'v2', 'v3'): continue
                rs = []; byfam = collections.defaultdict(list)
                for cfg, name in insts_for(comp):
                    x, y = med(a, cfg, name, m), med(b, cfg, name, m)
                    if x and y: rs.append(x / y); byfam[f"{cfg}/{fam.get((cfg, name), '?')}"].append(x / y)
                if not rs: continue
                fams = "  ".join(f"{f}={gm(v):.3f}" for f, v in sorted(byfam.items()))
                print(f"  {a:6s} n={len(rs):3d} gm={gm(rs):.4f} [{min(rs):.3f}, {max(rs):.3f}]   {fams}")
    print(f"-- clang vs gcc, same arm ({m})")
    for base in order:
        rs = [med(f"{base}-c", c, n, m) / med(f"{base}-g", c, n, m) for c, n in same if med(f"{base}-c", c, n, m) and med(f"{base}-g", c, n, m)]
        if rs: print(f"  {base}: clang/gcc n={len(rs)} gm={gm(rs):.4f} [{min(rs):.3f}, {max(rs):.3f}]")
    sp = [spread(a, c, n, m) for a in arms for c, n in insts_for(comp_of(a)) if spread(a, c, n, m)]
    if sp: print(f"-- within-arm rep spread (max/min over reps): median {statistics.median(sp):.4f}, 90th pct {sorted(sp)[int(.9*len(sp))]:.4f}, max {max(sp):.4f}")
