#!/usr/bin/env python3
# usage: make-comment.py RESULTS MACHINE META > comment.md
# Builds the markdown comment to post on the tracking issue: machine summary, the
# fixed-work check and geometric means (from analyse-campaign.py), and per-instance
# median ratios, all small enough for one GitHub comment.
import sys, os, re, subprocess, statistics, collections

res, machine, meta = sys.argv[1:4]
here = os.path.dirname(os.path.abspath(__file__))
mach = open(machine).read()

def field(pat):
    m = re.search(pat, mach, re.M)
    return m.group(1).strip() if m else None

cpu = field(r'^Model name:\s*(.*)$') or field(r'^(Apple .*)$') or '?'
flags = (field(r'^Flags:\s*(.*)$') or '').split()
interesting = [f for f in ('popcnt', 'avx2', 'avx512f', 'avx512bw', 'avx512vl', 'avx512_vpopcntdq', 'bmi2') if f in flags]

runs = collections.defaultdict(list)
for l in open(res):
    f = l.split()
    if len(f) != 11: continue
    arm, cfg, name, rep, pad, ins, cyc, rt, nodes, st, rc = f
    runs[(arm, cfg, name)].append(dict(ins=ins, cyc=cyc, rt=rt))
have = lambda m: any(r[m] not in ('NA', '<not') for rs in runs.values() for r in rs)
metrics = [m for m in ('cyc', 'ins', 'rt') if have(m)]
primary = metrics[0] if metrics else 'rt'

analysis = subprocess.run([sys.executable, os.path.join(here, 'analyse-campaign.py'), res, meta] + metrics,
                          capture_output=True, text=True).stdout

def med(arm, cfg, name, m):
    xs = [float(r[m]) for r in runs.get((arm, cfg, name), []) if r[m] not in ('NA', '<not')]
    return statistics.median(xs) if xs else None

arms = sorted({k[0] for k in runs})
order = ['main', 'pr1', 'pr2', 'pr3', 'v1', 'v2', 'v3']
insts = sorted({(k[1], k[2]) for k in runs})
reps = max(len(v) for v in runs.values()) if runs else 0

def table(m):
    out = []
    for comp, cname in (('g', 'gcc'), ('c', 'clang')):
        cols = [a for a in order if f"{a}-{comp}" in arms and a != 'main']
        if f"main-{comp}" not in arms: continue
        out.append(f"\n{cname}: `main` median {m}, then each arm's median divided by `main`'s\n")
        out.append("| cfg | instance | main | " + " | ".join(cols) + " |")
        out.append("|---|---|---:|" + "---:|" * len(cols))
        for cfg, name in insts:
            b = med(f"main-{comp}", cfg, name, m)
            if not b: continue
            cells = []
            for a in cols:
                x = med(f"{a}-{comp}", cfg, name, m)
                cells.append(f"{x / b:.3f}" if x else "–")
            out.append(f"| {cfg} | {name} | {b:.4g} | " + " | ".join(cells) + " |")
    return "\n".join(out)

print(f"## Results: {cpu}\n")
print(f"- **CPU:** {cpu}" + (f" (flags of interest: {', '.join(interesting)})" if interesting else ""))
for label, pat in (('OS', r'^((?:Linux|Darwin) \S+ \S+)$'), ('macOS', r'^ProductVersion:\s*(.*)$'),
                   ('boost', r'^boost (.*)$'), ('intel_pstate no_turbo', r'^intel_pstate no_turbo (.*)$'),
                   ('governor', r'^governor (.*)$'), ('pinning', r'^pinned core (.*)$'),
                   ('perf counters', r'^perf counters: (.*)$'), ('power', r"^Now drawing from '(.*)'$"),
                   ('gcc arm', r'^gcc arm compiler (.*)$'), ('clang arm', r'^clang arm compiler (.*)$'),
                   ('kit commit', r'^kit: (.*)$')):
    v = field(pat)
    if v: print(f"- **{label}:** {v}")
print(f"- **runs:** {sum(len(v) for v in runs.values())} ({len(insts)} instances, {len(arms)} arms, up to {reps} reps); metrics {', '.join(metrics)}")
print("\n### Analysis (`analyse-campaign.py`)\n\n```\n" + analysis.rstrip() + "\n```")
print(f"\n<details><summary>Per-instance medians ({primary})</summary>\n" + table(primary) + "\n\n</details>")
if primary != 'ins' and 'ins' in metrics:
    print("\n<details><summary>Per-instance medians (ins)</summary>\n" + table('ins') + "\n\n</details>")
print("\n<details><summary>machine file</summary>\n\n```\n" + mach.rstrip() + "\n```\n\n</details>")
