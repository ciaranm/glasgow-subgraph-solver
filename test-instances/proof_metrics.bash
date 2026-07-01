#!/bin/bash
# Measure proof-emission size across a fixed, representative matrix of instances and
# option combinations. This is the deterministic baseline that the solve-pipeline
# refactor (see dev_docs/preprocessor-refactor.md) diffs against, phase by phase:
# Phases 1-3 must leave these numbers unchanged (pure refactors), and Phase 4 onward
# should shrink the supplemental-graph count and PBP line count without changing the
# solution counts.
#
# For each configuration it records: the number of supplemental ("shape") graphs built,
# the OPB (model) line count, the PBP (proof) line count, the solution count / status,
# and the node count. It deliberately does NOT run VeriPB -- it is pure measurement, so
# it is fast and runs in every lane, doubling as a proof-emission smoke test over many
# feature combinations (under the sanitizers it exercises proof.cc with no VeriPB needed).
#
# The option combinations mirror ones the proof ctests already verify, so each is known
# to emit a checkable proof; here we only measure what it emits.
#
# Usage:
#   proof_metrics.bash <solver> <instances_dir> <workdir> [out.tsv]
#
# Prints a TSV table to stdout (and to out.tsv if given). Exits 0 iff every solver
# invocation succeeded.

set -u

if [ "$#" -lt 3 ]; then
    echo "usage: $0 <solver> <instances_dir> <workdir> [out.tsv]" 1>&2
    exit 2
fi

solver="$1"
inst="$2"
workdir="$3"
out="${4:-}"

# Proof logging always needs --no-clique-detection. SUB_MIN additionally turns off the
# supplemental graphs and NDS (the leanest proof); the other rows turn individual
# techniques back on so their emission is measured in isolation.
common="--no-clique-detection"
min="${common} --no-supplementals --no-nds"

# label | solver options | pattern | target
configs=(
    "decision|${min} --format csv|trident.csv|longtrident.csv"
    "count|${min} --count-solutions --format csv|trident.csv|longtrident.csv"
    "induced_unsat|${min} --induced --format csv|c3.csv|trident.csv"
    "supplementals|${common} --no-nds --format csv|trident.csv|longtrident.csv"
    "distance3|${common} --no-nds --distance3 --format csv|trident.csv|longtrident.csv"
    "nds|${common} --no-supplementals --format csv|trident.csv|longtrident.csv"
    "cliques|${common} --no-supplementals --no-nds --cliques --format csv|trident.csv|longtrident.csv"
    "locally_injective|${min} --locally-injective --count-solutions --format csv|trident.csv|longtrident.csv"
    "loopy_supplementals|${common} --no-nds --format lad|small|large"
    "staged|${common} --staged --format csv|trident.csv|longtrident.csv"
)

header=$'config\tshape_graphs\topb_lines\tpbp_lines\tsolutions\tnodes'
printf '%s\n' "${header}"
[ -n "${out}" ] && printf '%s\n' "${header}" > "${out}"

rc=0
for row in "${configs[@]}"; do
    IFS='|' read -r label opts pat tgt <<< "${row}"
    proof="${workdir}/pm_${label}"
    log="${proof}.log"

    # shellcheck disable=SC2086
    if ! "${solver}" ${opts} --prove "${proof}" "${inst}/${pat}" "${inst}/${tgt}" > "${log}" 2>&1; then
        echo "${label}: solver failed" 1>&2
        cat "${log}" 1>&2
        rc=1
        continue
    fi

    opb_lines=$(wc -l < "${proof}.opb")
    pbp_lines=$(wc -l < "${proof}.pbp")
    shape=$(sed -n 's/^shape_graphs = //p' "${log}"); shape="${shape:-NA}"
    sols=$(sed -n 's/^solution_count = //p' "${log}")
    [ -z "${sols}" ] && sols=$(sed -n 's/^status = //p' "${log}")
    sols="${sols:-NA}"
    nodes=$(sed -n 's/^nodes = //p' "${log}"); nodes="${nodes:-NA}"

    line=$(printf '%s\t%s\t%s\t%s\t%s\t%s' "${label}" "${shape}" "${opb_lines}" "${pbp_lines}" "${sols}" "${nodes}")
    printf '%s\n' "${line}"
    [ -n "${out}" ] && printf '%s\n' "${line}" >> "${out}"
done

exit "${rc}"
