#!/bin/bash
# Generate random disconnected patterns, count them with --decomposition while writing a
# counting statement and its proof, and check that:
#
#   - the count agrees with the plain solver's;
#   - count_statement agrees with both;
#   - VeriPB accepts the proof as showing that the solver's encoding is equienumerable
#     with count_statement's counting formula (phases 4 and 5 of
#     dev_docs/decomposition-counting.md), or, when there is no statement, as concluding
#     that there are no solutions.
#
# Patterns have three or four components of one to three vertices, at most six in all,
# and targets at least as many vertices. The second half of the instances have loops.
# The instances come from awk's random numbers, so they may differ between awk
# implementations; every one of them is checked the same way.
#
# Usage:
#   decomposition_proof_sweep.bash <solver> <veripb> <count_statement> <workdir>
#
# Exits 0 iff every instance passes.

set -u

if [ "$#" -ne 4 ]; then
    echo "usage: $0 <solver> <veripb> <count_statement> <workdir>" 1>&2
    exit 2
fi

solver="$1"
veripb="$2"
count_statement="$3"
workdir="$4"
mkdir -p "${workdir}" || exit 1

# awk -v seed=N -v loops=0|1 -v which=pattern|target: one LAD graph
generate='
BEGIN {
    srand(seed);
    ncomp = 3 + int(rand() * 2);
    do {
        m = 0;
        for (c = 0; c < ncomp; ++c) { size[c] = 1 + int(rand() * 3); first[c] = m; m += size[c]; }
    } while (m > 6);
    for (c = 0; c < ncomp; ++c) {
        for (v = first[c] + 1; v < first[c] + size[c]; ++v) { w = first[c] + int(rand() * (v - first[c])); pe[w, v] = 1; }
        for (v = first[c]; v < first[c] + size[c]; ++v)
            for (w = v + 1; w < first[c] + size[c]; ++w)
                if (rand() < 0.3) pe[v, w] = 1;
        if (loops && rand() < 0.3) pe[first[c], first[c]] = 1;
    }
    n = m + int(rand() * 3);
    if (n < 4) n = 4;
    for (v = 0; v < n; ++v)
        for (w = v; w < n; ++w)
            if ((v == w) ? (loops && rand() < 0.3) : (rand() < 0.5)) te[v, w] = 1;
    k = (which == "pattern") ? m : n;
    print k;
    for (v = 0; v < k; ++v) {
        line = ""; d = 0;
        for (w = v; w < k; ++w)
            if ((which == "pattern" && ((v, w) in pe)) || (which == "target" && ((v, w) in te))) { line = line " " w; ++d }
        print d line;
    }
}'

fails=0
checked=0
for loops in 0 1; do
    for seed in $(seq 1 15); do
        tag="dps_${loops}_${seed}"
        pat="${workdir}/${tag}_pattern.lad"
        tgt="${workdir}/${tag}_target.lad"
        awk -v seed="${seed}" -v loops="${loops}" -v which=pattern "${generate}" > "${pat}"
        awk -v seed="${seed}" -v loops="${loops}" -v which=target "${generate}" > "${tgt}"
        stmt="${workdir}/${tag}.stmt"
        proof="${workdir}/${tag}"
        rm -f "${stmt}" "${proof}.opb" "${proof}.pbp" "${proof}.counting.opb"
        checked=$((checked + 1))

        want=$("${solver}" --format lad --count-solutions "${pat}" "${tgt}" | grep '^solution_count')
        got=$("${solver}" --format lad --count-solutions --decomposition --decomposition-statement "${stmt}" \
            --prove "${proof}" "${pat}" "${tgt}" 2>&1 | grep '^solution_count')
        if [ -z "${want}" ] || [ "${want}" != "${got}" ]; then
            echo "${tag}: counts differ: '${want}' and '${got}'"
            fails=$((fails + 1))
            continue
        fi

        if [ -s "${stmt}" ]; then
            counted=$("${count_statement}" "${stmt}" --opb "${proof}.counting.opb" | grep '^count' | sed 's/count/solution_count/')
            if [ "${counted}" != "${want}" ]; then
                echo "${tag}: count_statement says '${counted}', not '${want}'"
                fails=$((fails + 1))
                continue
            fi
            verified=$("${veripb}" --force-checked-deletion "${proof}.opb" "${proof}.pbp" "${proof}.counting.opb" 2>&1 | grep '^s ')
            if ! echo "${verified}" | grep -q '^s VERIFIED OUTPUT EQUIENUMERABLE'; then
                echo "${tag}: the proof does not verify as equienumerable with the counting formula"
                fails=$((fails + 1))
            fi
        else
            verified=$("${veripb}" --force-checked-deletion "${proof}.opb" "${proof}.pbp" 2>&1 | grep '^s ')
            if ! echo "${verified}" | grep -q '^s VERIFIED COMPLETE ENUMERATION OF 0 SOLUTIONS'; then
                echo "${tag}: no statement, but the proof does not conclude that there are no solutions"
                fails=$((fails + 1))
            fi
        fi
    done
done

echo "${checked} instances checked, ${fails} failed"
[ "${fails}" -eq 0 ]
