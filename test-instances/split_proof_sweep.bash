#!/bin/bash
# Generate random instances with sparse patterns, count them with --split-components while
# writing a counting statement and its proof, and check that:
#
#   - the count agrees with the plain solver's;
#   - count_statement agrees with both;
#   - VeriPB accepts the proof as showing that the solver's encoding is equienumerable
#     with count_statement's counting formula ("Proofs", under "Splitting inside the
#     search", in dev_docs/decomposition-counting.md), or, when there is no statement, as
#     concluding that there are no solutions.
#
# Patterns have four to seven vertices and are sparse, so that the search often splits
# them, and targets have five to seven. A third of the instances have loops, and a third
# are non-injective. The signature-space bound cycles through 1, 4, 64 and a million, so
# that some splits are refused, some accepted with a few shared images, and some nested.
# Then two fixed instances from test-instances/split, described below.
# The instances come from awk's random numbers, so they may differ between awk
# implementations; every one of them is checked the same way.
#
# Usage:
#   split_proof_sweep.bash <solver> <veripb> <count_statement> <workdir>
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

# awk -v seed=N -v which=pattern|target: one LAD graph
generate='
function graph(n, density, loops,    v, w) {
    for (v = 0; v < n; ++v) {
        for (w = v + 1; w < n; ++w)
            if (rand() < density) e[v, w] = 1;
        if (loops && rand() < 0.2) e[v, v] = 1;
    }
}
BEGIN {
    srand(seed);
    loops = (seed % 3 == 2);
    m = 4 + int(rand() * 4);
    pd = 0.1 + rand() * 0.25;
    n = 5 + int(rand() * 3);
    td = 0.3 + rand() * 0.55;
    if (which == "pattern") graph(m, pd, loops);
    else { graph(m, pd, 0); delete e; graph(n, td, loops); m = n }
    print m;
    for (v = 0; v < m; ++v) {
        line = ""; d = 0;
        for (w = 0; w < m; ++w)
            if (((v < w) && ((v, w) in e)) || ((w < v) && ((w, v) in e)) || ((v == w) && ((v, v) in e))) { line = line " " w; ++d }
        print d line;
    }
}'

fails=0
checked=0
splits=0
nested=0

# check <tag> <pattern> <target> <bound> [extra options...]
check() {
    local tag="$1" pat="$2" tgt="$3" bound="$4"
    shift 4
    local stmt="${workdir}/${tag}.stmt"
    local proof="${workdir}/${tag}"
    rm -f "${stmt}" "${proof}.opb" "${proof}.pbp" "${proof}.counting.opb"
    checked=$((checked + 1))

    local want got out counted verified
    want=$("${solver}" --format lad --count-solutions "$@" "${pat}" "${tgt}" | grep '^solution_count')
    out=$("${solver}" --format lad --count-solutions --no-clique-detection "$@" --split-components "${bound}" --split-statement "${stmt}" \
        --prove "${proof}" "${pat}" "${tgt}" 2>&1)
    got=$(echo "${out}" | grep '^solution_count')
    if [ -z "${want}" ] || [ "${want}" != "${got}" ]; then
        echo "${tag}: counts differ: '${want}' and '${got}'"
        fails=$((fails + 1))
        return
    fi
    echo "${out}" | grep -q '^splits = [1-9]' && splits=$((splits + 1))
    echo "${out}" | grep -q '^nested_splits = [1-9]' && nested=$((nested + 1))

    if [ -s "${stmt}" ]; then
        counted=$("${count_statement}" "${stmt}" --opb "${proof}.counting.opb" | grep '^count' | sed 's/count/solution_count/')
        if [ "${counted}" != "${want}" ]; then
            echo "${tag}: count_statement says '${counted}', not '${want}'"
            fails=$((fails + 1))
            return
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
}

bounds=(1 4 64 1000000)
for seed in $(seq 1 150); do
    tag="sps_${seed}"
    pat="${workdir}/${tag}_pattern.lad"
    tgt="${workdir}/${tag}_target.lad"
    awk -v seed="${seed}" -v which=pattern "${generate}" > "${pat}"
    awk -v seed="${seed}" -v which=target "${generate}" > "${tgt}"
    if [ $((seed % 3)) -eq 1 ]; then
        check "${tag}" "${pat}" "${tgt}" "${bounds[$((seed % 4))]}" --noninjective
    else
        check "${tag}" "${pat}" "${tgt}" "${bounds[$((seed % 4))]}"
    fi
done

# Fixed instances, from a larger random search, where a split's combining step rejects a
# pair of entries that are both disjunctions. Unit propagation cannot see from either
# entry's node alone where it puts the vertices they collide on, so these are the
# instances whose proofs need the signature lemmas.
fixed="$(dirname "$0")/split"
for name in or_or_collision_1 or_or_collision_2; do
    check "${name}" "${fixed}/${name}_pattern.lad" "${fixed}/${name}_target.lad" 1000000
done

echo "${checked} instances checked, ${splits} of them split, ${nested} with nested splits, ${fails} failed"
[ "${fails}" -eq 0 ]
