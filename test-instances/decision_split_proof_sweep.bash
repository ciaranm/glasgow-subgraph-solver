#!/bin/bash
# Generate random instances whose patterns fall apart into components during search, solve
# them as decision problems with --split-decision while writing a proof, and check that:
#
#   - the answer agrees with the plain solver's;
#   - VeriPB accepts the proof as showing that answer.
#
# Half the patterns are small dense blocks joined through one or two hub vertices, whose
# components split apart and join again as the hubs are assigned and unassigned; the rest
# are sparse. A third have loops. The options cycle through injective and non-injective,
# with default restarts and with none, and staged. Then the fixed instances in
# test-instances/split-decision, described below.
# The instances come from awk's random numbers, so they may differ between awk
# implementations; every one of them is checked the same way.
#
# Usage:
#   decision_split_proof_sweep.bash <solver> <veripb> <workdir>
#
# Exits 0 iff every instance passes.

set -u

if [ "$#" -ne 3 ]; then
    echo "usage: $0 <solver> <veripb> <workdir>" 1>&2
    exit 2
fi

solver="$1"
veripb="$2"
workdir="$3"
mkdir -p "${workdir}" || exit 1

# awk -v seed=N -v which=pattern|target: one LAD graph
generate='
function edge(v, w) { if (v < w) e[v, w] = 1; else e[w, v] = 1 }
BEGIN {
    srand(seed);
    loops = (seed % 3 == 2);
    if (seed % 2) {
        # dense blocks of two to four vertices, joined through one or two hubs
        nh = 1 + int(rand() * 2); nb = 2 + int(rand() * 4); m = nh;
        for (b = 0; b < nb; ++b) {
            first[b] = m; size[b] = 2 + int(rand() * 3); m += size[b];
            for (v = first[b]; v < m; ++v)
                for (w = v + 1; w < m; ++w)
                    if (rand() < 0.7) edge(v, w);
        }
        for (h = 0; h < nh; ++h)
            for (b = 0; b < nb; ++b)
                for (k = 1 + int(rand() * 2); k > 0; --k)
                    edge(h, first[b] + int(rand() * size[b]));
        if (nh == 2 && rand() < 0.5) edge(0, 1);
    }
    else {
        m = 6 + int(rand() * 9);
        pd = 0.08 + rand() * 0.17;
        for (v = 0; v < m; ++v)
            for (w = v + 1; w < m; ++w)
                if (rand() < pd) edge(v, w);
    }
    if (loops)
        for (v = 0; v < m; ++v)
            if (rand() < 0.2) e[v, v] = 1;
    if (which == "target") {
        delete e;
        n = m + int(rand() * 10);
        td = 0.15 + rand() * 0.45;
        for (v = 0; v < n; ++v) {
            for (w = v + 1; w < n; ++w)
                if (rand() < td) e[v, w] = 1;
            if (loops && rand() < 0.3) e[v, v] = 1;
        }
        m = n;
    }
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
split=0
refuted=0
combined=0

# check <tag> <pattern> <target> [extra options...]
check() {
    local tag="$1" pat="$2" tgt="$3"
    shift 3
    local proof="${workdir}/${tag}"
    rm -f "${proof}.opb" "${proof}.pbp"
    checked=$((checked + 1))

    local want got out verified
    want=$("${solver}" --format lad --no-clique-detection "$@" "${pat}" "${tgt}" | grep '^status')
    out=$("${solver}" --format lad --no-clique-detection --split-decision "$@" --prove "${proof}" "${pat}" "${tgt}" 2>&1)
    got=$(echo "${out}" | grep '^status')
    if [ -z "${want}" ] || [ "${want}" != "${got}" ]; then
        echo "${tag} ($*): answers differ: '${want}' and '${got}'"
        fails=$((fails + 1))
        return
    fi
    echo "${out}" | grep -q '^decision_split_searches = [1-9]' && split=$((split + 1))
    echo "${out}" | grep -q '^decision_split_refutations = [1-9]' && refuted=$((refuted + 1))
    echo "${out}" | grep -q '^decision_split_solutions = [1-9]' && combined=$((combined + 1))

    verified=$("${veripb}" --force-checked-deletion "${proof}.opb" "${proof}.pbp" 2>&1 | grep '^s ')
    if [ "${want}" = "status = true" ]; then
        echo "${verified}" | grep -q '^s VERIFIED SATISFIABLE' && return
    else
        echo "${verified}" | grep -q '^s VERIFIED UNSATISFIABLE' && return
    fi
    echo "${tag} ($*): the proof does not verify ('${verified}')"
    fails=$((fails + 1))
}

for seed in $(seq 1 360); do
    tag="dsps_${seed}"
    pat="${workdir}/${tag}_pattern.lad"
    tgt="${workdir}/${tag}_target.lad"
    awk -v seed="${seed}" -v which=pattern "${generate}" > "${pat}"
    awk -v seed="${seed}" -v which=target "${generate}" > "${tgt}"
    case $((seed % 6)) in
        0) check "${tag}" "${pat}" "${tgt}" ;;
        1) check "${tag}" "${pat}" "${tgt}" --restarts none ;;
        2) check "${tag}" "${pat}" "${tgt}" --noninjective ;;
        3) check "${tag}" "${pat}" "${tgt}" --noninjective --restarts none ;;
        4) check "${tag}" "${pat}" "${tgt}" --staged ;;
        5) check "${tag}" "${pat}" "${tgt}" --staged --noninjective ;;
    esac
done

# Fixed instances, from a larger random search, where a component of a split has a solution
# and its search leaves core nogoods registered at levels below the split's node. Unless
# those are moved up to the node's level (Proof::keep_deeper_lines_at_active_level), a
# later backtrack deletes them while the node's own nogood does not yet subsume them, and
# the checked deletion fails.
fixed="$(dirname "$0")/split-decision"
check level_1 "${fixed}/level_1_pattern.lad" "${fixed}/level_1_target.lad"
check level_2 "${fixed}/level_2_pattern.lad" "${fixed}/level_2_target.lad"
check level_3 "${fixed}/level_3_pattern.lad" "${fixed}/level_3_target.lad" --restarts none

echo "${checked} instances checked, ${split} searched a component, ${refuted} refuted by one, ${combined} solved from components, ${fails} failed"

# otherwise the above has tested nothing
if [ "${split}" -lt 40 ] || [ "${refuted}" -lt 25 ] || [ "${combined}" -lt 12 ]; then
    echo "too few instances split, were refuted by a component, or were solved from components"
    fails=$((fails + 1))
fi
[ "${fails}" -eq 0 ]
