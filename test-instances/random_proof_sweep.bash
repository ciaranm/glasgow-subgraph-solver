#!/bin/bash
# Generate random instances and check that the solver's proof verifies on each,
# across a spread of graph structures and option combinations. This catches
# proof-logging bugs that the fixed instances miss -- both loopless and loopy
# instances are swept, since loops + supplemental graphs (issue #56) only break on
# certain structures. Correctness of the *counts* is checked separately and
# in-process by random_homomorphism_test.
#
# Deterministic: create_random_graph is seeded and the solver is run without
# restarts, so every proof is reproducible.
#
# Usage:
#   random_proof_sweep.bash <solver> <veripb> <create_random_graph> <workdir>
#
# Exits 0 iff every generated proof verifies.

set -u

if [ "$#" -ne 4 ]; then
    echo "usage: $0 <solver> <veripb> <create_random_graph> <workdir>" 1>&2
    exit 2
fi

solver="$1"
veripb="$2"
crg="$3"
workdir="$4"

fails=0
checked=0

# Two families: loopless (no --loops) and loopy (--loops). Supplemental graphs are on
# by default in both, so the exact-path / distance-3 / degree / NDS derivations are
# exercised, and the loopy family additionally drives the loop-cancellation path
# (issue #56). Each uses the same option spread.
for loops in "" "--loops 0.3"; do
    for seed in $(seq 1 8); do
        pat="${workdir}/rps_pattern_${seed}.csv"
        tgt="${workdir}/rps_target_${seed}.csv"
        # shellcheck disable=SC2086
        "${crg}" --seed "${seed}" ${loops} 5 0.5 > "${pat}"
        # shellcheck disable=SC2086
        "${crg}" --seed "$((seed + 100))" ${loops} 8 0.45 > "${tgt}"

        # Local injectivity uses the neighbourhood-injectivity encoding and the degree
        # and NDS eliminations built on it; it does not yet support supplemental graphs in
        # proofs, so its arms carry --no-supplementals. (On loopy instances the degree/NDS
        # filters are disabled anyway, issue #58, so those arms mainly exercise the encoding
        # and counting.)
        for opts in "" "--induced" "--count-solutions" "--induced --count-solutions" "--distance3" \
                "--locally-injective --no-supplementals" "--locally-injective --no-supplementals --count-solutions"; do
            checked=$((checked + 1))
            proof="${workdir}/rps_proof_${seed}_${checked}"

            # shellcheck disable=SC2086
            if ! "${solver}" --format csv --no-clique-detection ${opts} --prove "${proof}" \
                    "${pat}" "${tgt}" > "${proof}.log" 2>&1; then
                echo "loops='${loops}' seed ${seed}, opts '${opts}': solver failed" 1>&2
                cat "${proof}.log" 1>&2
                fails=$((fails + 1))
                continue
            fi

            if ! "${veripb}" "${proof}.opb" "${proof}.pbp" 2>&1 | grep -q 'VERIFIED'; then
                echo "loops='${loops}' seed ${seed}, opts '${opts}': proof did NOT verify" 1>&2
                "${veripb}" "${proof}.opb" "${proof}.pbp" 2>&1 | tail -4 1>&2
                fails=$((fails + 1))
            fi
        done
    done
done

echo "checked ${checked} random proofs, ${fails} failure(s)"
[ "${fails}" -eq 0 ]
