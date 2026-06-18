#!/bin/bash
# Generate random instances and check that the solver's proof verifies on each,
# across a spread of graph structures and option combinations. This catches
# proof-logging bugs that the fixed instances miss -- a derivation that only breaks
# on certain structures (e.g. issue #56 would have shown up here). Correctness of
# the *counts* is checked separately and in-process by random_homomorphism_test.
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

# loopless random graphs (no --loops): supplemental graphs are on by default, so the
# exact-path / distance derivations get exercised; loops + supplementals is issue #56.
for seed in $(seq 1 8); do
    pat="${workdir}/rps_pattern_${seed}.csv"
    tgt="${workdir}/rps_target_${seed}.csv"
    "${crg}" --seed "${seed}" 5 0.5 > "${pat}"
    "${crg}" --seed "$((seed + 100))" 8 0.45 > "${tgt}"

    for opts in "" "--induced" "--count-solutions" "--induced --count-solutions"; do
        checked=$((checked + 1))
        proof="${workdir}/rps_proof_${seed}_${checked}"

        # shellcheck disable=SC2086
        if ! "${solver}" --format csv --no-clique-detection ${opts} --prove "${proof}" \
                "${pat}" "${tgt}" > "${proof}.log" 2>&1; then
            echo "seed ${seed}, opts '${opts}': solver failed" 1>&2
            cat "${proof}.log" 1>&2
            fails=$((fails + 1))
            continue
        fi

        if ! "${veripb}" "${proof}.opb" "${proof}.pbp" 2>&1 | grep -q 'VERIFIED'; then
            echo "seed ${seed}, opts '${opts}': proof did NOT verify" 1>&2
            "${veripb}" "${proof}.opb" "${proof}.pbp" 2>&1 | tail -4 1>&2
            fails=$((fails + 1))
        fi
    done
done

echo "checked ${checked} random proofs, ${fails} failure(s)"
[ "${fails}" -eq 0 ]
