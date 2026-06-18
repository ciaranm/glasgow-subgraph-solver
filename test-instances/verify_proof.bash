#!/bin/bash
# Run a solver with proof logging on an instance, then check the proof with VeriPB.
# Used by the ctest proof-verification suite (see src/CMakeLists.txt).
#
# Usage:
#   verify_proof.bash <solver> <veripb> <workdir> <tag> -- <solver-args...>
#
# Exits 0 iff VeriPB prints "VERIFIED". No solver flags are added here: every flag
# (e.g. --no-clique-detection, --no-supplementals, --induced, the instance files)
# is passed explicitly after the "--", so each test declares exactly the feature
# combination it exercises.

set -u

if [ "$#" -lt 5 ]; then
    echo "usage: $0 <solver> <veripb> <workdir> <tag> -- <solver-args...>" 1>&2
    exit 2
fi

solver="$1"
veripb="$2"
workdir="$3"
tag="$4"
shift 4

if [ "${1:-}" = "--" ]; then
    shift
fi

proof="${workdir}/proof_${tag}"

"${solver}" --prove "${proof}" "$@" > "${proof}.solver.log" 2>&1
solver_status=$?
if [ "${solver_status}" -ne 0 ]; then
    echo "solver exited with status ${solver_status}:" 1>&2
    cat "${proof}.solver.log" 1>&2
    exit 1
fi

veripb_out=$("${veripb}" "${proof}.opb" "${proof}.pbp" 2>&1)
echo "${veripb_out}"

if echo "${veripb_out}" | grep -q 'VERIFIED'; then
    exit 0
else
    exit 1
fi
