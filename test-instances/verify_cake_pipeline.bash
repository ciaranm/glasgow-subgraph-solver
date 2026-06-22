#!/bin/bash
# Verify a subgraph-isomorphism proof end to end through the formally verified
# CakePB toolchain. Used by the ctest "cake_*" suite (see src/CMakeLists.txt),
# which is only registered when cake_pb_iso is found at configure time.
#
# The pipeline (mirroring the constraint solver's run_cake_pb_cp.bash) is:
#   1. solver --prove NAME           -> NAME.opb, NAME.pbp
#   2. cake_pb_iso pat tgt           -> NAME.cakeopb   (cake's own trusted OPB)
#   3. veripb NAME.cakeopb NAME.pbp --elaborate NAME.core.pbp
#                                    (elaborate the user-friendly proof, against
#                                     cake's OPB, down to the kernel subset cake reads)
#   4. cake_pb_iso pat tgt NAME.core.pbp  (cake checks the elaborated core proof)
#
# We check cake's OPB rather than the solver's, so loop->loop mappings verify even
# though they do not under plain VeriPB against the solver's own OPB (issue #49):
# cake derives the loop-correct encoding directly from the LAD files.
#
# Usage:
#   verify_cake_pipeline.bash <solver> <veripb> <cake> <workdir> <tag> <pat.lad> <tgt.lad> -- <extra solver args...>
#
# Exits 0 iff cake_pb_iso prints "VERIFIED" on the elaborated core proof.

set -u

if [ "$#" -lt 7 ]; then
    echo "usage: $0 <solver> <veripb> <cake> <workdir> <tag> <pat.lad> <tgt.lad> -- <extra solver args...>" 1>&2
    exit 2
fi

solver="$1"
veripb="$2"
cake="$3"
workdir="$4"
tag="$5"
pat="$6"
tgt="$7"
shift 7

if [ "${1:-}" = "--" ]; then
    shift
fi

proof="${workdir}/cake_${tag}"

# 1. solve with proof logging (proof logging requires --no-clique-detection)
"${solver}" --format lad --no-clique-detection --prove "${proof}" "$@" "${pat}" "${tgt}" > "${proof}.solver.log" 2>&1
solver_status=$?
if [ "${solver_status}" -ne 0 ]; then
    echo "solver exited with status ${solver_status}:" 1>&2
    cat "${proof}.solver.log" 1>&2
    exit 1
fi

# 2. cake's own trusted OPB encoding, derived from the LAD files
if ! "${cake}" "${pat}" "${tgt}" > "${proof}.cakeopb" 2>"${proof}.cake-opb.log"; then
    echo "cake_pb_iso failed to produce an OPB:" 1>&2
    cat "${proof}.cake-opb.log" 1>&2
    exit 1
fi

# 3. elaborate the solver's proof, against cake's OPB, to the kernel subset
elaborate_out=$("${veripb}" "${proof}.cakeopb" "${proof}.pbp" --elaborate "${proof}.core.pbp" 2>&1)
if ! echo "${elaborate_out}" | grep -q 'VERIFIED'; then
    echo "veripb elaboration did not verify:" 1>&2
    echo "${elaborate_out}" 1>&2
    exit 1
fi

# 4. cake checks the elaborated core proof
cake_out=$("${cake}" "${pat}" "${tgt}" "${proof}.core.pbp" 2>&1)
echo "${cake_out}"

if echo "${cake_out}" | grep -q 'VERIFIED'; then
    exit 0
else
    exit 1
fi
