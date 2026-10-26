#!/bin/bash
# Check a hand-written counting certificate (see dev_docs/decomposition-counting.md,
# phase 2): that the translator still turns the statement into exactly the counting
# formula and side-condition instance the proofs were written against, that the count is
# the expected one, and, if VeriPB is given, that both proofs still verify.
#
# Usage:
#   check_certificate.bash <count_statement> <veripb or ""> <certificate dir> <workdir> <expected count>
#
# The certificate directory holds statement, input.opb (the solver's encoding),
# output.opb and side.opb (the translator's output), main.pbp and side.pbp.

set -u

if [ "$#" -ne 5 ]; then
    echo "usage: $0 <count_statement> <veripb or \"\"> <certificate dir> <workdir> <expected count>" 1>&2
    exit 2
fi

count_statement="$1"
veripb="$2"
dir="$3"
work="$4"
expected="$5"

mkdir -p "${work}" || exit 1

out=$("${count_statement}" "${dir}/statement" --opb "${work}/output.opb" --side-conditions "${work}/side.opb") || {
    echo "count_statement failed: ${out}" 1>&2
    exit 1
}
echo "${out}"
if ! echo "${out}" | grep -qx "count = ${expected}"; then
    echo "expected count = ${expected}" 1>&2
    exit 1
fi

for f in output.opb side.opb; do
    if ! cmp "${dir}/${f}" "${work}/${f}"; then
        echo "the translator's ${f} no longer matches the one the proofs were written against" 1>&2
        exit 1
    fi
done

if [ -z "${veripb}" ]; then
    echo "no VeriPB, so the proofs are not checked"
    exit 0
fi

main_out=$("${veripb}" --force-checked-deletion "${dir}/input.opb" "${dir}/main.pbp" "${work}/output.opb" 2>&1)
echo "${main_out}"
if ! echo "${main_out}" | grep -q '^s VERIFIED OUTPUT EQUIENUMERABLE'; then
    echo "the main proof does not verify" 1>&2
    exit 1
fi

side_out=$("${veripb}" "${work}/side.opb" "${dir}/side.pbp" 2>&1)
echo "${side_out}"
if ! echo "${side_out}" | grep -q '^s VERIFIED UNSATISFIABLE'; then
    echo "the side-condition proof does not verify" 1>&2
    exit 1
fi
