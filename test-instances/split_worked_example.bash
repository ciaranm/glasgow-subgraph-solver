#!/bin/bash
# The worked example in dev_docs/split-counting.md: the subproblem with G at 1 and H at 4,
# pinned by vertex labels, has 28 solutions, and its statement renders exactly as the
# document shows it (test-instances/split/worked_example_pinned.render). The value ordering
# is fixed, because the default one is random, and standard libraries turn the same random
# numbers into different choices, so the statement's alternatives would come out in a
# different order on different platforms.
#
# Usage:
#   split_worked_example.bash <solver> <count_statement> <workdir>

set -u

if [ "$#" -ne 3 ]; then
    echo "usage: $0 <solver> <count_statement> <workdir>" 1>&2
    exit 2
fi

solver="$1"
count_statement="$2"
workdir="$3"
mkdir -p "${workdir}" || exit 1
dir="$(dirname "$0")/split"

out=$("${solver}" --count-solutions --value-ordering none --split-components 1000 --split-statement "${workdir}/pinned.stmt" \
    --format csv "${dir}/worked_example_pinned_pattern.csv" "${dir}/worked_example_pinned_target.csv") || exit 1
if ! echo "${out}" | grep -q '^solution_count = 28$'; then
    echo "expected 28 solutions:"
    echo "${out}"
    exit 1
fi

"${count_statement}" --render --format csv --pattern "${dir}/worked_example_pinned_pattern.csv" \
    --target "${dir}/worked_example_pinned_target.csv" "${workdir}/pinned.stmt" > "${workdir}/pinned.render" || exit 1
diff "${dir}/worked_example_pinned.render" "${workdir}/pinned.render"
