#!/bin/bash
# Exercise the gss-graph JSON format on the real test corpus, end to end through the
# binaries rather than in process.
#
# Two things are checked for each instance pair:
#
#   1. Converting to JSON and solving gives exactly the answer the original format
#      gives -- same solution count, same node count. A format that loses or alters
#      a graph property would show up here as a different search.
#   2. Converting the JSON again is byte-identical. The writer claims to emit a
#      canonical form, and that claim is what makes the round-trip test in
#      json_graph_test meaningful; this checks it on real instances too.
#
# Usage:
#   json_round_trip.bash <solver> <convert_to_json> <instances_dir> <workdir>

set -u

if [ "$#" -ne 4 ]; then
    echo "usage: $0 <solver> <convert_to_json> <instances_dir> <workdir>" 1>&2
    exit 2
fi

solver="$1"
convert="$2"
instances="$3"
workdir="$4"

mkdir -p "$workdir" || exit 1

status=0

# Extract one "key = value" line from the solver's output.
stat_of() {
    sed -n "s/^$2 = //p" "$1"
}

# convert <output.json> <format> <input>
convert_one() {
    if ! "$convert" --format "$2" "$3" > "$1"; then
        echo "FAIL: could not convert $3 to JSON" 1>&2
        status=1
        return 1
    fi

    # Converting the JSON again must reproduce it exactly.
    if ! "$convert" --format json "$1" > "$1.again"; then
        echo "FAIL: could not re-convert $1" 1>&2
        status=1
        return 1
    fi
    if ! cmp -s "$1" "$1.again"; then
        echo "FAIL: $1 is not a fixed point of the writer" 1>&2
        diff -u "$1" "$1.again" 1>&2
        status=1
        return 1
    fi

    return 0
}

# compare <name> <format> <pattern> <target> [solver args...]
compare() {
    local name="$1" format="$2" pattern="$3" target="$4"
    shift 4

    local p="$workdir/$name.pattern.json" t="$workdir/$name.target.json"
    convert_one "$p" "$format" "$pattern" || return
    convert_one "$t" "$format" "$target" || return

    local original="$workdir/$name.original.out" via_json="$workdir/$name.json.out"
    "$solver" --count-solutions --format "$format" "$@" "$pattern" "$target" > "$original" 2>&1
    "$solver" --count-solutions --format json "$@" "$p" "$t" > "$via_json" 2>&1

    for key in solution_count nodes; do
        local a b
        a="$(stat_of "$original" "$key")"
        b="$(stat_of "$via_json" "$key")"
        if [ -z "$a" ]; then
            echo "FAIL: $name: no $key in the original-format run" 1>&2
            cat "$original" 1>&2
            status=1
        elif [ "$a" != "$b" ]; then
            echo "FAIL: $name: $key is $a from $format but $b via JSON" 1>&2
            status=1
        else
            echo "ok: $name: $key = $a via both $format and JSON"
        fi
    done
}

compare trident csv "$instances/trident.csv" "$instances/longtrident.csv"
compare c3_induced csv "$instances/c3.csv" "$instances/c3c2.csv" --induced
compare small_large lad "$instances/small" "$instances/large"
compare li_exact_path csv "$instances/li_exact_path_pattern.csv" "$instances/li_exact_path_target.csv"

exit $status
