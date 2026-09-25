#!/bin/sh
# Encode one pattern with the local-polytope PB encoding and run RoundingSat on it, with a 60s
# limit. Prints: tag, status, objective (scaled by 10^4), seconds.
#
#   GRAPH3=~/graph3 ROUNDINGSAT=/path/to/roundingsat ./rs_one.sh pattern.csv tag
#
# Leaves opb/<tag>.opb and opb/<tag>.out behind.

set -e
here=$(dirname "$0")
f=$1
t=$2
mkdir -p opb
"$here/wsip" opb-strong "$f" "$GRAPH3/graph3.csv" > "opb/$t.raw"
python3 "$here/renum.py" "opb/$t.raw" > "opb/$t.opb"
rm "opb/$t.raw"
s=$(perl -MTime::HiRes=time -e 'print time')
perl -e 'alarm 60; exec @ARGV' "$ROUNDINGSAT" "opb/$t.opb" > "opb/$t.out" 2>&1 || true
e=$(perl -MTime::HiRes=time -e 'print time')
st=$(grep -E '^s ' "opb/$t.out" | cut -c3-)
o=$(grep -E '^o ' "opb/$t.out" | tail -1 | cut -c3-)
printf "%s\t%s\t%s\t%.2f\n" "$t" "${st:-TIMEOUT}" "$o" "$(echo "$e - $s" | bc)"
