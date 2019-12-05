#!/bin/bash
input="tests.txt"
while IFS=' ' read -r p t
do
  lol=$(../glasgow_bigraph_solver big/$t.big big/$p.big --print-all-solutions )
  readarray -t y <<<"$lol"
  echo "$p $t"
  n=5
  a=(${y[n]})
  d="${a[0]}"
  while [ $d != "status" ]
  do
    echo ${y[n]}
    (( n++ ))
    a=(${y[n]})
    d="${a[0]}"
  done
  echo "---"
done < $input


