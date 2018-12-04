#!/bin/bash -x

if ! grep '^status = true$' <(./glasgow_subgraph_solver --format lad test-instances/{small,large} ) ; then
    echo "non-induced test failed" 1>&1
    exit 1
fi

if ! grep '^status = false$' <(./glasgow_subgraph_solver --induced --format lad test-instances/{small,large} ) ; then
    echo "induced test failed" 1>&1
    exit 1
fi

if ! grep '^solution_count = 6$' <(./glasgow_subgraph_solver --enumerate --format lad test-instances/{small,large} ) ; then
    echo "non-induced enumerate test failed" 1>&1
    exit 1
fi

true

