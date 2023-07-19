#!/bin/bash

set -x

if ! grep '^status = true$' <(./build/glasgow_subgraph_solver --format lad test-instances/small test-instances/large ) ; then
    echo "non-induced test failed" 1>&1
    exit 1
fi

if ! grep '^status = false$' <(./build/glasgow_subgraph_solver --induced --format lad test-instances/small test-instances/large ) ; then
    echo "induced test failed" 1>&1
    exit 1
fi

if ! grep '^solution_count = 6$' <(./build/glasgow_subgraph_solver --count-solutions --format lad test-instances/small test-instances/large ) ; then
    echo "non-induced enumerate test failed" 1>&1
    exit 1
fi

if ! grep '^solution_count = 12$' <(./build/glasgow_subgraph_solver --count-solutions test-instances/trident.csv test-instances/longtrident.csv ) ; then
    echo "trident enumerate test failed" 1>&1
    exit 1
fi

if ! grep '^solution_count = 6$' <(./build/glasgow_subgraph_solver --count-solutions --induced test-instances/c3.csv test-instances/c3c2.csv ) ; then
    echo "induced cyclic enumerate test failed" 1>&1
    exit 1
fi

true

