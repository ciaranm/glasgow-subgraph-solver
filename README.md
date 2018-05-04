To build, type 'make'. You will need a C++14 compiler (we use GCC 7.3) and
Boost (built with threads enabled).

To run:

    ./solve_subgraph_isomorphism [ --induced ] algorithm-name pattern-file target-file

You may need to increase the stack space, for larger graphs. In bash this is
done as follows:

    ulimit -s 1048576

The algorithm variations are:

    simple
    restarting

There is also a customisable-sequential algorithm which has lots of parameters
for you to tune.

The pattern and target files should be in the format used by LAD (see the
instances/ directory). We can also read DIMACS format, using --dimacs.
