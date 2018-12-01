The Glasgow Subgraph Solver
===========================

This is a solver for subgraph isomorphism (induced and non-induced) problems,
based upon a series of papers by subsets of Ciaran McCreesh, Patrick Prosser,
and James Trimble, at the University of Glasgow.

If you use this software for a research paper, please consider citing the
following paper:

    https://dblp.org/rec/html/conf/cp/McCreeshP15

If you use this solver in a non-research setting, please get in touch if you
can. This software is an output of taxpayer funded research, and it is very
helpful for us if we can demonstrate real-world impact when we write grant
applications.

If you have queries, please contact:

    ciaran.mccreesh@glasgow.ac.uk

Compiling
---------

To build, type 'make'. You will need a C++17 compiler (we use GCC 7.3) and
Boost (built with threads enabled).

Running
-------

To run:

    ./glasgow_subgraph_solver [ --induced ] pattern-file target-file

You may need to increase the stack space, for larger graphs. In bash this is
done as follows:

    ulimit -s 1048576

We try to auto-detect the input format. We can read LAD, Labelled LAD, and
DIMACS 2 formatted graphs. If in doubt, it's best to use the format described here:

    https://perso.liris.cnrs.fr/christine.solnon/SIP.html

