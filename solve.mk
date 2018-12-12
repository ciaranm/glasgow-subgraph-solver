TARGET := glasgow_subgraph_solver

SOURCES := \
    formats/input_graph.cc \
    formats/graph_file_error.cc \
    formats/lad.cc \
    formats/dimacs.cc \
    formats/csv.cc \
    formats/read_file_format.cc \
    fixed_bit_set.cc \
    glasgow_subgraph_solver.cc \
    params.cc \
    restarts.cc \
    result.cc \
    solver.cc \
    verify.cc

TGT_LDLIBS := $(boost_ldlibs)
TGT_PREREQS := run-tests.bash

TGT_POSTMAKE := bash ./run-tests.bash

