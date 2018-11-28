TARGET := solve_subgraph_isomorphism

SOURCES := \
    formats/input_graph.cc \
    formats/graph_file_error.cc \
    formats/lad.cc \
    formats/dimacs.cc \
    formats/csv.cc \
    formats/read_file_format.cc \
    fixed_bit_set.cc \
    params.cc \
    result.cc \
    solve_subgraph_isomorphism.cc \
    solver.cc

TGT_LDLIBS := $(boost_ldlibs)

