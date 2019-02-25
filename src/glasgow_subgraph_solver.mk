TARGET := glasgow_subgraph_solver

SOURCES := \
    glasgow_subgraph_solver.cc

TGT_LDLIBS := libcommon.a $(boost_ldlibs)
TGT_PREREQS := run-tests.bash libcommon.a

TGT_POSTMAKE := bash ./run-tests.bash

