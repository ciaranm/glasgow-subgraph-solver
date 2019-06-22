TARGET := glasgow_subgraph_solver

SOURCES := \
    glasgow_subgraph_solver.cc

TGT_PREREQS := run-tests.bash libcommon.a
ifeq ($(shell uname -s), Linux)
TGT_LDLIBS := libcommon.a $(boost_ldlibs) -lstdc++fs
else
TGT_LDLIBS := libcommon.a $(boost_ldlibs)
endif

TGT_POSTMAKE := bash ./run-tests.bash

