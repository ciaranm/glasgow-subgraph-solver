TARGET := glasgow_common_bigraph_solver

SOURCES := \
    glasgow_common_bigraph_solver.cc

TGT_PREREQS := run-tests.bash libcommon.a
ifeq ($(shell uname -s), Linux)
TGT_LDLIBS := libcommon.a $(boost_ldlibs) -lstdc++fs
else
TGT_LDLIBS := libcommon.a $(boost_ldlibs)
endif


