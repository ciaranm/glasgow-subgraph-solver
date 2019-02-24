TARGET := glasgow_clique_solver

SOURCES := \
    glasgow_clique_solver.cc

TGT_LDLIBS := $(boost_ldlibs) libcommon.a
TGT_PREREQS := libcommon.a

