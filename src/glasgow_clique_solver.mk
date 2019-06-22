TARGET := glasgow_clique_solver

SOURCES := \
    glasgow_clique_solver.cc

TGT_LDLIBS := libcommon.a $(boost_ldlibs) -lstdc++fs
TGT_PREREQS := libcommon.a

