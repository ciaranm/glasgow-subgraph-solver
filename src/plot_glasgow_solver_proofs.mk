TARGET := plot_glasgow_solver_proofs

SOURCES := \
    plot_glasgow_solver_proofs.cc

ifeq ($(shell uname -s), Linux)
TGT_LDLIBS := $(boost_ldlibs) -lstdc++fs
else
TGT_LDLIBS := $(boost_ldlibs)
endif

