TARGET := plot_glasgow_solver_outputs

SOURCES := \
    plot_glasgow_solver_outputs.cc

ifeq ($(shell uname -s), Linux)
TGT_LDLIBS := $(boost_ldlibs) -lstdc++fs
else
TGT_LDLIBS := $(boost_ldlibs)
endif

