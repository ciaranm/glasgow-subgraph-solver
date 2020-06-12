BUILD_DIR := intermediate
TARGET_DIR := ./

SUBMAKEFILES := \
    src/common.mk \
    src/glasgow_subgraph_solver.mk \
    src/glasgow_clique_solver.mk \
    src/glasgow_common_subgraph_solver.mk \
    src/sip_to_opb.mk \
    src/sip_to_lad.mk \
    src/plot_glasgow_solver_outputs.mk \
    src/plot_glasgow_solver_proofs.mk \
    src/create_random_graph.mk

override CXXFLAGS += -O3 -march=native -std=c++17 -Isrc/ -W -Wall -g -ggdb3 -pthread

ifeq ($(shell uname -s), Linux)
override LDFLAGS += -pthread -lstdc++fs
boost_ldlibs := -lboost_thread -lboost_system -lboost_program_options -lboost_iostreams
else
override LDFLAGS += -pthread
boost_ldlibs := -lboost_thread-mt -lboost_system-mt -lboost_program_options-mt -lboost_filesystem-mt -lboost_iostreams-mt
endif
