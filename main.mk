BUILD_DIR := intermediate
TARGET_DIR := ./

SUBMAKEFILES := common.mk glasgow_subgraph_solver.mk glasgow_clique_solver.mk

boost_ldlibs := -lboost_thread -lboost_system -lboost_program_options

override CXXFLAGS += -O3 -march=native -std=c++17 -I./ -W -Wall -g -ggdb3 -pthread
override LDFLAGS += -pthread
