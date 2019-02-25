BUILD_DIR := intermediate
TARGET_DIR := ./

SUBMAKEFILES := src/common.mk src/glasgow_subgraph_solver.mk src/glasgow_clique_solver.mk

boost_ldlibs := -lboost_thread -lboost_system -lboost_program_options

override CXXFLAGS += -O3 -march=native -std=c++17 -Isrc/ -W -Wall -g -ggdb3 -pthread
override LDFLAGS += -pthread
