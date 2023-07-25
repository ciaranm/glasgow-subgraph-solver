set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

find_package(Boost 1.74 COMPONENTS iostreams container thread REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

add_library(glasgow_subgraphs
        clique.cc
        common_subgraph.cc
        configuration.cc
        homomorphism.cc
        restarts.cc
        sip_decomposer.cc
        timeout.cc
        innards/cheap_all_different.cc
        innards/graph_traits.cc
        innards/homomorphism_domain.cc
        innards/homomorphism_model.cc
        innards/homomorphism_searcher.cc
        innards/homomorphism_traits.cc
        innards/lackey.cc
        innards/proof.cc
        innards/svo_bitset.cc
        innards/symmetries.cc
        innards/thread_utils.cc
        innards/verify.cc
        innards/watches.cc
        formats/csv.cc
        formats/dimacs.cc
        formats/graph_file_error.cc
        formats/input_graph.cc
        formats/lad.cc
        formats/read_file_format.cc
        formats/vfmcs.cc)

target_link_libraries(glasgow_subgraphs LINK_PUBLIC ${Boost_LIBRARIES})
