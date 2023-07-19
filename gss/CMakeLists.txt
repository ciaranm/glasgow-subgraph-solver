set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

find_package(Boost 1.74 COMPONENTS iostreams container thread REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

add_library(glasgow_subgraphs
        cheap_all_different.cc
        clique.cc
        common_subgraph.cc
        configuration.cc
        create_random_graph.cc
        graph_traits.cc
        homomorphism.cc
        homomorphism_domain.cc
        homomorphism_model.cc
        homomorphism_searcher.cc
        homomorphism_traits.cc
        lackey.cc
        proof.cc
        restarts.cc
        sip_decomposer.cc
        svo_bitset.cc
        symmetries.cc
        thread_utils.cc
        timeout.cc
        verify.cc
        watches.cc
        formats/csv.cc
        formats/dimacs.cc
        formats/graph_file_error.cc
        formats/input_graph.cc
        formats/lad.cc
        formats/read_file_format.cc
        formats/vfmcs.cc)

target_link_libraries(glasgow_subgraphs LINK_PUBLIC ${Boost_LIBRARIES})
