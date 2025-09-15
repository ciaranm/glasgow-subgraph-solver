set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

find_path(GMP_INCLUDE_DIR gmp.h)
find_library(GMP_LIBRARY gmp)
find_library(GMPXX_LIBRARY gmpxx)

if (NOT GMP_INCLUDE_DIR OR NOT GMP_LIBRARY OR NOT GMPXX_LIBRARY)
    message(FATAL_ERROR "Could not find GMP/GMPXX libraries. Please install GMP and try again.")
endif()

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
        formats/vfmcs.cc
        utils/hashing_utils.hh
        utils/vertex_name_map.hh
        utils/cout_formatting.hh
        utils/cout_formatting.cc
)

target_include_directories(glasgow_subgraphs PUBLIC ${GMP_INCLUDE_DIR})
target_link_libraries(glasgow_subgraphs PUBLIC
        ${cxxopts_LIBRARIES}
        ${GMP_LIBRARY}
        ${GMPXX_LIBRARY}
        nlohmann_json::nlohmann_json
)

add_executable(subgraph_isomorphism_test subgraph_isomorphism_test.cc)
target_link_libraries(subgraph_isomorphism_test PRIVATE glasgow_subgraphs Catch2::Catch2WithMain)
add_test(NAME subgraph_isomorphism_test COMMAND $<TARGET_FILE:subgraph_isomorphism_test>)

add_executable(homomorphism_test homomorphism_test.cc)
target_link_libraries(homomorphism_test PRIVATE glasgow_subgraphs Catch2::Catch2WithMain)
add_test(NAME homomorphism_test COMMAND $<TARGET_FILE:homomorphism_test>)
