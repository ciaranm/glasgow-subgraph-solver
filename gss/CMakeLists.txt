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
        innards/clique_size_constraints.cc
        innards/cost_bound.cc
        innards/filter_activations.cc
        innards/graph_traits.cc
        innards/homomorphism_domain.cc
        innards/homomorphism_model.cc
        innards/homomorphism_proofs.cc
        innards/homomorphism_searcher.cc
        innards/homomorphism_traits.cc
        innards/proof.cc
        innards/reification.cc
        innards/supplemental_graphs.cc
        innards/svo_bitset.cc
        innards/thread_utils.cc
        innards/verify.cc
        innards/watches.cc
        formats/csv.cc
        formats/dimacs.cc
        formats/graph_file_error.cc
        formats/input_graph.cc
        formats/json_graph.cc
        formats/lad.cc
        formats/read_file_format.cc
        formats/vfmcs.cc
        utils/hashing_utils.hh
        utils/vertex_name_map.hh)

target_include_directories(glasgow_subgraphs PUBLIC ${GMP_INCLUDE_DIR})
target_link_libraries(glasgow_subgraphs PUBLIC
        ${cxxopts_LIBRARIES}
        ${GMP_LIBRARY}
        ${GMPXX_LIBRARY}
)
# PRIVATE: nlohmann/json is slow to compile and is confined to json_graph.cc, so it
# must not leak into anything that merely links the library.
target_link_libraries(glasgow_subgraphs PRIVATE nlohmann_json::nlohmann_json)

add_executable(subgraph_isomorphism_test subgraph_isomorphism_test.cc)
target_link_libraries(subgraph_isomorphism_test PRIVATE glasgow_subgraphs Catch2::Catch2WithMain)
add_test(NAME subgraph_isomorphism_test COMMAND $<TARGET_FILE:subgraph_isomorphism_test>)

add_executable(weighted_homomorphism_test weighted_homomorphism_test.cc)
target_link_libraries(weighted_homomorphism_test PRIVATE glasgow_subgraphs Catch2::Catch2WithMain)
add_test(NAME weighted_homomorphism_test COMMAND $<TARGET_FILE:weighted_homomorphism_test>)

add_executable(homomorphism_test homomorphism_test.cc)
target_link_libraries(homomorphism_test PRIVATE glasgow_subgraphs Catch2::Catch2WithMain)
add_test(NAME homomorphism_test COMMAND $<TARGET_FILE:homomorphism_test>)

add_executable(restarts_test restarts_test.cc)
target_link_libraries(restarts_test PRIVATE glasgow_subgraphs Catch2::Catch2WithMain)
add_test(NAME restarts_test COMMAND $<TARGET_FILE:restarts_test>)

add_executable(loooong_test loooong_test.cc)
target_link_libraries(loooong_test PRIVATE glasgow_subgraphs Catch2::Catch2WithMain)
add_test(NAME loooong_test COMMAND $<TARGET_FILE:loooong_test>)

add_executable(svo_bitset_test innards/svo_bitset_test.cc)
target_link_libraries(svo_bitset_test PRIVATE glasgow_subgraphs Catch2::Catch2WithMain)
add_test(NAME svo_bitset_test COMMAND $<TARGET_FILE:svo_bitset_test>)

add_executable(json_graph_test formats/json_graph_test.cc)
target_link_libraries(json_graph_test PRIVATE glasgow_subgraphs Catch2::Catch2WithMain)
add_test(NAME json_graph_test COMMAND $<TARGET_FILE:json_graph_test>)

add_executable(formats_test formats/formats_test.cc)
target_link_libraries(formats_test PRIVATE glasgow_subgraphs Catch2::Catch2WithMain)
add_test(NAME formats_test COMMAND $<TARGET_FILE:formats_test>)

add_executable(clique_test clique_test.cc)
target_link_libraries(clique_test PRIVATE glasgow_subgraphs Catch2::Catch2WithMain)
add_test(NAME clique_test COMMAND $<TARGET_FILE:clique_test>)

add_executable(common_subgraph_test common_subgraph_test.cc)
target_link_libraries(common_subgraph_test PRIVATE glasgow_subgraphs Catch2::Catch2WithMain)
add_test(NAME common_subgraph_test COMMAND $<TARGET_FILE:common_subgraph_test>)

add_executable(homomorphism_matrix_test homomorphism_matrix_test.cc)
target_link_libraries(homomorphism_matrix_test PRIVATE glasgow_subgraphs Catch2::Catch2WithMain)
add_test(NAME homomorphism_matrix_test COMMAND $<TARGET_FILE:homomorphism_matrix_test>)

add_executable(random_homomorphism_test random_homomorphism_test.cc)
target_link_libraries(random_homomorphism_test PRIVATE glasgow_subgraphs Catch2::Catch2WithMain)
add_test(NAME random_homomorphism_test COMMAND $<TARGET_FILE:random_homomorphism_test>)

add_executable(option_sweep_test option_sweep_test.cc)
target_link_libraries(option_sweep_test PRIVATE glasgow_subgraphs Catch2::Catch2WithMain)
# The sweep checks its per-cell outcomes against a table checked in next to it; see
# dev_docs/option-compatibility.md for what the table is for and how to regenerate it.
target_compile_definitions(option_sweep_test PRIVATE
        GSS_OPTION_SWEEP_GOLDEN="${CMAKE_CURRENT_SOURCE_DIR}/option_sweep_golden.tsv")
add_test(NAME option_sweep_test COMMAND $<TARGET_FILE:option_sweep_test>)

add_executable(extra_shapes_test extra_shapes_test.cc)
target_link_libraries(extra_shapes_test PRIVATE glasgow_subgraphs Catch2::Catch2WithMain)
add_test(NAME extra_shapes_test COMMAND $<TARGET_FILE:extra_shapes_test>)
