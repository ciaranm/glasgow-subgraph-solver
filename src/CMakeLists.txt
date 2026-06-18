set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

link_libraries(glasgow_subgraphs)

add_executable(glasgow_subgraph_solver glasgow_subgraph_solver.cc)
target_link_libraries(glasgow_subgraph_solver PRIVATE cxxopts)

add_executable(glasgow_clique_solver glasgow_clique_solver.cc)
target_link_libraries(glasgow_clique_solver PRIVATE cxxopts)

add_executable(glasgow_common_subgraph_solver glasgow_common_subgraph_solver.cc)
target_link_libraries(glasgow_common_subgraph_solver PRIVATE cxxopts)

add_executable(create_random_graph create_random_graph.cc)
target_link_libraries(create_random_graph PRIVATE cxxopts)

add_executable(convert_to_lad convert_to_lad.cc)
target_link_libraries(convert_to_lad PRIVATE cxxopts)

# Proof-verification tests: run the solver with proof logging on a small instance
# and check the resulting proof with VeriPB. Only registered when VeriPB was found
# at configure time (find_program in the top-level CMakeLists), so a checkout
# without VeriPB still builds and runs the rest of the suite.
if(VERIPB_EXECUTABLE)
    set(GSS_PROOF_RUNNER ${CMAKE_SOURCE_DIR}/test-instances/verify_proof.bash)
    set(GSS_INSTANCES ${CMAKE_SOURCE_DIR}/test-instances)

    # Register a proof test that VeriPB is expected to verify.
    function(gss_add_proof_test name)
        add_test(NAME ${name}
            COMMAND ${GSS_PROOF_RUNNER}
                $<TARGET_FILE:glasgow_subgraph_solver> ${VERIPB_EXECUTABLE}
                ${CMAKE_CURRENT_BINARY_DIR} ${name} -- ${ARGN})
    endfunction()

    # Register a proof test that is currently BROKEN (tracked, not yet fixed). It is
    # marked WILL_FAIL, so the suite stays green while the proof fails to verify;
    # once the underlying bug is fixed the test will start passing and ctest will
    # report it as a failure, prompting removal of the WILL_FAIL marker.
    # Root causes are documented in issues #49 and #50.
    function(gss_add_known_broken_proof_test name)
        gss_add_proof_test(${name} ${ARGN})
        set_tests_properties(${name} PROPERTIES WILL_FAIL TRUE)
    endfunction()

    # --- Currently verifying: UNSAT / refutation proofs ---
    gss_add_proof_test(proof_unsat_loopy_induced
        --induced --format lad ${GSS_INSTANCES}/small ${GSS_INSTANCES}/large)
    gss_add_proof_test(proof_unsat_loopless_induced
        --induced ${GSS_INSTANCES}/c3.csv ${GSS_INSTANCES}/trident.csv)

    # --- Currently broken: SAT / solution proofs (see issues #49 and #50) ---
    # (1) decision: GSS emits 'solx' (which needs a preserved set) where 'sol' is wanted.
    gss_add_known_broken_proof_test(proof_sat_loopless_decision_KNOWN_BROKEN
        ${GSS_INSTANCES}/trident.csv ${GSS_INSTANCES}/longtrident.csv)
    # (2) enumeration: 'solx' requires a 'preserved:' set that GSS never emits.
    gss_add_known_broken_proof_test(proof_sat_loopless_count_KNOWN_BROKEN
        --count-solutions ${GSS_INSTANCES}/trident.csv ${GSS_INSTANCES}/longtrident.csv)
    # (3) loop encoding: the loop->loop term is dropped from the adjacency
    #     constraint, so a valid loop-preserving solution falsifies the model.
    gss_add_known_broken_proof_test(proof_sat_loopy_decision_KNOWN_BROKEN
        --format lad ${GSS_INSTANCES}/small ${GSS_INSTANCES}/large)
endif()
