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

    # --- Refutation (UNSAT) proofs ---
    gss_add_proof_test(proof_unsat_loopy_induced
        --induced --format lad ${GSS_INSTANCES}/small ${GSS_INSTANCES}/large)
    gss_add_proof_test(proof_unsat_loopless_induced
        --induced ${GSS_INSTANCES}/c3.csv ${GSS_INSTANCES}/trident.csv)

    # --- Solution and enumeration proofs ---
    # A preserved set is declared in the OPB and every solution is logged with the
    # solx rule at the top proof level (so it survives backtracking deletions).
    # decision: a single solution, concluded SATISFIABLE.
    gss_add_proof_test(proof_sat_loopless_decision
        ${GSS_INSTANCES}/trident.csv ${GSS_INSTANCES}/longtrident.csv)
    # complete enumeration: every solution logged, concluded ENUMERATION_COMPLETE n.
    gss_add_proof_test(proof_count_loopless
        --count-solutions ${GSS_INSTANCES}/trident.csv ${GSS_INSTANCES}/longtrident.csv)
    # partial enumeration: stopped at a solution limit, concluded ENUMERATION_PARTIAL n.
    gss_add_proof_test(proof_partial_loopless
        --print-all-solutions --solution-limit 5 ${GSS_INSTANCES}/trident.csv ${GSS_INSTANCES}/longtrident.csv)

    # loop->loop mappings: the target self-loop term is kept in the adjacency
    # constraint (issue #49), so loop-preserving solutions verify.
    gss_add_proof_test(proof_sat_loopy_decision
        --format lad ${GSS_INSTANCES}/small ${GSS_INSTANCES}/large)
    gss_add_proof_test(proof_count_loopy
        --count-solutions --format lad ${GSS_INSTANCES}/small ${GSS_INSTANCES}/large)

    # End-to-end verified-pipeline tests through the CakePB checker. Only registered
    # when cake_pb_iso was found (see CAKE_PB_ISO_EXECUTABLE in the top-level
    # CMakeLists); a checkout without it just skips these. The pipeline checks the
    # solver's proof against cake's own (formally verified) OPB encoding.
    if(CAKE_PB_ISO_EXECUTABLE)
        set(GSS_CAKE_RUNNER ${CMAKE_SOURCE_DIR}/test-instances/verify_cake_pipeline.bash)
        function(gss_add_cake_test name pat tgt)
            add_test(NAME ${name}
                COMMAND ${GSS_CAKE_RUNNER}
                    $<TARGET_FILE:glasgow_subgraph_solver> ${VERIPB_EXECUTABLE} ${CAKE_PB_ISO_EXECUTABLE}
                    ${CMAKE_CURRENT_BINARY_DIR} ${name}
                    ${GSS_INSTANCES}/${pat} ${GSS_INSTANCES}/${tgt} -- ${ARGN})
        endfunction()

        # decision and enumeration of a triangle into a triangle (6 isomorphisms)
        gss_add_cake_test(cake_decision triangle.lad triangle.lad)
        gss_add_cake_test(cake_count triangle.lad triangle.lad --count-solutions)
        # complete enumeration of zero solutions (no triangle in a path)
        gss_add_cake_test(cake_unsat triangle.lad path4.lad --count-solutions)
        # loop->loop mapping, checked end to end by the verified pipeline
        gss_add_cake_test(cake_loop_decision loop_vertex.lad loop_target.lad)
        gss_add_cake_test(cake_loop_count loop_vertex.lad loop_target.lad --count-solutions)
    endif()
endif()
