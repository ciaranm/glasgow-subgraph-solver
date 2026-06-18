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

    # Register a proof test: run <target> --prove on an instance, then check with
    # VeriPB. The second argument is the solver target; everything after it is passed
    # to the solver verbatim, so each test spells out the feature combination it
    # exercises (which is what the coverage build measures).
    function(gss_add_proof_test name target)
        add_test(NAME ${name}
            COMMAND ${GSS_PROOF_RUNNER}
                $<TARGET_FILE:${target}> ${VERIPB_EXECUTABLE}
                ${CMAKE_CURRENT_BINARY_DIR} ${name} -- ${ARGN})
    endfunction()

    # The subgraph solver always needs --no-clique-detection for proof logging. SUB_MIN
    # is the minimal proof-compatible combination (also no supplemental graphs / NDS),
    # used where a test isn't specifically exercising those features.
    set(SUB glasgow_subgraph_solver --no-clique-detection)
    set(SUB_MIN ${SUB} --no-supplementals --no-nds)

    # --- Refutation (UNSAT) proofs ---
    gss_add_proof_test(proof_unsat_loopy_induced ${SUB_MIN}
        --induced --format lad ${GSS_INSTANCES}/small ${GSS_INSTANCES}/large)
    gss_add_proof_test(proof_unsat_loopless_induced ${SUB_MIN}
        --induced ${GSS_INSTANCES}/c3.csv ${GSS_INSTANCES}/trident.csv)

    # --- Solution and enumeration proofs ---
    # A preserved set is declared in the OPB and every solution is logged with the
    # solx rule at the top proof level (so it survives backtracking deletions).
    gss_add_proof_test(proof_sat_loopless_decision ${SUB_MIN}
        ${GSS_INSTANCES}/trident.csv ${GSS_INSTANCES}/longtrident.csv)
    gss_add_proof_test(proof_count_loopless ${SUB_MIN}
        --count-solutions ${GSS_INSTANCES}/trident.csv ${GSS_INSTANCES}/longtrident.csv)
    gss_add_proof_test(proof_partial_loopless ${SUB_MIN}
        --print-all-solutions --solution-limit 5 ${GSS_INSTANCES}/trident.csv ${GSS_INSTANCES}/longtrident.csv)
    gss_add_proof_test(proof_partial_count ${SUB_MIN}
        --count-solutions --solution-limit 5 ${GSS_INSTANCES}/trident.csv ${GSS_INSTANCES}/longtrident.csv)

    # loop->loop mappings verify under plain VeriPB now the self-loop term is kept in
    # the adjacency constraint (issue #49). (Supplemental graphs with loops are issue
    # #56, hence --no-supplementals here.)
    gss_add_proof_test(proof_sat_loopy_decision ${SUB_MIN}
        --format lad ${GSS_INSTANCES}/small ${GSS_INSTANCES}/large)
    gss_add_proof_test(proof_count_loopy ${SUB_MIN}
        --count-solutions --format lad ${GSS_INSTANCES}/small ${GSS_INSTANCES}/large)

    # --- Preprocessing / supplemental-graph proofs (loopless; loops are issue #56) ---
    # supplemental (exact-path) graphs: create_exact_path_graphs derivations.
    gss_add_proof_test(proof_supplementals ${SUB} --no-nds
        ${GSS_INSTANCES}/trident.csv ${GSS_INSTANCES}/longtrident.csv)
    # distance-3 filtering: create_distance3_graphs derivations.
    gss_add_proof_test(proof_distance3 ${SUB} --no-nds --distance3
        ${GSS_INSTANCES}/trident.csv ${GSS_INSTANCES}/longtrident.csv)
    # neighbourhood degree sequences: incompatible_by_nds derivations.
    gss_add_proof_test(proof_nds ${SUB} --no-supplementals
        ${GSS_INSTANCES}/trident.csv ${GSS_INSTANCES}/longtrident.csv)
    # clique-size constraints, then the clique-size argument that eliminates an
    # assignment (K4 has a triangle in each neighbourhood; K3,3 is triangle-free, so
    # degree passes but the clique bound fires): the hom-clique proof.
    gss_add_proof_test(proof_cliques ${SUB} --no-supplementals --no-nds --cliques
        ${GSS_INSTANCES}/trident.csv ${GSS_INSTANCES}/longtrident.csv)
    gss_add_proof_test(proof_hom_clique ${SUB} --no-supplementals --no-nds --cliques
        --format lad ${GSS_INSTANCES}/k4.lad ${GSS_INSTANCES}/k33.lad)
    # a pattern vertex whose degree exceeds every target's: domain empties at init.
    gss_add_proof_test(proof_empty_domain ${SUB_MIN}
        --format lad ${GSS_INSTANCES}/star5.lad ${GSS_INSTANCES}/path5.lad)
    # verbose proof comments (show_domains / propagated).
    gss_add_proof_test(proof_verbose ${SUB_MIN} --verbose-proofs
        ${GSS_INSTANCES}/trident.csv ${GSS_INSTANCES}/longtrident.csv)

    # --- Clique and common-subgraph solver proofs (optimisation, BOUNDS) ---
    gss_add_proof_test(proof_clique glasgow_clique_solver ${GSS_INSTANCES}/clique.clq)
    gss_add_proof_test(proof_mcs glasgow_common_subgraph_solver
        ${GSS_INSTANCES}/small ${GSS_INSTANCES}/large)
    gss_add_proof_test(proof_mcs_connected glasgow_common_subgraph_solver --connected
        ${GSS_INSTANCES}/small ${GSS_INSTANCES}/large)
    gss_add_proof_test(proof_mcs_clique glasgow_common_subgraph_solver --clique
        ${GSS_INSTANCES}/small ${GSS_INSTANCES}/large)

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
        # (--no-supplementals: supplemental graphs with loops are issue #56)
        gss_add_cake_test(cake_loop_decision loop_vertex.lad loop_target.lad --no-supplementals)
        gss_add_cake_test(cake_loop_count loop_vertex.lad loop_target.lad --no-supplementals --count-solutions)
    endif()
endif()
