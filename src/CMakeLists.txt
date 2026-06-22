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

# Proof-size metrics over a fixed, representative matrix (see
# dev_docs/preprocessor-refactor.md). Pure measurement: it emits a TSV of the
# supplemental-graph count and OPB/PBP line counts per configuration, the baseline that
# the solve-pipeline refactor's later phases diff against. It needs no VeriPB, so it is
# registered unconditionally and runs in every lane, doubling as a proof-emission smoke
# test (under the sanitizers it exercises proof.cc across many feature combinations).
add_test(NAME proof_metrics
    COMMAND ${CMAKE_SOURCE_DIR}/test-instances/proof_metrics.bash
        $<TARGET_FILE:glasgow_subgraph_solver>
        ${CMAKE_SOURCE_DIR}/test-instances
        ${CMAKE_CURRENT_BINARY_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}/proof_metrics.tsv)

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
    # the adjacency constraint (issue #49).
    gss_add_proof_test(proof_sat_loopy_decision ${SUB_MIN}
        --format lad ${GSS_INSTANCES}/small ${GSS_INSTANCES}/large)
    gss_add_proof_test(proof_count_loopy ${SUB_MIN}
        --count-solutions --format lad ${GSS_INSTANCES}/small ${GSS_INSTANCES}/large)

    # Loops together with supplemental graphs (issue #56): the loop-bearing adjacency
    # constraints are reduced to their loop-cancelled form before being summed into the
    # supplemental-graph / degree / distance-3 pols (Proof::loop_fix_adjacencies), and
    # the induced self-loop constraint completes the model. These exercise that path.
    gss_add_proof_test(proof_supplementals_loopy ${SUB} --no-nds
        --format lad ${GSS_INSTANCES}/small ${GSS_INSTANCES}/large)
    gss_add_proof_test(proof_supplementals_loopy_induced ${SUB} --no-nds
        --induced --format lad ${GSS_INSTANCES}/small ${GSS_INSTANCES}/large)
    gss_add_proof_test(proof_count_loopy_supplementals ${SUB} --no-nds
        --count-solutions --format lad ${GSS_INSTANCES}/small ${GSS_INSTANCES}/large)
    gss_add_proof_test(proof_distance3_loopy ${SUB} --no-nds --distance3
        --format lad ${GSS_INSTANCES}/small ${GSS_INSTANCES}/large)

    # Local injectivity: the neighbourhood-injectivity encoding + the degree and NDS
    # pigeonholes using it. Currently requires --no-supplementals (the supplemental-graph
    # derivations still need locally-injective variants). Loopless and loopy (small/large
    # have a self-loop, exercising the closed-neighbourhood form).
    gss_add_proof_test(proof_locally_injective_unsat ${SUB_MIN}
        --locally-injective --format lad ${GSS_INSTANCES}/star5.lad ${GSS_INSTANCES}/path5.lad)
    gss_add_proof_test(proof_locally_injective_count ${SUB_MIN}
        --locally-injective --count-solutions ${GSS_INSTANCES}/trident.csv ${GSS_INSTANCES}/longtrident.csv)
    gss_add_proof_test(proof_locally_injective_loopy ${SUB_MIN}
        --locally-injective --count-solutions --format lad ${GSS_INSTANCES}/small ${GSS_INSTANCES}/large)
    # NDS pigeonhole under local injectivity: --no-supplementals keeps NDS on, so the
    # "injectivity in the square" step sums the neighbourhood-injectivity (@linj)
    # constraints in place of the global ones (trident -> longtrident fires it for the
    # degree-2 vertex).
    gss_add_proof_test(proof_locally_injective_nds ${SUB} --no-supplementals
        --locally-injective --count-solutions ${GSS_INSTANCES}/trident.csv ${GSS_INSTANCES}/longtrident.csv)
    # Exact-path (distance-2) supplemental graphs under local injectivity: the
    # create_exact_path_graphs derivation sums @linj in place of the global injectivity,
    # both to drop the "q maps to t" term (site 1) and in the insufficient-paths pigeonhole
    # (site 2). This instance is dense enough to fire site 2 (a g>=2 exact-path graph with a
    # target vertex sharing too few common neighbours).
    gss_add_proof_test(proof_locally_injective_supplementals ${SUB} --no-nds --format csv
        --locally-injective --count-solutions
        ${GSS_INSTANCES}/li_exact_path_pattern.csv ${GSS_INSTANCES}/li_exact_path_target.csv)
    # Induced non-edge under local injectivity: two non-adjacent pattern vertices may both
    # map to the same loopless target (it is not adjacent to itself, so the non-edge is
    # preserved), which needs the target kept in the induced non-edge's permitted set. Two
    # edgeless vertices into two edgeless vertices -> all 4 mappings, including the two that
    # share a target.
    gss_add_proof_test(proof_locally_injective_induced ${SUB_MIN} --format lad
        --locally-injective --induced --count-solutions
        ${GSS_INSTANCES}/two_isolated.lad ${GSS_INSTANCES}/two_isolated_target.lad)

    # --- Preprocessing / supplemental-graph proofs (loopless) ---
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

    # Randomised proof sweep: generate random instances and check the proof verifies
    # on each, across decision/induced/counting combinations. Catches proof-logging
    # bugs that fixed instances miss; complements the in-process correctness oracle in
    # random_homomorphism_test.
    add_test(NAME proof_random_sweep
        COMMAND ${CMAKE_SOURCE_DIR}/test-instances/random_proof_sweep.bash
            $<TARGET_FILE:glasgow_subgraph_solver> ${VERIPB_EXECUTABLE}
            $<TARGET_FILE:create_random_graph> ${CMAKE_CURRENT_BINARY_DIR})

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
        gss_add_cake_test(cake_loop_decision loop_vertex.lad loop_target.lad --no-supplementals)
        gss_add_cake_test(cake_loop_count loop_vertex.lad loop_target.lad --no-supplementals --count-solutions)
        # loops together with supplemental graphs (issue #56): the loop-cancelled
        # adjacency derivations are checked against cake's own OPB encoding too.
        gss_add_cake_test(cake_loop_supplementals small large --no-nds)
        gss_add_cake_test(cake_loop_supplementals_count small large --no-nds --count-solutions)
    endif()
endif()
