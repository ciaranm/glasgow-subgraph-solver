set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

link_libraries(glasgow_subgraphs)

find_package(Boost 1.74 COMPONENTS program_options container REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

add_executable(glasgow_subgraph_solver glasgow_subgraph_solver.cc)
target_link_libraries(glasgow_subgraph_solver LINK_PUBLIC ${Boost_LIBRARIES})

add_executable(glasgow_clique_solver glasgow_clique_solver.cc)
target_link_libraries(glasgow_clique_solver LINK_PUBLIC ${Boost_LIBRARIES})

add_executable(glasgow_common_subgraph_solver glasgow_common_subgraph_solver.cc)
target_link_libraries(glasgow_common_subgraph_solver LINK_PUBLIC ${Boost_LIBRARIES})

add_executable(create_random_graph create_random_graph.cc)
target_link_libraries(create_random_graph LINK_PUBLIC ${Boost_LIBRARIES})

add_executable(convert_to_lad convert_to_lad.cc)
target_link_libraries(convert_to_lad LINK_PUBLIC ${Boost_LIBRARIES})
