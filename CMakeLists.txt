cmake_minimum_required(VERSION 3.18)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

project(
        Glasgow-Subgraph-Solver
        VERSION 1.0
        DESCRIPTION "Solvers for subgraph-finding problems like clique, subgraph isomorphism, and maximum common (connected) subgraph"
        HOMEPAGE_URL "https://github.com/ciaranm/glasgow-subgraph-solver"
        LANGUAGES CXX)

add_compile_options(-W)
add_compile_options(-Wall)

include(CheckCXXCompilerFlag)
unset(COMPILER_SUPPORTS_MARCH_NATIVE CACHE)
CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
if (COMPILER_SUPPORTS_MARCH_NATIVE)
    add_compile_options(-march=native)
endif (COMPILER_SUPPORTS_MARCH_NATIVE)

unset(COMPILER_SUPPORTS_NO_RESTRICT CACHE)
CHECK_CXX_COMPILER_FLAG(-Wno-restrict COMPILER_SUPPORTS_NO_RESTRICT)
if (COMPILER_SUPPORTS_NO_RESTRICT)
    add_compile_options(-Wno-restrict) # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=104336
endif (COMPILER_SUPPORTS_NO_RESTRICT)

unset(COMPILER_SUPPORTS_NO_STRINGOP_OVERREAD CACHE)
CHECK_CXX_COMPILER_FLAG(-Wno-stringop-overread COMPILER_SUPPORTS_NO_STRINGOP_OVERREAD)
if (COMPILER_SUPPORTS_NO_STRINGOP_OVERREAD)
    add_compile_options(-Wno-stringop-overread) # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=104336
endif (COMPILER_SUPPORTS_NO_STRINGOP_OVERREAD)

add_compile_options(-g)
add_compile_options(-g3)
add_compile_options(-gdwarf-3)
add_compile_options(-pthread)

add_link_options(-pthread)

if (NOT GCS_DEBUG_MODE)
    add_compile_options(-O3)
endif (NOT GCS_DEBUG_MODE)

Include(FetchContent)
Set(FETCHCONTENT_QUIET FALSE)
FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG 6e79e682b726f524310d55dec8ddac4e9c52fb5f # v3.4.0
)
FetchContent_MakeAvailable(Catch2)

SET(CATCH_CONFIG_ENABLE_ALL_STRINGMAKERS ON)

enable_testing()

include_directories(.)
add_subdirectory(gss)
add_subdirectory(src)
