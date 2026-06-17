cmake_minimum_required(VERSION 3.21)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

project(
        Glasgow-Subgraph-Solver
        VERSION 1.0
        DESCRIPTION "Solvers for subgraph-finding problems like clique, subgraph isomorphism, and maximum common (connected) subgraph"
        HOMEPAGE_URL "https://github.com/ciaranm/glasgow-subgraph-solver"
        LANGUAGES CXX)

# Build type configuration.
# Three types are supported: Release (default), Debug, and Sanitize.
# See CMakePresets.json and the README for usage.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release Debug Sanitize)

# Per-type flags, applied automatically by CMake for the active build type on top
# of the always-on flags set below.
# Release: full optimisation. -march=native is added later (Release only) after the
#   compiler capability check.
set(CMAKE_CXX_FLAGS_RELEASE "-O3"
    CACHE STRING "Flags for release builds")
# Debug: no optimisation so the debugger can follow every step.
set(CMAKE_CXX_FLAGS_DEBUG "-O0"
    CACHE STRING "Flags for debug builds")
# Sanitize: light optimisation keeps stack traces readable while ASan + UBSan catch
#   memory errors and undefined behaviour.
set(CMAKE_CXX_FLAGS_SANITIZE "-O1 -fsanitize=address,undefined -fno-omit-frame-pointer"
    CACHE STRING "Flags for sanitizer builds (ASan + UBSan)")
set(CMAKE_EXE_LINKER_FLAGS_SANITIZE "-fsanitize=address,undefined"
    CACHE STRING "Linker flags for sanitizer builds")
set(CMAKE_SHARED_LINKER_FLAGS_SANITIZE "-fsanitize=address,undefined"
    CACHE STRING "Linker flags for sanitizer builds")

# Compiler flag checks always run so results are available below.
include(CheckCXXCompilerFlag)
unset(COMPILER_SUPPORTS_MARCH_NATIVE CACHE)
CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)

unset(COMPILER_SUPPORTS_NO_RESTRICT CACHE)
CHECK_CXX_COMPILER_FLAG(-Wno-restrict COMPILER_SUPPORTS_NO_RESTRICT)

unset(COMPILER_SUPPORTS_NO_STRINGOP_OVERREAD CACHE)
CHECK_CXX_COMPILER_FLAG(-Wno-stringop-overread COMPILER_SUPPORTS_NO_STRINGOP_OVERREAD)

# Global compile/link options only apply when building as the top-level project.
# When embedded as a subproject via FetchContent the parent controls these settings.
if(PROJECT_IS_TOP_LEVEL)
    add_compile_options(-W)
    add_compile_options(-Wall)

    if(COMPILER_SUPPORTS_NO_RESTRICT)
        add_compile_options(-Wno-restrict) # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=104336
    endif()
    if(COMPILER_SUPPORTS_NO_STRINGOP_OVERREAD)
        add_compile_options(-Wno-stringop-overread) # https://gcc.gnu.org/bugzilla/show_bug.cgi?id=104336
    endif()

    # Always include full debug info, even in release builds, so crash reports are useful.
    add_compile_options(-g)
    add_compile_options(-g3)
    add_compile_options(-gdwarf-3)
    add_compile_options(-pthread)
    add_link_options(-pthread)

    # Native CPU tuning is worthwhile for release performance but unhelpful for
    # debugging or for portable (e.g. sanitizer) builds.
    if(COMPILER_SUPPORTS_MARCH_NATIVE)
        add_compile_options($<$<CONFIG:Release>:-march=native>)
    endif()
endif()

enable_testing()

Include(FetchContent)
Set(FETCHCONTENT_QUIET FALSE)
FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG 6e79e682b726f524310d55dec8ddac4e9c52fb5f # v3.4.0
)
FetchContent_Declare(
        cxxopts
        GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
        GIT_TAG 44380e5a44706ab7347f400698c703eb2a196202 # v3.3.1
)

FetchContent_MakeAvailable(Catch2)
FetchContent_MakeAvailable(cxxopts)

SET(CATCH_CONFIG_ENABLE_ALL_STRINGMAKERS ON)

# VeriPB, if present, is used by the proof-verification tests (see gss/CMakeLists.txt).
# Tests that need it are only registered when it is found, so a checkout without
# VeriPB still builds and runs the rest of the suite.
find_program(VERIPB_EXECUTABLE veripb)

include_directories(.)
add_subdirectory(gss)
add_subdirectory(src)
