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
# Four types are supported: Release (default), Debug, Sanitize, and Coverage.
# See CMakePresets.json and the README for usage.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Release Debug Sanitize Coverage)

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
# Coverage: no optimisation + gcov instrumentation, for measuring which proof-writing
#   branches the test suite actually triggers (see the `coverage` target below).
#
# FORCE is required: CMake pre-creates CMAKE_<LANG>_FLAGS_<CONFIG> as empty cache
# entries for the active build type during compiler detection (which runs in the
# project() call above), so a plain `set(... CACHE ...)` would be a silent no-op and
# the build type would get none of its flags.
set(CMAKE_CXX_FLAGS_SANITIZE "-O1 -fsanitize=address,undefined -fno-omit-frame-pointer"
    CACHE STRING "Flags for sanitizer builds (ASan + UBSan)" FORCE)
set(CMAKE_EXE_LINKER_FLAGS_SANITIZE "-fsanitize=address,undefined"
    CACHE STRING "Linker flags for sanitizer builds" FORCE)
set(CMAKE_SHARED_LINKER_FLAGS_SANITIZE "-fsanitize=address,undefined"
    CACHE STRING "Linker flags for sanitizer builds" FORCE)
set(CMAKE_CXX_FLAGS_COVERAGE "-O0 --coverage -fprofile-update=atomic"
    CACHE STRING "Flags for coverage builds" FORCE)
set(CMAKE_EXE_LINKER_FLAGS_COVERAGE "--coverage"
    CACHE STRING "Linker flags for coverage builds" FORCE)
set(CMAKE_SHARED_LINKER_FLAGS_COVERAGE "--coverage"
    CACHE STRING "Linker flags for coverage builds" FORCE)

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

# Prefer a prebuilt Catch2 (e.g. from Homebrew or apt) so it isn't compiled from
# source on every build; fall back to fetching a pinned version otherwise. Catch2
# v3 is a compiled library, so this is a worthwhile build-time saving (especially
# on CI). Any v3 is fine -- the test macros are stable across minor versions.
find_package(Catch2 3 QUIET)
if(NOT Catch2_FOUND)
    FetchContent_Declare(
            Catch2
            GIT_REPOSITORY https://github.com/catchorg/Catch2.git
            GIT_TAG 6e79e682b726f524310d55dec8ddac4e9c52fb5f # v3.4.0
    )
    FetchContent_MakeAvailable(Catch2)
endif()

FetchContent_Declare(
        cxxopts
        GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
        GIT_TAG 44380e5a44706ab7347f400698c703eb2a196202 # v3.3.1
)
FetchContent_MakeAvailable(cxxopts)

SET(CATCH_CONFIG_ENABLE_ALL_STRINGMAKERS ON)

# VeriPB, if present, is used by the proof-verification tests (see gss/CMakeLists.txt).
# Tests that need it are only registered when it is found, so a checkout without
# VeriPB still builds and runs the rest of the suite.
find_program(VERIPB_EXECUTABLE veripb)

# cake_pb_iso (the formally verified CakePB checker for subgraph isomorphism), if
# present, is used by the end-to-end verified-pipeline tests: VeriPB elaborates the
# proof against cake's OPB, then cake_pb_iso checks the elaborated core proof. It is
# not on the PATH by default, so point CMake at it with
#   -DCAKE_PB_ISO_EXECUTABLE=/path/to/cake_pb_iso
# when you want these tests; otherwise they are silently skipped.
find_program(CAKE_PB_ISO_EXECUTABLE cake_pb_iso)

include_directories(.)
add_subdirectory(gss)
add_subdirectory(src)

# Coverage reporting. With a Coverage build, `cmake --build <dir> --target coverage`
# runs the test suite and reports branch coverage of the proof-writing code, so we
# can see which conditional proof features no test triggers. Requires gcovr.
if(CMAKE_BUILD_TYPE STREQUAL "Coverage")
    find_program(GCOVR_EXECUTABLE gcovr)
    if(GCOVR_EXECUTABLE)
        add_custom_target(coverage
            COMMAND ${CMAKE_COMMAND} -E echo "Running tests under coverage instrumentation..."
            COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure || ${CMAKE_COMMAND} -E true
            COMMAND ${GCOVR_EXECUTABLE} --root ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR}
                    --filter "${CMAKE_SOURCE_DIR}/gss/innards/proof.cc"
                    --exclude ".*/_deps/.*" --gcov-ignore-errors=no_working_dir_found
                    --branches --decisions --txt --print-summary
                    --txt "${CMAKE_BINARY_DIR}/proof-coverage.txt"
                    --html-details "${CMAKE_BINARY_DIR}/proof-coverage.html"
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            COMMENT "Generating proof.cc branch-coverage report"
            USES_TERMINAL VERBATIM)
    else()
        message(WARNING "gcovr not found: the `coverage` target will be unavailable")
    endif()
endif()
