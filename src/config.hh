/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_CONFIG_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_CONFIG_HH 1

#ifdef __APPLE__
#  define STD_FS_IS_BOOST 1
#else
#  ifdef __GNUC__
#    if __GNUC_MAJOR < 8
#      define STD_FS_IS_EXPERIMENTAL 1
#    else
#      define STD_FS_IS_STD 1
#    endif
#  else
#    define STD_FS_IS_STD 1
#  endif
#endif

#endif
