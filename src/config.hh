/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_CONFIG_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_SRC_CONFIG_HH 1

#ifdef __GNUC__
#  if __GNUC_MAJOR < 8
#    define STD_FS_IS_EXPERIMENTAL 1
#  else
#    undef STD_FS_IS_EXPERIMENTAL
#  endif
#else
#  undef STD_FS_IS_EXPERIMENTAL
#endif

#endif
