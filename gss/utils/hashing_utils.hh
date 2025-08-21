#ifndef GLASGOW_SUBGRAPH_SOLVER_HASHING_UTILS_HH
#define GLASGOW_SUBGRAPH_SOLVER_HASHING_UTILS_HH

#include <functional>

// direct implementation of boost::hash_combine - it is the "best implementation" but this avoids the dependency
template <class T>
inline void hash_combine(std::size_t& seed, T const& v)
{
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

#endif // GLASGOW_SUBGRAPH_SOLVER_HASHING_UTILS_HH
