#ifndef GLASGOW_GUARD_VALUE_ORDERING_HH
#define GLASGOW_GUARD_VALUE_ORDERING_HH 1

namespace gss
{
    enum class ValueOrdering
    {
        None,
        Biased,
        Degree,
        AntiDegree,
        Random
    };
}

#endif
