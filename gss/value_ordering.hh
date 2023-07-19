/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_GUARD_VALUE_ORDERING_HH
#define GLASGOW_GUARD_VALUE_ORDERING_HH 1

enum class ValueOrdering
{
    None,
    Biased,
    Degree,
    AntiDegree,
    Random
};

#endif
