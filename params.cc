/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "params.hh"

using std::string;

UnsupportedConfiguration::UnsupportedConfiguration(const string & message) throw () :
    _what(message)
{
}

auto UnsupportedConfiguration::what() const throw () -> const char *
{
    return _what.c_str();
}

