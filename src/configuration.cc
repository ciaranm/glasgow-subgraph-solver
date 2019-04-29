/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "configuration.hh"

using std::string;

UnsupportedConfiguration::UnsupportedConfiguration(const string & message) noexcept :
    _what(message)
{
}

auto UnsupportedConfiguration::what() const noexcept -> const char *
{
    return _what.c_str();
}

