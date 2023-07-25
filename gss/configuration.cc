#include <gss/configuration.hh>

using namespace gss;

using std::string;

UnsupportedConfiguration::UnsupportedConfiguration(const string & message) noexcept :
    _what(message)
{
}

auto UnsupportedConfiguration::what() const noexcept -> const char *
{
    return _what.c_str();
}
