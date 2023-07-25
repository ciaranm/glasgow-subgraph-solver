#include <gss/formats/graph_file_error.hh>
#include <gss/formats/input_graph.hh>

#include <fstream>

using std::string;

GraphFileError::GraphFileError(const string & filename, const string & message, bool x) noexcept :
    _what("Error reading graph file '" + filename + "': " + message),
    _exists(x)
{
}

GraphFileError::GraphFileError(const string & message) noexcept :
    _what("Error creating graph: " + message),
    _exists(true)
{
}

auto GraphFileError::what() const noexcept -> const char *
{
    return _what.c_str();
}

auto GraphFileError::file_at_least_existed() const noexcept -> bool
{
    return _exists;
}
