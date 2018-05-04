/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "graph_file_error.hh"
#include "graph.hh"

#include <fstream>

using std::ifstream;
using std::string;

GraphFileError::GraphFileError(const string & filename, const string & message) throw () :
    _what("Error reading graph file '" + filename + "': " + message)
{
}

auto GraphFileError::what() const throw () -> const char *
{
    return _what.c_str();
}


