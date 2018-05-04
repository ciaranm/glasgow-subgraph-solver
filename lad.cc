/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "lad.hh"
#include "graph.hh"

#include <fstream>

using std::ifstream;
using std::string;

namespace
{
    auto read_word(ifstream & infile) -> int
    {
        int x;
        infile >> x;
        return x;
    }
}

auto read_lad(const string & filename) -> Graph
{
    Graph result(0);

    ifstream infile{ filename };
    if (! infile)
        throw GraphFileError{ filename, "unable to open file" };

    result.resize(read_word(infile));
    if (! infile)
        throw GraphFileError{ filename, "error reading size" };

    for (int r = 0 ; r < result.size() ; ++r) {
        int c_end = read_word(infile);
        if (! infile)
            throw GraphFileError{ filename, "error reading edges count" };

        for (int c = 0 ; c < c_end ; ++c) {
            int e = read_word(infile);

            if (e < 0 || e >= result.size())
                throw GraphFileError{ filename, "edge index out of bounds" };

            result.add_edge(r, e);
        }
    }

    string rest;
    if (infile >> rest)
        throw GraphFileError{ filename, "EOF not reached, next text is \"" + rest + "\"" };
    if (! infile.eof())
        throw GraphFileError{ filename, "EOF not reached" };

    return result;
}

