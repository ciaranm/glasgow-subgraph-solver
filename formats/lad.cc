/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "formats/lad.hh"
#include "formats/input_graph.hh"

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

auto read_lad(ifstream && infile, const string & filename) -> InputGraph
{
    InputGraph result(0);

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

auto read_labelled_lad(ifstream && infile, const string & filename, LabelsMap & labels_map) -> InputGraph
{
    InputGraph result(0);

    result.resize(read_word(infile));
    if (! infile)
        throw GraphFileError{ filename, "error reading size" };

    for (int r = 0 ; r < result.size() ; ++r) {
        int l = read_word(infile);
        if (! infile)
            throw GraphFileError{ filename, "error reading label" };

        int c_end = read_word(infile);
        if (! infile)
            throw GraphFileError{ filename, "error reading edges count" };

        result.set_vertex_label(r, labels_map.remap_vertex_label(l));

        for (int c = 0 ; c < c_end ; ++c) {
            int e = read_word(infile);

            if (e < 0 || e >= result.size())
                throw GraphFileError{ filename, "edge index out of bounds" };

            int l = read_word(infile);
            if (l < 0)
                throw GraphFileError{ filename, "edge label invalid" };

            result.add_edge(r, e, labels_map.remap_edge_label(l));
        }
    }

    string rest;
    if (infile >> rest)
        throw GraphFileError{ filename, "EOF not reached, next text is \"" + rest + "\"" };
    if (! infile.eof())
        throw GraphFileError{ filename, "EOF not reached" };

    return result;
}

