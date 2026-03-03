#include <gss/formats/flat.hh>
#include <gss/formats/input_graph.hh>

#include <fstream>
#include <map>
#include <string>

using std::istream;
using std::map;
using std::pair;
using std::stoi;
using std::string;
using std::to_string;

auto read_flat(istream && infile, const string & filename) -> InputGraph
{
    unsigned v;
    infile >> v;
    if (! infile)
        throw GraphFileError{filename, "error reading size", true};

    InputGraph result(v, false, false);

    for (int r = 0; r < result.size(); ++r)
        result.set_vertex_name(r, to_string(r));

    unsigned e;
    infile >> e;

    for (unsigned i = 0; i < e; ++i) {
        int f, t;
        if (! (infile >> f >> t))
            break;
        result.add_edge(f, t);
    }

    string rest;
    if (((infile >> rest) && (! rest.empty())) || ! infile.eof())
        throw GraphFileError{filename, "EOF not reached", true};

    return result;
}
