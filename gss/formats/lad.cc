#include <gss/formats/input_graph.hh>
#include <gss/formats/lad.hh>

#include <fstream>
#include <map>
#include <string>

using std::ifstream;
using std::map;
using std::pair;
using std::stoi;
using std::string;
using std::to_string;

namespace
{
    auto read_word(ifstream & infile) -> int
    {
        int x;
        infile >> x;
        return x;
    }

    auto read_any_lad(ifstream && infile, const string & filename,
        bool directed,
        bool vertex_labels,
        bool edge_labels) -> InputGraph
    {
        InputGraph result{0, vertex_labels, edge_labels};

        result.resize(read_word(infile));
        if (! infile)
            throw GraphFileError{filename, "error reading size", true};

        for (int r = 0; r < result.size(); ++r) {
            result.set_vertex_name(r, to_string(r));

            if (vertex_labels) {
                int l = read_word(infile);
                if (! infile)
                    throw GraphFileError{filename, "error reading label", true};

                result.set_vertex_label(r, to_string(l));
            }

            int c_end = read_word(infile);
            if (! infile)
                throw GraphFileError{filename, "error reading edges count", true};

            for (int c = 0; c < c_end; ++c) {
                int e = read_word(infile);

                if (e < 0 || e >= result.size())
                    throw GraphFileError{filename, "edge index out of bounds", true};

                if (edge_labels) {
                    int l = read_word(infile);
                    if (l < 0)
                        throw GraphFileError{filename, "edge label invalid", true};

                    result.add_directed_edge(r, e, to_string(l));
                }
                else if (directed)
                    result.add_directed_edge(r, e, "");
                else
                    result.add_edge(r, e);
            }
        }

        string rest;
        while (infile >> rest) {
            auto equals_pos = rest.find('=');
            if (string::npos == equals_pos)
                throw GraphFileError{filename, "EOF not reached, next text is \"" + rest + "\"", true};
            else {
                string before = rest.substr(0, equals_pos), after = rest.substr(equals_pos + 1);
                result.set_vertex_name(stoi(before), after);
            }
        }
        if (! infile.eof())
            throw GraphFileError{filename, "EOF not reached", true};

        return result;
    }
}

auto read_lad(ifstream && infile, const string & filename) -> InputGraph
{
    return read_any_lad(move(infile), filename, false, false, false);
}

auto read_directed_lad(ifstream && infile, const string & filename) -> InputGraph
{
    return read_any_lad(move(infile), filename, true, false, false);
}

auto read_labelled_lad(ifstream && infile, const string & filename) -> InputGraph
{
    return read_any_lad(move(infile), filename, true, true, true);
}

auto read_vertex_labelled_lad(ifstream && infile, const string & filename) -> InputGraph
{
    return read_any_lad(move(infile), filename, false, true, false);
}
