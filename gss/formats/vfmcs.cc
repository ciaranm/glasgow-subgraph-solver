#include <gss/formats/graph_file_error.hh>
#include <gss/formats/vfmcs.hh>

#include <fstream>
#include <string>

using std::ifstream;
using std::move;
using std::string;
using std::to_string;

namespace
{
    auto read_word(ifstream & infile) -> unsigned
    {
        unsigned char a, b;
        a = static_cast<unsigned char>(infile.get());
        b = static_cast<unsigned char>(infile.get());
        return unsigned(a) | (unsigned(b) << 8);
    }

    auto read_vfmcs(ifstream && infile, const string & filename, bool vertex_labels, bool directed) -> InputGraph
    {
        int size = read_word(infile);
        if (! infile)
            throw GraphFileError{filename, "error reading size", true};

        InputGraph result{size, vertex_labels, false};

        // to be like the CP 2011 labelling scheme...
        int m = result.size() * 33 / 100;
        int p = 1, k1 = 0, k2 = 0;
        while (p < m && k1 < 16) {
            p *= 2;
            k1 = k2;
            k2++;
        }

        for (int r = 0; r < result.size(); ++r) {
            unsigned l = read_word(infile) >> (16 - k1);
            if (vertex_labels)
                result.set_vertex_label(r, to_string(l));
        }

        if (! infile)
            throw GraphFileError{filename, "error reading attributes", true};

        for (int r = 0; r < result.size(); ++r) {
            int c_end = read_word(infile);
            if (! infile)
                throw GraphFileError{filename, "error reading edges count", true};

            for (int c = 0; c < c_end; ++c) {
                unsigned e = read_word(infile);

                if (e >= unsigned(result.size()))
                    throw GraphFileError{filename, "edge index " + to_string(e) + " out of bounds", true};

                if (directed)
                    result.add_directed_edge(r, e, "directed");
                else
                    result.add_edge(r, e);
                read_word(infile);
            }
        }

        infile.peek();
        if (! infile.eof())
            throw GraphFileError{filename, "EOF not reached", true};

        return result;
    }
}

auto read_unlabelled_undirected_vfmcs(ifstream && infile, const string & filename) -> InputGraph
{
    return read_vfmcs(move(infile), filename, false, false);
}

auto read_vertex_labelled_undirected_vfmcs(ifstream && infile, const string & filename) -> InputGraph
{
    return read_vfmcs(move(infile), filename, true, false);
}

auto read_vertex_labelled_directed_vfmcs(std::ifstream && infile, const std::string & filename) -> InputGraph
{
    return read_vfmcs(move(infile), filename, true, true);
}
