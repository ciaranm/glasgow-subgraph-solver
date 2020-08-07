/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "formats/bigraph.hh"
#include "formats/input_graph.hh"

#include <fstream>
#include <iostream>
#include <map>

using std::ifstream;
using std::string;
using std::to_string;
using std::cout;

namespace
{
    auto read_num(ifstream & infile) -> int
    {
        int x;
        infile >> x;
        return x;
    }

    auto read_str(ifstream & infile) -> string
    {
        string s;
        infile >> s;
        return s;
    }

    auto read_char(ifstream & infile) -> char
    {
        char c;
        infile >> c;
        return c;
    }
}

auto read_target_bigraph(ifstream && infile, const string & filename) -> InputGraph
{
    InputGraph result{ 0, true, true, true };

    int r = read_num(infile);
    int n = read_num(infile);
    int s = read_num(infile);
    int h = read_num(infile);
    result.resize(r + n + s);

    for (int i = 0 ; i != r ; ++i) {
        result.set_vertex_label(i, "ROOT");
        result.set_vertex_name(i, "ROOT" + to_string(i));
    }

    for (int i = (r + n) ; i < (r + n + s) ; ++i) {
        result.set_vertex_label(i, "SITE");
        result.set_vertex_name(i, "SITE" + to_string(i));
    }

    for (int i = r ; i < (r + n) ; ++i) {
        result.set_vertex_label(i, read_str(infile));
        result.set_vertex_name(i, to_string(i - r));
    }

    for (int i = 0 ; i != (r + n) ; ++i)
        for (int j = r ; j != (n + s + r) ; ++j) {
            char x = read_char(infile);
            if (x == '1')
                result.add_directed_edge(i, j, "dir");
        }

    for (int i=0; i<h; i++) {
        std::pair<bool, std::vector<int> > he;
        he.second.resize(r+n+s);
        string x = read_str(infile);

        while(x != "f" && x != "t") {
            int index = std::stoi(x);
            he.second[index-1+r]++;
            x = read_str(infile);
        }
        if(x == "t") he.first = true;
        else he.first = false;

        result.add_hyperedge(he);
    }

    return result;
}

auto read_pattern_bigraph(ifstream && infile, const string & filename) -> InputGraph
{
    InputGraph result{ 0, true, true, true };

    int r = read_num(infile);
    int n = read_num(infile);
    int s = read_num(infile);
    int h = read_num(infile);
    result.resize(n);

    for (int i=0; i<n; i++) result.set_vertex_label(i, read_str(infile));

    for (int i=0; i<(r+n); i++) for(int j=0; j<(n+s); j++) {
        char x = read_char(infile);
        if (i >= r && j < n && x == '1') result.add_directed_edge(i-r, j, "dir");
        if (i < r && j < n && x == '1') {
            result.set_child_of_root(j);
            result.add_pattern_root_edge(i, j);
        }
        if (j >= n && i >= r && x == '1') {
            result.set_parent_of_site(i-r);
            result.add_pattern_site_edge(j-n, i-r);
        }
    }

    for (int i=0; i<h; i++) {
        std::pair<bool, std::vector<int> > he;
        he.second.resize(n);
        string x = read_str(infile);

        while(x != "f" && x != "t") {
            int index = std::stoi(x);
            he.second[index-1]++;
            x = read_str(infile);
        }
        if(x == "t") he.first = true;
        else he.first = false;

        result.add_hyperedge(he);
    }

    return result;
}

// GLasgowBigraphSolver.cc -> Rewrite stuff
// Unary constraints to check in-degree and out-degree in the place graph to deal with that one problem Michele had
// Consider sites, roots can be ignored for now?
// Pass in a function


