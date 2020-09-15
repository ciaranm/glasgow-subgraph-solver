/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "formats/bigraph.hh"
#include "formats/input_graph.hh"

#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <utility>

using std::ifstream;
using std::pair;
using std::string;
using std::stoi;
using std::to_string;
using std::vector;

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

auto read_target_bigraph(ifstream && infile, const string &) -> InputGraph
{
   
    InputGraph result{ 0, true, true, true };

    std::vector<std::string> labels;
    string t = read_str(infile);
    while(t.find('}') == -1){
        t = read_str(infile);
        labels.push_back(t.substr(0, t.find(':')));
    }

    int r = read_num(infile);
    int n = read_num(infile);
    int s = read_num(infile);
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
        result.set_vertex_label(i, labels.at(i - r));
        result.set_vertex_name(i, to_string(i - r));
    }

    for (int i = 0 ; i != (r + n) ; ++i)
        for (int j = r ; j != (n + s + r) ; ++j) {
            char x = read_char(infile);
            if (x == '1')
                result.add_directed_edge(i, j, "dir");
        }

    string h = read_str(infile);
    while (h == "({},") {
        pair<bool, vector<int> > he;
        he.second.resize(r + n + s);
        he.first = (read_str(infile) != "{},");

        read_char(infile);

        while(true){
            string e = read_str(infile);
            string c = read_str(infile);
            he.second[stoi(e.substr(1, e.find(',')-1)) + r] = stoi(c.substr(0, c.find(')')));

            if (c.find('}') != -1) break;
        }
        
        result.add_hyperedge(move(he));
        h = read_str(infile);
    }

    return result;
}

auto read_pattern_bigraph(ifstream && infile, const string &) -> InputGraph
{
    InputGraph result{ 0, true, true, true };

    std::vector<std::string> labels;
    string t = read_str(infile);
    while(t.find('}') == -1){
        t = read_str(infile);
        labels.push_back(t.substr(0, t.find(':')));
    }

    int r = read_num(infile);
    int n = read_num(infile);
    int s = read_num(infile);
    result.resize(n);

    for (int i = 0 ; i != n ; ++i)
        result.set_vertex_label(i, labels.at(i));

    for (int i = 0 ; i != (r + n) ; ++i)
        for (int j = 0 ; j != (n + s) ; ++j) {
            char x = read_char(infile);
            if (i >= r && j < n && x == '1') result.add_directed_edge(i - r, j, "dir");
            else if (i < r && j < n && x == '1') {
                result.set_child_of_root(j);
                result.add_pattern_root_edge(i, j);
            }
            else if (j >= n && i >= r && x == '1') {
                result.set_parent_of_site(i - r);
                result.add_pattern_site_edge(j - n, i - r);
            }
        }

    string h = read_str(infile);
    while (h == "({},") {
        pair<bool, vector<int> > he;
        he.second.resize(n);
        he.first = (read_str(infile) != "{},");

        read_char(infile);

        while(true){
            string e = read_str(infile);
            string c = read_str(infile);
            he.second[stoi(e.substr(1, e.find(',')-1))] = stoi(c.substr(0, c.find(')')));

            if (c.find('}') != -1) break;
        }
    
        result.add_hyperedge(move(he));
        h = read_str(infile);        
    }

    return result;
}

