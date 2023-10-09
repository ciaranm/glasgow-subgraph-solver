/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "formats/bigraph.hh"
#include "formats/input_graph.hh"
#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <iostream>

using std::istream;
using std::pair;
using std::string;
using std::stoi;
using std::to_string;
using std::vector;
using std::cout;

namespace
{
    auto read_num(istream & infile) -> int
    {
        int x;
        infile >> x;
        return x;
    }

    auto read_str(istream & infile) -> string
    {
        string s;
        infile >> s;
        return s;
    }

    auto read_char(istream & infile) -> char
    {
        char c;
        infile >> c;
        return c;
    }
}

auto read_target_bigraph(istream && infile, const string &, bool directed) -> InputGraph
{
    InputGraph result{ 0, true, true, true };

    std::vector<std::string> labels;
    string t = read_str(infile);
    while (t.find('}') == string::npos) {
        t = read_str(infile);
        labels.push_back(t.substr(0, t.find(':')));
    }

    int r = read_num(infile);
    int n = read_num(infile);
    int s = read_num(infile);
    result.resize(r + n + s);

    // Add dummy root nodes
    for (int i = 0 ; i != r ; ++i) {
        result.set_vertex_label(i, "ROOT");
        result.set_vertex_name(i, "ROOT" + to_string(i));
    }

    // Add dummy site nodes
    for (int i = (r + n) ; i < (r + n + s) ; ++i) {
        result.set_vertex_label(i, "SITE");
        result.set_vertex_name(i, "SITE" + to_string(i));
    }

    // Add standard place graph nodes
    for (int i = r ; i < (r + n) ; ++i) {
        result.set_vertex_label(i, labels.at(i - r));
        result.set_vertex_name(i, to_string(i - r));
    }

    // Read adjacency matrix and add directed edges
    for (int i = 0 ; i != (r + n) ; ++i)
        for (int j = r ; j != (n + s + r) ; ++j) {
            char x = read_char(infile);
            if (x == '1')
                result.add_directed_edge(i, j, "dir");
        }


    // Add link graph nodes (this is where all the scary stuff happens...)
    string h = read_str(infile);
    int closed_link_count = 0;
    int open_link_count = 0;
    while (h == "({},") {
        pair<bool, vector<int> > he;
        vector<int> outgoing_he;
        he.second.resize(r + n + s);
        if(directed)
            outgoing_he.resize(r + n + s);

        string link_name = read_str(infile);

        he.first = (link_name == "{},");
        string open_link_dir;
        if (!he.first) {
            if(directed) {
                link_name = link_name.substr(2, link_name.length()-3);
                open_link_dir = (read_char(infile) == 'i') ? "i" : "o";
                read_str(infile);
            }
            else
                link_name = link_name.substr(1, link_name.length()-3);
        }
        read_char(infile);

        // Deal with hanging link case
        string e = read_str(infile);
        if (e == "})") {
            h = read_str(infile);
            continue;           
        }
        
        if(directed) {
            string ci = read_str(infile);
            string co = read_str(infile);

            he.second[stoi(e.substr(1, e.find(',') - 1)) + r] = stoi(ci.substr(0, ci.length()-1));
            outgoing_he[stoi(e.substr(1, e.find(',') - 1)) + r] = stoi(co.substr(0, co.find(')')));
            while (co.find('}') == string::npos) {
                e = read_str(infile);
                ci = read_str(infile);
                co = read_str(infile);
                he.second[stoi(e.substr(1, e.find(',') - 1)) + r] = stoi(ci.substr(0, ci.length()-1));
                outgoing_he[stoi(e.substr(1, e.find(',') - 1)) + r] = stoi(co.substr(0, co.find(')')));
            }
        }
        else {
            string c = read_str(infile);
            he.second[stoi(e.substr(1, e.find(',') - 1)) + r] = stoi(c.substr(0, c.find(')')));
            while (c.find('}') == string::npos) {
                e = read_str(infile);
                c = read_str(infile);
                he.second[stoi(e.substr(1, e.find(',') - 1)) + r] = stoi(c.substr(0, c.find(')')));     
            }
        }

        string port_id;
        if(he.first){
            port_id = ":CLX:" + to_string(closed_link_count);
            closed_link_count++;
        }
        else{ 
            port_id = ":OPX:" + link_name;
            open_link_count++;
        }

        int ports_connected = 0;
        for(int i=r; i<(r+n); i++) {
            for(int j=0;j<he.second[i];j++){
                result.add_link_node();
                result.set_vertex_label(result.size()-1, "LINK");
                result.add_directed_edge(i, result.size()-1, "dir");
                if(!directed || he.first)
                    result.set_vertex_name(result.size()-1, port_id + ":" + to_string(ports_connected)); 
                else
                    result.set_vertex_name(result.size()-1, port_id + ":" + to_string(ports_connected) + ":" + open_link_dir); 
                ports_connected++;
            }
            if(directed)
                for(int j=0;j<outgoing_he[i];j++){
                    result.add_link_node();
                    result.set_vertex_label(result.size()-1, "LINK");
                    result.add_directed_edge(result.size()-1, i, "dir");
                    if(he.first)
                        result.set_vertex_name(result.size()-1, port_id + ":" + to_string(ports_connected)); 
                    else
                        result.set_vertex_name(result.size()-1, port_id + ":" + to_string(ports_connected) + ":" + open_link_dir); 
                    ports_connected++;
                }
        }

        for(int i=(result.size()-ports_connected);i<result.size();i++)
            for(int j=(result.size()-ports_connected);j<result.size();j++)
                if(i != j)
                    result.add_directed_edge(i, j, "dir");

        if(he.first){
            result.add_link_node();
            result.set_vertex_name(result.size()-1, "C_LINK_" + to_string(closed_link_count-1)); 
            result.set_vertex_label(result.size()-1, "ANCHOR");
            for(int i=(result.size()-ports_connected-1);i<result.size()-1;i++)
                result.add_directed_edge(i, result.size()-1, "dir");           
        }                   

        h = read_str(infile);
    }
   
    return result;
}

auto read_pattern_bigraph(istream && infile, const string &, bool directed) -> InputGraph
{
    InputGraph result{ 0, true, true, true };

    vector<string> labels;
    string t = read_str(infile);
    while (t.find('}') == string::npos) {
        t = read_str(infile);
        labels.push_back(t.substr(0, t.find(':')));
    }

    int r = read_num(infile);
    int n = read_num(infile);
    int s = read_num(infile);


    // Add place graph nodes
    
    result.resize(n);
    for (int i = 0 ; i != n ; ++i)
        result.set_vertex_label(i, labels.at(i));

    // Add adjacency matrix, and identify children of root nodes and parents of site nodes for later constraints
    for (int i = 0 ; i != (r + n) ; ++i)
        for (int j = 0 ; j != (n + s) ; ++j) {
            char x = read_char(infile);
            if (i >= r && j < n && x == '1')
                result.add_directed_edge(i - r, j, "dir");
            else if (i < r && j < n && x == '1') {
                result.set_child_of_root(j);
                result.add_pattern_root_edge(i, j);
            }
            else if (j >= n && i >= r && x == '1') {
                result.set_parent_of_site(i - r);
                result.add_pattern_site_edge(j - n, i - r);
            }
        }

    // Add link graph nodes
    string h = read_str(infile);
    int closed_link_count = 0;
    int open_link_count = 0;
    while (h == "({},") {
        pair<bool, vector<int> > he;
        vector<int> outgoing_he;
        he.second.resize(n);
        if(directed)
            outgoing_he.resize(n);

        string link_name = read_str(infile);
        he.first = (link_name == "{},");
        string open_link_dir;
        if (!he.first) {
            if(directed) {
                link_name = link_name.substr(2, link_name.length()-3);
                open_link_dir = (read_char(infile) == 'i') ? "i" : "o";
                read_str(infile);
            }
            else
                link_name = link_name.substr(1, link_name.length()-3);
        }
        read_char(infile);

        // Deal with hanging link case
        string e = read_str(infile);
        if (e == "})") {
            h = read_str(infile);
            continue;           
        }

        if(directed) {
            string ci = read_str(infile);
            string co = read_str(infile);
            he.second[stoi(e.substr(1, e.find(',') - 1))] = stoi(ci.substr(0, ci.length()-1));
            outgoing_he[stoi(e.substr(1, e.find(',') - 1))] = stoi(co.substr(0, co.find(')')));
            while (co.find('}') == string::npos) {
                e = read_str(infile);
                ci = read_str(infile);
                co = read_str(infile);
                he.second[stoi(e.substr(1, e.find(',') - 1))] = stoi(ci.substr(0, ci.length()-1));
                outgoing_he[stoi(e.substr(1, e.find(',') - 1))] = stoi(co.substr(0, co.find(')')));
            }
        }
        else {
            string c = read_str(infile);
            he.second[stoi(e.substr(1, e.find(',') - 1))] = stoi(c.substr(0, c.find(')')));
            while (c.find('}') == string::npos) {
                e = read_str(infile);
                c = read_str(infile);
                he.second[stoi(e.substr(1, e.find(',') - 1)) ] = stoi(c.substr(0, c.find(')')));     
            }
        }

        string port_id;
        if(he.first){
            port_id = ":CLX:" + to_string(closed_link_count);
            closed_link_count++;
        }
        else{ 
            port_id = ":OPX:" + link_name;
            open_link_count++;
        }

        int ports_connected = 0;
        for(int i=0; i<n; i++) {
            for(int j=0;j<he.second[i];j++){
                result.add_link_node();                
                result.set_vertex_label(result.size()-1, "LINK");
                result.add_directed_edge(i, result.size()-1, "dir");
                if(!directed || he.first)
                    result.set_vertex_name(result.size()-1, port_id + ":" + to_string(ports_connected)); 
                else
                    result.set_vertex_name(result.size()-1, port_id + ":" + to_string(ports_connected) + ":" + open_link_dir); 
                ports_connected++;          
            }
            if(directed)
                for(int j=0;j<outgoing_he[i];j++){
                    result.add_link_node();                
                    result.set_vertex_label(result.size()-1, "LINK");
                    result.add_directed_edge(result.size()-1, i, "dir");
                    if(he.first)
                        result.set_vertex_name(result.size()-1, port_id + ":" + to_string(ports_connected)); 
                    else
                        result.set_vertex_name(result.size()-1, port_id + ":" + to_string(ports_connected) + ":" + open_link_dir); 
                    ports_connected++;          
                }
        }
                
        for(int i=(result.size()-ports_connected);i<result.size();i++)
            for(int j=(result.size()-ports_connected);j<result.size();j++)
                if(i != j)
                    result.add_directed_edge(i, j, "dir");
        
        if(he.first){
            result.add_link_node();            
            result.set_vertex_name(result.size()-1, "C_LINK_" + to_string(closed_link_count-1)); 
            result.set_vertex_label(result.size()-1, "ANCHOR");
            for(int i=(result.size()-ports_connected-1);i<result.size()-1;i++)
                result.add_directed_edge(i, result.size()-1, "dir");
        }            

            
        h = read_str(infile);
    }
    
    return result;
}

