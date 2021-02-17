/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "formats/bigraph.hh"
#include "formats/input_graph.hh"
#include <fstream>
#include <map>
#include <string>
#include <utility>
#include <iostream>

using std::ifstream;
using std::pair;
using std::string;
using std::stoi;
using std::to_string;
using std::vector;
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

auto read_target_bigraph(ifstream && infile, const string &) -> InputGraph
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
        he.second.resize(r + n + s);
        he.first = (read_str(infile) == "{},");

        read_char(infile);

        while (true) {
            string e = read_str(infile);
            string c = read_str(infile);
            he.second[stoi(e.substr(1, e.find(',') - 1)) + r] = stoi(c.substr(0, c.find(')')));

            if (c.find('}') != string::npos)
                break;
        }

        string port_id;
        if(he.first){
            port_id = "C" + to_string(closed_link_count);
            closed_link_count++;
        }
        else{ 
            port_id = "L" + to_string(open_link_count);    
            open_link_count++;
        }

        int ports_connected = 0;
        for(int i=r; i<(r+n); i++)
            for(int j=0;j<he.second[i];j++){
                result.add_link_node();
                result.add_directed_edge(i, result.size()-1, "dir");
                result.set_vertex_name(result.size()-1, port_id + "_" + to_string(ports_connected)); 
                ports_connected++;
            }

        for(int i=(result.size()-ports_connected);i<result.size();i++)
            for(int j=(result.size()-ports_connected);j<result.size();j++)
                if(i != j)
                    result.add_directed_edge(i, j, "dir");

        if(he.first){
            result.add_link_node();
            result.set_vertex_name(result.size()-1, "C_LINK_" + to_string(closed_link_count-1)); 
            for(int i=(result.size()-ports_connected);i<result.size();i++)
                result.add_directed_edge(i, result.size()-1, "dir");
        }                   

        h = read_str(infile);
              
        // Create x duplicates of each hyperedge (represented by a node labelled "LINK" which is adjacent to all nodes that the hyperedge
        // is connected to) where x is the number of pattern hyperedges, to allow many-to-one assignment
        //if(!he.first) {       
            //for(int x=0;x<pattern_graph.get_no_link_nodes();x++){     
            //    result.add_link_node();       
             //   result.set_vertex_name(result.size()-1, "C_LINK_" + to_string(closed_link_count)); 

            //   for(int i=r; i<(r+n); ++i){
            //        result.add_link_adjacency(i, result.size()-1, he.second[i]);
             //       for(int j=0;j<he.second[i];j++)
             //           result.add_directed_edge(i, result.size()-1, "dir");
             //   }
            //}
        //}
        //else {
           // for(int x=0;x<pattern_graph.get_no_link_nodes();x++){  
            //    result.add_link_node();
            //    result.set_vertex_name(result.size()-1, "O_LINK_" + to_string(open_link_count));
            //    for(int i=r; i<(r+n); ++i) {
            //        result.add_link_adjacency(i, result.size()-1, he.second[i]);
           //         for(int j=0;j<he.second[i];j++)
           //             result.add_directed_edge(i, result.size()-1, "dir"); 
           ///     }
           // }
        //}

    }

    // (Ciaran's idea) Naive method of connecting edges between all link nodes which are not clones of one another, i.e. share different names
    //for(int i=0;i<result.size();i++)
    //    for(int j=0;j<result.size();j++)
    //        if(result.vertex_label(i) == "LINK" && result.vertex_label(j) == "LINK" && result.vertex_name(i) != result.vertex_name(j) && i != j)
    //            result.add_directed_edge(i, j, "dir");
   
    return result;

}

auto read_pattern_bigraph(ifstream && infile, const string &) -> InputGraph
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
        he.second.resize(n);
        he.first = (read_str(infile) == "{},");

        read_char(infile);

        while (true) {
            string e = read_str(infile);
            string c = read_str(infile);
            he.second[stoi(e.substr(1, e.find(',') - 1))] = stoi(c.substr(0, c.find(')')));

            if (c.find('}') != string::npos)
                break;
        }

        string port_id;
        if(he.first){
            port_id = "C" + to_string(closed_link_count);
            closed_link_count++;
        }
        else{ 
            port_id = "L" + to_string(open_link_count);    
            open_link_count++;
        }

        int ports_connected = 0;
        for(int i=0; i<n; i++)
            for(int j=0;j<he.second[i];j++){
                result.add_link_node();                
                result.add_directed_edge(i, result.size()-1, "dir");                
                result.set_vertex_name(result.size()-1, port_id + "_" + to_string(ports_connected)); 
                ports_connected++;          
            }
                
        for(int i=(result.size()-ports_connected);i<result.size();i++)
            for(int j=(result.size()-ports_connected);j<result.size();j++)
                if(i != j)
                    result.add_directed_edge(i, j, "dir");
        
        if(he.first){
            result.add_link_node();
            result.set_vertex_name(result.size()-1, "C_LINK_" + to_string(closed_link_count-1)); 
            for(int i=(result.size()-ports_connected);i<result.size();i++)
                result.add_directed_edge(i, result.size()-1, "dir");
        }            

            
        h = read_str(infile);
            
        
        // Mark as either open or closed link in the vertex name for identifying in constraints later
        //if(!he.first)
        //    result.set_vertex_name(result.size()-1, "C_LINK_" + to_string(closed_link_count++)); 
       // else
       //     result.set_vertex_name(result.size()-1, "O_LINK_" + to_string(open_link_count++)); 

        // Represent the hyperedge by connecting all adjacencies to the link node
       // for(int i=0; i<n; ++i) {
        //    result.add_link_adjacency(i, result.size()-1, he.second[i]);
        //    for(int j=0;j<he.second[i];j++)
        //        result.add_directed_edge(i, result.size()-1, "dir");
        //}
    }

    return result;
}
