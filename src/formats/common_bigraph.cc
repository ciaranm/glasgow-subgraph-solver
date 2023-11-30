/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <algorithm>
#include "homomorphism.hh"
#include "formats/common_bigraph.hh"
#include "formats/input_graph.hh"
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <set>

using std::istream;
using std::string;
using std::cout;
using std::to_string;
using std::set;

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

Entity::Entity(int i, string ctrl, int ar){
    id = i;
    control = ctrl;
    arity = ar;
    is_leaf = false;
}
Entity::Entity(){}

bool operator==(const Entity &e1, const Entity &e2){
    return e1.id == e2.id;
}

auto Entity::copy() -> Entity
{
    Entity e = Entity(id, control, arity);
    e.sites.insert(sites.begin(), sites.end());
    e.regions.insert(regions.begin(), regions.end());
    return e;
}

Bigraph::Bigraph(){}

auto Bigraph::copy() -> Bigraph
{
    Bigraph big = Bigraph();
    big.reachability = reachability;

    big.sites.insert(sites.begin(), sites.end());
    big.regions.insert(regions.begin(), regions.end());
    big.original_size = original_size;

    for(unsigned int i=0;i<entities.size();i++){
        big.entities.push_back(entities[i].copy());
    }
    return big;
}

auto Bigraph::toString() const -> string
{
    string out = "{";
    for(unsigned int i=0;i<entities.size();i++) {
        out += '(' + to_string(entities[i].id) + ", " + entities[i].control + ':' + to_string(entities[i].arity) + ')';
        if(i < entities.size()-1) out += ',';
    }
    out += "}\n" + to_string(regions.size()) + ' ' + to_string(entities.size()) + ' ' + to_string(sites.size()) + '\n';
    
    std::vector<std::vector<bool> > matrix(
        regions.size() + entities.size(),
        std::vector<bool>(sites.size() + entities.size()));
        
    
    for(unsigned int i=0;i<entities.size();i++) {
        Entity e = entities[i];
        if(e.parent_index != -1){
            auto it = find(entities.begin(), entities.end(), entities[e.parent_index]);
            int k = it - entities.begin();
            matrix[regions.size()+k][i] = 1;
        }
    }

    int index = 0;
    for(auto r : regions) {
        for(unsigned int i=0;i<entities.size();i++) {
            if(find(entities[i].regions.begin(), entities[i].regions.end(), r) != entities[i].regions.end())
                matrix[index][i] = 1;
        }
        index++;
    }
    index = 0;
    for(auto s : sites) {
        for(unsigned int i=0;i<entities.size();i++) {
            if(find(entities[i].sites.begin(), entities[i].sites.end(), s) != entities[i].sites.end())
                matrix[regions.size()+i][entities.size()+index] = 1;
        }
        index++;
    }

    for(unsigned int i=0;i<matrix.size();i++) {
        for(unsigned int j=0;j<matrix[i].size();j++)
            out += to_string(matrix[i][j]);
        out += '\n';
    }

    for(unsigned int i=0;i<hyperedges.size();i++) {
        bool lazy_flag = false;
        out += "({}, {" + hyperedges[i].first + "}, {";
        for(unsigned int j=0;j<hyperedges[i].second.size();j++) {
            if(hyperedges[i].second[j] > 0) {
                if(lazy_flag)
                    out += ", ";
                out += "(" + to_string(j) + ", " + to_string(hyperedges[i].second[j]) + ")";
                lazy_flag = true;
            }
        } 
        out += "})\n";
    }
    for(unsigned int i=0;i<closures.size();i++) {
        bool lazy_flag = false;
        out += "({}, {}, {";
        for(unsigned int j=0;j<closures[i].adjacencies.size();j++) {
            if(closures[i].adjacencies[j] > 0) {
                if(lazy_flag)
                    out += ", ";
                out += "(" + to_string(j) + ", " + to_string(closures[i].adjacencies[j]) + ")";
                lazy_flag = true;
            }
        }
        out += "})\n";
    }

    return out;
}

auto Bigraph::encode(bool target, bool special_lts_case) const -> InputGraph
{
    InputGraph result{ 0, true, true, true };

    if(target) result.resize(entities.size() + regions.size() + sites.size());
    else if (special_lts_case) result.resize(entities.size() + regions.size());
    else result.resize(entities.size());

    std::vector<int> id_map;
    id_map.resize(original_size);

    if(target || special_lts_case) {
        for(unsigned int i=0; i<regions.size();i++) {
            result.set_vertex_label(i, "ROOT");
            result.set_vertex_name(i, "ROOT" + to_string(i));
        }
        for(unsigned int i=regions.size(); i<regions.size()+entities.size();i++) {
            result.set_vertex_label(i, entities[i-regions.size()].control);
            result.set_vertex_name(i, to_string(i-regions.size()));
            id_map[entities[i-regions.size()].id] = i-regions.size()+1;
        }
        if(target)
            for(unsigned int i=regions.size()+entities.size(); i<regions.size()+entities.size()+sites.size();i++) {
                result.set_vertex_label(i, "SITE");
                result.set_vertex_name(i, "SITE" + to_string(i));
            }
    }
    else {
        for(unsigned int i=0; i<entities.size();i++) {
            result.set_vertex_label(i, entities[i].control);
            id_map[entities[i].id] = i+1;        
        }
    }
    
    int index = 0;    
    for(auto r : regions) {
        for(unsigned int i=0;i<entities.size();i++) {
            if(find(entities[i].regions.begin(), entities[i].regions.end(), r) != entities[i].regions.end()) {
                if(target || special_lts_case)
                    result.add_directed_edge(index, regions.size() + i, "dir");
                else {
                    result.set_child_of_root(i);
                    result.add_pattern_root_edge(index, i);
                }
            }
        }
        index++;
    }
    
    if(special_lts_case) {
        for(unsigned int i=0;i<entities.size();i++) {
            if(entities[i].parent_index == -1 && entities[i].regions.size() == 0) {
                result.set_child_of_root(i+regions.size());
                result.add_pattern_root_edge(index, i+regions.size());
                index++;
            }
        }
    }

    index = 0;
    for(auto s : sites) {
        for(unsigned int i=0;i<entities.size();i++) {
            if(find(entities[i].sites.begin(), entities[i].sites.end(), s) != entities[i].sites.end()) {
                if(target)
                    result.add_directed_edge(regions.size() + i, regions.size() + entities.size() + index, "dir");
                else if(special_lts_case) {
                    result.set_parent_of_site(regions.size() + i);
                    result.add_pattern_site_edge(index, regions.size() + i);                    
                }
                else {
                    result.set_parent_of_site(i);
                    result.add_pattern_site_edge(index, i);
                }
            }
        }
        index++;
    }
    for(unsigned int i=0;i<entities.size();i++) {
        Entity e = entities[i];
        if(e.parent_index != -1){
            auto it = find(entities.begin(), entities.end(), entities[e.parent_index]);
            int k = it - entities.begin();
            if(target || special_lts_case)
                result.add_directed_edge(regions.size() + k, regions.size() + i, "dir");
            else
                result.add_directed_edge(k, i, "dir");
        }
    }

    // Encode hyperedges
    for(unsigned int i=0;i<hyperedges.size();i++) {
        int ports_connected = 0;
        string port_id = ":OPX:" + hyperedges[i].first;
        for(unsigned int j=0;j<hyperedges[i].second.size();j++){
            for(unsigned int k=0;k<hyperedges[i].second[j];k++){
                result.add_link_node();                
                result.set_vertex_label(result.size()-1, "LINK");
                if(target || special_lts_case)
                    result.add_directed_edge(regions.size()+id_map[j]-1, result.size()-1, "dir");     
                else
                    result.add_directed_edge(id_map[j]-1, result.size()-1, "dir");     
                result.set_vertex_name(result.size()-1, port_id + ":" + to_string(ports_connected));     
                ports_connected++;       
            }
        }
        for(int j=(result.size()-ports_connected);j<result.size();j++)
            for(int k=(result.size()-ports_connected);k<result.size();k++)
                if(j != k)
                    result.add_directed_edge(j, k, "dir");        
    }

    // Encode closures
    for(unsigned int i=0;i<closures.size();i++) {
        int ports_connected = 0;
        string port_id = ":CLX:" + to_string(i);
        for(unsigned int j=0;j<closures[i].adjacencies.size();j++){
            for(unsigned int k=0;k<closures[i].adjacencies[j];k++){
                result.add_link_node();                
                result.set_vertex_label(result.size()-1, "LINK");
                if(target || special_lts_case)
                    result.add_directed_edge(regions.size()+id_map[j]-1, result.size()-1, "dir");     
                else {
                    result.add_directed_edge(id_map[j]-1, result.size()-1, "dir");}
                result.set_vertex_name(result.size()-1, port_id + ":" + to_string(ports_connected)); 
                ports_connected++;           
            }
        }
        for(int j=(result.size()-ports_connected);j<result.size();j++)
            for(int k=(result.size()-ports_connected);k<result.size();k++)
                if(j != k)
                    result.add_directed_edge(j, k, "dir");  

        result.add_link_node();            
        result.set_vertex_name(result.size()-1, "C_LINK_" + to_string(i)); 
        result.set_vertex_label(result.size()-1, "ANCHOR");
        for(int j=(result.size()-ports_connected-1);j<result.size()-1;j++)
            result.add_directed_edge(j, result.size()-1, "dir");
    }

    return result;
}

auto split_regions(Bigraph a) -> Bigraph
{
    if(a.regions.size() > 0) {
        int max_region = *a.regions.rbegin();
        for(int i=0;i<a.entities.size();i++) {
            if(a.entities[i].regions.size() > 0) {
                for(int j=i+1;j<a.entities.size();j++) {
                    if(a.entities[j].regions.size() > 0 && *(a.entities[i].regions.begin()) == *(a.entities[j].regions.begin())) {
                        max_region++;
                        a.entities[j].regions.clear();
                        a.entities[j].regions.insert(max_region);
                        a.regions.insert(max_region);
                    }
                }
            }
        }
    }
    return a;
}

auto free_hyperedges(Bigraph a) -> Bigraph
{
    std::vector<std::pair<string, std::vector<int>>> new_he_set;
    for(int i=0;i<a.hyperedges.size();i++) {
        for(int j=0;j<a.hyperedges[i].second.size();j++){
            for(int k=0;k<a.hyperedges[i].second[j];k++){
                std::pair<string, std::vector<int>> he;
                he.first = a.hyperedges[i].first;
                he.second.resize(a.hyperedges[i].second.size());
                he.second[j] = 1;
                new_he_set.push_back(he);
            }
        }
    }
    a.hyperedges = new_he_set;
    return a;
}

auto make_place_RPO(Bigraph big1, Bigraph big2, Bigraph solution, std::vector<int> mapping) -> Bigraph
{

    // close all shared sites
    for(int i=0;i<solution.entities.size();i++) {
        if(solution.entities[i].sites.size() > 1) {
            int start = *solution.entities[i].sites.begin();
            for(int s : solution.entities[i].sites){
                if(s != start) {
                    solution.sites.erase(s); 
                }
            }
            solution.entities[i].sites.erase(++solution.entities[i].sites.begin(), solution.entities[i].sites.end());
        }
    }

    // open all necessary regions
    int max_region;
    if(solution.regions.size() > 0)
        max_region = *solution.regions.rbegin();
    else max_region = 0;
    for(int i=0;i<solution.entities.size();i++) {
        if(solution.entities[i].regions.size() == 0 && solution.entities[i].parent_index == -1 &&
        (big1.entities[solution.entities[i].id].regions.size() > 0 || big2.entities[mapping[solution.entities[i].id]].regions.size() > 0 ||
        big1.entities[solution.entities[i].id].parent_index != -1 || big2.entities[mapping[solution.entities[i].id]].parent_index != -1)) {
            max_region++;
            solution.entities[i].regions.insert(max_region);
            solution.regions.insert(max_region);
        }
    }

    //open all necessary sites
    int max_site = solution.sites.size();
    if(solution.sites.size() > 0)
        max_site = *solution.sites.rbegin();
    else max_site = 0;
    for(int i=0;i<solution.entities.size();i++) {
        if(solution.entities[i].sites.size() == 0 &&
        (big1.entities[solution.entities[i].id].sites.size() > 0 || big2.entities[mapping[solution.entities[i].id]].sites.size() > 0 ||
        solution.entities[i].child_indices.size() != big1.entities[solution.entities[i].id].child_indices.size() ||
        solution.entities[i].child_indices.size() != big2.entities[solution.entities[i].id].child_indices.size())) {
            max_site++;
            solution.entities[i].sites.insert(max_site);
            solution.sites.insert(max_site);            
        }
    }    

    // close all shared roots
    for(int i=0;i<solution.entities.size();i++) {
        for(int j=i+1;j<solution.entities.size();j++) {
            if(solution.entities[i].regions.size() > 0 && solution.entities[j].regions.size() > 0 &&
                *solution.entities[i].regions.begin() != *solution.entities[j].regions.begin() &&
                big1.entities[solution.entities[i].id].regions.size() > 0 && big2.entities[mapping[solution.entities[i].id]].regions.size() > 0 &&
                big1.entities[solution.entities[j].id].regions.size() > 0 && big2.entities[mapping[solution.entities[j].id]].regions.size() > 0 &&
                *big1.entities[solution.entities[i].id].regions.begin() == *big1.entities[solution.entities[j].id].regions.begin() &&
                *big2.entities[mapping[solution.entities[i].id]].regions.begin() == *big2.entities[mapping[solution.entities[j].id]].regions.begin()) {
                    solution.regions.erase(*solution.entities[j].regions.begin());
                    solution.entities[j].regions.clear();
                    solution.entities[j].regions.insert(*solution.entities[i].regions.begin());
            }
        }
    }
      
    return solution;
}

auto read_bigraph(istream && infile, const string &) -> Bigraph
{
    Bigraph big = Bigraph();
    int i = 0;
    string t = read_str(infile);

    while (t.find('}') == string::npos) {
        t = read_str(infile);
        big.entities.push_back(Entity(i, t.substr(0, t.find(':')), stoi(t.substr(t.find(':')+1, t.find(')')))));
        i++;
    }

    int r = read_num(infile);
    int n = read_num(infile);
    int s = read_num(infile);    
    big.largest_component_index = n-1;

    for(int i=0;i<r;i++) big.regions.insert(i);
    for(int i=0;i<s;i++) big.sites.insert(i);

    for (int i = 0 ; i != (r + n) ; ++i) {
        for (int j = 0 ; j != (n + s) ; ++j) {       
            char x = read_char(infile);
            if (x == '1') {
                if(i < r) {
                    big.entities[j].regions.insert(i);
                }
                else if (j >= n) {
                    big.entities[i-r].sites.insert(j-n);
                }
                else {
                    big.entities[j].parent_index = i-r;
                    big.entities[i-r].child_indices.push_back(j);
                }
            }
        } 
    }

    string h = read_str(infile);
    int no_closures = 0;
    while (h == "({},") {
        std::pair<string, std::vector<int>> he;
        he.second.resize(n);
        he.first = read_str(infile);

        bool is_closed = (he.first == "{},");
        if(!is_closed)
            he.first = he.first.substr(1, he.first.length()-3);
        else {
            he.first = "closure_e" + to_string(no_closures);
            no_closures++;
        }
        read_char(infile);

        // Deal with hanging link case
        string e = read_str(infile);
        if (e == "})") {
            h = read_str(infile);
            continue;           
        }

        string c = read_str(infile);
        int ports = 0;
        he.second[stoi(e.substr(1, e.find(',') - 1))] = stoi(c.substr(0, c.find(')')));
        while (c.find('}') == string::npos) {
            e = read_str(infile);
            c = read_str(infile);
            he.second[stoi(e.substr(1, e.find(',') - 1)) ] = stoi(c.substr(0, c.find(')')));     
            ports += stoi(c.substr(0, c.find(')')));
        }

        if(is_closed) {
            Closure c;
            c.name = he.first;
            c.adjacencies = he.second;
            c.port_count = ports;
            c.id = no_closures - 1;
            big.closures.push_back(c);
        }
        else
            big.hyperedges.push_back(he);

        h = read_str(infile);
    }

    big.reachability.resize(n, std::vector<bool>(n));

    for(int i=0;i<n;i++) {
        big.reachability[i][i] = 1;
        if(big.entities[i].child_indices.size() == 0) {
            Entity t = big.entities[i];
            big.entities[i].is_leaf = true;
            while(t.parent_index != -1) {
                for(int j=0;j<n;j++)
                    if(big.reachability[big.entities[t.parent_index].id][j] == 0 && big.reachability[t.id][j] == 1) {
                        big.reachability[big.entities[t.parent_index].id][j] = 1;
                    }
                t = big.entities[t.parent_index];
            }
        }
    }

    big.original_size = n;
    return big;
}

auto full_decomp(Bigraph big) -> std::vector<Bigraph>
{
    std::vector<Bigraph> components;
    for(unsigned int i=0;i<big.entities.size();i++){
        components.push_back(Bigraph());
        components[i].original_size = big.entities.size();
        components[i].largest_component_index = i;
        components[i].reachability = big.reachability;
        components[i].nogood_id = std::pow(2, big.entities.size() + big.closures.size() - (i + 1));
        components[i].entities.push_back(Entity(big.entities[i].id, big.entities[i].control, big.entities[i].arity));
        components[i].entities[0].is_leaf = big.entities[i].is_leaf;

        components[i].regions.insert(big.entities[i].regions.begin(), big.entities[i].regions.end());
        components[i].entities[0].regions.insert(big.entities[i].regions.begin(), big.entities[i].regions.end());

        components[i].sites.insert(big.entities[i].sites.begin(), big.entities[i].sites.end());
        components[i].entities[0].sites.insert(big.entities[i].sites.begin(), big.entities[i].sites.end());

        if(big.entities[i].parent_index != -1) {
            components[i].regions.insert((big.entities[i].id * -1) - 1);
            components[i].entities[0].regions.insert((big.entities[i].id * -1) - 1);
        }

        for(unsigned int j=0;j<big.entities[i].child_indices.size();j++) {
            components[i].sites.insert((big.entities[big.entities[i].child_indices[j]].id * -1) - 1);
            components[i].entities[0].sites.insert((big.entities[big.entities[i].child_indices[j]].id * -1) - 1);
        }

        for(int j=0;j<big.hyperedges.size();j++) {
            for(int k=0;k<big.hyperedges[j].second[i];k++) {
                std::pair<string, std::vector<int>> he;
                he.first = big.hyperedges[j].first;
                he.second.resize(big.entities.size());
                he.second[i] = 1;
                components[i].hyperedges.push_back(he);
            }
        }
        for(int j=0;j<big.closures.size();j++) {
            for(int k=0;k<big.closures[j].adjacencies[i];k++) {
                std::pair<string, std::vector<int>> he;
                he.first = big.closures[j].name;
                he.second.resize(big.entities.size());
                he.second[i] = 1;
                components[i].hyperedges.push_back(he);
            }
        }
    }

    for(int i=0;i<big.closures.size();i++) {
        components.push_back(Bigraph());
        components[i+big.entities.size()].nogood_id = std::pow(2, big.entities.size() + big.closures.size() - (i + 1));
        components[i+big.entities.size()].closures.push_back(big.closures[i]);
        components[i+big.entities.size()].largest_component_index = big.entities.size()+i;
    }

    return components;
}

// Status codes: 0 is successful composition, -1 is domain wipeout, -2 is invalid composition, -3 is something that should never happen
auto element_compose(Bigraph a, Bigraph b, Bigraph c, bool lts) -> std::pair<int, Bigraph>
{
    // Add closure only if all adjacent ports exist
    if(b.entities.size() == 0 && b.closures.size() == 1) {

        // check all and add closure to the thingy, you know the drill
        std::vector<int> copy;
        copy.resize(b.closures[0].adjacencies.size());
        for(int i=0; i<a.hyperedges.size();i++) {
            if(a.hyperedges[i].first == b.closures[0].name) {
                for(int j=0;j<a.hyperedges[i].second.size();j++) {
                    copy[j] += a.hyperedges[i].second[j];
                }
                a.hyperedges.erase(a.hyperedges.begin() + i);
                i -= 1;
            }
        }

        if(copy == b.closures[0].adjacencies) {

            std::vector<std::pair<bool, std::vector<int>>> new_map_set;
            for(auto m1: a.mappings) {
                for(auto m2: b.mappings) {
                    if(std::find(m1.second.begin(), m1.second.end(), m2.second[b.closures[0].id + c.entities.size()]) == m1.second.end())
                        continue;
                    auto mapped_he = c.closures[m2.second[b.closures[0].id + c.entities.size()] - c.entities.size()];
                    bool matchable = true;
                    for(int i=0;i<b.closures[0].adjacencies.size();i++) {
                        if(b.closures[0].adjacencies[i] > 0 && (m1.second[i] == -1 || b.closures[0].adjacencies[i] != mapped_he.adjacencies[m1.second[i]])) {
                            matchable = false;
                            break;
                        }
                    }
                    if(matchable) {
                        std::vector<int> new_map = m1.second;
                        new_map[b.closures[0].id + c.entities.size()] = m2.second[b.closures[0].id + c.entities.size()];
                    }
                }
            }
            if(new_map_set.size() == 0)
                return std::pair<int, Bigraph>(-1, a);
            a.mappings = new_map_set;

            a.closures.push_back(b.closures[0]);
            a.largest_component_index = std::max(a.largest_component_index, b.largest_component_index);
            a.nogood_id = a.nogood_id | b.nogood_id;
            return std::pair<int, Bigraph>(0, a);
        }
        return std::pair<int, Bigraph>(-2, a);
    }

    if(b.entities.size() != 1)
        return std::pair<int, Bigraph>(-3, a);

    for(unsigned int i=0;i<a.entities.size();i++)
        if(a.entities[i].id == b.entities[0].id)
            return std::pair<int, Bigraph>(-3, a);

    bool is_tensor_possible = true;
    bool has_true_solution = false;
    int below_index = -1;
    int above_index = -1;

    // Look for compatible region/site pairs (have the same id)
    for(unsigned int i=0; i<a.entities.size();i++) {
        if(!lts && below_index == -1 && (b.entities[0].regions.size() > 0 && *(b.entities[0].regions.begin()) < 0)) {
            auto below_comp = find(a.entities[i].sites.begin(), a.entities[i].sites.end(), *(b.entities[0].regions.begin()));
            if(below_comp != a.entities[i].sites.end()) below_index = i;
        }
        if(above_index == -1 && (a.entities[i].regions.size() > 0 && *(a.entities[i].regions.begin()) < 0)) {
            auto above_comp = find(b.entities[0].sites.begin(), b.entities[0].sites.end(), *(a.entities[i].regions.begin()));
            if(above_comp != b.entities[0].sites.end()) above_index = i;
        }
        if(above_index > -1 && below_index > -1) break;

        if(a.reachability[a.entities[i].id][b.entities[0].id] == 1 || a.reachability[b.entities[0].id][a.entities[i].id] == 1)
            is_tensor_possible = false;
    }

    // If no compatible region/site pairs, check if tensor product is allowed, return null if not
    if(below_index == -1 && above_index == -1){
        if(! is_tensor_possible || (lts && !b.entities[0].is_leaf))
            return std::pair<int, Bigraph>(-2, a);

        // check no overlapping mapping
        std::vector<std::pair<bool, std::vector<int>>> new_map_set;
        for(auto m1 : a.mappings) {
            for(auto m2 : b.mappings) {
                if(std::find(m1.second.begin(), m1.second.end(), m2.second[b.entities[0].id]) == m1.second.end()) {
                    std::vector<int> new_map = m1.second;
                    new_map[b.entities[0].id] = m2.second[b.entities[0].id];
                    if(m1.first && m2.first)
                        has_true_solution = true;
                    new_map_set.push_back(make_pair(m1.first && m2.first, new_map));
                }
            }
        }
        if(new_map_set.size() == 0)
            return std::pair<int, Bigraph>(-1, a);
            
        a.mappings = new_map_set;

        a.entities.push_back(b.entities[0].copy());
        a.regions.insert(a.entities[a.entities.size()-1].regions.begin(), a.entities[a.entities.size()-1].regions.end());
        a.sites.insert(a.entities[a.entities.size()-1].sites.begin(), a.entities[a.entities.size()-1].sites.end());
        a.largest_component_index = std::max(a.largest_component_index, b.largest_component_index);
        a.nogood_id = a.nogood_id | b.nogood_id;

        for(int i=0;i<b.hyperedges.size();i++)
            a.hyperedges.push_back(b.hyperedges[i]);
        return std::pair<int, Bigraph>(1 - has_true_solution, a);
    }

    // for each candidate mapping, check if the parent-child relation still exists
    a.entities.push_back(b.entities[0].copy());
    std::vector<std::pair<bool, std::vector<int>>> new_map_set;
    for(auto m1 : a.mappings) {
        for(auto m2 : b.mappings) {
            if(std::find(m1.second.begin(), m1.second.end(), m2.second[b.entities[0].id]) != m1.second.end())
                continue;
            if(below_index > -1 && c.entities[m2.second[b.entities[0].id]].parent_index != m1.second[a.entities[below_index].id])
                continue;
            if(above_index > -1 && c.entities[m1.second[a.entities[above_index].id]].parent_index != m2.second[b.entities[0].id])
                continue;
            std::vector<int> new_map = m1.second;
            new_map[b.entities[0].id] = m2.second[b.entities[0].id];

            if (!lts || above_index == -1) {
                if(m1.first)
                    has_true_solution = true;
                new_map_set.push_back(make_pair(m1.first, new_map));
            }
            else {                
                if(! m2.first)
                    new_map_set.push_back(make_pair(false, new_map));
                bool valid_flag = true;
                for(unsigned int i=0;i<a.entities.size();i++) {
                    if (i == above_index)
                        continue;
                    if (m1.second[a.entities[i].id] != -1 && a.entities[i].parent_index == -1 && c.entities[m1.second[a.entities[i].id]].parent_index != -1){
                        valid_flag = false;
                        break;
                    }     
                }    
                if (valid_flag)
                    has_true_solution = true;

                new_map_set.push_back(make_pair(valid_flag, new_map));
            }
        }
    }

    if(new_map_set.size() == 0)
        return std::pair<int, Bigraph>(-1, a);
    a.mappings = new_map_set;            

    // Connect the compatible region/site pairs accordingly and return the new structure
    if(below_index > -1){
        a.entities[below_index].child_indices.push_back(a.entities.size()-1);
        a.entities[a.entities.size()-1].parent_index = below_index;

        a.entities[below_index].sites.erase(*b.entities[0].regions.begin());
        a.sites.erase(*b.entities[0].regions.begin());
        a.entities[a.entities.size()-1].regions.clear();

        a.sites.insert(a.entities[a.entities.size()-1].sites.begin(), a.entities[a.entities.size()-1].sites.end());
        a.largest_component_index = std::max(a.largest_component_index, b.largest_component_index);        
    }
    if(above_index > -1){
        a.entities[a.entities.size()-1].child_indices.push_back(above_index);
        a.entities[above_index].parent_index = a.entities.size()-1;

        a.entities[a.entities.size()-1].sites.erase(*a.entities[above_index].regions.begin());
        a.regions.erase(*a.entities[above_index].regions.begin()); 
        a.entities[above_index].regions.clear();

        a.regions.insert(a.entities[a.entities.size()-1].regions.begin(), a.entities[a.entities.size()-1].regions.end());
        a.sites.insert(a.entities[a.entities.size()-1].sites.begin(), a.entities[a.entities.size()-1].sites.end());
        a.largest_component_index = std::max(a.largest_component_index, b.largest_component_index);        
    }
    for(int i=0;i<b.hyperedges.size();i++)
        a.hyperedges.push_back(b.hyperedges[i]);
    a.nogood_id = a.nogood_id | b.nogood_id;

    return std::pair<int, Bigraph>(1 - has_true_solution, a);
}