/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include <algorithm>
#include "formats/common_bigraph.hh"
#include "formats/input_graph.hh"
#include <fstream>
#include <iostream>
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

    for(unsigned int i=0;i<entities.size();i++){
        big.entities.push_back(entities[i].copy());
    }
    return big;
}

auto Bigraph::toString() const -> string
{
    string out = "{";
    for(unsigned int i=0;i<entities.size();i++) {
        out += '(' + to_string(i) + ", " + entities[i].control + ':' + to_string(entities[i].arity) + ')';
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

    return out;
}

auto Bigraph::encode(bool target) const -> InputGraph
{
    InputGraph result{ 0, true, true, true };

    if(target) result.resize(entities.size() + regions.size() + sites.size());
    else result.resize(entities.size());

    if(target) {
        for(unsigned int i=0; i<regions.size();i++) {
            result.set_vertex_label(i, "ROOT");
            result.set_vertex_name(i, "ROOT" + to_string(i));
        }
        for(unsigned int i=regions.size(); i<regions.size()+entities.size();i++) {
            result.set_vertex_label(i, entities[i-regions.size()].control);
            result.set_vertex_name(i, to_string(i-regions.size()));
        }
        for(unsigned int i=regions.size()+entities.size(); i<regions.size()+entities.size()+sites.size();i++) {
            result.set_vertex_label(i, "SITE");
            result.set_vertex_name(i, "SITE" + to_string(i));
        }
    }
    else {
        for(unsigned int i=0; i<entities.size();i++)
            result.set_vertex_label(i, entities[i].control);
    }

    int index = 0;
    for(auto r : regions) {
        for(unsigned int i=0;i<entities.size();i++) {
            if(find(entities[i].regions.begin(), entities[i].regions.end(), r) != entities[i].regions.end()) {
                if(target)
                    result.add_directed_edge(index, regions.size() + i, "dir");
                else {
                    result.set_child_of_root(i);
                    result.add_pattern_root_edge(index, i);
                }
            }
        }
        index++;
    }
    index = 0;
    for(auto s : sites) {
        for(unsigned int i=0;i<entities.size();i++) {
            if(find(entities[i].sites.begin(), entities[i].sites.end(), s) != entities[i].sites.end()) {
                if(target)
                    result.add_directed_edge(regions.size() + i, regions.size() + entities.size() + index, "dir");
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
            if(target)
                result.add_directed_edge(regions.size() + k, regions.size() + i, "dir");
            else
                result.add_directed_edge(k, i, "dir");
        }
    }
    return result;
}

auto free_all_entities(Bigraph a) -> Bigraph
{
    if(a.regions.size() > 0) {
        int max_region = *a.regions.rbegin();
        for(int i=0;i<a.entities.size();i++) {
            if(a.entities[i].parent_index == -1 && a.entities[i].regions.size() == 0) {
                max_region++;
                a.entities[i].regions.insert(max_region);
                a.regions.insert(max_region);
            }
            else if(a.entities[i].regions.size() > 0) {
                for(int j=i;j<a.entities.size();j++) {
                    if(*a.entities[i].regions.begin() == *a.entities[j].regions.begin()) {
                        max_region++;
                        a.entities[j].regions.clear();
                        a.entities[j].regions.insert(max_region);
                        a.regions.insert(max_region);
                    }
                }
            }
        }
    }
    if(a.sites.size() > 0) {
        int max_site = *a.sites.rbegin();
        for(int i=0;i<a.entities.size();i++) {
            if(a.entities[i].sites.size() == 0) {
                max_site++;
                a.entities[i].sites.insert(max_site);
                a.sites.insert(max_site);
            }
        }
    }
    return a;
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

    big.reachability.resize(n, std::vector<bool>(n));

    for(int i=0;i<n;i++) {
        big.reachability[i][i] = 1;
        if(big.entities[i].child_indices.size() == 0) {
            Entity t = big.entities[i];
            while(t.parent_index != -1) {
                for(int j=0;j<n;j++)
                    if(big.reachability[big.entities[t.parent_index].id][j] == 0 && big.reachability[t.id][j] == 1) {
                        big.reachability[big.entities[t.parent_index].id][j] = 1;
                    }
                t = big.entities[t.parent_index];
            }
        }
    }

    return big;
}

auto full_decomp(Bigraph big) -> std::vector<Bigraph>
{
    std::vector<Bigraph> components;
    for(unsigned int i=0;i<big.entities.size();i++){
        components.push_back(Bigraph());
        components[i].largest_component_index = i;
        components[i].reachability = big.reachability;
        components[i].entities.push_back(Entity(big.entities[i].id, big.entities[i].control, big.entities[i].arity));

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

    }
    return components;
}

auto element_compose(Bigraph a, Bigraph b) -> std::optional<Bigraph>
{
    if(b.entities.size() != 1)
        return std::nullopt;

    for(unsigned int i=0;i<a.entities.size();i++)
        if(a.entities[i].id == b.entities[0].id)
            return std::nullopt;

    bool is_tensor_possible = true;
    for(unsigned int i=0; i<a.entities.size();i++) {
        if(b.entities[0].regions.size() > 0 && *(b.entities[0].regions.begin()) < 0) {
            auto below_comp = find(a.entities[i].sites.begin(), a.entities[i].sites.end(), *(b.entities[0].regions.begin()));
            if(below_comp != a.entities[i].sites.end()) {
                a.entities.push_back(b.entities[0].copy());

                a.entities[i].child_indices.push_back(a.entities.size()-1);
                a.entities[a.entities.size()-1].parent_index = i;

                a.entities[i].sites.erase(*b.entities[0].regions.begin());
                a.sites.erase(*b.entities[0].regions.begin());
                a.entities[a.entities.size()-1].regions.clear();

                a.sites.insert(a.entities[a.entities.size()-1].sites.begin(), a.entities[a.entities.size()-1].sites.end());
                a.largest_component_index = std::max(a.largest_component_index, b.largest_component_index);
                return a;
            }
        }
        else if(a.entities[i].regions.size() > 0 && *(a.entities[i].regions.begin()) < 0) {
            auto above_comp = find(b.entities[0].sites.begin(), b.entities[0].sites.end(), *(a.entities[i].regions.begin()));
            if(above_comp != b.entities[0].sites.end()) {
                a.entities.push_back(b.entities[0].copy());

                a.entities[a.entities.size()-1].child_indices.push_back(i);
                a.entities[i].parent_index = a.entities.size()-1;

                a.entities[a.entities.size()-1].sites.erase(*a.entities[i].regions.begin());
                a.regions.erase(*a.entities[i].regions.begin()); 
                a.entities[i].regions.clear();

                a.regions.insert(a.entities[a.entities.size()-1].regions.begin(), a.entities[a.entities.size()-1].regions.end());
                a.sites.insert(a.entities[a.entities.size()-1].sites.begin(), a.entities[a.entities.size()-1].sites.end());
                a.largest_component_index = std::max(a.largest_component_index, b.largest_component_index);
                return a;
            }
        }

        if(a.reachability[a.entities[i].id][b.entities[0].id] == 1 || a.reachability[b.entities[0].id][a.entities[i].id] == 1)
            is_tensor_possible = false;
    }

    if(! is_tensor_possible) 
        return std::nullopt;

    a.entities.push_back(b.entities[0].copy());
    a.regions.insert(a.entities[a.entities.size()-1].regions.begin(), a.entities[a.entities.size()-1].regions.end());
    a.sites.insert(a.entities[a.entities.size()-1].sites.begin(), a.entities[a.entities.size()-1].sites.end());
    a.largest_component_index = std::max(a.largest_component_index, b.largest_component_index);
    return a;
}