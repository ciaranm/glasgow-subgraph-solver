#include <gss/innards/automorphisms.hh>

// #include <dejavu.h>

#include <vector>
#include "automorphisms.hh"

using std::string_view;
using std::vector;

using namespace gss::innards;

auto gss::innards::automorphisms_as_order_constraints(const InputGraph & i, const bool with_generators) -> OrderConstraints {
    std::vector<int> base;

    return with_generators ? automorphisms_as_order_constraints_with_generators(i, base) : automorphisms_as_order_constraints(i, base);
}

/**
 * Uses orbits
*/
auto gss::innards::automorphisms_as_order_constraints(const InputGraph & i, std::vector<int> base) -> OrderConstraints
{
    dejavu::static_graph g;     // Declare static graph
    unsigned long n_simple_edges = 0;       // Number of undirected edges
    unsigned long nv = i.size();
    if (!i.directed()) {
        i.for_each_edge([&](int f, int t, string_view) {
            if (f < t)
                ++n_simple_edges;
        });
    }
    else {
        for (int v = 0, v_end = i.size(); v != v_end; ++v) {
            for (int u = v+1, u_end = i.size(); u != u_end; ++u) {
                if (i.adjacent(v,u) && i.adjacent(u,v)) {
                    ++n_simple_edges;
                }
                else if (i.adjacent(u, v) || i.adjacent(v,u)) {
                    n_simple_edges += 3;
                    nv += 2;
                }
            }
        }
    }
    
    g.initialize_graph(nv, n_simple_edges);       // Initialise g with nv, ne
    vector<int> vertices;       // List of vertices
    for (int v = 0, v_end = i.size(); v != v_end; ++v) {        // For each vertex
        bool loop = i.adjacent(v, v);       // Check if v is adjacent to itself
        int deg = i.degree(v);
        if (i.directed()) {
            for (int u = 0, u_end = i.size(); u != u_end; u++) {
                if (u == v) continue;
                deg += i.adjacent(u, v);
            }
        }
        vertices.push_back(g.add_vertex(loop, deg - loop));   // Add colour and non-looping degree
    }
    i.for_each_edge([&](int f, int t, string_view) {
        if (!i.directed()) {
            if (f < t) 
                g.add_edge(vertices[f], vertices[t]);       // Add edges (undirected)
        }
        else {
            if (i.adjacent(t,f) && f < t) {
                g.add_edge(vertices[f], vertices[t]);
            }
            else {
                vertices.push_back(g.add_vertex(2, 2));
                vertices.push_back(g.add_vertex(3, 2));
                g.add_edge(vertices[f], vertices[vertices.size() - 2]);
                g.add_edge(vertices[vertices.size() - 2], vertices[vertices.size() - 1]);
                g.add_edge(vertices[t], vertices[vertices.size() - 1]);
            }
        }
    });

    dejavu::groups::random_schreier rschreier{nv};        // Random Schreier structure with domain size nv
    // std::vector<int> base;

    rschreier.set_base(base);       // Empty base to begin with

    dejavu::hooks::schreier_hook hook(rschreier);
    dejavu::solver s;
    s.automorphisms(&g, hook.get_hook());       // Compute automorphisms of g

    OrderConstraints result;        // List of constraint pairs
    std::set<std::pair<int,int>> unique_list;

    bool stab_trivial = false;
    dejavu::groups::orbit o{nv};      // Orbit structure of size nv
    while (! stab_trivial) {        // While stabiliser is non-trivial
        rschreier.get_stabilizer_orbit(base.size(), o);     // Get orbit partition
        stab_trivial = true;
        for (int v = 0, v_end = i.size(); v != v_end; ++v) {        // For each vertex
            if (o.orbit_size(v) > 1) {      // If stabiliser orbit(v) is non-trivial
                stab_trivial = false;
                base.push_back(v);          // Add v to base
            }
        }

        if (! stab_trivial)
            rschreier.set_base(base);       // Update base
    }

    std::cout << "aut_group_sz = " << s.get_automorphism_group_size() << "\n";

    std::cout << "Base:\n";
    for (unsigned x = 0; x < base.size(); ++x) {        // For each base element
        auto v = base.at(x);
        std::cout << v << " : ";
        
        for (auto & o : rschreier.get_fixed_orbit(x)) {     // For each orbit point of v
            std::cout << o << " ";
            if (o != v)     // If v is mapped to somewhere else
                unique_list.emplace(v, o);        // Enforce v<o
        }
        std::cout << "\n";
    }

    for (std::pair<int, int> p: unique_list) {
        result.emplace_back(i.vertex_name(p.first), i.vertex_name(p.second));
    }

    return result;
}


/**
 * Uses generators
*/
auto gss::innards::automorphisms_as_order_constraints_with_generators(const InputGraph & i, std::vector<int> base) -> OrderConstraints
{
    dejavu::static_graph g;     // Declare static graph
    unsigned long n_simple_edges = 0;       // Number of undirected edges
    i.for_each_edge([&](int f, int t, string_view) {
        if (f < t)
            ++n_simple_edges;
    });
    g.initialize_graph(i.size(), n_simple_edges);       // Initialise g with nv, ne
    vector<int> vertices;       // List of vertices
    for (int v = 0, v_end = i.size(); v != v_end; ++v) {        // For each vertex
        bool loop = i.adjacent(v, v);       // Check if v is adjacent to itself
        vertices.push_back(g.add_vertex(loop ? 1 : 0, loop ? i.degree(v) - 1 : i.degree(v)));   // Add colour and non-looping degree
    }
    i.for_each_edge([&](int f, int t, string_view) {
        if (f < t)
            g.add_edge(vertices[f], vertices[t]);       // Add edges (undirected)
    });

    dejavu::groups::random_schreier rschreier{i.size()};        // Random Schreier structure with domain size nv
    // std::vector<int> base;

    rschreier.set_base(base);       // Empty base to begin with

    dejavu::hooks::schreier_hook hook(rschreier);
    dejavu::solver s;
    s.automorphisms(&g, hook.get_hook());       // Compute automorphisms of g

    OrderConstraints result;        // List of constraint pairs
    std::set<std::pair<int,int>> unique_list;

    bool stab_trivial = false;
    dejavu::groups::orbit o{i.size()};      // Orbit structure of size nv
    while (! stab_trivial) {        // While stabiliser is non-trivial
        rschreier.get_stabilizer_orbit(base.size(), o);     // Get orbit partition
        stab_trivial = true;
        for (int v = 0, v_end = i.size(); v != v_end; ++v) {        // For each vertex
            if (o.orbit_size(v) > 1) {      // If stabiliser orbit(v) is non-trivial
                stab_trivial = false;
                base.push_back(v);          // Add v to base
            }
        }

        if (! stab_trivial)
            rschreier.set_base(base);       // Update base
    }

    std::cout << "Generators:\n";
    dejavu::groups::automorphism_workspace a(i.size());
    std::set<std::set<std::vector<int>>> gens;      // Generators are sets of lists (cycles)
    for (int x = 0; x < rschreier.get_number_of_generators(); x++) {
        rschreier.get_generator(x, a);
        std::set<vector<int>> cycle_notation;
        std::set<int> counted;
        for (int y = 0; y < i.size(); y++) {
            if (counted.find(y) == counted.end()) {
                int goes_to = a.p()[y];
                std::vector<int> cycle;
                cycle.push_back(y);
                counted.emplace(y);
                while (goes_to != y) {
                    cycle.push_back(goes_to);
                    counted.emplace(goes_to);
                    goes_to = a.p()[goes_to];
                }
                if (cycle.size() > 1) {
                    cycle_notation.emplace(cycle);
                }
            }
        }
        if (cycle_notation.size() > 0) {
            gens.emplace(cycle_notation);
        }
    }

    for (auto gen : gens) {
        bool first = true;
        for (auto cycle : gen) {
            std::cout << "( ";
            std::cout << cycle[0] << " ";
            for (unsigned int v_index = 1; v_index < cycle.size(); v_index++) {
                    std::cout << cycle[v_index] << " ";
                    if (first) {
                        unique_list.emplace(cycle[0], cycle[v_index]);
                        first = false;      // Only use the first cycle of the generator for correctness - see proof (too big for the margin)
                    }
            }
            std::cout << ")";
        }
        std::cout << "\n";
    }

    for (std::pair<int, int> p: unique_list) {
        result.emplace_back(i.vertex_name(p.first), i.vertex_name(p.second));
    }
    
    return result;
}

auto gss::innards::initialise_dynamic_structure(dejavu::groups::random_schreier &r, std::vector<innards::SVOBitset> m, const bool directed) -> void {
    dejavu::static_graph g;     // Declare static graph

    if (!directed) {
        unsigned long n_simple_edges = 0;       // Number of undirected edges
        for (int i = 0, i_end = m.size(); i != i_end; ++i) {
            for (int j = i + 1; j != i_end; ++j) {
                if (m[i].test(j)) ++n_simple_edges;
            }
        }
        g.initialize_graph(m.size(), n_simple_edges);       // Initialise g with nv, ne    
        vector<int> vertices;       // List of vertices
        for (int v = 0, v_end = m.size(); v != v_end; ++v) {        // For each vertex
            vertices.push_back(g.add_vertex(m[v].test(v), m[v].count() - m[v].test(v)));   // Add colour and non-looping degree
        }
        for (int i = 0, i_end = m.size(); i != i_end; ++i) {
            for (int j = i + 1; j != i_end; ++j) {
                if (m[i].test(j)) g.add_edge(vertices[i], vertices[j]);
            }
        }
    }
    else {
        unsigned long n_simple_edges = 0;   // Number of directed edges
        unsigned long n_vertices = m.size();    // Number of vertices
        for (int i = 0, i_end = m.size(); i != i_end; ++i) {
            for (int j = i + 1, j_end = m.size(); j != j_end; ++j) {
                if (m[i].test(j) && m[j].test(i)) {
                    ++n_simple_edges;
                }
                else if (m[i].test(j) || m[j].test(i)) {
                    n_simple_edges += 3;    // Undirect that edge!
                    n_vertices += 2;        // Colours to represent direction
                }
                // if (m[i].test(j)) n_simple_edges += 3;
                // if (m[j].test(i)) n_simple_edges += 3;
            }
        }
        g.initialize_graph(n_vertices, n_simple_edges);       // Initialise g
        vector<int> vertices;       // List of vertices
        for (int v = 0, v_end = m.size(); v != v_end; ++v) {        // For each vertex
            int deg = 0;
            for (int u = 0, u_end = m.size(); u != u_end; ++u) {
                deg += m[v].test(u) || m[u].test(v);
            }
            deg -= 2 * m[v].test(v);
            vertices.push_back(g.add_vertex(m[v].test(v), deg));   // Add colour and non-looping degree
        }
        for (int i = 0, i_end = m.size(); i != i_end; ++i) {
            for (int j = i + 1; j != i_end; ++j) {
                if (m[i].test(j) && m[j].test(i))  {
                    g.add_edge(vertices[i], vertices[j]);
                }
                else if (m[i].test(j)) {
                    vertices.push_back(g.add_vertex(2, 2));
                    vertices.push_back(g.add_vertex(3, 2));
                    g.add_edge(vertices[i], vertices[vertices.size() - 2]);
                    g.add_edge(vertices[vertices.size() - 2], vertices[vertices.size() - 1]);
                    g.add_edge(vertices[j], vertices[vertices.size() - 1]);
                }
                else if (m[j].test(i)) {
                    vertices.push_back(g.add_vertex(2, 2));
                    vertices.push_back(g.add_vertex(3, 2));
                    g.add_edge(vertices[j], vertices[vertices.size() - 2]);
                    g.add_edge(vertices[vertices.size() - 2], vertices[vertices.size() - 1]);
                    g.add_edge(vertices[i], vertices[vertices.size() - 1]);
                }
            }
        }
    }    

    vector<int> empty_base;

    r.set_base(empty_base);       // Let dejavu start with an empty base

    dejavu::hooks::schreier_hook hook(r);
    dejavu::solver s;
    s.set_print(false);

    s.automorphisms(&g, hook.get_hook());       // Compute automorphisms of g

    // std::cout << "aut_group_sz = " << s.get_automorphism_group_size() << "\n";

    // std::cout << "Initialised rschreier!\n";
}


auto gss::innards::dynamic_order_constraints(int sz, vector<int> base, dejavu::groups::random_schreier &rschreier, std::vector<std::pair<unsigned int, unsigned int>> &constraints) -> void {

    //TODO we only actually care about adding one extra layer of the chain here, maybe there is a way to limit depth of layer computation in Dejavu?
    rschreier.set_base(base);           // Recompute the stabiliser chain with the provided (partial) base
    std::vector<std::pair<unsigned int, unsigned int>> cons;

    dejavu::groups::orbit o{sz};      // Orbit structure of size sz
    // std::cout << "Base:\n";
    // int x = base.size() - 1;            // Only interested in the most recently added base point
    for (int x = 0; x < base.size(); x++) {
        rschreier.get_stabilizer_orbit(x, o);       // Retrieve stabiliser orbit at this point
        int y = base.at(x);
        // std::cout << y << " : ";
        for (int z = 0; z != sz; ++z) {     // For each vertex z in the graph
            if (o.are_in_same_orbit(y,z) && y != z) {       // If y,z are symmetric (but not equal) at this stabiliser chain layer 
                // std::cout << z << " ";
                cons.push_back(std::make_pair(y,z));     // Add constraint y<z
            }
        }
        // std::cout << "\n";
    }

    constraints = cons;
}
