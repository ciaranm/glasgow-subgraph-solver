#include <gss/innards/automorphisms.hh>

// #include <dejavu.h>

#include <vector>
#include "automorphisms.hh"
#include "math.h"

using std::string_view;
using std::vector;

using namespace gss::innards;

auto gss::innards::automorphisms_as_order_constraints(const InputGraph & i, const bool with_generators, const bool degree_sequence, vector<int> &orbit_sizes) -> OrderConstraints {
    std::vector<int> base;

    if (degree_sequence) {
        base.resize(i.size());
        std::iota(base.begin(), base.end(), 0);
        stable_sort(base.begin(), base.end(), [&](int a, int b) -> bool {
            return -i.degree(a) < -i.degree(b);
        });
    }

    return with_generators ? automorphisms_as_order_constraints_with_generators(i, base, orbit_sizes) : automorphisms_as_order_constraints(i, base, orbit_sizes);
}

/**
 * Uses orbits
*/
auto gss::innards::automorphisms_as_order_constraints(const InputGraph & i, std::vector<int> base, std::vector<int> &orbit_sizes) -> OrderConstraints
{
    dejavu::static_graph g;     // Declare static graph
    unsigned long n_simple_edges = 0;       // Number of undirected edges
    int nv = i.size();
    if (!i.directed()) {
        i.for_each_edge([&](int f, int t, string_view) {
            if (f < t)
                ++n_simple_edges;
        });
        g.initialize_graph(nv, n_simple_edges);       // Initialise g with nv, ne
        vector<int> vertices;       // List of vertices
        for (int v = 0, v_end = i.size(); v != v_end; ++v) {        // For each vertex
            bool loop = i.adjacent(v, v);       // Check if v is adjacent to itself
            int deg = i.degree(v);
            vertices.push_back(g.add_vertex(loop, deg - loop));   // Add colour and non-looping degree
        }
        i.for_each_edge([&](int f, int t, string_view) {
            if (f < t) 
                g.add_edge(vertices[f], vertices[t]);       // Add edges (undirected)
        });
    }
    else {
        vector<int> degrees;
        for (int v = 0; v < i.size(); v++) {
            degrees.push_back(i.degree(v));
        }
        i.for_each_edge([&](int f, int t, string_view) {
            if (i.adjacent(t, f)) {
                if (f < t) {
                    ++n_simple_edges;           // Bidirectional edges are undirected already
                }
            }
            else {
                nv += 2;     // Add extra vertices to indicate direction
                n_simple_edges += 3;
                ++degrees[t];
            }
        });
        g.initialize_graph(nv, n_simple_edges);
        vector<int> vertices;       // List of vertices
        for (int v = 0, v_end = i.size(); v != v_end; ++v) {        // For each vertex in the orignal graph
            bool loop = i.adjacent(v, v);       // Check if v is adjacent to itself
            vertices.push_back(g.add_vertex(loop ? 1 : 0, loop ? degrees[v] - 1 : degrees[v]));   // Add colour and non-looping degree
        }
        vector<int> direction_vertices;
        for (int v = i.size(), v_end = nv; v != v_end; v += 2) { // For each additional vertex (direction indicators)
            direction_vertices.push_back(g.add_vertex(2, 2));   // From
            direction_vertices.push_back(g.add_vertex(3, 2));   // To
        }

        i.for_each_edge([&](int f, int t, string_view) {
            if (i.adjacent(t, f)) {
                if (f < t) {
                    g.add_edge(vertices[f], vertices[t]);
                }
            }
            else {
                g.add_edge(vertices[f], direction_vertices[0]);
                g.add_edge(direction_vertices[0], direction_vertices[1]);
                g.add_edge(vertices[t], direction_vertices[1]);
                direction_vertices.erase(direction_vertices.begin());
                direction_vertices.erase(direction_vertices.begin());
            }
        });
    }

    dejavu::groups::random_schreier rschreier{nv};        // Random Schreier structure with domain size nv

    rschreier.set_base(base);       // Empty base to begin with

    dejavu::hooks::schreier_hook hook(rschreier);
    dejavu::solver s;
    s.set_print(false);
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

    for (unsigned x = 0; x < base.size(); ++x) {        // For each base element
        const int v = base.at(x);
        orbit_sizes.push_back(1);
        for (auto & o : rschreier.get_fixed_orbit(x)) {     // For each orbit point of v
            if (o != v) {     // If v is mapped to somewhere else
                unique_list.emplace(v, o);        // Enforce v<o
                ++orbit_sizes[x];
            }
        }
    }

    for (std::pair<int, int> p: unique_list) {
        result.emplace_back(i.vertex_name(p.first), i.vertex_name(p.second));
    }

    return result;
}


/**
 * Uses generators
*/
auto gss::innards::automorphisms_as_order_constraints_with_generators(const InputGraph & i, std::vector<int> base, std::vector<int> &orbit_sizes) -> OrderConstraints
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

    rschreier.set_base(base);       // Empty base to begin with

    dejavu::hooks::schreier_hook hook(rschreier);
    dejavu::solver s;
    s.set_print(false);
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

    for (size_t x = 0; x < base.size(); x++) {
        orbit_sizes.push_back(o.orbit_size(base.at(x)));
    }

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
            for (unsigned int v_index = 1; v_index < cycle.size(); v_index++) {
                    if (first) {
                        unique_list.emplace(cycle[0], cycle[v_index]);
                        first = false;      // Only use the first cycle of the generator for correctness - see proof (too big for the margin)
                    }
            }
        }
    }

    for (std::pair<int, int> p: unique_list) {
        result.emplace_back(i.vertex_name(p.first), i.vertex_name(p.second));
    }
    
    return result;
}

/**
 * Build a schreier structure for symmetry calculations during search
 */
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
}

/**
 * Calculate any non-trivial orbits with the new (partial) base given, and make constraints from them
 */
auto gss::innards::dynamic_order_constraints(int sz, vector<int> base, vector<int> &orbit_sizes, dejavu::groups::random_schreier &rschreier, std::vector<std::pair<unsigned int, unsigned int>> &constraints) -> void {

    //TODO we only actually care about adding one extra layer of the chain here, maybe there is a way to limit depth of layer computation in Dejavu?
    rschreier.set_base(base);           // Recompute the stabiliser chain with the provided (partial) base
    std::vector<std::pair<unsigned int, unsigned int>> cons;
    std::vector<std::vector<int>> orbs;

    dejavu::groups::orbit o{sz};      // Orbit structure of size sz
    // int x = base.size() - 1;            // Only interested in the most recently added base point
    for (unsigned int x = 0; x < base.size(); x++) {
        rschreier.get_stabilizer_orbit(x, o);       // Retrieve stabiliser orbit at this point
        int y = base.at(x);
        orbit_sizes[x] = o.orbit_size(y);
        for (int z = 0; z != sz; ++z) {     // For each vertex z in the graph
            if (o.are_in_same_orbit(y,z) && y != z) {       // If y,z are symmetric (but not equal) at this stabiliser chain layer 
                cons.push_back(std::make_pair(y,z));     // Add constraint y<z
            }
        }
    }

    constraints = cons;
}

/**
 * Calculate the generating set for the automorphism group of an input graph (picking a random base)
 */
auto::gss::innards::generating_set(const InputGraph &i) -> std::vector<std::vector<unsigned int>> {
    std::vector<int> base;
    return generating_set(i, base);
}

/**
 * Calulcate the generating set for the automorphism group of an input graph (given a [partial] base)
 */
auto gss::innards::generating_set(const InputGraph &i, std::vector<int> base) -> std::vector<std::vector<unsigned int>> {
    dejavu::static_graph g;     // Declare static graph
    unsigned long n_simple_edges = 0;       // Number of undirected edges
    int n_undirected_vertices = i.size();
    if (i.directed()) {
        vector<int> degrees;
        for (int v = 0; v < i.size(); v++) {
            degrees.push_back(i.degree(v));
        }
        i.for_each_edge([&](int f, int t, string_view) {
            if (i.adjacent(t, f)) {
                if (f < t) {
                    ++n_simple_edges;           // Bidirectional edges are undirected already
                }
            }
            else {
                n_undirected_vertices += 2;     // Add extra vertices to indicate direction
                n_simple_edges += 3;
                ++degrees[t];
            }
        });
        g.initialize_graph(n_undirected_vertices, n_simple_edges);
        vector<int> vertices;       // List of vertices
        for (int v = 0, v_end = i.size(); v != v_end; ++v) {        // For each vertex in the orignal graph
            bool loop = i.adjacent(v, v);       // Check if v is adjacent to itself
            vertices.push_back(g.add_vertex(loop ? 1 : 0, loop ? degrees[v] - 1 : degrees[v]));   // Add colour and non-looping degree
        }
        vector<int> direction_vertices;
        for (int v = i.size(), v_end = n_undirected_vertices; v != v_end; v += 2) { // For each additional vertex (direction indicators)
            direction_vertices.push_back(g.add_vertex(2, 2));   // From
            direction_vertices.push_back(g.add_vertex(3, 2));   // To
        }
        i.for_each_edge([&](int f, int t, string_view) {
            if (i.adjacent(t, f)) {
                if (f < t) {
                    g.add_edge(vertices[f], vertices[t]);
                }
            }
            else {
                int lower = f < t ? f : t;
                int higher = f < t ? t : f;
                g.add_edge(vertices[lower], direction_vertices[0]);
                g.add_edge(direction_vertices[0], direction_vertices[1]);
                g.add_edge(vertices[higher], direction_vertices[1]);
                direction_vertices.erase(direction_vertices.begin());
                direction_vertices.erase(direction_vertices.begin());
            }
        });
    }
    else {
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
    }

    dejavu::groups::random_schreier rschreier{n_undirected_vertices};        // Random Schreier structure with domain size nv

    rschreier.set_base(base);       // Empty base to begin with

    dejavu::hooks::schreier_hook hook(rschreier);
    dejavu::solver s;
    s.set_print(false);
    s.automorphisms(&g, hook.get_hook());       // Compute automorphisms of g

    dejavu::groups::automorphism_workspace a(n_undirected_vertices);
    std::vector<std::vector<unsigned int>> mappings;            // Store generators (and identity) in a sensible notation
    std::vector<unsigned int> identity;                      // We need an identity since we're composing elements from two sets
    for (int y = 0; y < i.size(); y++) {
        identity.push_back(y);                      // Every element maps to itself
    }
    mappings.push_back(identity);
    for (int x = 0; x < rschreier.get_number_of_generators(); x++) {
        rschreier.get_generator(x, a);              // Get the xth generator of Aut(G), p
        std::vector<unsigned int> mapping;
        for (int y = 0; y < i.size(); y++) {
            mapping.push_back(a.p()[y]);            // mapping[i] maps to p*i
        }
        mappings.push_back(mapping);
    }

    return mappings;
}

auto gss::innards::dynamic_generating_set(std::vector<int> & base, int sz, dejavu::groups::random_schreier &rschreier, std::vector<std::vector<int>> & generators) -> void {
    rschreier.set_base(base);       // Empty base to begin with

    dejavu::groups::automorphism_workspace a(sz);
    generators.erase(generators.begin() + 1, generators.end());     // Leave the identity
    std::vector<int> mapping;
    for (int x = 0; x < rschreier.get_number_of_generators(); x++) {
        rschreier.get_generator(x, a);              // Get the xth generator of Aut(G), p
        mapping.clear();
        for (size_t y = 0; y < base.size(); y++) {
            mapping.push_back(a.p()[y]);            // mapping[i] maps to p*i
            // std::cout << y << "->" << mapping[y] << " ";
        }
        // std::cout << "\n";
        generators.push_back(mapping);
    }
}

auto gss::innards::invert_automorphism(std::vector<unsigned int> aut) -> std::vector<unsigned int> {
    std::vector<unsigned int> res(aut.size());
    for (unsigned int i = 0; i < aut.size(); i++) {
        res[aut[i]] = i; 
    }
    return res;
}

auto gss::innards::invert_list(std::vector<std::vector<unsigned int>> ls) -> std::vector<std::vector<unsigned int>> {
    std::vector<std::vector<unsigned int>> res;
    for (unsigned int i = 0; i < ls.size(); i++) {
        res.push_back(invert_automorphism(ls[i]));
    }
    return res;
    
}