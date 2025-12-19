#include <gss/innards/automorphisms.hh>

#include <vector>
#include "automorphisms.hh"
#include "math.h"

using std::string_view;
using std::vector;

using namespace gss::innards;

/**
 * Uses orbits
*/
auto gss::innards::automorphisms_as_order_constraints(const InputGraph & i, std::vector<int> & base, std::vector<int> &orbit_sizes, const bool degree_sequence) -> OrderConstraints
{
    dejavu::static_graph g = build_static_graph(i);
    int nv = g.get_sgraph()->v_size;

    dejavu::groups::random_schreier rschreier{nv};        // Random Schreier structure with domain size nv

    rschreier.set_base(base);       // Empty base to begin with

    dejavu::hooks::schreier_hook hook(rschreier);
    dejavu::solver s;
    s.set_print(false);
    s.automorphisms(&g, hook.get_hook());       // Compute automorphisms of g

    OrderConstraints result;        // List of constraint pairs
    std::set<std::pair<int,int>> unique_list;

    std::vector<int> vertex_order(i.size());
    std::iota(vertex_order.begin(), vertex_order.end(), 0);
    if (degree_sequence) {
        stable_sort(vertex_order.begin(), vertex_order.end(), [&](int a, int b) -> bool {
            return -i.degree(a) < -i.degree(b);
        });
    }

    bool stab_trivial = false;
    dejavu::groups::orbit o{nv};      // Orbit structure of size nv
    while (! stab_trivial) {        // While stabiliser is non-trivial
        rschreier.get_stabilizer_orbit(base.size(), o);     // Get orbit partition
        stab_trivial = true;
        for (auto & v : vertex_order) {        // For each vertex
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
        int o_sz = 1;
        for (auto & o : rschreier.get_fixed_orbit(x)) {     // For each orbit point of v
            if (o != v) {     // If v is mapped to somewhere else
                unique_list.emplace(v, o);        // Enforce v<o
                o_sz++;
            }
        }
        orbit_sizes.push_back(o_sz);
    }

    for (std::pair<int, int> p: unique_list) {
        result.emplace_back(i.vertex_name(p.first), i.vertex_name(p.second));
    }

    return result;
}

/**
 * Build a schreier structure for symmetry calculations during search - return true if nontrivial
 */
auto gss::innards::initialise_dynamic_structure(dejavu::groups::random_schreier &r, std::vector<innards::SVOBitset> m, const bool directed) -> bool {
    dejavu::static_graph g = build_static_graph(m, directed);     // Declare static graph

    vector<int> empty_base;

    r.set_base(empty_base);       // Let dejavu start with an empty base

    dejavu::hooks::schreier_hook hook(r);
    dejavu::solver s;
    s.set_print(false);

    s.automorphisms(&g, hook.get_hook());       // Compute automorphisms of g

    dejavu::big_number sz = s.get_automorphism_group_size();
    if (sz.mantissa == 1 && sz.exponent == 0) return false;

    return true;
}

/**
 * Calculate any non-trivial orbits with the new (partial) base given, and make constraints from them
 */
auto gss::innards::dynamic_order_constraints(int sz, vector<int> &base, vector<int> &orbit_sizes, dejavu::groups::random_schreier &rschreier, std::vector<std::pair<unsigned int, unsigned int>> &constraints) -> void {
    // TODO rewrite with get_fixed_orbit
    //TODO we only actually care about adding one extra layer of the chain here, maybe there is a way to limit depth of layer computation in Dejavu?
    rschreier.set_base(base);           // Recompute the stabiliser chain with the provided (partial) base
    std::vector<std::pair<unsigned int, unsigned int>> cons;
    std::vector<std::vector<int>> orbs;

    dejavu::groups::orbit o{sz};      // Orbit structure of size sz
    // int x = base.size() - 1;            // Only interested in the most recently added base point
    for (unsigned int x = 0; x < base.size(); x++) {
        rschreier.get_stabilizer_orbit(x, o);       // Retrieve stabiliser orbit at this point
        int y = base.at(x);
        orbit_sizes[y] = o.orbit_size(y);
        for (int z = 0; z != sz; ++z) {     // For each vertex z in the graph
            if (o.are_in_same_orbit(y,z) && y != z) {       // If y,z are symmetric (but not equal) at this stabiliser chain layer 
                cons.push_back(std::make_pair(y,z));     // Add constraint y<z
            }
        }
    }

    constraints = std::move(cons);
}

/**
 * Calulcate the generating set for the automorphism group of an input graph (given a [partial] base)
 */
auto gss::innards::coset_reps(const InputGraph &i, std::vector<int> & orbit_sizes, std::vector<int> & base, const bool degree_sequence) -> std::vector<std::vector<unsigned int>> {
    dejavu::static_graph g = build_static_graph(i);
    int nv = g.get_sgraph()->v_size;

    orbit_sizes.resize(i.size(), 1);

    dejavu::groups::random_schreier rschreier{nv};        // Random Schreier structure with domain size nv

    rschreier.set_base(base);       // Empty base to begin with

    std::vector<int> vertex_order(i.size());                // Decide on a fixed vertex order
    std::iota(vertex_order.begin(), vertex_order.end(), 0);
    if (degree_sequence) {
        stable_sort(vertex_order.begin(), vertex_order.end(), [&](int a, int b) -> bool {
            return -i.degree(a) < -i.degree(b);
        });
    }

    dejavu::hooks::schreier_hook hook(rschreier);
    dejavu::solver s;
    s.set_print(false);
    s.automorphisms(&g, hook.get_hook());       // Compute automorphisms of g

    dejavu::groups::automorphism_workspace a(nv);
    std::vector<std::vector<unsigned int>> mappings;            // Store generators (and identity) in a sensible notation
    std::vector<unsigned int> identity;                      // We need an identity since we're composing elements from two sets
    for (int y = 0; y < i.size(); y++) {
        identity.push_back(y);                      // Every element maps to itself
    }
    mappings.push_back(identity);

    bool stab_trivial = false;
    dejavu::groups::orbit o{nv};      // Orbit structure of size nv
    while (! stab_trivial) {        // While stabiliser is non-trivial
        rschreier.get_stabilizer_orbit(base.size(), o);     // Get orbit partition
        stab_trivial = true;
        for (int v: vertex_order) {        // For each vertex
            if (o.orbit_size(v) > 1) {      // If stabiliser orbit(v) is non-trivial
                stab_trivial = false;
                base.push_back(v);          // Add v to base
            }
        }

        if (! stab_trivial)
            rschreier.set_base(base);       // Update base
    }

    for (int x = 0; x < base.size(); x++) {
        std::vector<int> orb = rschreier.get_fixed_orbit(x);              // Get orbit
        orbit_sizes[base.at(x)] = orb.size();
        std::vector<unsigned int> mapping;
        for (auto v : orb) {
            if (v == base.at(x)) continue;
            if (rschreier.get_transversal_element(x,v,a)) {
                mapping.clear();
                for (int y = 0; y < i.size(); y++) {
                    mapping.push_back(a.p()[y]);            // mapping[i] maps to p*i
                }
                mappings.push_back(mapping);
            }
            else {
                // warn
            }
        }

    }

    return mappings;
}

auto gss::innards::dynamic_coset_reps(std::vector<int> & base, int sz, dejavu::groups::random_schreier &rschreier, std::vector<std::vector<unsigned int>> & reps, std::vector<std::vector<unsigned int>> & invs, std::vector<int> & orbit_sizes) -> void {
    rschreier.set_base(base);       // Empty base to begin with

    dejavu::groups::automorphism_workspace a(sz);       // TODO maybe only need to look at the most recent element
    reps.erase(reps.begin() + 1, reps.end());     // Leave the identity
    invs.erase(invs.begin() + 1, invs.end());
    std::vector<unsigned int> mapping;

    // int x = base.size() - 1;            // Only interested in the most recently added base point
    for (unsigned int x = 0; x < base.size(); x++) {
        vector<int> orb = rschreier.get_fixed_orbit(x);       // Retrieve stabiliser orbit at this point
        orbit_sizes[base.at(x)] = orb.size();
        for (auto v : orb) {
            if (v == base.at(x)) continue;
            if (rschreier.get_transversal_element(x,v,a)) {
                mapping.clear();
                for (int y = 0; y < sz; y++) {
                    mapping.push_back(a.p()[y]);
                }
                invs.push_back(mapping);
                reps.push_back(invert_automorphism(mapping));
            }
            else {
                // warn
            }
        }
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

auto gss::innards::build_static_graph(const InputGraph &i) -> dejavu::static_graph {
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
                g.add_edge(vertices[f], direction_vertices[0]);
                g.add_edge(direction_vertices[0], direction_vertices[1]);
                g.add_edge(vertices[t], direction_vertices[1]);
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

    return g;
}

auto gss::innards::build_static_graph(std::vector<innards::SVOBitset> m, const bool directed) -> dejavu::static_graph {
    dejavu::static_graph g;
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
    return g;  
}