#include <gss/innards/automorphisms.hh>

#include <dejavu.h>

#include <vector>

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


auto gss::innards::dynamic_order_constraints(std::vector<innards::SVOBitset> m, vector<int> base) -> std::vector<std::pair<unsigned int, unsigned int>> {
    
    dejavu::static_graph g;     // Declare static graph
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

    int nv = static_cast<int>(m.size());        // Integer cast to suppress narrowing conversion warning
    
    dejavu::groups::random_schreier rschreier{nv};        // Random Schreier structure with domain size nv

    vector<int> empty_base;

    rschreier.set_base(empty_base);       // Let dejavu start with an empty base

    dejavu::hooks::schreier_hook hook(rschreier);
    dejavu::solver s;
    s.set_print(false);
    s.automorphisms(&g, hook.get_hook());       // Compute automorphisms of g

    std::vector<std::pair<unsigned int, unsigned int>> result;        // List of constraint pairs

    rschreier.set_base(base);           // Recompute the stabiliser chain with the provided (partial) base

    dejavu::groups::orbit o{nv};      // Orbit structure of size nv
    // std::cout << "Base:\n";
    for (int x = 0, x_end = base.size(); x != x_end; ++x) {
        rschreier.get_stabilizer_orbit(x, o);
        int y = base.at(x);
        // std::cout << y << " : ";
        for (int z = 0, z_end = m.size(); z != z_end; ++z) {
            if (o.are_in_same_orbit(y,z) && y != z) {
                // std::cout << z << " ";
                result.push_back(std::make_pair(y,z));
            }
        }
        // std::cout << "\n";
    }

    return result;
}
