#include <gss/innards/automorphisms.hh>

#include <dejavu.h>

#include <vector>

using std::string_view;
using std::vector;

using namespace gss::innards;

auto gss::innards::automorphisms_as_order_constraints(const InputGraph & i) -> OrderConstraints {
    std::vector<int> base;

    return automorphisms_as_order_constraints(i, base);
}

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

    std::cout << "Generators:\n";
    dejavu::groups::automorphism_workspace a(i.size());
    std::set<std::set<std::vector<int>>> gens;
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
        for (auto cycle : gen) {
            std::cout << "( ";
            std::cout << cycle[0] << " ";
            for (unsigned int v_index = 1; v_index < cycle.size(); v_index++) {
                std::cout << cycle[v_index] << " ";
                unique_list.emplace(cycle[0], cycle[v_index]);
            }
            std::cout << ")";
        }
        std::cout << "\n";
    }

    std::cout << "Base:\n";
    for (unsigned x = 0; x < base.size(); ++x) {        // For each base element
        auto v = base.at(x);
        std::cout << v << " : ";
        std::vector<int> orb = rschreier.get_fixed_orbit(x);
        for (unsigned int o = 0; o < orb.size(); o++) {     // For each orbit point of v
            std::cout << orb[o] << " ";
            for (unsigned int o2 = o + 1; o2 < orb.size(); o2++) {
                if (orb[o] != orb[o2])     // If o is mapped to somewhere else
                    unique_list.emplace(std::min(orb[o], orb[o2]), std::max(orb[o], orb[o2]));        // Enforce o<o2
            }
        }
        std::cout << "\n";
    }

    for (std::pair<int, int> p: unique_list) {
        result.emplace_back(i.vertex_name(p.first), i.vertex_name(p.second));
    }
    
    return result;
}


// auto gss::innards::dynamic_order_constraints(dejavu::static_graph g, vector<int> base, int size) -> std::vector<std::pair<unsigned int, unsigned int>> {
//     dejavu::groups::random_schreier rschreier{size};        // Random Schreier structure with domain size nv
//     // std::vector<int> base;

//     rschreier.set_base(base);       // Empty base to begin with

//     dejavu::hooks::schreier_hook hook(rschreier);
//     dejavu::solver s;
//     s.automorphisms(&g, hook.get_hook());       // Compute automorphisms of g

//     std::vector<std::pair<unsigned int, unsigned int>> result;        // List of constraint pairs

//     bool stab_trivial = false;
//     dejavu::groups::orbit o{size};      // Orbit structure of size nv
//     while (! stab_trivial) {        // While stabiliser is non-trivial
//         rschreier.get_stabilizer_orbit(base.size(), o);     // Get orbit partition
//         stab_trivial = true;
//         for (int v = 0, v_end = size; v != v_end; ++v) {        // For each vertex
//             if (o.orbit_size(v) > 1) {      // If stabiliser orbit(v) is non-trivial
//                 stab_trivial = false;
//                 base.push_back(v);          // Add v to base
//             }
//         }

//         if (! stab_trivial)
//             rschreier.set_base(base);       // Update base
//     }

//     std::cout << "Base:\n";
//     for (unsigned x = 0; x < base.size(); ++x) {        // For each base element
//         auto v = base.at(x);
//         std::cout << v << " : ";
//         for (auto & o : rschreier.get_fixed_orbit(x)) {     // For each orbit point of v
//             std::cout << o << " ";
//             if (o != v)     // If v is mapped to somewhere else
//                 result.push_back(std::make_pair(v, o));        // Enforce v<o
//         }
//         std::cout << "\n";
//     }

//     return result;
// }
