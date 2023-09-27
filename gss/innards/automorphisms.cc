#include <gss/innards/automorphisms.hh>

#include <dejavu.h>

#include <vector>

using std::string_view;
using std::vector;

using namespace gss::innards;

auto gss::innards::automorphisms_as_order_constraints(const InputGraph & i) -> OrderConstraints
{
    dejavu::static_graph g;
    g.initialize_graph(i.size(), i.number_of_directed_edges() / 2);
    vector<int> vertices;
    for (int v = 0, v_end = i.size(); v != v_end; ++v)
        vertices.push_back(g.add_vertex(0, i.degree(v)));
    i.for_each_edge([&](int f, int t, string_view) {
        if (f < t)
            g.add_edge(vertices[f], vertices[t]);
    });

    dejavu::groups::random_schreier rschreier{i.size()};
    std::vector<int> base;
    rschreier.set_base(base);

    dejavu::hooks::schreier_hook hook(rschreier);
    dejavu::solver s;
    s.automorphisms(&g, hook.get_hook());

    OrderConstraints result;

    bool stab_trivial = false;
    dejavu::groups::orbit o{i.size()};
    while (! stab_trivial) {
        rschreier.get_stabilizer_orbit(base.size(), o);
        stab_trivial = true;
        for (int v = 0, v_end = i.size(); v != v_end; ++v) {
            if (o.orbit_size(v) > 1) {
                stab_trivial = false;
                base.push_back(v);
            }
        }

        if (! stab_trivial)
            rschreier.set_base(base);
    }

    for (unsigned x = 0; x < base.size(); ++x) {
        auto v = base.at(x);
        for (auto & o : rschreier.get_fixed_orbit(x)) {
            if (o != v)
                result.emplace_back(i.vertex_name(v), i.vertex_name(o));
        }
    }

    return result;
}
