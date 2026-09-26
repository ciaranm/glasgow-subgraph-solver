#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_GSS_INNARDS_COST_BOUND_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_GSS_INNARDS_COST_BOUND_HH 1

#include <gss/innards/homomorphism_domain.hh>
#include <gss/innards/reification.hh>

#include <chrono>
#include <list>
#include <string>
#include <vector>

namespace gss::innards
{
    /**
     * What the cost bound needs to know about an instance: what each target vertex
     * costs, and which pattern and target vertices are edge-vertices standing for an
     * edge (see reification.hh). Vertices numbered below the original sizes are original
     * vertices; the edge-vertices come after them. A graph that was not reified has no
     * edge-vertices, and its original size is its size.
     */
    struct CostData
    {
        std::vector<long long> target_costs;
        int pattern_original_size, target_original_size;
        std::vector<EdgeVertexEndpoints> pattern_edge_vertices, target_edge_vertices;
        bool directed;
    };

    /**
     * A lower bound on the cost of any completion of a partial mapping, which prunes
     * the values of original pattern vertices that cannot be part of a mapping cheaper
     * than the incumbent.
     *
     * The objective is the sum of the costs of the target vertices used, which is a
     * unary cost on each pattern vertex's image. Taken at face value that gives a weak
     * bound, since an unassigned edge-vertex could take the cheapest edge anywhere. So
     * the bound puts the pairwise structure back: each edge-vertex whose endpoints are
     * both unassigned becomes a cost on the pair of their images, read off the values
     * still in its domain, and one with a single endpoint assigned becomes a cost on the
     * other endpoint. That is a problem with unary and pairwise costs over the original
     * pattern vertices. Its local-polytope relaxation is tightened by a few sweeps of dual
     * block-coordinate ascent (in the style of MPLP), and the reparametrised unary costs
     * are then combined by a minimum-cost assignment, which accounts for injectivity and
     * whose reduced costs say which values to remove.
     *
     * Everything is exact integer arithmetic. The dual messages are rounded down, which
     * keeps the reparametrisation exact (any messages give a valid bound, provided the
     * pairwise residuals are computed against the same messages), so no rounding can
     * make the bound exceed the true optimum and prune it. The price is a slightly weaker
     * ascent than with exact halves.
     *
     * The prototype this came from, and the measurements that chose this bound over the
     * alternatives, are on the prototype/weighted-matching branch.
     */
    class CostBound
    {
    private:
        const CostData & _data;
        unsigned _pattern_size, _target_size;
        int _dual_sweeps;
        bool _pruning;

        // From the last successful propagate(): per original pattern vertex, per original
        // target vertex, its reparametrised unary cost.
        std::vector<std::vector<long long>> _scores;

        unsigned long long _calls = 0, _removals = 0;
        std::chrono::nanoseconds _time{0};

        auto propagate_timed(const std::vector<int> & assigned, std::vector<HomomorphismDomain> & domains,
            long long upper_bound, bool & changed) -> bool;

    public:
        /**
         * \param pruning if false, propagate() does nothing, which is what proof logging
         *     needs until the bound's reasoning can be certified; the searcher still uses
         *     the object for the costs of mappings.
         * \throw UnsupportedConfiguration if the costs are large enough that summing them
         *     over a mapping could overflow.
         */
        CostBound(const CostData & data, unsigned pattern_size, unsigned target_size, int dual_sweeps = 5, bool pruning = true);

        /**
         * \param assigned the target vertex of each assigned pattern vertex, or -1.
         * \param domains the domains of the pattern vertices, of which unfixed ones are
         *     the unassigned vertices.
         * \param upper_bound the cost that a mapping must beat, or a value at least
         *     infinity() if there is none yet.
         * \param changed set to true if any value was removed.
         * \return false if no completion can beat upper_bound.
         */
        auto propagate(const std::vector<int> & assigned, std::vector<HomomorphismDomain> & domains,
            long long upper_bound, bool & changed) -> bool;

        /**
         * After a successful propagate(), how promising a value of an unassigned original
         * pattern vertex is: lower is better.
         */
        auto score(unsigned pattern_vertex, unsigned target_vertex) const -> long long;

        /**
         * The cost of a complete mapping.
         */
        auto cost_of(const std::vector<int> & assigned) const -> long long;

        auto is_original_pattern_vertex(unsigned pattern_vertex) const -> bool;

        auto pattern_original_size() const -> unsigned;

        static auto infinity() -> long long;

        /**
         * How often the bound ran, what it removed, and how long it took, for the
         * result's extra stats.
         */
        auto add_extra_stats(std::list<std::string> & stats) const -> void;
    };
}

#endif
