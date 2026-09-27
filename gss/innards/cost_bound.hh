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
     * Why the last call of CostBound::propagate() failed or removed values, as multipliers
     * on the constraints of the proof's model (HomomorphismProofs::emit_reified_model).
     * Adding up the objective-improving constraint and, with these multipliers, each
     * original vertex's exactly-one, each target vertex's at-most-one, and the linking
     * equalities, gives sum(-reduced cost * variable) >= bound - (upper bound - 1): at the
     * node that conflicts, or propagates every value the bound removed. Values removed for
     * lack of support need no certificate, since propagation over the linking equalities
     * finds them.
     */
    struct CostBoundCertificate
    {
        // Whether there is anything to certify, and if so, what: a bound reaching the
        // incumbent (or removing values), or an assignment step with no injective
        // assignment at all, which is a Hall violator.
        enum class Kind
        {
            None,
            Bound,
            HallViolator
        };
        Kind kind = Kind::None;

        long long bound = 0;

        // Multipliers, of any sign, on sum_t x(u, t) = 1 for original pattern vertex u.
        std::vector<std::pair<unsigned, long long>> exactly_one;
        // Non-negative multipliers on sum_u x(u, t) <= 1 for original target vertex t.
        std::vector<std::pair<unsigned, long long>> at_most_one;

        // Multipliers, of any sign, on sum_y z(a, b, x, y) = x(a, x) (on_a) or on
        // sum_x z(a, b, x, y) = x(b, y) (not on_a), for a < b.
        struct Link
        {
            unsigned a, b;
            bool on_a;
            unsigned value;
            long long multiplier;
        };
        std::vector<Link> links;

        // For a Hall violator: original pattern vertices whose candidates, together, are
        // fewer than they are.
        std::vector<unsigned> hall_rows, hall_columns;

        // Every (pattern vertex, target vertex) the call removed, whatever the reason.
        std::vector<std::pair<unsigned, unsigned>> removed;
    };

    /**
     * A lower bound on the cost of any completion of a partial mapping, which prunes
     * the values of original pattern vertices that cannot be part of a mapping cheaper
     * than the incumbent.
     *
     * The objective is the sum of the costs of the target vertices used, which is a
     * unary cost on each pattern vertex's image. Taken at face value that gives a weak
     * bound, since an unassigned edge-vertex could take the cheapest edge anywhere. So
     * the bound puts the pairwise structure back. Every original pattern vertex is a row,
     * with its current candidates (just its value, once it has one), and every pair of
     * adjacent original vertices is a pairwise cost on their candidates, the sum over the
     * edge-vertices between them of the cost of the one each would use, read off the
     * values still in its domain. A pair with a single candidate on one side folds its
     * costs into the other side exactly; the rest have their local-polytope relaxation
     * tightened by a few sweeps of dual block-coordinate ascent (in the style of MPLP).
     * The reparametrised unary costs are then combined by a minimum-cost assignment,
     * which accounts for injectivity, solved as a square problem padded with zero-cost
     * rows so that its dual is exact, and whose reduced costs say which values to remove.
     * The bound used is exactly the value of the dual, so that the multipliers it records
     * (see CostBoundCertificate) prove it.
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

        // From the last successful propagate(): per original pattern vertex, per original
        // target vertex, its reparametrised unary cost.
        std::vector<std::vector<long long>> _scores;

        unsigned long long _calls = 0, _removals = 0;
        std::chrono::nanoseconds _time{0};

        bool _want_certificates = false;
        CostBoundCertificate _certificate;

        auto propagate_timed(const std::vector<int> & assigned, std::vector<HomomorphismDomain> & domains,
            long long upper_bound, bool & changed) -> bool;

    public:
        /**
         * \throw UnsupportedConfiguration if the costs are large enough that summing them
         *     over a mapping could overflow.
         */
        CostBound(const CostData & data, unsigned pattern_size, unsigned target_size, int dual_sweeps = 5);

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

        /**
         * Record a CostBoundCertificate for each call from now on, for proof logging.
         */
        auto want_certificates() -> void;

        /**
         * After a call of propagate() that failed or removed values, why.
         */
        auto certificate() const -> const CostBoundCertificate &;

        static auto infinity() -> long long;

        /**
         * How often the bound ran, what it removed, and how long it took, for the
         * result's extra stats.
         */
        auto add_extra_stats(std::list<std::string> & stats) const -> void;
    };
}

#endif
