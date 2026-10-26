#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_GSS_COUNTING_STATEMENT_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_GSS_COUNTING_STATEMENT_HH 1

#include <gss/loooong.hh>

#include <exception>
#include <iosfwd>
#include <string>
#include <utility>
#include <vector>

namespace gss
{
    /// A statement in the counting language of dev_docs/decomposition-counting.md: a
    /// DAG of reified gates over multi-valued variables, whose count is easy to compute
    /// and which translates to a PB formula with that many solutions over the
    /// variables. Every variable is an exactly-one group, so a literal fixes one value.
    ///
    /// Node ids, variable ids and domain ids are numbered from zero in order of
    /// creation, and a node may only refer to nodes created before it.
    enum class CountingNodeKind
    {
        Cube, ///< one value for each variable in scope; counts 1
        False, ///< nothing, over an explicit scope; counts 0
        Product, ///< children over disjoint scopes; counts the product
        Sum, ///< children over equal scopes that never hold together; counts the sum
        Difference ///< a and not b, where b implies a; counts #a - #b
    };

    struct CountingNode
    {
        CountingNodeKind kind;

        /// For a cube, (variable, value) pairs.
        std::vector<std::pair<int, int>> literals;

        /// For false, the scope.
        std::vector<int> variables;

        /// For product, sum and difference (which has exactly two, a then b).
        std::vector<int> children;
    };

    struct CountingStatement
    {
        std::vector<std::vector<std::string>> domains;
        std::vector<std::pair<std::string, int>> variables;
        std::vector<CountingNode> nodes;
        int root = -1;

        auto add_domain(std::vector<std::string> values) -> int;
        auto add_variable(std::string name, int domain) -> int;
        auto add_cube(std::vector<std::pair<int, int>> literals) -> int;
        auto add_false(std::vector<int> variables) -> int;
        auto add_product(std::vector<int> children) -> int;
        auto add_sum(std::vector<int> children) -> int;
        auto add_difference(int a, int b) -> int;
    };

    /// Thrown for a malformed statement, or one that breaks a side condition that can
    /// be checked syntactically.
    class CountingStatementError : public std::exception
    {
    private:
        std::string _what;

    public:
        explicit CountingStatementError(const std::string & message) noexcept;

        auto what() const noexcept -> const char * override;
    };

    /// A side condition that cannot be checked syntactically, and so needs its own
    /// proof before the count can be believed.
    struct CountingObligation
    {
        enum class Kind
        {
            Disjoint, ///< no two of this sum's children hold together
            Implies ///< this difference's b implies its a
        };

        Kind kind;
        int node;
    };

    struct CountingStatementCount
    {
        loooong count;
        std::vector<CountingObligation> obligations;
    };

    auto read_counting_statement(std::istream &) -> CountingStatement;

    auto write_counting_statement(std::ostream &, const CountingStatement &) -> void;

    /// Checks the syntactic side conditions, throwing CountingStatementError if any
    /// fails, and returns the root's count together with the side conditions that it
    /// could not check. The count is only right if those hold.
    auto count_counting_statement(const CountingStatement &) -> CountingStatementCount;

    /// Writes the PB translation: an exactly-one group for each variable, using the
    /// certified encoding's names (x<variable>_<value>, with @al1 / @am1 labels), a
    /// definition for every node's gate in both directions, and the root as a unit. Its
    /// count projected onto the x variables is the root's count, if the side conditions
    /// hold.
    auto write_counting_statement_opb(std::ostream &, const CountingStatement &) -> void;

    auto counting_statement_gate_name(int node) -> std::string;
}

#endif
