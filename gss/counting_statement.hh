#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_GSS_COUNTING_STATEMENT_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_GSS_COUNTING_STATEMENT_HH 1

#include <gss/loooong.hh>

#include <exception>
#include <iosfwd>
#include <string>
#include <utility>
#include <vector>

// Counting statements: the file format in which --decomposition states its counting
// argument, and its translation into a pseudo-Boolean formula. The format, what it
// means, and the exact translation are specified in dev_docs/counting-statements.md,
// which is the reference; this is one implementation of it.

namespace gss
{
    class CountingStatementError : public std::exception
    {
    private:
        std::string _what;

    public:
        explicit CountingStatementError(const std::string & message) noexcept;

        auto what() const noexcept -> const char * override;
    };

    enum class CountingNodeKind
    {
        Cube,
        And,
        Or,
        Minus
    };

    struct CountingNode
    {
        CountingNodeKind kind;

        /// For a cube, its (variable, value) pairs, in increasing order of variable.
        std::vector<std::pair<int, int>> assignment;

        /// For anything else, its children's ids.
        std::vector<int> children;
    };

    struct CountingStatement
    {
        /// Variable u takes one of the values 0 to domain_sizes[u] - 1.
        std::vector<int> domain_sizes;
        std::vector<CountingNode> nodes;
        int root = 0;
    };

    /// Writes a statement as its nodes are created, so that it need not be held in
    /// memory. Each add_ function returns the new node's id. Nothing is checked here:
    /// read the result back with read_counting_statement() to check it.
    class CountingStatementWriter
    {
    private:
        std::ostream & _out;
        int _next_id = 0;

        auto start_node(const char * kind) -> int;

    public:
        CountingStatementWriter(std::ostream & out, const std::vector<int> & domain_sizes);

        auto add_cube(const std::vector<std::pair<int, int>> & assignment) -> int;
        auto add_and(const std::vector<int> & children) -> int;
        auto add_or(const std::vector<int> & children) -> int;
        auto add_minus(int a, int b) -> int;

        /// Closes the document. Must be called exactly once, last.
        auto finish(int root) -> void;

        auto number_of_nodes() const -> int;
    };

    /// Parses a statement and checks that it is well-formed, throwing
    /// CountingStatementError if it is not.
    auto read_counting_statement(std::istream &) -> CountingStatement;

    /// Checks that a statement is well-formed, throwing CountingStatementError if not.
    auto check_counting_statement(const CountingStatement &) -> void;

    /// The nodes whose side conditions cannot be discharged syntactically, so must be
    /// proved by refuting the side-condition instance. The statement must be well-formed.
    auto counting_side_conditions(const CountingStatement &) -> std::vector<int>;

    /// The count of every node, by id. The statement must be well-formed. Throws
    /// CountingStatementError if a minus would go negative, which means that one of its
    /// side conditions does not hold.
    auto count_counting_nodes(const CountingStatement &) -> std::vector<loooong>;

    /// The count of the root. Only right if the side conditions hold.
    auto count_counting_statement(const CountingStatement &) -> loooong;

    /// The counting formula, as an OPB file. It mentions only the statement's own
    /// variables, so a proof that it is equienumerable with some input has to move the
    /// input's preserved set over to them.
    auto write_counting_opb(std::ostream &, const CountingStatement &) -> void;

    /// The side-condition instance, as an OPB file, which must be unsatisfiable. Writes
    /// nothing and returns false if there are no side conditions left to prove.
    auto write_side_condition_opb(std::ostream &, const CountingStatement &) -> bool;
}

#endif
