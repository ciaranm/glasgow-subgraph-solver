#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_GSS_COUNTING_STATEMENT_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_GSS_COUNTING_STATEMENT_HH 1

#include <gss/loooong.hh>

#include <exception>
#include <functional>
#include <iosfwd>
#include <string>
#include <utility>
#include <vector>

// Counting statements: the file format in which a counting search states its argument, and
// its translation into a pseudo-Boolean formula. The format, what it means, and the exact
// translation are specified in dev_docs/counting-statements.md, which is the reference;
// this is one implementation of it.

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
        Or
    };

    struct CountingNode
    {
        CountingNodeKind kind;

        /// The number of the node's objects.
        loooong count;

        /// For a cube, its (variable, value) pairs, in increasing order of variable.
        std::vector<std::pair<int, int>> assignment;

        /// For anything else, its children's ids, in increasing order.
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
    /// memory, working out each node's count as it goes. Each add_ function returns the
    /// new node's id. The children of an and or an or must be given in strictly
    /// increasing order, and anything else is only checked when the result is read back
    /// with read_counting_statement().
    class CountingStatementWriter
    {
    private:
        std::ostream & _out;
        std::vector<loooong> _counts;
        std::function<void(int, const CountingNode &)> _observer;

        auto add(CountingNode node) -> int;

    public:
        CountingStatementWriter(std::ostream & out, const std::vector<int> & domain_sizes);

        auto add_cube(const std::vector<std::pair<int, int>> & assignment) -> int;
        auto add_and(const std::vector<int> & children) -> int;
        auto add_or(const std::vector<int> & children) -> int;

        /// Called with each node's id and contents as soon as it is written.
        auto set_observer(std::function<void(int, const CountingNode &)>) -> void;

        /// Closes the document. Must be called exactly once, last.
        auto finish(int root) -> void;

        auto number_of_nodes() const -> int;

        /// A node's count, as written.
        auto count(int node) const -> const loooong &;
    };

    /// Parses a statement and checks that it is well-formed, counts included, throwing
    /// CountingStatementError if it is not.
    auto read_counting_statement(std::istream &) -> CountingStatement;

    /// Checks that a statement is well-formed, counts included, throwing
    /// CountingStatementError if not.
    auto check_counting_statement(const CountingStatement &) -> void;

    /// The names of the tagged formula's variables: variable u takes value t; node i is
    /// active; or node i takes its child in position j.
    auto counting_variable_name(int u, int t) -> std::string;
    auto counting_activity_name(int i) -> std::string;
    auto counting_tag_name(int i, int j) -> std::string;

    /// For each node, its selectors (see dev_docs/counting-statements.md): one variable
    /// per parent, in order of parent.
    auto counting_selectors(const CountingStatement &) -> std::vector<std::vector<std::string>>;

    /// The tagged formula, as an OPB file. It mentions only the statement's own
    /// variables, so a proof that it is equienumerable with some input has to move the
    /// input's preserved set over to them. Its solutions, projected onto its preserved
    /// set, correspond one-to-one with the root's objects, so there are as many of them as
    /// the root's count.
    auto write_tagged_opb(std::ostream &, const CountingStatement &) -> void;
}

#endif
