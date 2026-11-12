#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_GSS_COUNTING_STATEMENT_RENDER_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_GSS_COUNTING_STATEMENT_RENDER_HH 1

#include <gss/counting_statement.hh>
#include <gss/loooong.hh>

#include <functional>
#include <iosfwd>
#include <string>
#include <vector>

// Writes a counting statement (dev_docs/counting-statements.md) for people to read, as a
// compressed form of the solutions it counts. None of this is needed to check a count, and
// none of it is part of what a verified checker would have to trust: it only reads a
// statement that read_counting_statement() has already accepted.
//
// A statement's root has count objects, and each object determines an assignment to every
// variable. When the statement was proved for an instance, the objects are its solutions,
// one each.

namespace gss
{
    /// Names for a statement's variables and values, when printing.
    struct CountingStatementNames
    {
        std::function<std::string(int)> variable = [](int u) { return std::to_string(u); };
        std::function<std::string(int)> value = [](int t) { return std::to_string(t); };
    };

    /// The statement, written for people to read, as a list of definitions starting with
    /// the root's. A definition is one or more alternatives, each a set of assignments
    /// (fixed values, and a table of rows for the rest) times some distinct nodes, written
    /// in place, times some named definitions. Its objects are, for each alternative, each
    /// row of its table together with an object of each distinct node and of each name.
    /// Nodes other than distinct nodes that are shared, or that are multiplied by something
    /// else with more than one object, get a name; everything else is multiplied out in
    /// place.
    auto render_counting_statement(std::ostream &, const CountingStatement &, const CountingStatementNames & = {}) -> void;
}

#endif
