#include <gss/counting_statement_views.hh>

#include <algorithm>
#include <deque>
#include <iterator>
#include <map>
#include <numeric>
#include <ostream>
#include <set>
#include <sstream>
#include <utility>

#include <gmp.h>

using namespace gss;

using std::deque;
using std::map;
using std::ostream;
using std::ostringstream;
using std::pair;
using std::string;
using std::to_string;
using std::unique_ptr;
using std::vector;

namespace
{
    auto unrank(const CountingStatement & statement, int node, loooong index, vector<int> & values) -> void
    {
        auto & n = statement.nodes.at(node);
        switch (n.kind) {
        case CountingNodeKind::Cube:
            for (auto & [u, t] : n.assignment) {
                if (-1 != values.at(u))
                    throw CountingStatementError{"an object assigns variable " + to_string(u) + " twice"};
                values.at(u) = t;
            }
            return;

        case CountingNodeKind::And:
            for (auto & c : n.children) {
                loooong quotient, remainder;
                mpz_fdiv_qr(quotient.raw(), remainder.raw(), index.raw(), statement.nodes.at(c).count.raw());
                unrank(statement, c, remainder, values);
                index = quotient;
            }
            return;

        case CountingNodeKind::Or:
            for (auto & c : n.children) {
                auto & count = statement.nodes.at(c).count;
                if (index < count) {
                    unrank(statement, c, index, values);
                    return;
                }
                index -= count;
            }
            throw CountingStatementError{"index out of range"};
        }
    }

    // One alternative of a definition: a cube, times some named nodes.
    struct Alternative
    {
        vector<pair<int, int>> cube;
        vector<int> names;
    };

    struct Renderer
    {
        const CountingStatement & statement;
        const CountingStatementNames & names;
        vector<bool> named;

        auto count(int node) const -> const loooong &
        {
            return statement.nodes.at(node).count;
        }

        auto name(int node) const -> string
        {
            return node == statement.root ? "root" : "n" + to_string(node);
        }

        // A node's objects, as alternatives: a named node other than the one being defined
        // is left as its name, and everything else is multiplied out. At most one child of
        // an and has more than one alternative, because any child with more than one
        // object beside another such child is named.
        auto expand(int node, bool defining) const -> vector<Alternative>
        {
            if (named[node] && ! defining)
                return {Alternative{{}, {node}}};

            auto & n = statement.nodes.at(node);
            switch (n.kind) {
            case CountingNodeKind::Cube:
                return {Alternative{n.assignment, {}}};

            case CountingNodeKind::Or: {
                vector<Alternative> result;
                for (auto & c : n.children)
                    for (auto & a : expand(c, false))
                        result.push_back(std::move(a));
                return result;
            }

            case CountingNodeKind::And: {
                vector<Alternative> result{Alternative{}};
                for (auto & c : n.children) {
                    auto child = expand(c, false);
                    vector<Alternative> next;
                    for (auto & a : result)
                        for (auto & b : child) {
                            Alternative both = a;
                            both.cube.insert(both.cube.end(), b.cube.begin(), b.cube.end());
                            both.names.insert(both.names.end(), b.names.begin(), b.names.end());
                            next.push_back(std::move(both));
                        }
                    result = std::move(next);
                }
                return result;
            }
            }
            return {};
        }

        // An alternative with the same variables and names as others goes in one table, with
        // a row of values each.
        struct Table
        {
            vector<int> variables;
            vector<int> names;
            vector<vector<int>> rows;
        };

        // The distinct values of the rows on these columns, in order of first appearance.
        static auto projection(const vector<vector<int>> & rows, const vector<unsigned> & columns) -> vector<vector<int>>
        {
            vector<vector<int>> result;
            std::set<vector<int>> seen;
            for (auto & row : rows) {
                vector<int> p;
                for (auto & c : columns)
                    p.push_back(row[c]);
                if (seen.insert(p).second)
                    result.push_back(std::move(p));
            }
            return result;
        }

        // Splits a table's columns into blocks such that the table is exactly the product of
        // its projections onto them. Two columns go in one block if their joint projection
        // is smaller than the product of theirs; the result is used only if the sizes then
        // multiply up to the number of rows, which, with distinct rows, makes it exact.
        static auto factor_table(const vector<vector<int>> & rows, const vector<unsigned> & columns) -> vector<vector<unsigned>>
        {
            if (columns.empty())
                return {};
            auto whole = vector<vector<unsigned>>{columns};
            if (columns.size() < 2 || projection(rows, columns).size() != rows.size())
                return whole;

            vector<unsigned> block_of(columns.size());
            std::iota(block_of.begin(), block_of.end(), 0);
            auto find = [&](unsigned i) {
                while (block_of[i] != i)
                    i = block_of[i] = block_of[block_of[i]];
                return i;
            };
            vector<std::size_t> sizes;
            for (auto & c : columns)
                sizes.push_back(projection(rows, {c}).size());
            for (unsigned i = 0; i < columns.size(); ++i)
                for (unsigned j = i + 1; j < columns.size(); ++j)
                    if (find(i) != find(j) && projection(rows, {columns[i], columns[j]}).size() < sizes[i] * sizes[j])
                        block_of[find(i)] = find(j);

            map<unsigned, vector<unsigned>> blocks;
            for (unsigned i = 0; i < columns.size(); ++i)
                blocks[find(i)].push_back(columns[i]);
            vector<vector<unsigned>> result;
            loooong product{1};
            for (auto & [_, b] : blocks) {
                product *= loooong{(unsigned long)projection(rows, b).size()};
                result.push_back(b);
            }
            if (product != loooong{(unsigned long)rows.size()})
                return whole;
            return result;
        }

        // One block of a table: "{u v} = {3 4}" if its rows are every way of putting its
        // variables at the same set of values, one each; "{u v} distinct in {3 4 5}" if they
        // are every way of putting them at distinct values of a set; and otherwise
        // "(u v) in {(3 4) (4 7)}".
        auto block_text(const Table & table, const vector<unsigned> & block) const -> string
        {
            auto rows = projection(table.rows, block);
            auto list = [&](const vector<int> & items, auto && f) {
                string s;
                for (auto & i : items)
                    s += (s.empty() ? "" : " ") + f(i);
                return s;
            };
            vector<int> variables;
            for (auto & c : block)
                variables.push_back(table.variables[c]);

            if (block.size() >= 2) {
                // the rows are distinct, so if each has distinct values from a set of d,
                // and there are d (d - 1) ... (d - k + 1) of them, they are all such
                std::set<int> values;
                bool distinct = true;
                for (auto & row : rows) {
                    std::set<int> in_row(row.begin(), row.end());
                    if (in_row.size() != row.size())
                        distinct = false;
                    values.insert(row.begin(), row.end());
                }
                loooong ways{1};
                for (unsigned long i = 0; i < block.size(); ++i)
                    ways *= loooong{(unsigned long)(values.size() - i)};
                if (distinct && ways == loooong{(unsigned long)rows.size()}) {
                    vector<int> sorted_values(values.begin(), values.end());
                    return "{" + list(variables, names.variable) + (values.size() == block.size() ? "} = {" : "} distinct in {") + list(sorted_values, names.value) + "}";
                }
            }

            auto tuple = [&](const vector<int> & items, auto && f) {
                return 1 == items.size() ? list(items, f) : "(" + list(items, f) + ")";
            };
            std::sort(rows.begin(), rows.end());
            string text = tuple(variables, names.variable) + " in {";
            for (unsigned r = 0; r < rows.size(); ++r)
                text += (0 == r ? "" : " ") + tuple(rows[r], names.value);
            return text + "}";
        }

        // One table of alternatives, on one line: its fixed values, a table of the rest, and its
        // names, with how many objects it has.
        auto table_text(const Table & table) const -> string
        {
            vector<bool> fixed(table.variables.size(), true);
            for (auto & row : table.rows)
                for (unsigned i = 0; i < row.size(); ++i)
                    if (row[i] != table.rows.front()[i])
                        fixed[i] = false;

            vector<string> parts;
            string fixed_part;
            for (unsigned i = 0; i < table.variables.size(); ++i)
                if (fixed[i])
                    fixed_part += (fixed_part.empty() ? "" : " ") + names.variable(table.variables[i]) + "->" + names.value(table.rows.front()[i]);
            if (! fixed_part.empty())
                parts.push_back("(" + fixed_part + ")");

            vector<unsigned> free;
            for (unsigned i = 0; i < table.variables.size(); ++i)
                if (! fixed[i])
                    free.push_back(i);
            for (auto & block : factor_table(table.rows, free))
                parts.push_back(block_text(table, block));

            loooong objects{(unsigned long)table.rows.size()};
            for (auto & n : table.names) {
                parts.push_back(name(n));
                objects *= count(n);
            }

            if (parts.empty())
                parts.push_back("()");
            string line;
            for (auto & p : parts)
                line += (line.empty() ? "" : " * ") + p;
            ostringstream objects_text;
            objects_text << objects;
            return line + "  [" + objects_text.str() + "]";
        }

        // Writes text, which starts at column start, breaking it at spaces so that lines are
        // at most about width columns, and indenting the continuations. It does not break
        // inside a tuple of values or a set.
        static auto write_wrapped(ostream & out, const string & text, std::size_t start, const string & indent) -> void
        {
            constexpr std::size_t width = 100;
            std::size_t column = start, line_start = 0, last_break = string::npos;
            // in a table's {...}, break only between its tuples, and never inside a set's
            bool in_table = false, in_set = false;
            int depth = 0;
            string result;
            for (std::size_t i = 0; i < text.size(); ++i) {
                char c = text[i];
                if ('{' == c) {
                    if (i >= 3 && 0 == text.compare(i - 3, 3, "in "))
                        in_table = true;
                    else
                        in_set = true;
                }
                else if ('}' == c)
                    in_table = in_set = false;
                else if ('(' == c && in_table)
                    ++depth;
                else if (')' == c && in_table)
                    --depth;
                else if (' ' == c && ! in_set && ! (in_table && depth > 0) && i + 1 < text.size() && ' ' != text[i + 1] && '=' != text[i + 1] && '[' != text[i + 1] && ! (i > 0 && '=' == text[i - 1]))
                    last_break = i;
                ++column;
                if (column > width && string::npos != last_break && last_break > line_start) {
                    result += text.substr(line_start, last_break - line_start) + "\n" + indent;
                    line_start = last_break + 1;
                    column = indent.size() + (i + 1 - line_start);
                    last_break = string::npos;
                }
            }
            result += text.substr(line_start);
            out << result << '\n';
        }

        auto write_definition(ostream & out, int node, deque<int> & queue, vector<bool> & queued) const -> void
        {
            auto alternatives = expand(node, true);
            for (auto & a : alternatives) {
                std::sort(a.cube.begin(), a.cube.end());
                for (auto & n : a.names)
                    if (! queued[n]) {
                        queued[n] = true;
                        queue.push_back(n);
                    }
            }

            // what every alternative assigns alike is written once, in front
            vector<pair<int, int>> common = alternatives.front().cube;
            for (auto & a : alternatives) {
                vector<pair<int, int>> both;
                std::set_intersection(common.begin(), common.end(), a.cube.begin(), a.cube.end(), std::back_inserter(both));
                common = std::move(both);
            }
            if (1 == alternatives.size())
                common.clear();

            vector<Table> tables;
            map<pair<vector<int>, vector<int>>, unsigned> table_of;
            for (auto & a : alternatives) {
                vector<int> variables, row;
                for (auto & assignment : a.cube)
                    if (! std::binary_search(common.begin(), common.end(), assignment)) {
                        variables.push_back(assignment.first);
                        row.push_back(assignment.second);
                    }
                auto key = pair{variables, a.names};
                auto [where, inserted] = table_of.emplace(key, unsigned(tables.size()));
                if (inserted)
                    tables.push_back(Table{std::move(variables), a.names, {}});
                tables[where->second].rows.push_back(std::move(row));
            }

            ostringstream head_stream;
            head_stream << name(node) << " [" << count(node) << "] = ";
            auto head = head_stream.str();
            if (1 == tables.size() && common.empty()) {
                write_wrapped(out, head + table_text(tables.front()), 0, "    ");
                return;
            }

            string common_text;
            for (auto & [u, t] : common)
                common_text += (common_text.empty() ? "" : " ") + names.variable(u) + "->" + names.value(t);
            if (1 == tables.size()) {
                write_wrapped(out, head + "(" + common_text + ") * " + table_text(tables.front()), 0, "    ");
                return;
            }
            write_wrapped(out, head + (common.empty() ? "" : "(" + common_text + ") * ") + "one of", 0, "    ");
            for (auto & table : tables)
                write_wrapped(out, "  | " + table_text(table), 0, "        ");
        }
    };
}

auto gss::unrank_counting_statement(const CountingStatement & statement, const loooong & index) -> vector<int>
{
    if (index < loooong{0} || ! (index < statement.nodes.at(statement.root).count))
        throw CountingStatementError{"index out of range"};
    vector<int> values(statement.domain_sizes.size(), -1);
    unrank(statement, statement.root, index, values);
    for (unsigned u = 0; u < values.size(); ++u)
        if (-1 == values[u])
            throw CountingStatementError{"an object leaves variable " + to_string(u) + " unassigned"};
    return values;
}

struct CountingStatementSampler::Imp
{
    const CountingStatement & statement;
    gmp_randstate_t state;

    explicit Imp(const CountingStatement & s) :
        statement(s)
    {
        gmp_randinit_mt(state);
    }

    ~Imp()
    {
        gmp_randclear(state);
    }
};

CountingStatementSampler::CountingStatementSampler(const CountingStatement & statement, unsigned long seed) :
    _imp(std::make_unique<Imp>(statement))
{
    gmp_randseed_ui(_imp->state, seed);
}

CountingStatementSampler::~CountingStatementSampler() = default;

auto CountingStatementSampler::operator()() -> vector<int>
{
    loooong index;
    mpz_urandomm(index.raw(), _imp->state, _imp->statement.nodes.at(_imp->statement.root).count.raw());
    return unrank_counting_statement(_imp->statement, index);
}

auto gss::write_counting_object(ostream & out, const vector<int> & values, const CountingStatementNames & names) -> void
{
    out << "mapping =";
    for (unsigned u = 0; u < values.size(); ++u)
        out << " (" << names.variable(int(u)) << " -> " << names.value(values[u]) << ")";
    out << '\n';
}

auto gss::render_counting_statement(ostream & out, const CountingStatement & statement, const CountingStatementNames & names) -> void
{
    auto & nodes = statement.nodes;
    auto n = nodes.size();

    // Only what the root reaches is rendered. Children have smaller ids than their parents.
    vector<bool> reached(n, false);
    vector<unsigned> parents(n, 0);
    reached.at(statement.root) = true;
    for (int i = statement.root; i >= 0; --i)
        if (reached[i])
            for (auto & c : nodes[i].children) {
                reached[c] = true;
                ++parents[c];
            }

    // A node with more than one object gets a name if it is shared, or if it is multiplied
    // by another child of an and with more than one object. A node with one object is a
    // single assignment, and is always written out.
    Renderer renderer{statement, names, vector<bool>(n, false)};
    loooong one{1};
    for (unsigned i = 0; i < n; ++i) {
        if (! reached[i] || nodes[i].kind != CountingNodeKind::And)
            continue;
        unsigned many = 0;
        for (auto & c : nodes[i].children)
            if (one < nodes[c].count)
                ++many;
        if (many >= 2)
            for (auto & c : nodes[i].children)
                if (one < nodes[c].count)
                    renderer.named[c] = true;
    }
    for (unsigned i = 0; i < n; ++i)
        if (reached[i] && one < nodes[i].count && parents[i] >= 2)
            renderer.named[i] = true;
    renderer.named.at(statement.root) = true;

    deque<int> queue{statement.root};
    vector<bool> queued(n, false);
    queued.at(statement.root) = true;
    while (! queue.empty()) {
        auto node = queue.front();
        queue.pop_front();
        renderer.write_definition(out, node, queue, queued);
    }
}
