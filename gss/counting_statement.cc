#include <gss/counting_statement.hh>

#include <algorithm>
#include <climits>
#include <istream>
#include <iterator>
#include <optional>
#include <ostream>
#include <string>

using namespace gss;

using std::istream;
using std::istreambuf_iterator;
using std::optional;
using std::ostream;
using std::pair;
using std::string;
using std::to_string;
using std::vector;

CountingStatementError::CountingStatementError(const string & message) noexcept :
    _what(message)
{
}

auto CountingStatementError::what() const noexcept -> const char *
{
    return _what.c_str();
}

CountingStatementWriter::CountingStatementWriter(ostream & out, const vector<int> & domain_sizes) :
    _out(out)
{
    _out << "(\n(version 1)\n(domain";
    for (auto & n : domain_sizes)
        _out << " " << n;
    _out << ")\n(nodes\n";
}

auto CountingStatementWriter::start_node(const char * kind) -> int
{
    _out << "(" << _next_id << " " << kind;
    return _next_id++;
}

auto CountingStatementWriter::add_cube(const vector<pair<int, int>> & assignment) -> int
{
    auto id = start_node("cube");
    for (auto & [u, t] : assignment)
        _out << " (" << u << " " << t << ")";
    _out << ")\n";
    return id;
}

auto CountingStatementWriter::add_and(const vector<int> & children) -> int
{
    auto id = start_node("and");
    for (auto & c : children)
        _out << " " << c;
    _out << ")\n";
    return id;
}

auto CountingStatementWriter::add_or(const vector<int> & children) -> int
{
    auto id = start_node("or");
    for (auto & c : children)
        _out << " " << c;
    _out << ")\n";
    return id;
}

auto CountingStatementWriter::add_minus(int a, int b) -> int
{
    auto id = start_node("minus");
    _out << " " << a << " " << b << ")\n";
    return id;
}

auto CountingStatementWriter::finish(int root) -> void
{
    _out << ")\n(root " << root << ")\n)\n";
}

auto CountingStatementWriter::number_of_nodes() const -> int
{
    return _next_id;
}

namespace
{
    auto is_space(int c) -> bool
    {
        return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\f' || c == '\v';
    }

    // Reads the grammar in dev_docs/counting-statements.md directly off a stream of
    // tokens, without building an s-expression tree first, so that a large statement
    // is held in memory only once.
    class Reader
    {
    private:
        istreambuf_iterator<char> _pos, _end;
        long _line = 1;

        enum class TokenKind
        {
            Open,
            Close,
            Atom,
            End
        };

        TokenKind _kind = TokenKind::End;
        string _atom;

        auto advance() -> void
        {
            while (_pos != _end && is_space(*_pos)) {
                if (*_pos == '\n')
                    ++_line;
                ++_pos;
            }
            _atom.clear();
            if (_pos == _end)
                _kind = TokenKind::End;
            else if (*_pos == '(') {
                _kind = TokenKind::Open;
                ++_pos;
            }
            else if (*_pos == ')') {
                _kind = TokenKind::Close;
                ++_pos;
            }
            else {
                _kind = TokenKind::Atom;
                while (_pos != _end && ! is_space(*_pos) && *_pos != '(' && *_pos != ')')
                    _atom.push_back(*_pos++);
            }
        }

        auto describe() const -> string
        {
            switch (_kind) {
            case TokenKind::Open: return "'('";
            case TokenKind::Close: return "')'";
            case TokenKind::Atom: return "'" + _atom + "'";
            case TokenKind::End: return "end of input";
            }
            return "?";
        }

        [[noreturn]] auto fail(const string & expected) const -> void
        {
            throw CountingStatementError{"line " + to_string(_line) + ": expected " + expected + ", found " + describe()};
        }

    public:
        explicit Reader(istream & in) :
            _pos(in),
            _end()
        {
            advance();
        }

        auto at_open() const -> bool
        {
            return _kind == TokenKind::Open;
        }

        auto at_atom() const -> bool
        {
            return _kind == TokenKind::Atom;
        }

        auto expect_open() -> void
        {
            if (_kind != TokenKind::Open)
                fail("'('");
            advance();
        }

        auto expect_close() -> void
        {
            if (_kind != TokenKind::Close)
                fail("')'");
            advance();
        }

        auto expect_end() -> void
        {
            if (_kind != TokenKind::End)
                fail("end of input");
        }

        auto expect_keyword(const string & keyword) -> void
        {
            if (_kind != TokenKind::Atom || _atom != keyword)
                fail("'" + keyword + "'");
            advance();
        }

        auto read_atom() -> string
        {
            if (_kind != TokenKind::Atom)
                fail("an atom");
            auto result = _atom;
            advance();
            return result;
        }

        auto read_number() -> int
        {
            if (_kind != TokenKind::Atom)
                fail("a number");
            bool ok = (! _atom.empty()) && (_atom == "0" || _atom[0] != '0');
            long long value = 0;
            for (char c : _atom) {
                if (c < '0' || c > '9' || value > INT_MAX) {
                    ok = false;
                    break;
                }
                value = value * 10 + (c - '0');
            }
            if ((! ok) || value > INT_MAX)
                fail("a number");
            advance();
            return int(value);
        }

        auto line() const -> long
        {
            return _line;
        }
    };

    [[noreturn]] auto fail_node(int id, const string & message) -> void
    {
        throw CountingStatementError{"node " + to_string(id) + ": " + message};
    }

    // The scope of every node, as a sorted list of variables, checking well-formedness
    // along the way.
    auto check_and_compute_scopes(const CountingStatement & statement) -> vector<vector<int>>
    {
        for (auto & n : statement.domain_sizes)
            if (n < 0)
                throw CountingStatementError{"negative domain size"};
        int number_of_variables = statement.domain_sizes.size();

        vector<vector<int>> scopes;
        scopes.reserve(statement.nodes.size());
        for (int id = 0; id < int(statement.nodes.size()); ++id) {
            auto & node = statement.nodes[id];
            for (auto & c : node.children)
                if (c < 0 || c >= id)
                    fail_node(id, "child " + to_string(c) + " is not an earlier node");

            vector<int> scope;
            switch (node.kind) {
            case CountingNodeKind::Cube:
                if (! node.children.empty())
                    fail_node(id, "a cube has no children");
                for (auto & [u, t] : node.assignment) {
                    if (u < 0 || u >= number_of_variables)
                        fail_node(id, "variable " + to_string(u) + " out of range");
                    if (t < 0 || t >= statement.domain_sizes[u])
                        fail_node(id, "value " + to_string(t) + " out of range for variable " + to_string(u));
                    if ((! scope.empty()) && scope.back() >= u)
                        fail_node(id, "cube variables are not strictly increasing");
                    scope.push_back(u);
                }
                break;

            case CountingNodeKind::And:
                for (auto & c : node.children) {
                    vector<int> merged;
                    std::merge(scope.begin(), scope.end(), scopes[c].begin(), scopes[c].end(), std::back_inserter(merged));
                    if (std::adjacent_find(merged.begin(), merged.end()) != merged.end())
                        fail_node(id, "children of an and have overlapping scopes");
                    scope = std::move(merged);
                }
                break;

            case CountingNodeKind::Or:
                if (node.children.empty())
                    fail_node(id, "an or has no children");
                for (auto & c : node.children)
                    if (scopes[c] != scopes[node.children.front()])
                        fail_node(id, "children of an or have different scopes");
                scope = scopes[node.children.front()];
                break;

            case CountingNodeKind::Minus:
                if (node.children.size() != 2)
                    fail_node(id, "a minus has exactly two children");
                if (scopes[node.children[0]] != scopes[node.children[1]])
                    fail_node(id, "children of a minus have different scopes");
                scope = scopes[node.children[0]];
                break;
            }

            if (node.kind != CountingNodeKind::Cube && ! node.assignment.empty())
                fail_node(id, "only a cube has an assignment");

            scopes.push_back(std::move(scope));
        }

        if (statement.root < 0 || statement.root >= int(statement.nodes.size()))
            throw CountingStatementError{"root " + to_string(statement.root) + " is not a node"};
        auto & root_scope = scopes[statement.root];
        bool covers_everything = int(root_scope.size()) == number_of_variables;
        for (int u = 0; covers_everything && u < number_of_variables; ++u)
            covers_everything = root_scope[u] == u;
        if (! covers_everything)
            throw CountingStatementError{"the root's scope is not every variable"};

        return scopes;
    }

    // An or of cubes is disjoint if no two of its children are the same cube: they have
    // equal scopes, so different cubes give some variable different values.
    auto syntactically_disjoint(const CountingStatement & statement, const CountingNode & node) -> bool
    {
        vector<const vector<pair<int, int>> *> cubes;
        for (auto & c : node.children) {
            if (statement.nodes[c].kind != CountingNodeKind::Cube)
                return false;
            cubes.push_back(&statement.nodes[c].assignment);
        }
        std::sort(cubes.begin(), cubes.end(), [](auto a, auto b) { return *a < *b; });
        return cubes.end() == std::adjacent_find(cubes.begin(), cubes.end(), [](auto a, auto b) { return *a == *b; });
    }

    auto y(int u, int t) -> string
    {
        return "y" + to_string(u) + "_" + to_string(t);
    }

    auto g(int i) -> string
    {
        return "g" + to_string(i);
    }

    auto write_definitions(ostream & out, const CountingStatement & statement) -> void
    {
        for (int u = 0; u < int(statement.domain_sizes.size()); ++u) {
            for (int t = 0; t < statement.domain_sizes[u]; ++t)
                out << "1 " << y(u, t) << " ";
            out << ">= 1 ;\n";
            for (int t = 0; t < statement.domain_sizes[u]; ++t)
                out << "1 " << y(u, t) << " ";
            out << "<= 1 ;\n";
        }

        for (int i = 0; i < int(statement.nodes.size()); ++i) {
            auto & node = statement.nodes[i];
            switch (node.kind) {
            case CountingNodeKind::Cube:
            case CountingNodeKind::And: {
                vector<string> literals;
                if (node.kind == CountingNodeKind::Cube)
                    for (auto & [u, t] : node.assignment)
                        literals.push_back(y(u, t));
                else
                    for (auto & c : node.children)
                        literals.push_back(g(c));

                auto k = literals.size();
                if (0 == k)
                    out << "1 " << g(i) << " >= 1 ;\n";
                else {
                    out << k << " ~" << g(i);
                    for (auto & l : literals)
                        out << " 1 " << l;
                    out << " >= " << k << " ;\n";
                    out << "1 " << g(i);
                    for (auto & l : literals)
                        out << " 1 ~" << l;
                    out << " >= 1 ;\n";
                }
            } break;

            case CountingNodeKind::Or: {
                auto k = node.children.size();
                out << "1 ~" << g(i);
                for (auto & c : node.children)
                    out << " 1 " << g(c);
                out << " >= 1 ;\n";
                out << k << " " << g(i);
                for (auto & c : node.children)
                    out << " 1 ~" << g(c);
                out << " >= " << k << " ;\n";
            } break;

            case CountingNodeKind::Minus: {
                auto a = node.children[0], b = node.children[1];
                out << "1 ~" << g(i) << " 1 " << g(a) << " >= 1 ;\n";
                out << "1 ~" << g(i) << " 1 ~" << g(b) << " >= 1 ;\n";
                out << "1 " << g(i) << " 1 ~" << g(a) << " 1 " << g(b) << " >= 1 ;\n";
            } break;
            }
        }
    }
}

auto gss::read_counting_statement(istream & in) -> CountingStatement
{
    Reader reader{in};
    CountingStatement result;

    reader.expect_open();

    reader.expect_open();
    reader.expect_keyword("version");
    reader.expect_keyword("1");
    reader.expect_close();

    reader.expect_open();
    reader.expect_keyword("domain");
    while (reader.at_atom())
        result.domain_sizes.push_back(reader.read_number());
    reader.expect_close();

    reader.expect_open();
    reader.expect_keyword("nodes");
    while (reader.at_open()) {
        reader.expect_open();
        auto line = reader.line();
        auto id = reader.read_number();
        if (id != int(result.nodes.size()))
            throw CountingStatementError{"line " + to_string(line) + ": expected node " + to_string(result.nodes.size()) + ", found node " + to_string(id)};

        auto kind = reader.read_atom();
        CountingNode node;
        if (kind == "cube") {
            node.kind = CountingNodeKind::Cube;
            while (reader.at_open()) {
                reader.expect_open();
                auto u = reader.read_number();
                auto t = reader.read_number();
                reader.expect_close();
                node.assignment.emplace_back(u, t);
            }
        }
        else if (kind == "and" || kind == "or" || kind == "minus") {
            node.kind = kind == "and" ? CountingNodeKind::And : kind == "or" ? CountingNodeKind::Or
                                                                             : CountingNodeKind::Minus;
            while (reader.at_atom())
                node.children.push_back(reader.read_number());
        }
        else
            throw CountingStatementError{"line " + to_string(line) + ": unknown node kind '" + kind + "'"};

        reader.expect_close();
        result.nodes.push_back(std::move(node));
    }
    reader.expect_close();

    reader.expect_open();
    reader.expect_keyword("root");
    result.root = reader.read_number();
    reader.expect_close();

    reader.expect_close();
    reader.expect_end();

    check_counting_statement(result);
    return result;
}

auto gss::check_counting_statement(const CountingStatement & statement) -> void
{
    check_and_compute_scopes(statement);
}

auto gss::counting_side_conditions(const CountingStatement & statement) -> vector<int>
{
    vector<int> result;
    for (int i = 0; i < int(statement.nodes.size()); ++i) {
        auto & node = statement.nodes[i];
        if (node.kind == CountingNodeKind::Minus || (node.kind == CountingNodeKind::Or && ! syntactically_disjoint(statement, node)))
            result.push_back(i);
    }
    return result;
}

auto gss::count_counting_nodes(const CountingStatement & statement) -> vector<loooong>
{
    vector<loooong> counts;
    counts.reserve(statement.nodes.size());
    for (int i = 0; i < int(statement.nodes.size()); ++i) {
        auto & node = statement.nodes[i];
        switch (node.kind) {
        case CountingNodeKind::Cube:
            counts.emplace_back(1);
            break;
        case CountingNodeKind::And: {
            loooong r = 1;
            for (auto & c : node.children)
                r *= counts[c];
            counts.push_back(r);
        } break;
        case CountingNodeKind::Or: {
            loooong r = 0;
            for (auto & c : node.children)
                r += counts[c];
            counts.push_back(r);
        } break;
        case CountingNodeKind::Minus: {
            auto r = counts[node.children[0]] - counts[node.children[1]];
            if (r < loooong{0})
                fail_node(i, "minus is negative, so its second child cannot imply its first");
            counts.push_back(r);
        } break;
        }
    }
    return counts;
}

auto gss::count_counting_statement(const CountingStatement & statement) -> loooong
{
    return count_counting_nodes(statement).at(statement.root);
}

auto gss::write_counting_opb(ostream & out, const CountingStatement & statement) -> void
{
    out << "preserved:";
    for (int u = 0; u < int(statement.domain_sizes.size()); ++u)
        for (int t = 0; t < statement.domain_sizes[u]; ++t)
            out << " " << y(u, t);
    out << " ;\n";

    write_definitions(out, statement);

    out << "1 " << g(statement.root) << " >= 1 ;\n";
}

auto gss::write_side_condition_opb(ostream & out, const CountingStatement & statement) -> bool
{
    auto conditions = counting_side_conditions(statement);
    if (conditions.empty())
        return false;

    write_definitions(out, statement);

    for (auto & i : conditions) {
        auto & node = statement.nodes[i];
        auto v = "v" + to_string(i);
        if (node.kind == CountingNodeKind::Or) {
            out << "2 ~" << v;
            for (auto & c : node.children)
                out << " 1 " << g(c);
            out << " >= 2 ;\n";
        }
        else {
            out << "1 ~" << v << " 1 " << g(node.children[1]) << " >= 1 ;\n";
            out << "1 ~" << v << " 1 ~" << g(node.children[0]) << " >= 1 ;\n";
        }
    }

    for (auto & i : conditions)
        out << "1 v" << i << " ";
    out << ">= 1 ;\n";
    return true;
}
