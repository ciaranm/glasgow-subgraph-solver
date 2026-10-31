#include <gss/counting_statement.hh>

#include <algorithm>
#include <climits>
#include <functional>
#include <istream>
#include <iterator>
#include <ostream>
#include <sstream>
#include <string>

using namespace gss;

using std::function;
using std::istream;
using std::istreambuf_iterator;
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
    _out << "(\n(version 2)\n(domain";
    for (auto & n : domain_sizes)
        _out << " " << n;
    _out << ")\n(nodes\n";
}

auto CountingStatementWriter::add(CountingNode node) -> int
{
    int id = _counts.size();
    for (unsigned i = 0; i < node.children.size(); ++i)
        if (node.children[i] < 0 || node.children[i] >= id || (i > 0 && node.children[i - 1] >= node.children[i]))
            throw CountingStatementError{"node " + to_string(id) + ": children are not earlier nodes in strictly increasing order"};

    switch (node.kind) {
    case CountingNodeKind::Cube:
        node.count = loooong{1};
        break;
    case CountingNodeKind::And:
        node.count = loooong{1};
        for (auto & c : node.children)
            node.count *= _counts[c];
        break;
    case CountingNodeKind::Or:
        node.count = loooong{0};
        for (auto & c : node.children)
            node.count += _counts[c];
        break;
    }

    _out << "(" << id << " " << node.count << " ";
    switch (node.kind) {
    case CountingNodeKind::Cube:
        _out << "cube";
        for (auto & [u, t] : node.assignment)
            _out << " (" << u << " " << t << ")";
        break;
    case CountingNodeKind::And:
    case CountingNodeKind::Or:
        _out << (node.kind == CountingNodeKind::And ? "and" : "or");
        for (auto & c : node.children)
            _out << " " << c;
        break;
    }
    _out << ")\n";

    _counts.push_back(node.count);
    if (_observer)
        _observer(id, node);
    return id;
}

auto CountingStatementWriter::add_cube(const vector<pair<int, int>> & assignment) -> int
{
    return add(CountingNode{CountingNodeKind::Cube, 0, assignment, {}});
}

auto CountingStatementWriter::add_and(const vector<int> & children) -> int
{
    return add(CountingNode{CountingNodeKind::And, 0, {}, children});
}

auto CountingStatementWriter::add_or(const vector<int> & children) -> int
{
    return add(CountingNode{CountingNodeKind::Or, 0, {}, children});
}

auto CountingStatementWriter::set_observer(function<void(int, const CountingNode &)> observer) -> void
{
    _observer = std::move(observer);
}

auto CountingStatementWriter::finish(int root) -> void
{
    _out << ")\n(root " << root << ")\n)\n";
}

auto CountingStatementWriter::number_of_nodes() const -> int
{
    return _counts.size();
}

auto CountingStatementWriter::count(int node) const -> const loooong &
{
    return _counts.at(node);
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

        // A number of any size: a count.
        auto read_count() -> loooong
        {
            if (_kind != TokenKind::Atom)
                fail("a count");
            bool ok = (! _atom.empty()) && (_atom == "0" || _atom[0] != '0');
            for (char c : _atom)
                if (c < '0' || c > '9')
                    ok = false;
            if (! ok)
                fail("a count");
            loooong result;
            mpz_set_str(result.raw(), _atom.c_str(), 10);
            advance();
            return result;
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

    auto y(int u, int t) -> string
    {
        return "y" + to_string(u) + "_" + to_string(t);
    }

    auto a(int i) -> string
    {
        return "a" + to_string(i);
    }

    auto t(int i, int j) -> string
    {
        return "t" + to_string(i) + "_" + to_string(j);
    }
}

auto gss::counting_variable_name(int u, int value) -> string
{
    return y(u, value);
}

auto gss::counting_activity_name(int i) -> string
{
    return a(i);
}

auto gss::counting_tag_name(int i, int j) -> string
{
    return t(i, j);
}

auto gss::check_counting_statement(const CountingStatement & statement) -> void
{
    for (auto & n : statement.domain_sizes)
        if (n < 0)
            throw CountingStatementError{"negative domain size"};
    int number_of_variables = statement.domain_sizes.size();

    // the scope of every node, as a sorted list of variables
    vector<vector<int>> scopes;
    scopes.reserve(statement.nodes.size());
    for (int id = 0; id < int(statement.nodes.size()); ++id) {
        auto & node = statement.nodes[id];
        for (unsigned i = 0; i < node.children.size(); ++i) {
            auto c = node.children[i];
            if (c < 0 || c >= id)
                fail_node(id, "child " + to_string(c) + " is not an earlier node");
            if (i > 0 && node.children[i - 1] >= c)
                fail_node(id, "children are not strictly increasing");
        }

        vector<int> scope;
        loooong count;
        switch (node.kind) {
        case CountingNodeKind::Cube:
            if (! node.children.empty())
                fail_node(id, "a cube has no children");
            for (auto & [u, value] : node.assignment) {
                if (u < 0 || u >= number_of_variables)
                    fail_node(id, "variable " + to_string(u) + " out of range");
                if (value < 0 || value >= statement.domain_sizes[u])
                    fail_node(id, "value " + to_string(value) + " out of range for variable " + to_string(u));
                if ((! scope.empty()) && scope.back() >= u)
                    fail_node(id, "cube variables are not strictly increasing");
                scope.push_back(u);
            }
            count = loooong{1};
            break;

        case CountingNodeKind::And:
            count = loooong{1};
            for (auto & c : node.children) {
                vector<int> merged;
                std::merge(scope.begin(), scope.end(), scopes[c].begin(), scopes[c].end(), std::back_inserter(merged));
                if (std::adjacent_find(merged.begin(), merged.end()) != merged.end())
                    fail_node(id, "children of an and have overlapping scopes");
                scope = std::move(merged);
                count *= statement.nodes[c].count;
            }
            break;

        case CountingNodeKind::Or:
            if (node.children.empty())
                fail_node(id, "an or has no children");
            count = loooong{0};
            for (auto & c : node.children) {
                if (scopes[c] != scopes[node.children.front()])
                    fail_node(id, "children of an or have different scopes");
                count += statement.nodes[c].count;
            }
            scope = scopes[node.children.front()];
            if (scope.empty())
                fail_node(id, "an or has an empty scope");
            break;
        }

        if (node.kind != CountingNodeKind::Cube && ! node.assignment.empty())
            fail_node(id, "only a cube has an assignment");
        if (count != node.count) {
            std::ostringstream message;
            message << "its count is " << count << ", not " << node.count;
            fail_node(id, message.str());
        }

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
}

auto gss::read_counting_statement(istream & in) -> CountingStatement
{
    Reader reader{in};
    CountingStatement result;

    reader.expect_open();

    reader.expect_open();
    reader.expect_keyword("version");
    reader.expect_keyword("2");
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

        CountingNode node{CountingNodeKind::Cube, reader.read_count(), {}, {}};
        auto kind = reader.read_atom();
        if (kind == "cube") {
            while (reader.at_open()) {
                reader.expect_open();
                auto u = reader.read_number();
                auto value = reader.read_number();
                reader.expect_close();
                node.assignment.emplace_back(u, value);
            }
        }
        else if (kind == "and" || kind == "or") {
            node.kind = kind == "and" ? CountingNodeKind::And : CountingNodeKind::Or;
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

auto gss::counting_selectors(const CountingStatement & statement) -> vector<vector<string>>
{
    vector<vector<string>> result(statement.nodes.size());
    for (int p = 0; p < int(statement.nodes.size()); ++p) {
        auto & node = statement.nodes[p];
        for (unsigned j = 0; j < node.children.size(); ++j)
            result[node.children[j]].push_back(node.kind == CountingNodeKind::And ? a(p) : t(p, j));
    }
    return result;
}

auto gss::write_tagged_opb(ostream & out, const CountingStatement & statement) -> void
{
    out << "preserved:";
    for (int u = 0; u < int(statement.domain_sizes.size()); ++u)
        for (int value = 0; value < statement.domain_sizes[u]; ++value)
            out << " " << y(u, value);
    for (int i = 0; i < int(statement.nodes.size()); ++i)
        if (statement.nodes[i].kind == CountingNodeKind::Or)
            for (unsigned j = 0; j < statement.nodes[i].children.size(); ++j)
                out << " " << t(i, j);
    out << " ;\n";

    for (int u = 0; u < int(statement.domain_sizes.size()); ++u) {
        for (int value = 0; value < statement.domain_sizes[u]; ++value)
            out << "1 " << y(u, value) << " ";
        out << ">= 1 ;\n";
        for (int value = 0; value < statement.domain_sizes[u]; ++value)
            out << "1 " << y(u, value) << " ";
        out << "<= 1 ;\n";
    }

    out << "1 " << a(statement.root) << " >= 1 ;\n";

    auto selectors = counting_selectors(statement);
    for (int i = 0; i < int(statement.nodes.size()); ++i) {
        auto & node = statement.nodes[i];
        if (i != statement.root) {
            out << "1 ~" << a(i);
            for (auto & s : selectors[i])
                out << " 1 " << s;
            out << " >= 1 ;\n";
            for (auto & s : selectors[i])
                out << "1 " << a(i) << " 1 ~" << s << " >= 1 ;\n";
        }

        switch (node.kind) {
        case CountingNodeKind::Cube:
            for (auto & [u, value] : node.assignment)
                out << "1 ~" << a(i) << " 1 " << y(u, value) << " >= 1 ;\n";
            break;

        case CountingNodeKind::And:
            break;

        case CountingNodeKind::Or: {
            auto k = node.children.size();
            out << "1 ~" << a(i);
            for (unsigned j = 0; j < k; ++j)
                out << " 1 " << t(i, j);
            out << " >= 1 ;\n";
            for (unsigned j = 0; j < k; ++j)
                out << "1 ~" << t(i, j) << " 1 " << a(i) << " >= 1 ;\n";
            if (k >= 2) {
                for (unsigned j = 0; j < k; ++j)
                    out << "1 ~" << t(i, j) << " ";
                out << ">= " << k - 1 << " ;\n";
            }
        } break;
        }
    }
}
