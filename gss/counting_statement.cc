#include <gss/counting_statement.hh>

#include <algorithm>
#include <istream>
#include <limits>
#include <ostream>
#include <set>
#include <sstream>

using namespace gss;

using std::istream;
using std::istringstream;
using std::ostream;
using std::pair;
using std::set;
using std::string;
using std::stringstream;
using std::to_string;
using std::vector;

namespace
{
    const string magic = "gss-counting-statement";
    const int version = 1;

    auto kind_name(CountingNodeKind kind) -> string
    {
        switch (kind) {
        case CountingNodeKind::Cube: return "cube";
        case CountingNodeKind::False: return "false";
        case CountingNodeKind::Product: return "product";
        case CountingNodeKind::Sum: return "sum";
        case CountingNodeKind::Difference: return "difference";
        }
        throw CountingStatementError{"unknown node kind"};
    }

    auto literal_name(const CountingStatement & s, int variable, int value) -> string
    {
        auto & [name, domain] = s.variables[variable];
        return "x" + name + "_" + s.domains[domain][value];
    }

    auto check_node_reference(int node, int child) -> void
    {
        if (child < 0 || child >= node)
            throw CountingStatementError{"node " + to_string(node) + " refers to node " + to_string(child) + ", which is not defined before it"};
    }

    auto check_variable(const CountingStatement & s, int node, int variable) -> void
    {
        if (variable < 0 || variable >= int(s.variables.size()))
            throw CountingStatementError{"node " + to_string(node) + " refers to undefined variable " + to_string(variable)};
    }
}

auto CountingStatement::add_domain(vector<string> values) -> int
{
    domains.push_back(std::move(values));
    return domains.size() - 1;
}

auto CountingStatement::add_variable(string name, int domain) -> int
{
    variables.emplace_back(std::move(name), domain);
    return variables.size() - 1;
}

auto CountingStatement::add_cube(vector<pair<int, int>> literals) -> int
{
    nodes.push_back(CountingNode{CountingNodeKind::Cube, std::move(literals), {}, {}});
    return nodes.size() - 1;
}

auto CountingStatement::add_false(vector<int> scope) -> int
{
    nodes.push_back(CountingNode{CountingNodeKind::False, {}, std::move(scope), {}});
    return nodes.size() - 1;
}

auto CountingStatement::add_product(vector<int> children) -> int
{
    nodes.push_back(CountingNode{CountingNodeKind::Product, {}, {}, std::move(children)});
    return nodes.size() - 1;
}

auto CountingStatement::add_sum(vector<int> children) -> int
{
    nodes.push_back(CountingNode{CountingNodeKind::Sum, {}, {}, std::move(children)});
    return nodes.size() - 1;
}

auto CountingStatement::add_difference(int a, int b) -> int
{
    nodes.push_back(CountingNode{CountingNodeKind::Difference, {}, {}, {a, b}});
    return nodes.size() - 1;
}

CountingStatementError::CountingStatementError(const string & message) noexcept :
    _what(message)
{
}

auto CountingStatementError::what() const noexcept -> const char *
{
    return _what.c_str();
}

auto gss::counting_statement_gate_name(int node) -> string
{
    // the certified encoding's own names all start with x, or with conn for connected
    // MCS, so this cannot collide with them
    return "cs" + to_string(node);
}

auto gss::write_counting_statement(ostream & out, const CountingStatement & s) -> void
{
    out << magic << " " << version << '\n';
    for (unsigned d = 0; d < s.domains.size(); ++d) {
        out << "domain " << d << " " << s.domains[d].size();
        for (auto & v : s.domains[d])
            out << " " << v;
        out << '\n';
    }
    for (unsigned v = 0; v < s.variables.size(); ++v)
        out << "variable " << v << " " << s.variables[v].first << " " << s.variables[v].second << '\n';
    for (unsigned n = 0; n < s.nodes.size(); ++n) {
        auto & node = s.nodes[n];
        out << kind_name(node.kind) << " " << n;
        switch (node.kind) {
        case CountingNodeKind::Cube:
            out << " " << node.literals.size();
            for (auto & [variable, value] : node.literals)
                out << " " << variable << " " << value;
            break;
        case CountingNodeKind::False:
            out << " " << node.variables.size();
            for (auto & v : node.variables)
                out << " " << v;
            break;
        case CountingNodeKind::Product:
        case CountingNodeKind::Sum:
            out << " " << node.children.size();
            for (auto & c : node.children)
                out << " " << c;
            break;
        case CountingNodeKind::Difference:
            out << " " << node.children.at(0) << " " << node.children.at(1);
            break;
        }
        out << '\n';
    }
    out << "root " << s.root << '\n';
}

auto gss::read_counting_statement(istream & in) -> CountingStatement
{
    CountingStatement s;
    string line;
    int line_number = 0;
    bool seen_header = false, seen_root = false;

    auto fail = [&](const string & message) -> CountingStatementError {
        return CountingStatementError{"line " + to_string(line_number) + ": " + message};
    };

    while (getline(in, line)) {
        ++line_number;
        if (line.empty() || line[0] == '*')
            continue;

        istringstream tokens{line};
        string word;
        tokens >> word;

        auto read_int = [&]() -> int {
            long long v;
            if (! (tokens >> v) || v < 0 || v > std::numeric_limits<int>::max())
                throw fail("expected a non-negative integer");
            return int(v);
        };
        auto read_name = [&]() -> string {
            string n;
            if (! (tokens >> n))
                throw fail("expected a name");
            return n;
        };
        auto read_id = [&](int expected) {
            if (read_int() != expected)
                throw fail("expected id " + to_string(expected) + ", ids must be consecutive from zero");
        };

        if (! seen_header) {
            if (word != magic)
                throw fail("expected '" + magic + "'");
            if (read_int() != version)
                throw fail("unsupported version");
            seen_header = true;
        }
        else if (seen_root)
            throw fail("nothing may follow the root");
        else if (word == "domain") {
            read_id(s.domains.size());
            vector<string> values(read_int());
            for (auto & v : values)
                v = read_name();
            s.add_domain(std::move(values));
        }
        else if (word == "variable") {
            read_id(s.variables.size());
            auto name = read_name();
            auto domain = read_int();
            if (domain >= int(s.domains.size()))
                throw fail("undefined domain " + to_string(domain));
            s.add_variable(std::move(name), domain);
        }
        else if (word == "cube") {
            read_id(s.nodes.size());
            vector<pair<int, int>> literals(read_int());
            for (auto & [variable, value] : literals) {
                variable = read_int();
                value = read_int();
            }
            s.add_cube(std::move(literals));
        }
        else if (word == "false") {
            read_id(s.nodes.size());
            vector<int> variables(read_int());
            for (auto & v : variables)
                v = read_int();
            s.add_false(std::move(variables));
        }
        else if (word == "product" || word == "sum") {
            read_id(s.nodes.size());
            vector<int> children(read_int());
            for (auto & c : children)
                c = read_int();
            if (word == "product")
                s.add_product(std::move(children));
            else
                s.add_sum(std::move(children));
        }
        else if (word == "difference") {
            read_id(s.nodes.size());
            auto a = read_int();
            auto b = read_int();
            s.add_difference(a, b);
        }
        else if (word == "root") {
            s.root = read_int();
            seen_root = true;
        }
        else
            throw fail("unknown statement '" + word + "'");

        string extra;
        if (tokens >> extra)
            throw fail("unexpected '" + extra + "'");
    }

    if (! seen_header)
        throw CountingStatementError{"missing '" + magic + "' header"};
    if (! seen_root)
        throw CountingStatementError{"missing root"};
    return s;
}

auto gss::count_counting_statement(const CountingStatement & s) -> CountingStatementCount
{
    CountingStatementCount result;
    vector<vector<int>> scopes(s.nodes.size());
    vector<loooong> counts(s.nodes.size());

    for (auto & [name, domain] : s.variables)
        if (domain < 0 || domain >= int(s.domains.size()))
            throw CountingStatementError{"variable " + name + " has an undefined domain"};

    // a scratch marking of variables, stamped with the node that set it
    vector<int> seen_by(s.variables.size(), -1);

    auto check_equal_scopes = [&](int n, const CountingNode & node) {
        for (auto & c : node.children)
            if (scopes[c] != scopes[node.children.front()])
                throw CountingStatementError{"node " + to_string(n) + " is a " + kind_name(node.kind) + " of children with different scopes"};
    };

    for (int n = 0; n < int(s.nodes.size()); ++n) {
        auto & node = s.nodes[n];
        for (auto & c : node.children)
            check_node_reference(n, c);

        switch (node.kind) {
        case CountingNodeKind::Cube:
            for (auto & [variable, value] : node.literals) {
                check_variable(s, n, variable);
                if (value < 0 || value >= int(s.domains[s.variables[variable].second].size()))
                    throw CountingStatementError{"node " + to_string(n) + " gives variable " + to_string(variable) + " an undefined value"};
                if (seen_by[variable] == n)
                    throw CountingStatementError{"cube " + to_string(n) + " mentions variable " + to_string(variable) + " twice"};
                seen_by[variable] = n;
                scopes[n].push_back(variable);
            }
            counts[n] = loooong{1};
            break;

        case CountingNodeKind::False:
            for (auto & variable : node.variables) {
                check_variable(s, n, variable);
                if (seen_by[variable] == n)
                    throw CountingStatementError{"false " + to_string(n) + " mentions variable " + to_string(variable) + " twice"};
                seen_by[variable] = n;
                scopes[n].push_back(variable);
            }
            counts[n] = loooong{0};
            break;

        case CountingNodeKind::Product:
            counts[n] = loooong{1};
            for (auto & c : node.children) {
                for (auto & variable : scopes[c]) {
                    if (seen_by[variable] == n)
                        throw CountingStatementError{"product " + to_string(n) + " has children whose scopes overlap on variable " + to_string(variable)};
                    seen_by[variable] = n;
                    scopes[n].push_back(variable);
                }
                counts[n] *= counts[c];
            }
            break;

        case CountingNodeKind::Sum: {
            if (node.children.empty())
                throw CountingStatementError{"sum " + to_string(n) + " has no children, so no scope; use false instead"};
            check_equal_scopes(n, node);
            scopes[n] = scopes[node.children.front()];

            // two cubes over the same scope hold together exactly when they are the same
            // cube, so disjointness among the cubes is checkable here, and only the rest
            // of the children need a proof
            set<vector<pair<int, int>>> cubes;
            bool all_cubes = true;
            counts[n] = loooong{0};
            for (auto & c : node.children) {
                if (s.nodes[c].kind == CountingNodeKind::Cube) {
                    auto literals = s.nodes[c].literals;
                    sort(literals.begin(), literals.end());
                    if (! cubes.insert(literals).second)
                        throw CountingStatementError{"sum " + to_string(n) + " has two children that are the same cube"};
                }
                else if (s.nodes[c].kind != CountingNodeKind::False)
                    all_cubes = false;
                counts[n] += counts[c];
            }
            if (! all_cubes && node.children.size() > 1)
                result.obligations.push_back(CountingObligation{CountingObligation::Kind::Disjoint, n});
            break;
        }

        case CountingNodeKind::Difference:
            if (node.children.size() != 2)
                throw CountingStatementError{"difference " + to_string(n) + " must have two children"};
            check_equal_scopes(n, node);
            scopes[n] = scopes[node.children[0]];
            counts[n] = counts[node.children[0]] - counts[node.children[1]];
            if (counts[n] < loooong{0})
                throw CountingStatementError{"difference " + to_string(n) + " has a negative count, so b cannot imply a"};
            if (node.children[0] != node.children[1])
                result.obligations.push_back(CountingObligation{CountingObligation::Kind::Implies, n});
            break;
        }

        sort(scopes[n].begin(), scopes[n].end());
    }

    if (s.root < 0 || s.root >= int(s.nodes.size()))
        throw CountingStatementError{"the root is not a defined node"};
    if (int(scopes[s.root].size()) != int(s.variables.size()))
        throw CountingStatementError{"the root's scope is not every variable"};

    result.count = counts[s.root];
    return result;
}

auto gss::write_counting_statement_opb(ostream & out, const CountingStatement & s) -> void
{
    stringstream constraints;
    long number_of_constraints = 0, number_of_variables = s.nodes.size();

    auto emit = [&](const string & c) {
        constraints << c << " ;\n";
        ++number_of_constraints;
    };

    // the exactly-one groups, in exactly the certified encoding's syntax, so that the
    // proof does not need to rewrite them (see Proof::create_cp_variable)
    for (unsigned v = 0; v < s.variables.size(); ++v) {
        auto & [name, domain] = s.variables[v];
        auto size = s.domains[domain].size();
        number_of_variables += size;

        constraints << "* vertex " << v << " domain\n";
        string al1 = "@al1" + name, am1 = "@am1" + name;
        for (unsigned i = 0; i < size; ++i) {
            al1 += " 1 " + literal_name(s, v, i);
            am1 += " -1 " + literal_name(s, v, i);
        }
        emit(al1 + " >= 1");
        emit(am1 + " >= -1");
    }

    // r <-> l_1 /\ ... /\ l_k, as k ~r + l_1 + ... + l_k >= k and r + ~l_1 + ... >= 1
    auto define_and = [&](const string & r, const vector<string> & literals) {
        string forwards = to_string(literals.size()) + " ~" + r, backwards = "1 " + r;
        for (auto & l : literals) {
            forwards += " 1 " + l;
            backwards += " 1 ~" + l;
        }
        if (! literals.empty())
            emit(forwards + " >= " + to_string(literals.size()));
        emit(backwards + " >= 1");
    };

    for (int n = 0; n < int(s.nodes.size()); ++n) {
        auto & node = s.nodes[n];
        auto r = counting_statement_gate_name(n);
        constraints << "* " << kind_name(node.kind) << " " << n << '\n';
        switch (node.kind) {
        case CountingNodeKind::Cube: {
            vector<string> literals;
            for (auto & [variable, value] : node.literals)
                literals.push_back(literal_name(s, variable, value));
            define_and(r, literals);
            break;
        }

        case CountingNodeKind::False:
            emit("1 ~" + r + " >= 1");
            break;

        case CountingNodeKind::Product: {
            vector<string> literals;
            for (auto & c : node.children)
                literals.push_back(counting_statement_gate_name(c));
            define_and(r, literals);
            break;
        }

        case CountingNodeKind::Sum: {
            string forwards = "1 ~" + r;
            for (auto & c : node.children) {
                forwards += " 1 " + counting_statement_gate_name(c);
                emit("1 " + r + " 1 ~" + counting_statement_gate_name(c) + " >= 1");
            }
            emit(forwards + " >= 1");
            break;
        }

        case CountingNodeKind::Difference: {
            auto a = counting_statement_gate_name(node.children.at(0)), b = counting_statement_gate_name(node.children.at(1));
            emit("1 ~" + r + " 1 " + a + " >= 1");
            emit("1 ~" + r + " 1 ~" + b + " >= 1");
            emit("1 " + r + " 1 ~" + a + " 1 " + b + " >= 1");
            break;
        }
        }
    }

    constraints << "* root\n";
    emit("1 " + counting_statement_gate_name(s.root) + " >= 1");

    out << "* #variable= " << number_of_variables << " #constraint= " << number_of_constraints << ";\n";
    out << constraints.str();
}
