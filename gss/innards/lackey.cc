#include <gss/innards/lackey.hh>

#include <fstream>
#include <map>
#include <mutex>

using namespace gss;
using namespace gss::innards;

using std::endl;
using std::function;
using std::ifstream;
using std::map;
using std::mutex;
using std::ofstream;
using std::string;
using std::to_string;
using std::unique_lock;

DisobedientLackeyError::DisobedientLackeyError(const std::string & m) noexcept :
    _what(m)
{
}

auto DisobedientLackeyError::what() const throw() -> const char *
{
    return _what.c_str();
}

struct Lackey::Imp
{
    mutex external_solver_mutex;

    ofstream send_to;
    ifstream read_from;
    const InputGraph & pattern_graph;
    const InputGraph & target_graph;

    long number_of_checks = 0, number_of_propagations = 0, number_of_deletions = 0, number_of_calls = 0;
};

Lackey::Lackey(const string & send_to_name, const string & read_from_name,
    const InputGraph & pattern_graph, const InputGraph & target_graph) :
    _imp(new Imp{{}, ofstream{send_to_name}, ifstream{read_from_name}, pattern_graph, target_graph})
{
    if ((! _imp->read_from) || (! _imp->send_to))
        throw DisobedientLackeyError{"error setting up lackey communication using " + send_to_name + " and " + read_from_name};
}

Lackey::~Lackey()
{
    if (_imp->send_to) {
        _imp->send_to << "Q 0" << endl;
    }
}

auto Lackey::check_solution(
    const VertexToVertexMapping & m,
    bool partial,
    bool all_solutions,
    const function<auto(int, int)->bool> & deletion) -> bool
{
    unique_lock<mutex> lock{_imp->external_solver_mutex};
    ++_imp->number_of_calls;

    string command;
    if (partial) {
        if (deletion) {
            ++_imp->number_of_propagations;
            command = "P";
        }
        else {
            ++_imp->number_of_checks;
            command = "C";
        }
    }
    else {
        ++_imp->number_of_checks;
        if (all_solutions)
            command = "A";
        else
            command = "F";
    }

    _imp->send_to << command << " " << m.size();
    for (auto & [p, t] : m)
        _imp->send_to << " " << _imp->pattern_graph.vertex_name(p) << " " << _imp->target_graph.vertex_name(t);
    _imp->send_to << endl;

    if (! _imp->send_to)
        throw DisobedientLackeyError{"error giving lackey its orders"};

    string operation;
    if (! (_imp->read_from >> operation) || operation != command)
        throw DisobedientLackeyError{"asked lackey to " + command + ", but it replied with '" + operation + "'"};

    bool result;
    string response;
    if (! (_imp->read_from >> response))
        throw DisobedientLackeyError{"asked lackey to " + command + ", but it gave no T/F"};
    else if (response == "T")
        result = true;
    else if (response == "F")
        result = false;
    else
        throw DisobedientLackeyError{"asked lackey to " + command + " but it replied with '" + operation + "' then '" + response + "'"};

    int n;
    if (! (_imp->read_from >> n))
        throw DisobedientLackeyError{"lackey replied with length '" + to_string(n) + "' to " + command + " query"};

    if (command == "S") {
        for (int i = 0; i < n; ++i) {
            string k, v;
            if (! (_imp->read_from >> k >> v))
                throw DisobedientLackeyError{"lackey gave bad response pair " + to_string(i) + " to " + command + " query"};
        }
    }
    else if (command == "C" || command == "P") {
        for (int i = 0; i < n; ++i) {
            string k, v;
            int m;
            if (! (_imp->read_from >> k >> m))
                throw DisobedientLackeyError{"lackey gave bad response pair " + k + " " + to_string(m) + " to " + command + " query"};
            auto p = _imp->pattern_graph.vertex_from_name(k);

            for (int j = 0; j < m; ++j) {
                if (! (_imp->read_from >> v))
                    throw DisobedientLackeyError{"lackey gave bad response pair " + k + " " + to_string(m) + " to " + command + " query"};

                if (deletion) {
                    auto t = _imp->target_graph.vertex_from_name(v);
                    if (p && t)
                        if (deletion(*p, *t))
                            ++_imp->number_of_deletions;
                }
            }
        }
    }

    return result;
}

auto Lackey::reduce_initial_bounds(
    const RestrictRangeFunction & restrict_range) -> bool
{
    unique_lock<mutex> lock{_imp->external_solver_mutex};
    ++_imp->number_of_calls;

    string command = "I";
    _imp->send_to << command << " " << 0 << endl;

    if (! _imp->send_to)
        throw DisobedientLackeyError{"error giving lackey its orders"};

    string operation;
    if (! (_imp->read_from >> operation) || operation != command)
        throw DisobedientLackeyError{"asked lackey to " + command + ", but it replied with '" + operation + "'"};

    string response;
    if (! (_imp->read_from >> response))
        throw DisobedientLackeyError{"asked lackey to " + command + ", but it gave no T/F"};
    else if (response == "T") {
        /* nothing */
    }
    else if (response == "F")
        return false;
    else
        throw DisobedientLackeyError{"asked lackey to " + command + " but it replied with '" + operation + "' then '" + response + "'"};

    int n;
    if (! (_imp->read_from >> n))
        throw DisobedientLackeyError{"lackey replied with length '" + to_string(n) + "' to " + command + " query"};

    for (int i = 0; i < n; ++i) {
        string k;
        int lower, upper;
        if (! (_imp->read_from >> k >> lower >> upper))
            throw DisobedientLackeyError{"lackey gave bad response triple " + to_string(i) + " to " + command + " query"};
        auto p = _imp->pattern_graph.vertex_from_name(k);
        if (p) {
            auto delete_one = [&](int v) {
                auto v_name = _imp->target_graph.vertex_from_name(to_string(v));
                if (v_name)
                    restrict_range(*p, *v_name);
            };

            for (int i = 1; i < lower; ++i)
                delete_one(i);
            for (int i = upper + 1; i < _imp->target_graph.size(); ++i)
                delete_one(i);
        }
    }

    return true;
}

auto Lackey::number_of_checks() const -> long
{
    return _imp->number_of_checks;
}

auto Lackey::number_of_propagations() const -> long
{
    return _imp->number_of_propagations;
}

auto Lackey::number_of_deletions() const -> long
{
    return _imp->number_of_deletions;
}

auto Lackey::number_of_calls() const -> long
{
    return _imp->number_of_calls;
}
