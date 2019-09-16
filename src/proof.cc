/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "proof.hh"

#include <algorithm>
#include <iterator>
#include <fstream>
#include <map>
#include <sstream>
#include <tuple>

using std::copy;
using std::endl;
using std::istreambuf_iterator;
using std::map;
using std::ofstream;
using std::ostreambuf_iterator;
using std::pair;
using std::string;
using std::stringstream;
using std::tuple;
using std::vector;

ProofError::ProofError(const string & m) noexcept :
    _message("Proof error: " + m)
{
}

auto ProofError::what() const noexcept -> const char *
{
    return _message.c_str();
}

struct Proof::Imp
{
    string opb_filename, log_filename;
    stringstream model_stream;
    ofstream proof_stream;
    bool levels;

    map<pair<int, int>, int> variable_mappings;
    map<int, int> at_least_one_value_constraints, at_most_one_value_constraints, injectivity_constraints;
    map<tuple<int, int, int, int>, int> adjacency_lines;

    int nb_constraints = 0;
    int proof_line = 0;
};

Proof::Proof(const std::string & opb_file, const std::string & log_file, bool l) :
    _imp(new Imp)
{
    _imp->opb_filename = opb_file;
    _imp->log_filename = log_file;
    _imp->levels = l;
}

Proof::Proof(Proof &&) = default;

Proof::~Proof() = default;

auto Proof::operator= (Proof &&) -> Proof & = default;

auto Proof::create_cp_variable(int pattern_vertex, int target_size) -> void
{
    for (int i = 0 ; i < target_size ; ++i)
        _imp->variable_mappings.emplace(pair{ pattern_vertex, i }, _imp->variable_mappings.size() + 1);

    _imp->model_stream << "* vertex " << pattern_vertex << " domain" << endl;
    for (int i = 0 ; i < target_size ; ++i)
        _imp->model_stream << "1 x" << _imp->variable_mappings[{ pattern_vertex, i }] << " ";
    _imp->model_stream << ">= 1 ;" << endl;
    _imp->at_least_one_value_constraints.emplace(pattern_vertex, ++_imp->nb_constraints);

    for (int i = 0 ; i < target_size ; ++i)
        _imp->model_stream << "-1 x" << _imp->variable_mappings[{ pattern_vertex, i }] << " ";
    _imp->model_stream << ">= -1 ;" << endl;
    _imp->at_most_one_value_constraints.emplace(pattern_vertex, ++_imp->nb_constraints);
}

auto Proof::create_injectivity_constraints(int pattern_size, int target_size) -> void
{
    for (int v = 0 ; v < target_size ; ++v) {
        _imp->model_stream << "* injectivity on value " << v << endl;

        for (int p = 0 ; p < pattern_size ; ++p) {
            auto x = _imp->variable_mappings.find(pair{ p, v });
            if (x != _imp->variable_mappings.end())
                _imp->model_stream << "-1 x" << x->second << " ";
        }
        _imp->model_stream << ">= -1 ;" << endl;
        _imp->injectivity_constraints.emplace(v, ++_imp->nb_constraints);
    }
}

auto Proof::create_forbidden_assignment_constraint(int p, int t) -> void
{
    _imp->model_stream << "* forbidden assignment" << endl;
    _imp->model_stream << "1 ~x" << _imp->variable_mappings[pair{ p, t }] << " >= 1 ;" << endl;
    ++_imp->nb_constraints;
}

auto Proof::create_adjacency_constraint(int p, int q, int t, const std::vector<int> & uu) -> void
{
    _imp->model_stream << "* adjacency for edge " << p << " -- " << q << " mapping to vertex " << t << endl;
    _imp->model_stream << "1 ~x" << _imp->variable_mappings[pair{ p, t }];
    for (auto & u : uu)
        _imp->model_stream << " 1 x" << _imp->variable_mappings[pair{ q, u }];
    _imp->model_stream << " >= 1 ;" << endl;
    _imp->adjacency_lines.emplace(tuple{ 0, p, q, t }, ++_imp->nb_constraints);
}

auto Proof::finalise_model() -> void
{
    ofstream f{ _imp->opb_filename };

    f << "* #variable= " << _imp->variable_mappings.size() << " #constraint= " << _imp->nb_constraints << endl;
    copy(istreambuf_iterator<char>{ _imp->model_stream }, istreambuf_iterator<char>{}, ostreambuf_iterator<char>{ f });
    _imp->model_stream.clear();

    if (! f)
        throw ProofError{ "Error writing opb file to '" + _imp->opb_filename + "'" };

    _imp->proof_stream = ofstream{ _imp->log_filename };
    _imp->proof_stream << "refutation using f l p u w c 0" << endl;
    _imp->proof_stream << "f " << _imp->nb_constraints << " 0" << endl;
    _imp->proof_line += _imp->nb_constraints;
    _imp->proof_stream << "l " << _imp->variable_mappings.size() << " 0" << endl;
    _imp->proof_line += _imp->variable_mappings.size() * 2;

    if (! _imp->proof_stream)
        throw ProofError{ "Error writing proof file to '" + _imp->log_filename + "'" };

}

auto Proof::finish_unsat_proof() -> void
{
    _imp->proof_stream << "* asserting that we've proved unsat" << endl;
    _imp->proof_stream << "u opb >= 1 ;" << endl;
    ++_imp->proof_line;
    _imp->proof_stream << "c " << _imp->proof_line << " 0" << endl;
}

auto Proof::failure_due_to_pattern_bigger_than_target() -> void
{
    _imp->proof_stream << "* failure due to the pattern being bigger than the target" << endl;

    // we get a hall violator by adding up all of the things
    _imp->proof_stream << "p 0";
    for (auto & [ _, line ] : _imp->at_least_one_value_constraints)
        _imp->proof_stream << " " << line << " +";
    for (auto & [ _, line ] : _imp->injectivity_constraints)
        _imp->proof_stream << " " << line << " +";
    _imp->proof_stream << " 0" << endl;
    ++_imp->proof_line;
}

auto Proof::incompatible_by_degrees(int g, int p, const vector<int> & n_p, int t, const vector<int> & n_t) -> void
{
    _imp->proof_stream << "* cannot map " << p << " to " << t << " due to degrees in graph pairs " << g << endl;

    _imp->proof_stream << "p 0";
    for (auto & n : n_p)
        _imp->proof_stream << " " << _imp->adjacency_lines[tuple{ g, p, n, t }] << " +";

    // if I map p to t, I have to map the neighbours of p to neighbours of t
    for (auto & n : n_t)
        _imp->proof_stream << " " << _imp->injectivity_constraints[n] << " +";

    _imp->proof_stream << " 0" << endl;
    ++_imp->proof_line;

    _imp->proof_stream << "u opb 1 ~x" << _imp->variable_mappings[pair{ p, t }] << " >= 1 ;" << endl;
    ++_imp->proof_line;

    _imp->proof_stream << "w " << (_imp->proof_line - 1) << " 0" << endl;
}

auto Proof::incompatible_by_nds(int g, int p, int t) -> void
{
    _imp->proof_stream << "* cannot map " << p << " to " << t << " due to nds in graph pairs " << g << endl;
}

auto Proof::initial_domain_is_empty(int p) -> void
{
    _imp->proof_stream << "* failure due to domain " << p << " being empty" << endl;
}

auto Proof::emit_hall_set_or_violator(const std::vector<int> & lhs, const std::vector<int> & rhs) -> void
{
    _imp->proof_stream << "* hall set or violator size " << lhs.size() << "/" << rhs.size() << endl;
    _imp->proof_stream << "p 0";
    for (auto & l : lhs)
        _imp->proof_stream << " " << _imp->at_least_one_value_constraints[l] << " +";
    for (auto & r : rhs)
        _imp->proof_stream << " " << _imp->injectivity_constraints[r] << " +";
    _imp->proof_stream << " 0" << endl;
    ++_imp->proof_line;
}

auto Proof::root_propagation_failed() -> void
{
    _imp->proof_stream << "* root node propagation failed" << endl;
}

auto Proof::guessing(int depth, int branch_v, int val) -> void
{
    _imp->proof_stream << "* [" << depth << "] guessing " << branch_v << "=" << val << endl;
}

auto Proof::propagation_failure(const std::vector<std::pair<int, int> > & decisions, int branch_v, int val) -> void
{
    _imp->proof_stream << "* [" << decisions.size() << "] propagation failure on " << branch_v << "=" << val << endl;
    _imp->proof_stream << "u opb";
    for (auto & [ var, val ] : decisions)
        _imp->proof_stream << " -1 x" << _imp->variable_mappings[pair{ var, val }];
    _imp->proof_stream << " >= -" << (decisions.size() - 1) << " ;" << endl;
    ++_imp->proof_line;
}

auto Proof::incorrect_guess(const std::vector<std::pair<int, int> > & decisions) -> void
{
    _imp->proof_stream << "* [" << decisions.size() << "] incorrect guess" << endl;
    _imp->proof_stream << "u opb";
    for (auto & [ var, val ] : decisions)
        _imp->proof_stream << " -1 x" << _imp->variable_mappings[pair{ var, val }];
    _imp->proof_stream << " >= -" << (decisions.size() - 1) << " ;" << endl;
    ++_imp->proof_line;
}

auto Proof::out_of_guesses(const std::vector<std::pair<int, int> > & decisions) -> void
{
    _imp->proof_stream << "* [" << decisions.size() << "] out of guesses" << endl;
}

auto Proof::unit_propagating(int var, int val) -> void
{
    _imp->proof_stream << "* unit propagating " << var << "=" << val << endl;
}

auto Proof::start_level(int level) -> void
{
    if (_imp->levels) {
        _imp->proof_stream << "lvlset " << level << endl;
        _imp->proof_stream << "lvlclear " << level << endl;
    }
}

auto Proof::back_up_to_level(int level) -> void
{
    if (_imp->levels)
        _imp->proof_stream << "lvlset " << level << endl;
}

auto Proof::back_up_to_top() -> void
{
    if (_imp->levels)
        _imp->proof_stream << "lvlset " << 0 << endl;
}

auto Proof::post_restart_nogood(const std::vector<std::pair<int, int> > & decisions) -> void
{
    _imp->proof_stream << "* [" << decisions.size() << "] restart nogood" << endl;
    _imp->proof_stream << "u opb";
    for (auto & [ var, val ] : decisions)
        _imp->proof_stream << " -1 x" << _imp->variable_mappings[pair{ var, val }];
    _imp->proof_stream << " >= -" << (decisions.size() - 1) << " ;" << endl;
    ++_imp->proof_line;
}
