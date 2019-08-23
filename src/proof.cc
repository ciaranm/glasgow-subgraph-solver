/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "proof.hh"

#include <algorithm>
#include <iterator>
#include <fstream>
#include <map>
#include <sstream>

using std::copy;
using std::endl;
using std::istreambuf_iterator;
using std::map;
using std::ofstream;
using std::ostreambuf_iterator;
using std::pair;
using std::string;
using std::stringstream;

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
    string opb_filename;

    stringstream model_stream;
    map<pair<int, int>, int> variable_mappings;
    int nb_constraints = 0;
};

Proof::Proof(const std::string & opb_file, const std::string & log_file) :
    _imp(new Imp{ })
{
    _imp->opb_filename = opb_file;
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
    for (int i = 0 ; i < target_size ; ++i)
        _imp->model_stream << "-1 x" << _imp->variable_mappings[{ pattern_vertex, i }] << " ";
    _imp->model_stream << ">= -1 ;" << endl;
    _imp->nb_constraints += 2;
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
        ++_imp->nb_constraints;
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
    ++_imp->nb_constraints;
}

auto Proof::finalise_model() -> void
{
    ofstream f{ _imp->opb_filename };

    f << "* #variable= " << _imp->variable_mappings.size() << " #constraint= " << _imp->nb_constraints << endl;
    copy(istreambuf_iterator<char>{ _imp->model_stream }, istreambuf_iterator<char>{}, ostreambuf_iterator<char>{ f });
    _imp->model_stream.clear();

    if (! f)
        throw ProofError{ "Error writing opb file to '" + _imp->opb_filename + "'" };
}

