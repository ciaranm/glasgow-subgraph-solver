#include <gss/innards/symmetry_model.hh>
#include <gss/innards/homomorphism_domain.hh>

#include <vector>
#include <numeric>

using namespace gss::innards;

struct SymmetryModel::Settings
{
    const symmetry_mode pattern_mode;
    const symmetry_mode target_mode;
    const flexibility pattern_flex;
    const flexibility target_flex;

    Settings(symmetry_mode p_mode, flexibility p_flex, symmetry_mode t_mode, flexibility t_flex) :
        pattern_mode(p_mode),
        target_mode(t_mode),
        pattern_flex(p_flex),
        target_flex(t_flex)
    {
    }
};

struct SymmetryModel::Order
{
    std::vector<int> values;
    int size;

    Order(int sz)
    {
        values.resize(sz);
        size = sz;
        std::iota(values.begin(), values.end(), 0);
    }
};

struct SymmetryModel::Base
{
    std::vector<int> values;                                            // values[i] := the i^th base point
    std::vector<std::vector<int>> orbits;                               // orbits[i] := the orbit of the value i
    std::vector<std::vector<std::vector<int>>> transversals;            // transversals[i] := the transversals of the value i
    std::vector<std::vector<std::vector<int>>> t_inverses;              // t_inverses[i] := the inverse of transversals[i]
    int aut_grp_sz;

    Base(int sz) {
        values.resize(sz);
        orbits.resize(sz);
        transversals.resize(sz);
    }
};

SymmetryModel::SymmetryModel(int p_size, int t_size, symmetry_mode p_mode, symmetry_mode t_mode, flexibility p_flex, flexibility t_flex) :
    pattern_base(new Base(p_size)),
    target_base(new Base(t_size)),
    pattern_order(new Order(p_size)),
    target_order(new Order(t_size)),
    sym_settings(new Settings(p_mode, p_flex, t_mode, t_flex))
{

}

SymmetryModel::~SymmetryModel() = default;

auto SymmetryModel::enforce_order(int lesser, int greater, bool pattern) -> void {
    if (pattern) {
        if (pattern_order->values[lesser] < pattern_order->values[greater]) return;
        pattern_order->values[lesser] = pattern_order->values[greater];
        for (int i = 0; i < pattern_order->values.size(); i++) {
            if (pattern_order->values[greater] <= pattern_order->values[i] && pattern_order->values[i] < pattern_order->values[lesser]) {
                pattern_order->values[i]++;
            }
        }
    }
    else {
        if (target_order->values[lesser] < target_order->values[greater]) return;
        target_order->values[lesser] = target_order->values[greater];
        for (int i = 0; i < target_order->values.size(); i++) {
            if (target_order->values[greater] <= target_order->values[i] && target_order->values[i] < target_order->values[lesser]) {
                target_order->values[i]++;
            }
        }
    }
}

auto SymmetryModel::is_less(int a, int b) -> bool {
    return target_order->values[a] < target_order->values[b];
}

auto SymmetryModel::num_symmetric_solutions(std::vector<int> &phi, int p_sz, int t_sz) -> int {
    if (sym_settings->pattern_mode == orbit) {
        if (sym_settings->target_mode == orbit) {
            // better hope we're doing dynamic...
            return 1;
        }
        else if (sym_settings->target_mode == transversal) {
            // probably fine
            return 1;
        }
        else {
            return pattern_base->aut_grp_sz;
        }
    }
    else if (sym_settings->pattern_mode == transversal) {
        if (sym_settings->target_mode == orbit) {
            // pain to reformulate
            return 1;
        }
        else if (sym_settings->target_mode == transversal) {
            // might be tricky with overlaps
            return 1;
        }
        else {
            return pattern_base->aut_grp_sz;
        }
    }
    else {
        if (sym_settings->target_mode == orbit) {
            int mult = 1;
            for (int i = 0; i < p_sz; i++) {
                mult *= target_base->orbits[phi[i]].size();
            }
            return mult;
        }
        else if (sym_settings->target_mode == transversal) {
            if (sym_settings->target_flex == dynamic) {
                int mult = 1;
                for (int i = 0; i < p_sz; i++) {
                    mult *= target_base->transversals[phi[i]].size();
                }
                return mult;
            } 
            else {      // may be clunky
                return 1;
            }
        }
        else {
            // Hmm... if we aren't breaking symmetries, how did we get here?
            return 1;
        }
    }
}

auto SymmetryModel::lex_least(std::vector<int> &phi, std::vector<HomomorphismDomain> &new_domains) -> bool {


    return true;
}


// TODO "overload" this for the different types
auto SymmetryModel::lex_less(std::vector<int> &phi_alpha, std::vector<int> &phi, std::vector<int> &p_inv, std::vector<int> &t_inv, std::vector<HomomorphismDomain> &new_domains) -> bool {
    for (int i = 0; i < pattern_order->size; i++) {
        if (phi[i] != -1 && phi_alpha[i] != -1) {
            if (phi_alpha[i] < phi[i]) {       // The permuted mapping is 'less than' the original
                return false;
            }
            else if (phi_alpha[i] == phi[i]) {     // The mapping is the same so far
                continue;
            }
            else if (phi_alpha[i] > phi[i]) {      // The original mapping is 'less than' the permutation
                break;                          // TODO we don't need to check this particular p_aut,t_aut combination again until we backtrack
            }
        }
        else if (phi[i] != -1) {
            for (auto &d : new_domains) {
                if (d.v == p_inv[i]) {               // Find the variable's domain
                    for (unsigned int x = 0; x < target_order->size; x++) {  // For each value...
                        if (t_inv[x] < phi[i]) {
                            d.values.reset(t_inv[x]);
                        }
                    }
                }
            }
            break;
        }
        else {                                      // TODO Could do arc consistency here
            break;
        }
    }

    return true;
}

auto SymmetryModel::push_back_base(bool pattern, int base_point) -> void {
    if (pattern) {
        pattern_base->values.push_back(base_point);
    }
    else {
        target_base->values.push_back(base_point);
    }
}

auto SymmetryModel::pop_back_base(bool pattern) -> void {
    if (pattern) {
        int last = pattern_base->values.back();
        if (sym_settings->pattern_mode == orbit) {      // Probably only do this if the orbit is non-trivial
            pattern_base->orbits[last].resize(0);
            pattern_base->orbits[last].push_back(last);
        }
        if (sym_settings->pattern_mode == transversal) {        // As above
            pattern_base->transversals[last].resize(0);
            pattern_base->t_inverses[last].resize(0);
        }
        pattern_base->values.pop_back();
    }
    else {
        int last = target_base->values.back();
        if (sym_settings->target_mode == orbit) {
            target_base->orbits[last].resize(0);
            target_base->orbits[last].push_back(last);
        }
        if (sym_settings->target_mode == transversal) {
            target_base->transversals[last].resize(0);
            target_base->t_inverses[last].resize(0);
        }
        target_base->values.pop_back();
    }
}

auto SymmetryModel::add_perm(bool pattern, std::vector<int> &perm, std::vector<int> &perm_inverse) -> void {

}


