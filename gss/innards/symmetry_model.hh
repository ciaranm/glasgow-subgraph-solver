#include <memory>
#include <vector>

#include "gss/innards/homomorphism_domain.hh"

namespace gss::innards
{
    enum symmetry_mode: int {
        none,
        orbit,
        transversal
    };

    enum flexibility: int {
        stat,
        flexi,
        dynamic
    };

    class SymmetryModel {
        private:
            struct Order;
            std::unique_ptr<Order> pattern_order;
            std::unique_ptr<Order> target_order;
            struct Base;
            std::unique_ptr<Base> pattern_base;
            std::unique_ptr<Base> target_base;
            struct Settings;
            std::unique_ptr<Settings> sym_settings;

            auto enforce_order(int lesser, int greater, bool pattern) -> void;
            auto is_less(int a, int b) -> bool;
        public:
            SymmetryModel(int p_size, int t_size, symmetry_mode p_mode, symmetry_mode t_mode, flexibility p_flex, flexibility t_flex);
            ~SymmetryModel();

            auto num_symmetric_solutions(std::vector<int> &phi, int p_size, int t_size) -> int;

            auto lex_least(std::vector<int> &phi, std::vector<HomomorphismDomain> &new_domains) -> bool;
            auto lex_less(std::vector<int> &phi_alpha, std::vector<int> &phi, std::vector<int> &p_inv, std::vector<int> &t_inv, std::vector<HomomorphismDomain> &new_domains) -> bool;

            auto push_back_base(bool pattern, int base_point) -> void;
            auto pop_back_base(bool pattern) -> void;

            auto add_perm(bool pattern, std::vector<int> &perm, std::vector<int> &perm_inverse) -> void;
    };

}