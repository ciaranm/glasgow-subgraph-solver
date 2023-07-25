#include <gss/innards/cheap_all_different.hh>

#include <tuple>
#include <type_traits>

using namespace gss;
using namespace gss::innards;

using std::conditional_t;
using std::shared_ptr;
using std::tuple;
using std::vector;

namespace
{
    template <bool proof_>
    auto cheap_all_different_with_optional_proofs(
        unsigned target_size,
        vector<HomomorphismDomain> & domains,
        const shared_ptr<Proof> & proof,
        const HomomorphismModel * const model) -> bool
    {
        // Pick domains smallest first; ties are broken by smallest .v first.
        // For each count p we have a linked list, whose first member is
        // first[p].  The element following x in one of these lists is next[x].
        // Any domain with a count greater than domains.size() is put
        // int the "count==domains.size()" bucket.
        // The "first" array is sized to be able to hold domains.size()+1
        // elements
        vector<int> first(target_size + 1, -1), next(target_size, -1);

        [[maybe_unused]] conditional_t<proof_, vector<NamedVertex>, tuple<>> lhs, hall_lhs, hall_rhs;

        // Iterate backwards, because we insert elements at the head of
        // lists and we want the sort to be stable
        for (int i = int(domains.size()) - 1; i >= 0; --i) {
            unsigned count = domains.at(i).count;
            if (count > domains.size())
                count = domains.size();
            next.at(i) = first.at(count);
            first.at(count) = i;
        }

        // counting all-different
        SVOBitset domains_so_far{target_size, 0}, hall{target_size, 0};
        unsigned neighbours_so_far = 0;

        [[maybe_unused]] conditional_t<proof_, unsigned, tuple<>> last_outputted_hall_size{};

        for (unsigned i = 0; i <= domains.size(); ++i) {
            // iterate over linked lists
            int domain_index = first[i];
            while (domain_index != -1) {
                auto & d = domains.at(domain_index);

                if constexpr (proof_)
                    lhs.push_back(model->pattern_vertex_for_proof(d.v));

                [[maybe_unused]] conditional_t<proof_, unsigned, tuple<>> old_d_values_count;
                if constexpr (proof_)
                    old_d_values_count = d.values.count();

                d.values.intersect_with_complement(hall);
                d.count = d.values.count();

                if constexpr (proof_)
                    if (last_outputted_hall_size != hall.count() && d.count != old_d_values_count) {
                        last_outputted_hall_size = hall.count();
                        proof->emit_hall_set_or_violator(hall_lhs, hall_rhs);
                    }

                if (0 == d.count)
                    return false;

                domains_so_far |= d.values;
                ++neighbours_so_far;

                unsigned domains_so_far_popcount = domains_so_far.count();

                if (domains_so_far_popcount < neighbours_so_far) {
                    // hall violator, so we fail (after outputting a proof)
                    if constexpr (proof_) {
                        vector<NamedVertex> rhs;
                        auto d = domains_so_far;
                        for (auto v = d.find_first(); v != decltype(d)::npos; v = d.find_first()) {
                            d.reset(v);
                            rhs.push_back(model->target_vertex_for_proof(v));
                        }
                        proof->emit_hall_set_or_violator(lhs, rhs);
                    }
                    return false;
                }
                else if (domains_so_far_popcount == neighbours_so_far) {
                    // equivalent to hall=domains_so_far
                    hall |= domains_so_far;
                    if constexpr (proof_) {
                        hall_lhs.clear();
                        for (auto & l : lhs)
                            hall_lhs.push_back(l);
                        hall_rhs.clear();
                        auto d = domains_so_far;
                        for (auto v = d.find_first(); v != decltype(d)::npos; v = d.find_first()) {
                            d.reset(v);
                            hall_rhs.push_back(model->target_vertex_for_proof(v));
                        }
                    }
                }
                domain_index = next[domain_index];
            }
        }

        return true;
    }
}

auto gss::innards::cheap_all_different(unsigned target_size, vector<HomomorphismDomain> & domains, const shared_ptr<Proof> & proof,
    const HomomorphismModel * const model) -> bool
{
    if (! proof.get())
        return cheap_all_different_with_optional_proofs<false>(target_size, domains, proof, model);
    else
        return cheap_all_different_with_optional_proofs<true>(target_size, domains, proof, model);
}
