#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_WATCHES_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_WATCHES_HH 1

#include <algorithm>
#include <list>
#include <utility>
#include <vector>

namespace gss::innards
{
    // A nogood, aways of the form (list of decisions) -> false, where the
    // last part is implicit. If there are at least two assignments, then the
    // first two assignments are the watches (and the literals are permuted
    // when the watches are updates).
    template <typename Decision_>
    struct Nogood
    {
        std::vector<Decision_> literals;
    };

    // Two watched literals for our nogoods store.
    template <typename Decision_, template <typename> typename WatchTable_>
    struct Watches
    {
        // nogoods stored here
        using NogoodStore = std::list<Nogood<Decision_>>;

        NogoodStore nogoods;

        // For each watched literal, we have a list of watched things, each of
        // which is an iterator into the global watch list (so we can reorder
        // the literal to keep the watches as the first two elements).
        using WatchList = std::list<typename NogoodStore::iterator>;

        WatchTable_<WatchList> table;

        // Rather than backjumping, we update the watch list on restarts (to make
        // parallel shenanigans easier).
        using NeedToWatch = std::list<typename NogoodStore::iterator>;

        NeedToWatch need_to_watch, gathered_need_to_watch;

        template <typename CanWatchFunction_, typename AssignmentIsNogoodFunction_>
        auto propagate(
            Decision_ current_assignment,
            const CanWatchFunction_ & can_watch,
            const AssignmentIsNogoodFunction_ & assignment_is_nogood) -> void
        {
            auto & watches_to_update = table[current_assignment];
            for (auto watch_to_update = watches_to_update.begin(); watch_to_update != watches_to_update.end();) {
                auto & nogood = **watch_to_update;

                // make the first watch the thing we just triggered
                if (nogood.literals[0] != current_assignment)
                    std::swap(nogood.literals[0], nogood.literals[1]);

                // can we find something else to watch?
                bool success = false;
                for (auto new_literal = next(nogood.literals.begin(), 2); new_literal != nogood.literals.end(); ++new_literal) {
                    if (can_watch(*new_literal)) {
                        // we can watch new_literal instead of current_assignment in this nogood
                        success = true;

                        // move the new watch to be the first item in the nogood
                        std::swap(nogood.literals[0], *new_literal);

                        // start watching it
                        table[nogood.literals[0]].push_back(*watch_to_update);

                        // remove the current watch, and update the loop iterator
                        watches_to_update.erase(watch_to_update++);

                        break;
                    }
                }

                // found something new? nothing to propagate (and we've already updated our loop iterator in the erase)
                if (success)
                    continue;

                // no new watch, this nogood will now propagate.
                assignment_is_nogood(nogood.literals[1]);

                // step the loop variable, only if we've not already erased it
                ++watch_to_update;
            }
        }

        // posts a nogood, which doesn't kick in until apply_new_nogoods() is
        // called.
        auto post_nogood(Nogood<Decision_> && nogood)
        {
            nogoods.emplace_back(std::move(nogood));
            need_to_watch.emplace_back(std::prev(nogoods.end()));
        }

        template <typename AssignmentIsNogoodFunction_>
        auto apply_new_nogoods(
            const AssignmentIsNogoodFunction_ & assignment_is_nogood) -> bool
        {
            for (auto & n : need_to_watch)
                if (apply_one_new_nogood(n, assignment_is_nogood))
                    return true;

            for (auto & n : gathered_need_to_watch)
                if (apply_one_new_nogood(n, assignment_is_nogood))
                    return true;

            return false;
        }

        template <typename AssignmentIsNogoodFunction_>
        auto apply_one_new_nogood(
            const typename NogoodStore::iterator & n,
            const AssignmentIsNogoodFunction_ & assignment_is_nogood) -> bool
        {
            if (n->literals.empty())
                return true;
            else if (1 == n->literals.size())
                assignment_is_nogood(n->literals[0]);
            else {
                table[n->literals[0]].push_back(n);
                table[n->literals[1]].push_back(n);
            }

            return false;
        }

        auto gather_nogoods_from(
            Watches & other)
        {
            for (auto & n : other.need_to_watch) {
                nogoods.emplace_back(*n);
                gathered_need_to_watch.emplace_back(std::prev(nogoods.end()));
            }
        }

        auto clear_new_nogoods() -> void
        {
            need_to_watch.clear();
            gathered_need_to_watch.clear();
        }
    };
}

#endif
