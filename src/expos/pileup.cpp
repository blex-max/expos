#include "pileup.hpp"

PileupColumn filter
(const PileupColumn &pc, std::function<bool(const bam_pileup1_t*)> predicate_fn)
{
    PileupColumn out;
    for (const auto p : pc) {
        if (predicate_fn (p)) {
            out.push_back (p);
        }
    }
    return out;
}
