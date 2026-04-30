#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <vector>

#include "expos/pileup.hpp"

TEST_CASE ("PileupMetrics merge") {
    PileupMetrics a;
    a.query_position = {10, 20, 30};
    a.normalised_as = {0.9, 0.8, 0.7};
    a.template_endpoints = {{5, 50}, {15, 60}};

    PileupMetrics b;
    b.query_position = {40, 50};
    b.normalised_as = {0.6, 0.5};
    b.template_endpoints = {{35, 80}};

    auto merged = merge_pileup_metrics (a, b);

    REQUIRE (merged.query_position == std::vector<uint64_t>{10, 20, 30, 40, 50});
    REQUIRE (merged.normalised_as == std::vector<double>{0.9, 0.8, 0.7, 0.6, 0.5});
    REQUIRE (merged.template_endpoints.size() == 3);
    REQUIRE (merged.template_endpoints[2].lmost == 35);
    REQUIRE (merged.template_endpoints[2].rmost == 80);

    // originals must be unmodified
    REQUIRE (a.query_position.size() == 3);
    REQUIRE (b.query_position.size() == 2);
}

TEST_CASE ("PileupMetrics merge empty") {
    PileupMetrics a;
    a.query_position = {5, 10};

    PileupMetrics empty;

    auto merged_ab = merge_pileup_metrics (a, empty);
    REQUIRE (merged_ab.query_position.size() == 2);

    auto merged_ba = merge_pileup_metrics (empty, a);
    REQUIRE (merged_ba.query_position.size() == 2);
}
