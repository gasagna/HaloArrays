#include "DArrays.hpp"
#include <algorithm>
#include <catch.hpp>
#include <iterator>
#include <iostream>
#include <array>
#include <set>

template <typename CONTAINER>
bool all_unique(const CONTAINER& container) {
    std::set<typename CONTAINER::value_type> seen;
    for (auto& el : container) {
        if ( seen.count(el) == 1 ) {
            return false;
        }
        seen.insert(el);
    }
    return true;
}

TEST_CASE("mpiwrapper", "mpiwrapper") {

    SECTION("1D") {
        auto bnds = DArrays::AllBoundaries<1>();
        std::vector<int> tags;
        std::transform(bnds.begin(), bnds.end(), 
                    std::back_inserter(tags), DArrays::MPI::message_tag<1>);
        REQUIRE( all_unique(tags) == true );
    }

    SECTION("2D") {
        auto bnds = DArrays::AllBoundaries<2>();
        std::vector<int> tags;
        std::transform(bnds.begin(), bnds.end(), 
                    std::back_inserter(tags), DArrays::MPI::message_tag<2>);
        REQUIRE( all_unique(tags) == true );
    }

    SECTION("3D") {
        auto bnds = DArrays::AllBoundaries<3>();
        std::vector<int> tags;
        std::transform(bnds.begin(), bnds.end(), 
                    std::back_inserter(tags), DArrays::MPI::message_tag<3>);
        REQUIRE( all_unique(tags) == true );
    }
}
