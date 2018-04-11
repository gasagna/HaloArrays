#include "DArrays.hpp"
#include <algorithm>
#include <catch.hpp>
#include <iterator>
#include <iostream>
#include <array>

// import all
using namespace DArrays;

TEST_CASE("mpiwrapper", "test_1") {

    SECTION("1D") {
        // use this grid layout for tests
        std::array<int, 1> layout_size = {27};

        SECTION("periodic = false") {
            std::array<int, 1> is_periodic = {true};

            // create layout
            DArrayLayout<1> layout(MPI_COMM_WORLD, layout_size, is_periodic);

            // create array 
            std::array<int, 1> array_size = {5*layout.nprocs()}; 
            std::array<int, 1> nhalo_out  = {3};
            std::array<int, 1> nhalo_in   = {2};
            DArray<double, 1> a(layout, array_size, nhalo_out, nhalo_in); 

            // fill array with rank
            std::fill(a.begin(), a.end(), layout.rank());

            // do it!
            a.swap_halo();

            // check in domain is not modified and left/right halo is filled with proc on left/right
            for (auto i : {0, 1, 2, 3, 4} )
                REQUIRE( a(i) == layout.rank() );

            for (auto i : {-2, -1} )
                REQUIRE( a(i) == layout.rank_of_neighbour_at(Boundary::LEFT, 0) );

            for (auto i : {5, 6} )
                REQUIRE( a(i) == layout.rank_of_neighbour_at(Boundary::RIGHT, 0) );
        }

        SECTION("periodic = true") {
            std::array<int, 1> is_periodic = {true};

            // create layout
            DArrayLayout<1> layout(MPI_COMM_WORLD, layout_size, is_periodic);

            // create array 
            std::array<int, 1> array_size = {5*layout.nprocs()}; 
            std::array<int, 1> nhalo_out  = {3};
            std::array<int, 1> nhalo_in   = {2};
            DArray<double, 1> a(layout, array_size, nhalo_out, nhalo_in); 

            // fill array with rank
            std::fill(a.begin(), a.end(), layout.rank());

            // do it!
            a.swap_halo();

            // check in domain is not modified and left/right halo is filled with proc on left/right
            for (auto i : {0, 1, 2, 3, 4} )
                REQUIRE( a(i) == layout.rank() );

            if (layout.has_neighbour_at(Boundary::LEFT, 0)) {
                for (auto i : {-2, -1} )
                    REQUIRE( a(i) == layout.rank_of_neighbour_at(Boundary::LEFT, 0) );
            } else {
                for (auto i : {-3, -2, -1} )
                    REQUIRE( a(i) == layout.rank() );
            }

            if (layout.has_neighbour_at(Boundary::RIGHT, 0)) {
                for (auto i : {5, 6} )
                    REQUIRE( a(i) == layout.rank_of_neighbour_at(Boundary::RIGHT, 0) );
            } else {
                for (auto i : {5, 6, 7} )
                    REQUIRE( a(i) == layout.rank() );
            }
        }
    }
}