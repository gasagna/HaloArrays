#include "DArrays.hpp"
#include <catch.hpp>
#include <array>
#include <iostream>

TEST_CASE("2D tests", "[2D tests]") {

    std::array<int, 2> proc_grid_size = {3, 4}; // to be run on 12 processors
    std::array<int, 2> array_size = {30, 40};
    std::array<int, 2> isperiodic = {false, false};
    std::array<int, 2> nhalo_out = {0, 0};
    int nhalo_in = 2;

    DArrays::DArray<double, 2> u(MPI_COMM_WORLD,
                                 proc_grid_size,
                                 array_size,
                                 isperiodic, nhalo_out, nhalo_in);

    SECTION("inbound indexing") {
        for (auto [i, j] : u.indices()) {
            u(i, j) = i*j;
            REQUIRE( u(i, j) == i*j );
        }
    }

    SECTION("out of bounds indexing") {
        if (u.is_on_left_boundary(0) == true )
            REQUIRE_THROWS( u(-1, 0) = 1 );

        if (u.is_on_right_boundary(0) == true )
            REQUIRE_THROWS( u(10, 0) = 1 );

        if (u.is_on_left_boundary(1) == true)
            REQUIRE_THROWS( u(0, -1) = 1 );

        if (u.is_on_right_boundary(1) == true)
            REQUIRE_THROWS( u(0, 10) = 1 );

        if (!u.is_on_boundary() == true) {
            REQUIRE_NOTHROW( u(-1,  0) = 1 );
            REQUIRE_NOTHROW( u(10,  0) = 1 );
            REQUIRE_NOTHROW( u( 0, -1) = 1 );
            REQUIRE_NOTHROW( u( 0, 10) = 1 );
            REQUIRE_NOTHROW( u(-1, -1) = 1 );
            REQUIRE_NOTHROW( u(10, 10) = 1 );
            REQUIRE_NOTHROW( u(10, -1) = 1 );
            REQUIRE_NOTHROW( u(-1, 10) = 1 );
        }
    }

}