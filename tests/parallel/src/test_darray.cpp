#include "DArrays.hpp"
#include <catch.hpp>
#include <array>
#include <iostream>

// import all
using namespace DArrays;

TEST_CASE("darray - 1D", "test_1]") {

    // use this grid layout for tests
    std::array<int, 1> layout_size = {27};

    SECTION("periodic = false") {
        std::array<int, 1> is_periodic = {false};

        // create layout
        DArrayLayout<1> layout(MPI_COMM_WORLD, layout_size, is_periodic);

        // create array 
        std::array<int, 1> array_size = {27*5}; 
        std::array<int, 1> nhalo_out  = {2};
        std::array<int, 1> nhalo_in   = {4};
        DArray<double, 1> a(layout, array_size, nhalo_out, nhalo_in); 

        // every processor in the domain interior indexes 
        // its local part plus nhalo_in boundary points 
        if (!layout.has_neighbour_at(Boundary::LEFT, 0) and 
            !layout.has_neighbour_at(Boundary::RIGHT, 0)) {
            for (int i : {-4, -3, -2, -1,  0, 1, 2, 3, 4,  5, 6, 7, 8}) {
                a(i) = i; 
                REQUIRE( a(i) == i );
            }
            REQUIRE_THROWS( a(-5) );
            REQUIRE_THROWS( a( 9) );

            // test nhalo_points
            REQUIRE( a.nhalo_points(Boundary::LEFT,   0) == 4 );
            REQUIRE( a.nhalo_points(Boundary::RIGHT,  0) == 4 );
        }

        // processors on the left and right domain parts might 
        // have a different number of points they can index
        if (!layout.has_neighbour_at(Boundary::LEFT, 0)) {
            for (int i : {-2, -1,  0, 1, 2, 3, 4,  5, 6, 7, 8}) {
                a(i) = i; 
                REQUIRE( a(i) == i );
            }
            REQUIRE_THROWS( a(-3) );
            REQUIRE_THROWS( a( 9) );

            // test raw_size
            std::array<int, 1> expected_2 = {11};
            REQUIRE( a.raw_size() == expected_2 );

            // test nelements
            REQUIRE( a.nelements() == 11 );   

            // test nhalo_points
            REQUIRE( a.nhalo_points(Boundary::LEFT,   0) == 2 );
            REQUIRE( a.nhalo_points(Boundary::RIGHT,  0) == 4 );
        }

        if (!layout.has_neighbour_at(Boundary::RIGHT, 0)) {
            for (int i : {-4, -3, -2, -1,  0, 1, 2, 3, 4,  5, 6}) {
                a(i) = i; 
                REQUIRE( a(i) == i );
            }
            REQUIRE_THROWS( a(-5) );
            REQUIRE_THROWS( a( 7) );

            // test raw_size
            std::array<int, 1> expected_2 = {11};
            REQUIRE( a.raw_size() == expected_2 );

            // test nelements
            REQUIRE( a.nelements() == 11 );   

            // test nhalo_points
            REQUIRE( a.nhalo_points(Boundary::LEFT,   0) == 4 );
            REQUIRE( a.nhalo_points(Boundary::RIGHT,  0) == 2 );
        }

        // test indices
        std::array<int, 5> expected_1 = {0, 1, 2, 3, 4};
        for ( auto [i] : a.indices() ) {
            REQUIRE( i == expected_1[i] );
        }

        // test size
        std::array<int, 1> expected_3 = {5};
        REQUIRE( a.size() == expected_3 );        
    }

    SECTION("periodic = true") {
        std::array<int, 1> is_periodic = {true};

        // create layout
        DArrayLayout<1> layout(MPI_COMM_WORLD, layout_size, is_periodic);

        // create array 
        std::array<int, 1> array_size = {27*5};
        std::array<int, 1> nhalo_out  = {2};
        std::array<int, 1> nhalo_in   = {4};
        DArray<double, 1> a(layout, array_size, nhalo_out, nhalo_in); 

        // all processors have nhalo_in
        for (int i : {-4, -3, -2, -1,  0, 1, 2, 3, 4,  5, 6, 7, 8}) {
            a(i) = i; 
            REQUIRE( a(i) == i );
        }
        REQUIRE_THROWS( a(-5) );
        REQUIRE_THROWS( a( 9) );

        // test nhalo_points
        REQUIRE( a.nhalo_points(Boundary::LEFT,  0) == 4 );
        REQUIRE( a.nhalo_points(Boundary::RIGHT, 0) == 4 );

        // test indices
        std::array<int, 5> expected_1 = {0, 1, 2, 3, 4};
        for ( auto [i] : a.indices() ) {
            REQUIRE( i == expected_1[i] );
        }

        // test raw_size
        std::array<int, 1> expected_2 = {13};
        REQUIRE( a.raw_size() == expected_2 );

        // test size
        std::array<int, 1> expected_3 = {5};
        REQUIRE( a.size() == expected_3 ); 
        REQUIRE( a.size(0) == expected_3[0] ); 
        REQUIRE_THROWS( a.size( 1) ); 
        REQUIRE_THROWS( a.size(-1) ); 

        // test nelements
        REQUIRE( a.nelements() == 13 );        
    }
}

TEST_CASE("darray - 2D", "test_2]") {

    // use this grid layout for tests
    std::array<int, 2> layout_size = {3, 9};

    SECTION("test grid size and array size compatibility") {
        std::array<int, 2> is_periodic = {false, false};

        // create layout
        DArrayLayout<2> layout(MPI_COMM_WORLD, layout_size, is_periodic);

        // create array 
        std::array<int, 2> array_size = {3*5 + 1, 9*5+1}; 
        std::array<int, 2> nhalo_out  = {2, 2};
        std::array<int, 2> nhalo_in   = {4, 2};

        auto fun = [&] () { 
            DArray<double, 2> a(layout, array_size, nhalo_out, nhalo_in); 
            return 1;
        };

        REQUIRE_THROWS( fun() );
    }

    SECTION("periodic = true, true") {
        std::array<int, 2> is_periodic = {true, true};

        // create layout
        DArrayLayout<2> layout(MPI_COMM_WORLD, layout_size, is_periodic);

        // create array 
        std::array<int, 2> array_size = {3*5, 9*6}; 
        std::array<int, 2> nhalo_out  = {2, 2};
        std::array<int, 2> nhalo_in   = {4, 4};
        DArray<double, 2> a(layout, array_size, nhalo_out, nhalo_in); 

        // test indexing
        a(0, 0) = 1; REQUIRE( a(0, 0) == 1 );
        a(1, 2) = 2; REQUIRE( a(1, 2) == 2 );
        REQUIRE_THROWS( a(-5,  -5) );
        REQUIRE_THROWS( a( 0,  14) );
        REQUIRE_THROWS( a(13,   0) );

        // test nhalo_points
        REQUIRE( a.nhalo_points(Boundary::LEFT,   0) == 4 );
        REQUIRE( a.nhalo_points(Boundary::RIGHT,  0) == 4 );
        REQUIRE( a.nhalo_points(Boundary::LEFT,   1) == 4 );
        REQUIRE( a.nhalo_points(Boundary::RIGHT,  1) == 4 );

        // test raw_size
        std::array<int, 2> expected_1 = {13, 14};
        REQUIRE( a.raw_size() == expected_1 );

        // test size
        std::array<int, 2> expected_2 = {5, 6};
        REQUIRE( a.size() == expected_2 );
        REQUIRE( a.size(0) == expected_2[0] );
        REQUIRE( a.size(1) == expected_2[1] );

        // test nelements
        REQUIRE( a.nelements() == 13*14 );   
    }
}