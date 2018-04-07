#include "DArrays.hpp"
#include <catch.hpp>
#include <array>
#include <iostream>

TEST_CASE("1D tests - array", "[1D-tests]") {

    // use this grid layout for tests
    std::array<int, 1> layout_size = {27};

    SECTION("case 1") {
        std::array<int, 1> is_periodic = {false};

        // create layout
        DArrays::DArrayLayout<1> layout(MPI_COMM_WORLD, 
                                        layout_size,
                                        is_periodic);

        // create array 
        std::array<int, 1> array_size = {27*5}; 
        std::array<int, 1> nhalo_out  = {2};
        int                nhalo_in   =  4;
        DArrays::DArray<double, 1> a(layout, array_size, nhalo_out, nhalo_in); 

        // every processor in the domain interior indexes 
        // its local part plus nhalo_in boundary points 
        if ( !layout.is_on_boundary() ) {
            for (int i : {-4, -3, -2, -1,  0, 1, 2, 3, 4,  5, 6, 7, 8}) {
                a(i) = i; 
                REQUIRE( a(i) == i );
            }
            REQUIRE_THROWS( a(-5) );
            REQUIRE_THROWS( a( 9) );

            // test nhalo
            REQUIRE( a.nhalo(DArrays::BoundaryTag::LEFT,   0) == 4 );
            REQUIRE( a.nhalo(DArrays::BoundaryTag::RIGHT,  0) == 4 );
        }

        // processors on the left and right domain parts might 
        // have a different number of points they can index
        if ( layout.is_on_boundary(DArrays::BoundaryTag::LEFT, 0) ) {
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

            // test nhalo
            REQUIRE( a.nhalo(DArrays::BoundaryTag::LEFT,   0) == 2 );
            REQUIRE( a.nhalo(DArrays::BoundaryTag::RIGHT,  0) == 4 );
        }

        if ( layout.is_on_boundary(DArrays::BoundaryTag::RIGHT, 0) ) {
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

            // test nhalo
            REQUIRE( a.nhalo(DArrays::BoundaryTag::LEFT,   0) == 4 );
            REQUIRE( a.nhalo(DArrays::BoundaryTag::RIGHT,  0) == 2 );
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

    SECTION("case 2") {
        std::array<int, 1> is_periodic = {true};

        // create layout
        DArrays::DArrayLayout<1> layout(MPI_COMM_WORLD, 
                                        layout_size,
                                        is_periodic);

        // create array 
        std::array<int, 1> array_size = {27*5};
        std::array<int, 1> nhalo_out  = {2};
        int                nhalo_in   =  4;
        DArrays::DArray<double, 1> a(layout, array_size, nhalo_out, nhalo_in); 

        // all processors have nhalo_in
        for (int i : {-4, -3, -2, -1,  0, 1, 2, 3, 4,  5, 6, 7, 8}) {
            a(i) = i; 
            REQUIRE( a(i) == i );
        }
        REQUIRE_THROWS( a(-5) );
        REQUIRE_THROWS( a( 9) );

        // test nhalo
        REQUIRE( a.nhalo(DArrays::BoundaryTag::LEFT,  0) == 4 );
        REQUIRE( a.nhalo(DArrays::BoundaryTag::RIGHT, 0) == 4 );

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