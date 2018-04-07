#include "DArrays.hpp"
#include <catch.hpp>
#include <array>
#include <iostream>

TEST_CASE("1D tests - subarray", "[1D-tests]") {

    // use this grid layout for tests
    std::array<int, 1> layout_size = {27};

    SECTION("periodic = true") {
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

        // raw index   0  1  2  3   4 5 6 7 8   9 10 11 12
        // element    -4 -3 -2 -1 | 0 1 2 3 4 | 5  6  7  8
        SECTION("left boundary") {
            DArrays::Boundary<1> bnd = {DArrays::BoundaryTag::LEFT};
            DArrays::SubArray<double, 1> sub_send(a, bnd, DArrays::BoundaryIntent::SEND);
            DArrays::SubArray<double, 1> sub_recv(a, bnd, DArrays::BoundaryIntent::RECV);
            
            // send
            std::array<int, 1> expected_1 = {4};
            REQUIRE( sub_send.size() == expected_1 );

            std::array<int, 1> expected_2 = {4};
            REQUIRE( sub_send.raw_origin() == expected_2 );
            
            // recv
            std::array<int, 1> expected_3 = {4};
            REQUIRE( sub_recv.size() == expected_3 );

            std::array<int, 1> expected_4 = {0};
            REQUIRE( sub_recv.raw_origin() == expected_4 );
        }

        SECTION("right boundary") {
            DArrays::Boundary<1> bnd = {DArrays::BoundaryTag::RIGHT};
            DArrays::SubArray<double, 1> sub_send(a, bnd, DArrays::BoundaryIntent::SEND);
            DArrays::SubArray<double, 1> sub_recv(a, bnd, DArrays::BoundaryIntent::RECV);

            // send
            std::array<int, 1> expected_1 = {4};
            REQUIRE( sub_send.size() == expected_1 );

            std::array<int, 1> expected_2 = {5};
            REQUIRE( sub_send.raw_origin() == expected_2 );

            std::array<int, 1> expected_3 = {4};
            REQUIRE( sub_recv.size() == expected_3 );

            std::array<int, 1> expected_4 = {9};
            REQUIRE( sub_recv.raw_origin() == expected_4 );
        }
    }

    SECTION("periodic = false") {
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

        // in domain
        // raw index   0  1  2  3   4 5 6 7 8   9 10 11 12
        // element    -4 -3 -2 -1 | 0 1 2 3 4 | 5  6  7  8
        // left boundary
        // raw index   0  1   2 3 4 5 6   7  8 9 10
        // element    -2 -1 | 0 1 2 3 4 | 5  6 7  8
        SECTION("left boundary") {
            DArrays::Boundary<1> bnd = {DArrays::BoundaryTag::LEFT};
            DArrays::SubArray<double, 1> sub_send(a, bnd, DArrays::BoundaryIntent::SEND);
            DArrays::SubArray<double, 1> sub_recv(a, bnd, DArrays::BoundaryIntent::RECV);
            
            // send
            if ( layout.has_neighbour_at(DArrays::BoundaryTag::LEFT, 0) ) {
                std::array<int, 1> expected_1 = {4};
                REQUIRE( sub_send.size() == expected_1 );

                std::array<int, 1> expected_2 = {4};
                REQUIRE( sub_send.raw_origin() == expected_2 );
            } else {
                std::array<int, 1> expected_1 = {2};
                REQUIRE( sub_send.size() == expected_1 );

                std::array<int, 1> expected_2 = {2};
                REQUIRE( sub_send.raw_origin() == expected_2 );
            }
            
            // recv
            if ( layout.has_neighbour_at(DArrays::BoundaryTag::LEFT, 0) ) {
                std::array<int, 1> expected_1 = {4};
                REQUIRE( sub_recv.size() == expected_1 );

                std::array<int, 1> expected_2 = {0};
                REQUIRE( sub_recv.raw_origin() == expected_2 );
            } else {
                std::array<int, 1> expected_1 = {2};
                REQUIRE( sub_recv.size() == expected_1 );

                std::array<int, 1> expected_2 = {0};
                REQUIRE( sub_recv.raw_origin() == expected_2 );
            }
        }
        
        // in domain
        // raw index   0  1  2  3   4 5 6 7 8   9 10 11 12
        // element    -4 -3 -2 -1 | 0 1 2 3 4 | 5  6  7  8
        // right boundary
        // raw index   0  1  2  3   4 5 6 7 8   9 10
        // element    -4 -3 -2 -1 | 0 1 2 3 4 | 5  6
        // left boundary
        // raw index   0  1   2 3 4 5 6   7  8 9 10
        // element    -2 -1 | 0 1 2 3 4 | 5  6 7  8
        SECTION("right boundary") {
            DArrays::Boundary<1> bnd = {DArrays::BoundaryTag::RIGHT};
            DArrays::SubArray<double, 1> sub_send(a, bnd, DArrays::BoundaryIntent::SEND);
            DArrays::SubArray<double, 1> sub_recv(a, bnd, DArrays::BoundaryIntent::RECV);
            
            // send
            if ( layout.has_neighbour_at(DArrays::BoundaryTag::RIGHT, 0) ) {
                std::array<int, 1> expected_1 = {4};
                REQUIRE( sub_send.size() == expected_1 );

                std::array<int, 1> expected_2;
                if ( layout.is_on_boundary(DArrays::BoundaryTag::LEFT, 0) ) {
                    expected_2 = {3};
                } else {
                    expected_2 = {5};
                }
                REQUIRE( sub_send.raw_origin() == expected_2 );
            } else {
                std::array<int, 1> expected_1 = {2};
                REQUIRE( sub_send.size() == expected_1 );

                std::array<int, 1> expected_2 = {7};
                REQUIRE( sub_send.raw_origin() == expected_2 );
            }
            
            // recv
            if ( layout.has_neighbour_at(DArrays::BoundaryTag::RIGHT, 0) ) {
                std::array<int, 1> expected_1 = {4};
                REQUIRE( sub_recv.size() == expected_1 );

                std::array<int, 1> expected_2;
                if ( layout.is_on_boundary(DArrays::BoundaryTag::LEFT, 0) ) {
                    expected_2 = {7};
                } else {
                    expected_2 = {9};
                }
                REQUIRE( sub_recv.raw_origin() == expected_2 );
            } else {
                std::array<int, 1> expected_1 = {2};
                REQUIRE( sub_recv.size() == expected_1 );

                std::array<int, 1> expected_2 = {9};
                REQUIRE( sub_recv.raw_origin() == expected_2 );
            }
        }

    }
}
