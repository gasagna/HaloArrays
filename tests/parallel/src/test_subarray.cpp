#include "DArrays.hpp"
#include <catch.hpp>
#include <array>
#include <iostream>

// import all
using namespace DArrays;

TEST_CASE("subarray - 1D", "test_1") {

    // use this grid layout for tests
    std::array<int, 1> layout_size = {27};

    SECTION("periodic = true") {
        std::array<int, 1> is_periodic = {true};

        // create layout
        DArrayLayout<1> layout(MPI_COMM_WORLD, 
                               layout_size,
                               is_periodic);

        // create array 
        std::array<int, 1> array_size = {27*5}; 
        std::array<int, 1> nhalo_out  = {2};
        std::array<int, 1> nhalo_in   = {4};
        DArray<double, 1> a(layout, array_size, nhalo_out, nhalo_in); 

        // raw index   0  1  2  3   4 5 6 7 8   9 10 11 12
        // element    -4 -3 -2 -1 | 0 1 2 3 4 | 5  6  7  8
        SECTION("left boundary") {
            HaloRegionSpec<1> bnd(Boundary::LEFT);
            SubArray<double, 1> sub_send(a, bnd, HaloIntent::SEND);
            SubArray<double, 1> sub_recv(a, bnd, HaloIntent::RECV);
            
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
            HaloRegionSpec<1> bnd(Boundary::RIGHT);
            SubArray<double, 1> sub_send(a, bnd, HaloIntent::SEND);
            SubArray<double, 1> sub_recv(a, bnd, HaloIntent::RECV);

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
        DArrayLayout<1> layout(MPI_COMM_WORLD, 
                               layout_size,
                               is_periodic);

        // create array 
        std::array<int, 1> array_size = {27*5}; 
        std::array<int, 1> nhalo_out  = {2};
        std::array<int, 1> nhalo_in   = {4};
        DArray<double, 1> a(layout, array_size, nhalo_out, nhalo_in); 

        // in domain
        // raw index   0  1  2  3   4 5 6 7 8   9 10 11 12
        // element    -4 -3 -2 -1 | 0 1 2 3 4 | 5  6  7  8
        // left boundary
        // raw index   0  1   2 3 4 5 6   7  8 9 10
        // element    -2 -1 | 0 1 2 3 4 | 5  6 7  8
        SECTION("left boundary") {
            HaloRegionSpec<1> bnd(Boundary::LEFT);
            SubArray<double, 1> sub_send(a, bnd, HaloIntent::SEND);
            SubArray<double, 1> sub_recv(a, bnd, HaloIntent::RECV);
            
            // send
            if ( layout.has_neighbour_at(Boundary::LEFT, 0) ) {
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
            if ( layout.has_neighbour_at(Boundary::LEFT, 0) ) {
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
            HaloRegionSpec<1> bnd(Boundary::RIGHT);
            SubArray<double, 1> sub_send(a, bnd, HaloIntent::SEND);
            SubArray<double, 1> sub_recv(a, bnd, HaloIntent::RECV);
            
            // send
            if ( layout.has_neighbour_at(Boundary::RIGHT, 0) ) {
                std::array<int, 1> expected_1 = {4};
                REQUIRE( sub_send.size() == expected_1 );

                std::array<int, 1> expected_2;
                if ( layout.has_neighbour_at(Boundary::LEFT, 0) ) {
                    expected_2 = {5};
                } else {
                    expected_2 = {3};
                }
                REQUIRE( sub_send.raw_origin() == expected_2 );
            } else {
                std::array<int, 1> expected_1 = {2};
                REQUIRE( sub_send.size() == expected_1 );

                std::array<int, 1> expected_2 = {7};
                REQUIRE( sub_send.raw_origin() == expected_2 );
            }
            
            // recv
            if ( layout.has_neighbour_at(Boundary::RIGHT, 0) ) {
                std::array<int, 1> expected_1 = {4};
                REQUIRE( sub_recv.size() == expected_1 );

                std::array<int, 1> expected_2;
                if ( layout.has_neighbour_at(Boundary::LEFT, 0) ) {
                    expected_2 = {9};
                } else {
                    expected_2 = {7};
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

TEST_CASE("subarray 2D", "test_2") {

    // use this grid layout for tests
    // 0   1  2  3  4  5  6  7  8
    // 9  10 11 12 13 14 15 16 17
    // 18 19 20 21 22 23 24 25 26
    std::array<int, 2> layout_size = {3, 9};

    SECTION("periodic = true, true") {
        std::array<int, 2> is_periodic = {true, true};

        // create layout
        DArrayLayout<2> layout(MPI_COMM_WORLD, layout_size, is_periodic);

        // create array 
        std::array<int, 2> array_size = {3*2, 9*2}; 
        std::array<int, 2> nhalo_out  = {1, 1};
        std::array<int, 2> nhalo_in   = {1, 1};
        DArray<double, 2> a(layout, array_size, nhalo_out, nhalo_in); 

        SECTION("L*") {
            HaloRegionSpec<2> bnd(Boundary::LEFT, Boundary::WILDCARD);
            SubArray<double, 2> sub_send(a, bnd, HaloIntent::SEND);
            SubArray<double, 2> sub_recv(a, bnd, HaloIntent::RECV);

            // send
            std::array<int, 2> expected_1 = {1, 4};
            REQUIRE( sub_send.size() == expected_1 );
            std::array<int, 2> expected_2 = {1, 0};
            REQUIRE( sub_send.raw_origin() == expected_2 );
            // recv
            std::array<int, 2> expected_3 = {1, 4};
            REQUIRE( sub_recv.size() == expected_3 );
            std::array<int, 2> expected_4 = {0, 0};
            REQUIRE( sub_recv.raw_origin() == expected_4 );
        }
        
        SECTION("R*") {
            HaloRegionSpec<2> bnd(Boundary::RIGHT, Boundary::WILDCARD);
            SubArray<double, 2> sub_send(a, bnd, HaloIntent::SEND);
            SubArray<double, 2> sub_recv(a, bnd, HaloIntent::RECV);

            // send
            std::array<int, 2> expected_1 = {1, 4};
            REQUIRE( sub_send.size() == expected_1 );
            std::array<int, 2> expected_2 = {2, 0};
            REQUIRE( sub_send.raw_origin() == expected_2 );
            // recv
            std::array<int, 2> expected_3 = {1, 4};
            REQUIRE( sub_recv.size() == expected_3 );
            std::array<int, 2> expected_4 = {3, 0};
            REQUIRE( sub_recv.raw_origin() == expected_4 );
        }
        
        SECTION("CL") {
            HaloRegionSpec<2> bnd(Boundary::CENTER, Boundary::LEFT);
            SubArray<double, 2> sub_send(a, bnd, HaloIntent::SEND);
            SubArray<double, 2> sub_recv(a, bnd, HaloIntent::RECV);

            // send
            std::array<int, 2> expected_1 = {2, 1};
            REQUIRE( sub_send.size() == expected_1 );
            std::array<int, 2> expected_2 = {1, 1};
            REQUIRE( sub_send.raw_origin() == expected_2 );
            // recv
            std::array<int, 2> expected_3 = {2, 1};
            REQUIRE( sub_recv.size() == expected_3 );
            std::array<int, 2> expected_4 = {1, 0};
            REQUIRE( sub_recv.raw_origin() == expected_4 );
        }

        SECTION("CR") {
            HaloRegionSpec<2> bnd(Boundary::CENTER, Boundary::RIGHT);
            SubArray<double, 2> sub_send(a, bnd, HaloIntent::SEND);
            SubArray<double, 2> sub_recv(a, bnd, HaloIntent::RECV);

            // send
            std::array<int, 2> expected_1 = {2, 1};
            REQUIRE( sub_send.size() == expected_1 );
            std::array<int, 2> expected_2 = {1, 2};
            REQUIRE( sub_send.raw_origin() == expected_2 );
            // recv
            std::array<int, 2> expected_3 = {2, 1};
            REQUIRE( sub_recv.size() == expected_3 );
            std::array<int, 2> expected_4 = {1, 3};
            REQUIRE( sub_recv.raw_origin() == expected_4 );
        }
    }    

    SECTION("periodic = true, false") {
        std::array<int, 2> is_periodic = {true, false};

        // create layout
        DArrayLayout<2> layout(MPI_COMM_WORLD, layout_size, is_periodic);

        // create array 
        std::array<int, 2> array_size = {3*5, 9*5}; 
        std::array<int, 2> nhalo_out  = {1, 2};
        std::array<int, 2> nhalo_in   = {3, 4};
        DArray<double, 2> a(layout, array_size, nhalo_out, nhalo_in); 

        SECTION("L*") {
            HaloRegionSpec<2> bnd(Boundary::LEFT, Boundary::WILDCARD);
            SubArray<double, 2> sub_send(a, bnd, HaloIntent::SEND);
            SubArray<double, 2> sub_recv(a, bnd, HaloIntent::RECV);

            SECTION("send") {
                std::array<int,                6> ranks    = {0, 1, 8, 18, 19, 26}; 
                std::array<std::array<int, 4>, 6> expected = {{{3, 11,  3, 0},   // 0
                                                               {3, 13,  3, 0},   // 1
                                                               {3, 11,  3, 0},   // 8
                                                               {3, 11,  3, 0},   // 18
                                                               {3, 13,  3, 0},   // 19
                                                               {3, 11,  3, 0}}}; // 26
                // for every test
                for (auto j : LinRange(6)) {
                    // if it is actually the processor we want to check
                    if (layout.rank() == ranks[j]) {
                        REQUIRE( sub_send.size(0)       == expected[j][0] );
                        REQUIRE( sub_send.size(1)       == expected[j][1] );
                        REQUIRE( sub_send.raw_origin(0) == expected[j][2] );
                        REQUIRE( sub_send.raw_origin(1) == expected[j][3] );
                    }
                }
            }
            
            SECTION("recv") {
                std::array<int,                6> ranks    = {0, 1, 8, 18, 19, 26}; 
                std::array<std::array<int, 4>, 6> expected = {{{3, 11,  0, 0},   // 0
                                                               {3, 13,  0, 0},   // 1
                                                               {3, 11,  0, 0},   // 8
                                                               {3, 11,  0, 0},   // 18
                                                               {3, 13,  0, 0},   // 19
                                                               {3, 11,  0, 0}}}; // 16
                // for every test
                for (auto j : LinRange(6)) {
                    // if it is actually the processor we want to check
                    if (layout.rank() == ranks[j]) {
                        REQUIRE( sub_recv.size(0)       == expected[j][0] );
                        REQUIRE( sub_recv.size(1)       == expected[j][1] );
                        REQUIRE( sub_recv.raw_origin(0) == expected[j][2] );
                        REQUIRE( sub_recv.raw_origin(1) == expected[j][3] );
                    }
                }
            }

        }

        SECTION("R*") {
            HaloRegionSpec<2> bnd(Boundary::RIGHT, Boundary::WILDCARD);
            SubArray<double, 2> sub_send(a, bnd, HaloIntent::SEND);
            SubArray<double, 2> sub_recv(a, bnd, HaloIntent::RECV);

            SECTION("send") {
                std::array<int,                6> ranks    = {0, 1, 8, 18, 19, 26}; 
                std::array<std::array<int, 4>, 6> expected = {{{3, 11,  5, 0},   // 0
                                                               {3, 13,  5, 0},   // 1
                                                               {3, 11,  5, 0},   // 8
                                                               {3, 11,  5, 0},   // 18
                                                               {3, 13,  5, 0},   // 19
                                                               {3, 11,  5, 0}}}; // 16
                // for every test
                for (auto j : LinRange(6)) {
                    // if it is actually the processor we want to check
                    if (layout.rank() == ranks[j]) {
                        REQUIRE( sub_send.size(0)       == expected[j][0] );
                        REQUIRE( sub_send.size(1)       == expected[j][1] );
                        REQUIRE( sub_send.raw_origin(0) == expected[j][2] );
                        REQUIRE( sub_send.raw_origin(1) == expected[j][3] );
                    }
                }
            }
            
            SECTION("recv") {
                std::array<int,                6> ranks    = {0, 1, 8, 18, 19, 26}; 
                std::array<std::array<int, 4>, 6> expected = {{{3, 11,  8, 0},   // 0
                                                               {3, 13,  8, 0},   // 1
                                                               {3, 11,  8, 0},   // 8
                                                               {3, 11,  8, 0},   // 18
                                                               {3, 13,  8, 0},   // 19
                                                               {3, 11,  8, 0}}}; // 16
                // for every test
                for (auto j : LinRange(6)) {
                    // if it is actually the processor we want to check
                    if (layout.rank() == ranks[j]) {
                        REQUIRE( sub_recv.size(0)       == expected[j][0] );
                        REQUIRE( sub_recv.size(1)       == expected[j][1] );
                        REQUIRE( sub_recv.raw_origin(0) == expected[j][2] );
                        REQUIRE( sub_recv.raw_origin(1) == expected[j][3] );
                    }
                }
            }
        }        

        SECTION("CL") {
            HaloRegionSpec<2> bnd(Boundary::CENTER, Boundary::LEFT);
            SubArray<double, 2> sub_send(a, bnd, HaloIntent::SEND);
            SubArray<double, 2> sub_recv(a, bnd, HaloIntent::RECV);

            SECTION("send") {
                std::array<int,                6> ranks    = {0, 1, 8, 18, 19, 26}; 
                std::array<std::array<int, 4>, 6> expected = {{{5, 2,  3, 2},   // 0
                                                               {5, 4,  3, 4},   // 1
                                                               {5, 4,  3, 4},   // 8
                                                               {5, 2,  3, 2},   // 18
                                                               {5, 4,  3, 4},   // 19
                                                               {5, 4,  3, 4}}}; // 26
                // for every test
                for (auto j : LinRange(6)) {
                    // if it is actually the processor we want to check
                    if (layout.rank() == ranks[j]) {
                        REQUIRE( sub_send.size(0)       == expected[j][0] );
                        REQUIRE( sub_send.size(1)       == expected[j][1] );
                        REQUIRE( sub_send.raw_origin(0) == expected[j][2] );
                        REQUIRE( sub_send.raw_origin(1) == expected[j][3] );
                    }
                }
            }
            
            SECTION("recv") {
                std::array<int,                6> ranks    = {0, 1, 8, 18, 19, 26}; 
                std::array<std::array<int, 4>, 6> expected = {{{5, 2,  3, 0},   // 0
                                                               {5, 4,  3, 0},   // 1
                                                               {5, 4,  3, 0},   // 8
                                                               {5, 2,  3, 0},   // 18
                                                               {5, 4,  3, 0},   // 19
                                                               {5, 4,  3, 0}}}; // 26
                // for every test
                for (auto j : LinRange(6)) {
                    // if it is actually the processor we want to check
                    if (layout.rank() == ranks[j]) {
                        REQUIRE( sub_recv.size(0)       == expected[j][0] );
                        REQUIRE( sub_recv.size(1)       == expected[j][1] );
                        REQUIRE( sub_recv.raw_origin(0) == expected[j][2] );
                        REQUIRE( sub_recv.raw_origin(1) == expected[j][3] );
                    }
                }
            }
        }       

        SECTION("CR") {
            HaloRegionSpec<2> bnd(Boundary::CENTER, Boundary::RIGHT);
            SubArray<double, 2> sub_send(a, bnd, HaloIntent::SEND);
            SubArray<double, 2> sub_recv(a, bnd, HaloIntent::RECV);

            SECTION("send") {
                std::array<int,                6> ranks    = {0, 1, 8, 18, 19, 26}; 
                std::array<std::array<int, 4>, 6> expected = {{{5, 4,  3, 5+2-4},   // 0
                                                               {5, 4,  3, 5+4-4},   // 1
                                                               {5, 2,  3, 5+4-2},   // 8
                                                               {5, 4,  3, 5+2-4},   // 18
                                                               {5, 4,  3, 5+4-4},   // 19
                                                               {5, 2,  3, 5+4-2}}}; // 26
                // for every test
                for (auto j : LinRange(6)) {
                    // if it is actually the processor we want to check
                    if (layout.rank() == ranks[j]) {
                        REQUIRE( sub_send.size(0)       == expected[j][0] );
                        REQUIRE( sub_send.size(1)       == expected[j][1] );
                        REQUIRE( sub_send.raw_origin(0) == expected[j][2] );
                        REQUIRE( sub_send.raw_origin(1) == expected[j][3] );
                    }
                }
            }
            
            SECTION("recv") {
                std::array<int,                6> ranks    = {0, 1, 8, 18, 19, 26}; 
                std::array<std::array<int, 4>, 6> expected = {{{5, 4, 3, 5+2},   // 0
                                                               {5, 4, 3, 5+4},   // 1
                                                               {5, 2, 3, 5+4},   // 8
                                                               {5, 4, 3, 5+2},   // 18
                                                               {5, 4, 3, 5+4},   // 19
                                                               {5, 2, 3, 5+4}}}; // 26
                // for every test
                for (auto j : LinRange(6)) {
                    // if it is actually the processor we want to check
                    if (layout.rank() == ranks[j]) {
                        REQUIRE( sub_recv.size(0)       == expected[j][0] );
                        REQUIRE( sub_recv.size(1)       == expected[j][1] );
                        REQUIRE( sub_recv.raw_origin(0) == expected[j][2] );
                        REQUIRE( sub_recv.raw_origin(1) == expected[j][3] );
                    }
                }
            }
        }    
    }          
}