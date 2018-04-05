#include "DArrays.hpp"
#include <catch.hpp>
#include <iostream>
#include <array>
#include <map>

TEST_CASE("1D tests", "[1D-tests]") {

    // use this layout size for tests
    std::array<int, 1> layout_size = {27};

    // get halo boundaries for 1d layout
    auto boundaries = DArrays::AllBoundaries<1>();

    // check numeric value for a few processors at key locations
    // this is the processor layout. Note MPI use row major layout

    //  26 |  0  1  ... 26 |  0

    // array of processors to be tested
    std::array<int,                3> procs    = {0, 10, 26}; 
    // this follows the order defined in haloregion.hpp                     
    std::array<std::array<int, 2>, 3> expected = {{{26,  1},   //  0
                                                   { 9, 11},   // 10
                                                   {25,  0}}}; // 26

    // initialize expected value                          rank is_periodic  is_on_boundary
    std::map<std::pair<int, bool>, bool> expected_1 = { {{ 0,   true},      false},
                                                        {{ 0,  false},      true},
                                                        {{10,   true},      false},
                                                        {{10,  false},      false},
                                                        {{26,   true},      false},
                                                        {{26,  false},      true}};

    // initialize expected value                          rank is_periodic  is_on_boundary(LEFT, 0)
    std::map<std::pair<int, bool>, bool> expected_2 = { {{ 0,   true},      false},
                                                        {{ 0,  false},      true},
                                                        {{10,   true},      false},
                                                        {{10,  false},      false},
                                                        {{26,   true},      false},
                                                        {{26,  false},      false}};    

    // initialize expected value                          rank is_periodic  is_on_boundary(RIGHT, 0)
    std::map<std::pair<int, bool>, bool> expected_3 = { {{ 0,   true},      false},
                                                        {{ 0,  false},      false},
                                                        {{10,   true},      false},
                                                        {{10,  false},      false},
                                                        {{26,   true},      false},
                                                        {{26,  false},      true}};                                                                                                            

    for (auto periodic1 : {true, false}) {
        std::array<int, 1>  is_periodic = {periodic1};
        DArrays::DArrayLayout<1> layout(MPI_COMM_WORLD, layout_size, is_periodic);

        // check is_periodic
        for (auto dim : LinRange(1))
            REQUIRE( layout.is_periodic(dim) == is_periodic[dim] ); 

        REQUIRE_THROWS( layout.is_periodic(-1) );
        REQUIRE_THROWS( layout.is_periodic( 1) );

        // check grid size
        for (auto dim : LinRange(1))
            REQUIRE( layout.size(dim) == layout_size[dim] );
            
        REQUIRE_THROWS( layout.size(-1) );
        REQUIRE_THROWS( layout.size( 1) );    

        // check is_on_boundary functions, only for the ranks we have defined
        std::pair<int, bool> key = {layout.rank(), layout.is_periodic(0)};

        if (expected_1.count(key) == 1)
            REQUIRE( layout.is_on_boundary() == expected_1[key] );

        if (expected_2.count(key) == 1)
            REQUIRE( layout.is_on_boundary(DArrays::BoundaryTag::LEFT, 0) == expected_2[key] );

        if (expected_3.count(key) == 1)
            REQUIRE( layout.is_on_boundary(DArrays::BoundaryTag::RIGHT, 0) == expected_3[key] );

        // for every processors we want to check
        for (auto j : LinRange(procs.size())) {
            // if it is actually the processor we want to check
            if (layout.rank() == procs[j]) {
                // for every neighbour
                for (auto i : LinRange(boundaries.size())) {
                    // we only check if there is a neighbour, otherwise it should throw
                    if ( layout.has_neighbour_at(boundaries[i]) ) {
                        REQUIRE( layout.rank_of_neighbour_at(boundaries[i]) == expected[j][i] );
                    } else {
                        REQUIRE_THROWS( layout.rank_of_neighbour_at(boundaries[i]) );
                    }
                }
            }
        }                                                      
        // for every neighbour
    }
}

TEST_CASE("2D tests", "[2D-tests]") {

    // use this layout size for tests
    std::array<int, 2> layout_size = {3, 9};

    // get halo boundaries for 2d layout
    auto boundaries = DArrays::AllBoundaries<2>();

    // check numeric value for a few processors at key locations
    // this is the processor layout. Note MPI use row major layout

    //  26 | 18 19 20 21 22 23 24 25 26 | 18    
    //  -  +  -  -  -  -  -  -  -  -  - +  -
    //   8 |  0  1  2  3  4  5  6  7  8 |  0
    //  17 |  9 10 11 12 13 14 15 16 17 |  9
    //  26 | 18 19 20 21 22 23 24 25 26 | 18
    //   - +  -  -  -  -  -  -  -  -  - +  -
    //   8 |  0  1  2  3  4  5  6  7  8 |  0

    // array of processors to be tested
    std::array<int,                9> procs    = {0, 18, 8, 13, 23, 26, 17, 4, 9}; 
    // this follows the order defined in haloregion.hpp                     
    std::array<std::array<int, 8>, 9> expected = {{{26, 18, 19,  8,  1, 17,  9, 10},   //  0
                                                   {17,  9, 10, 26, 19,  8,  0,  1},   // 18
                                                   {25, 26, 18,  7,  0, 16, 17,  9},   //  8 
                                                   { 3,  4,  5, 12, 14, 21, 22, 23},   // 13
                                                   {13, 14, 15, 22, 24,  4,  5,  6},   // 23
                                                   {16, 17,  9, 25, 18,  7,  8,  0},   // 26
                                                   { 7,  8,  0, 16,  9, 25, 26, 18},   // 17
                                                   {21, 22, 23,  3,  5, 12, 13, 14},   //  4
                                                   { 8,  0,  1, 17, 10, 26, 18, 19}}}; //  9

    // we do not text everything because it would be too tedious, we instead check at random
    // initialize expected value                                rank is_periodic(0) is_periodic(1) is_on_boundary
    std::map<std::tuple<int, bool, bool>, bool> expected_1 = { {{ 0,   true,          true},        false},
                                                               {{ 0,  false,         false},         true},
                                                               {{ 0,  false,          true},         true},
                                                               {{ 0,   true,         false},         true},                                                
                                                               {{ 8,   true,          true},        false},
                                                               {{ 8,  false,         false},         true},
                                                               {{ 8,  false,          true},         true},
                                                               {{ 8,   true,         false},         true},                                                
                                                               {{26,   true,          true},        false},
                                                               {{26,  false,         false},         true},
                                                               {{26,  false,          true},         true},
                                                               {{26,   true,         false},         true},                                    
                                                               {{18,   true,          true},        false},
                                                               {{18,  false,         false},         true},
                                                               {{18,  false,          true},         true},
                                                               {{18,   true,         false},         true},
                                                               {{ 4,   true,          true},        false},
                                                               {{ 4,  false,         false},         true},
                                                               {{ 4,  false,          true},         true},
                                                               {{ 4,   true,         false},        false},
                                                               {{14,   true,          true},        false},
                                                               {{14,  false,         false},        false},
                                                               {{14,  false,          true},        false},
                                                               {{14,   true,         false},        false},
                                                               {{17,   true,          true},        false},
                                                               {{17,  false,         false},         true},
                                                               {{17,  false,          true},        false},
                                                               {{17,   true,         false},         true}}; 

    // initialize expected value                                rank is_periodic(0) is_periodic(1) is_on_boundary(LEFT, 0)
    std::map<std::tuple<int, bool, bool>, bool> expected_2 = { {{ 0,  false,         false},        true},
                                                               {{ 1,   true,         false},       false},
                                                               {{ 1,   true,          true},       false},
                                                               {{ 1,   false,         true},        true},
                                                               {{ 1,   false,        false},        true},
                                                               {{11,  false,         false},       false},
                                                               {{18,  false,         false},       false} };

   // initialize expected value                                rank is_periodic(0) is_periodic(1) is_on_boundary(RIGHT, 0)
    std::map<std::tuple<int, bool, bool>, bool> expected_3 = { {{18,  false,         false},        true},
                                                               {{22,   true,         false},       false},
                                                               {{22,   true,          true},       false},
                                                               {{22,  false,          true},        true},
                                                               {{22,  false,         false},        true},
                                                               {{11,  false,         false},       false},
                                                               {{8,   false,         false},       false} };       

    // initialize expected value                                rank is_periodic(0) is_periodic(1) is_on_boundary(LEFT, 1)
    std::map<std::tuple<int, bool, bool>, bool> expected_4 = { {{ 0,  false,         false},        true},
                                                               {{ 9,  false,         false},        true},
                                                               {{ 9,  false,          true},       false},
                                                               {{ 9,   true,         false},        true},
                                                               {{ 9,   true,          true},       false},
                                                               {{18,  false,         false},        true},
                                                               {{18,  false,          true},       false},
                                                               {{18,   true,         false},        true},
                                                               {{18,   true,          true},       false},
                                                               {{ 3,  false,         false},       false},
                                                               {{ 3,  false,          true},       false},
                                                               {{ 3,   true,         false},       false},
                                                               {{ 3,   true,          true},       false},
                                                               {{12,  false,         false},       false},
                                                               {{12,  false,          true},       false},
                                                               {{12,   true,         false},       false},
                                                               {{12,   true,          true},       false},   
                                                               {{21,  false,         false},       false},
                                                               {{21,  false,          true},       false},
                                                               {{21,   true,         false},       false},
                                                               {{21,   true,          true},       false},
                                                               {{ 8,  false,         false},       false},
                                                               {{ 8,  false,          true},       false},
                                                               {{ 8,   true,         false},       false},
                                                               {{ 8,   true,          true},       false},
                                                               {{17,  false,         false},       false},
                                                               {{17,  false,          true},       false},
                                                               {{17,   true,         false},       false},
                                                               {{17,   true,          true},       false},
                                                               {{26,  false,         false},       false},
                                                               {{26,  false,          true},       false},
                                                               {{26,   true,         false},       false},
                                                               {{26,   true,          true},       false}};

    // initialize expected value                                rank is_periodic(0) is_periodic(1) is_on_boundary(RIGHT, 1)
    std::map<std::tuple<int, bool, bool>, bool> expected_5 = { {{ 0,  false,         false},       false},
                                                               {{ 8,  false,         false},        true},
                                                               {{ 8,  false,          true},       false},
                                                               {{ 8,   true,         false},        true},
                                                               {{ 8,   true,          true},       false}};                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             

    for (auto periodic1 : {true, false}) {
        for (auto periodic2 : {true, false}) {
            std::array<int, 2> is_periodic = {periodic1, periodic2};
            DArrays::DArrayLayout<2> layout(MPI_COMM_WORLD, layout_size, is_periodic);

            // check is_periodic
            for (auto dim : LinRange(2))
                REQUIRE( layout.is_periodic(dim) == is_periodic[dim] ); 

            REQUIRE_THROWS( layout.is_periodic(-1) );
            REQUIRE_THROWS( layout.is_periodic( 2) );

            // check grid size
            for (auto dim : LinRange(2))
                REQUIRE( layout.size(dim) == layout_size[dim] );
            
            REQUIRE_THROWS( layout.size(-1) );
            REQUIRE_THROWS( layout.size( 2) );

            // check is_on_boundary functions, only for the ranks we have defined
            std::tuple<int, bool, bool> key = {layout.rank(), layout.is_periodic(0), layout.is_periodic(1)};
            
            if (expected_1.count(key) == 1)
                REQUIRE( layout.is_on_boundary() == expected_1[key] );      

            if (expected_2.count(key) == 1)
                REQUIRE( layout.is_on_boundary(DArrays::BoundaryTag::LEFT,  0) == expected_2[key] );        

            if (expected_3.count(key) == 1)
                REQUIRE( layout.is_on_boundary(DArrays::BoundaryTag::RIGHT, 0) == expected_3[key] );          

            if (expected_4.count(key) == 1)
                REQUIRE( layout.is_on_boundary(DArrays::BoundaryTag::LEFT, 1) == expected_4[key] );               

            if (expected_5.count(key) == 1)
                REQUIRE( layout.is_on_boundary(DArrays::BoundaryTag::RIGHT, 1) == expected_5[key] );     

            // for every processors we want to check
            for (auto j : LinRange(procs.size())) {
                // if it is actually the processor we want to check
                if (layout.rank() == procs[j]) {
                    // for every neighbour
                    for (auto i : LinRange(boundaries.size())) {
                        // we only check if there is a neighbour, otherwise it should throw
                        if ( layout.has_neighbour_at(boundaries[i]) ) {
                            REQUIRE( layout.rank_of_neighbour_at(boundaries[i]) == expected[j][i] );
                        } else {
                            REQUIRE_THROWS( layout.rank_of_neighbour_at(boundaries[i]) );
                        }
                    }
                }
            }
        }
    }
}

TEST_CASE("3D tests", "[3D-tests]") {

    // use this layout size for tests
    std::array<int, 3> layout_size = {3, 3, 3};

    // get halo boundaries for 3d layout
    auto boundaries = DArrays::AllBoundaries<3>();

    // check numeric value for a few processors at key locations
    // this is the processor layout. Note MPI use row major layout

    //   24 | 18 21 24 | 18   25 | 19 22 25 | 19  26 | 20 23 26 | 20    
    //    - +  -  -  - + -     - +  -  -  - + -    - +  -  -  - + -  
    //    6 |  0  3  6 |  0    7 |  1  4  7 |  1   8 |  2  5  8 |  2    
    //   15 |  9 12 15 |  9   16 | 10 13 16 | 10  17 | 11 14 17 | 11    
    //   24 | 18 21 24 | 18   25 | 19 22 25 | 19  26 | 20 23 26 | 20    
    //    - +  -  -  - + -     - +  -  -  - + -    - +  -  -  - + -  
    //    6 |  0  3  6 |  0    7 |  1  4  7 |  1   8 |  2  5  8 |  2  

    // array of processors to be tested
    std::array<int,                 3> procs    = {13, 0, 24}; 
    // this follows the order defined in haloregion.hpp                     
    std::array<std::array<int, 26>, 3> expected = {{{ 0,  9, 18,  3, 12, 21,  6, 15, 24,  1, 10, 19,  4, 22,  7, 16, 25,  2, 11, 20,  5, 14, 23,  8, 17, 26},   // 13
                                                    {26,  8, 17, 20,  2, 11, 23,  5, 14, 24,  6, 15, 18,  9, 21,  3, 12, 25,  7, 16, 19,  1, 10, 22,  4, 13},   // 0
                                                    {14, 23,  5, 17, 26,  8, 11, 20,  2, 12, 21,  3, 15,  6,  9, 18,  0, 13, 22,  4, 16, 25,  7, 10, 19,  1}}}; // 24

    for (auto periodic1 : {true, false}) {
        for (auto periodic2 : {true, false}) {
            for (auto periodic3 : {true, false}) {
                std::array<int, 3> is_periodic = {periodic1, periodic2, periodic3};
                DArrays::DArrayLayout<3> layout(MPI_COMM_WORLD, layout_size, is_periodic);

                // check is_periodic
                for (auto dim : LinRange(3))
                    REQUIRE( layout.is_periodic(dim) == is_periodic[dim] ); 

                REQUIRE_THROWS( layout.is_periodic(-1) );
                REQUIRE_THROWS( layout.is_periodic( 3) );

                // check grid size
                for (auto dim : LinRange(3))
                    REQUIRE( layout.size(dim) == layout_size[dim] );
                
                REQUIRE_THROWS( layout.size(-1) );
                REQUIRE_THROWS( layout.size( 3) );

                // for 3d we do not check the boundary functions, because the code 
                // is generic and it's been tested for 2D

                // for every processors we want to check
                for (auto j : LinRange(procs.size())) {
                    // if it is actually the processor we want to check
                    if (layout.rank() == procs[j]) {
                        // for every neighbour
                        for (auto i : LinRange(boundaries.size())) {
                            // we only check if there is a neighbour, otherwise it should throw
                            if ( layout.has_neighbour_at(boundaries[i]) ) {
                                REQUIRE( layout.rank_of_neighbour_at(boundaries[i]) == expected[j][i] );
                            } else {
                                REQUIRE_THROWS( layout.rank_of_neighbour_at(boundaries[i]) );
                            }
                        }
                    }
                }
            }
        }
    }
}
