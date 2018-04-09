#include "DArrays.hpp"
#include <catch.hpp>
#include <iostream>
#include <array>
#include <map>

// import all
using namespace DArrays;

TEST_CASE("dlayout", "test_1") {

    // use this layout size for tests
    std::array<int, 1> layout_size = {27};

    // get all halo boundaries for 1d layout
    auto boundaries = std::get<1>(DArrays::_halospeclist);

    // check numeric value for a few processors at key locations
    // this is the processor layout. Note MPI use row major layout

    //  26 |  0  1  ... 26 |  0

    // array of processors to be tested
    std::array<int,                3> procs    = {0, 10, 26}; 

    // expected neighbours
    std::array<std::array<int, 2>, 3> expected = {{{26,  1},   //  0
                                                   { 9, 11},   // 10
                                                   {25,  0}}}; // 26
                                                                                                      
    for (auto periodic1 : {true, false}) {
        // construct layout
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

        // for every processors we want to check the neighbours of
        for (auto j : LinRange(procs.size())) {
            // if it is actually the processor we want to check
            if (layout.rank() == procs[j]) {
                // for every neighbour
                for (auto i : LinRange(boundaries.size())) {
                    // we only check if there is a neighbour, otherwise it should be null
                    if ( layout.has_neighbour_at(boundaries[i]) ) {
                        REQUIRE( layout.rank_of_neighbour_at(boundaries[i]) == expected[j][i] );
                    } else {
                        REQUIRE( layout.rank_of_neighbour_at(boundaries[i]) == MPI_PROC_NULL );
                    }
                }
            }
        }                                                      
    }
}

TEST_CASE("2D tests", "[2D-tests]") {

    // use this layout size for tests
    std::array<int, 2> layout_size = {3, 9};

    // get halo boundaries for 2d layout
    auto boundaries = std::get<2>(DArrays::_halospeclist);

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
    // these N, S, W, E neighbours
    std::array<std::array<int, 4>, 9> expected = {{{18,  9,  8,  1},   //  0
                                                   { 9,  0, 26, 19},   // 18
                                                   {26, 17,  7,  0},   //  8 
                                                   { 4, 22, 12, 14},   // 13
                                                   {14,  5, 22, 24},   // 23
                                                   {17,  8, 25, 18},   // 26
                                                   { 8, 26, 16,  9},   // 17
                                                   {22, 13,  3,  5},   //  4
                                                   { 0, 18, 17, 10}}}; //  9


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

            // for every processors we want to check
            for (auto j : LinRange(procs.size())) {
                // if it is actually the processor we want to check
                if (layout.rank() == procs[j]) {
                    // for every neighbour
                    for (auto i : LinRange(boundaries.size())) {
                        // we only check if there is a neighbour, otherwise it should be null
                        if ( layout.has_neighbour_at(boundaries[i]) ) {
                            REQUIRE( layout.rank_of_neighbour_at(boundaries[i]) == expected[j][i] );
                        } else {
                            REQUIRE( layout.rank_of_neighbour_at(boundaries[i]) == MPI_PROC_NULL );
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
    auto boundaries = std::get<3>(DArrays::_halospeclist);

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
    std::array<int,                 4> procs    = {13, 0, 24, 8}; 
    // this follows the order defined in haloregion.hpp                     
    std::array<std::array<int, 26>, 4> expected = {{{ 4, 22, 10, 16, 12, 14},   // 13
                                                    {18,  9,  6,  3,  2,  1},   // 0
                                                    {15,  6, 21, 18, 26, 25},
                                                    {26, 17,  5,  2,  7,  6}}}; // 24

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
                            // we only check if there is a neighbour, otherwise it should be null
                            if ( layout.has_neighbour_at(boundaries[i]) ) {
                                REQUIRE( layout.rank_of_neighbour_at(boundaries[i]) == expected[j][i] );
                            } else {
                                REQUIRE( layout.rank_of_neighbour_at(boundaries[i]) == MPI_PROC_NULL );
                            }
                        }
                    }
                }
            }
        }
    }
}
