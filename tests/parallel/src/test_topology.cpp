#include "DArrays.hpp"
#include <catch.hpp>
#include <array>
#include <iostream>

TEST_CASE("1D tests", "[1D-tests]") {

    // use this grid size for tests
    std::array<int, 1> proc_grid_size = {27};

    // get halo regions for 1d topology
    auto regions = DArrays::HaloRegions<1>();

    // check numeric value for a few processors at key locations
    // this is the processor layout. Note MPI use row major layout

    //  26 |  0  1  ... 26 |  0

    // array of processors to be tested
    std::array<int,                3> procs    = {0,  10, 26}; 
    // this follows the order defined in haloregion.hpp                     
    std::array<std::array<int, 2>, 3> expected = {{{26,  1},   //  0
                                                   { 9, 11},   // 10
                                                   {25,  0}}}; // 26

    for (auto periodic1 : {true, false}) {
        SECTION("is_periodic = xxx") {
            std::array<int, 1>  is_periodic = {periodic1};
            DArrays::Topology::DArrayTopology<1> topo(MPI_COMM_WORLD, 
                                                      proc_grid_size,
                                                      is_periodic);
            // for every processors we want to check
            for (auto j : LinRange(procs.size())) {
                // if it is actually the processor we want to check
                if (topo._comm_rank == procs[j]) {
                    // for every neighbour
                    for (auto i : LinRange(regions.size())) {
                        // we only check if there is a neighbour, otherwise it should throw
                        if ( !topo.has_grid_boundary_at(regions[i]) ) {
                            REQUIRE( topo.neighbour_proc_rank(regions[i]) == expected[j][i] );
                        } else {
                            REQUIRE_THROWS( topo.neighbour_proc_rank(regions[i]) );
                        }
                    }
                }
            }
        }
    }
}

TEST_CASE("2D tests", "[2D-tests]") {

    // use this grid size for tests
    std::array<int, 2> proc_grid_size = {3, 9};

    // get halo regions for 2d topology
    auto regions = DArrays::HaloRegions<2>();

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

    for (auto periodic1 : {true, false}) {
        for (auto periodic2 : {true, false}) {
            SECTION("is_periodic = xxx ") {
                std::array<int, 2> is_periodic = {periodic1, periodic2};
                DArrays::Topology::DArrayTopology<2> topo(MPI_COMM_WORLD, 
                                                          proc_grid_size,
                                                          is_periodic);

                // for every processors we want to check
                for (auto j : LinRange(procs.size())) {
                    // if it is actually the processor we want to check
                    if (topo._comm_rank == procs[j]) {
                        // for every neighbour
                        for (auto i : LinRange(regions.size())) {
                            // we only check if there is a neighbour, otherwise it should throw
                            if ( !topo.has_grid_boundary_at(regions[i]) ) {
                                REQUIRE( topo.neighbour_proc_rank(regions[i]) == expected[j][i] );
                            } else {
                                REQUIRE_THROWS( topo.neighbour_proc_rank(regions[i]) );
                            }
                        }
                    }
                }
            }
        }
    }
}

TEST_CASE("3D tests", "[3D-tests]") {

    // use this grid size for tests
    std::array<int, 3> proc_grid_size = {3, 3, 3};

    // get halo regions for 3d topology
    auto regions = DArrays::HaloRegions<3>();

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
                SECTION("is_periodic = xxx ") {
                    std::array<int, 3> is_periodic = {periodic1, periodic2, periodic3};
                    DArrays::Topology::DArrayTopology<3> topo(MPI_COMM_WORLD, 
                                                              proc_grid_size,
                                                              is_periodic);

                    // for every processors we want to check
                    for (auto j : LinRange(procs.size())) {
                        // if it is actually the processor we want to check
                        if (topo._comm_rank == procs[j]) {
                            // for every neighbour
                            for (auto i : LinRange(regions.size())) {
                                // we only check if there is a neighbour, otherwise it should throw
                                if ( !topo.has_grid_boundary_at(regions[i]) ) {
                                    REQUIRE( topo.neighbour_proc_rank(regions[i]) == expected[j][i] );
                                } else {
                                    REQUIRE_THROWS( topo.neighbour_proc_rank(regions[i]) );
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
