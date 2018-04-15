#include "DArrays.hpp"
#include <catch.hpp>
#include <array>
#include <iostream>

// import all
using namespace DArrays;

TEST_CASE("memoryorder - 2D", "test_1]") {

    // use this grid layout for tests
    std::array<int, 2> layout_size = {3, 9};
    std::array<int, 2> is_periodic = {true, true};

    // create layout
    DArrayLayout<2> layout(MPI_COMM_WORLD, layout_size, is_periodic);

    // create array 
    std::array<int, 2> array_size = {27*5};
    std::array<int, 2> nhalo_out  = {0, 0};
    std::array<int, 2> nhalo_in   = {0, 0};
    FortranMemoryOrder<2> m;
    CMemoryOrder<2> m;
    MemoryOrder<0, 1, 2> m;
    DArray<double, 2, m> a(layout, array_size, nhalo_out, nhalo_in); 
}