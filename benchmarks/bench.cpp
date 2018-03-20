#include "DArrays.hpp"
#include <iostream>
#include <numeric>
#include <vector>
#include <chrono>
#include <cmath>
#include <mpi.h>

int main () {

    DArrays::MPISession::Initialize();
    

    DArrays::Topology::DArrayTopology<3> topo(MPI_COMM_WORLD,          // comm
                                              {1,     1,     1},       // grid size
                                              {false, false, false});  // is periodic

    DArrays::DArray<double, 3> A(topo, {100, 100, 100}, {1, 1, 1}, 1);
    DArrays::DArray<double, 3> B(topo, {100, 100, 100}, {1, 1, 1}, 1);

    // fill with garbage
    for (auto [i, j, k] : A.indices())
        A(i, j, k) = i*j*k;

    // compute sum of array pieces
    double sum = 0;
    auto t0 = std::chrono::steady_clock::now();    

    for (auto [i, j, k] : A.indices())
        B(i, j, k) = A(i, j, k) - B(i, j, k + 1) - B(i, j, k - 1) -
                                  B(i, j + 1, k) - B(i, j - 1, k) -
                                  B(i + 1, j, k) - B(i - 1, j, k);
    
    auto t1 = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();
    std::cout << "us\n";
    // trick the compilar not to elide the loop
    if (B[0])
        std::cout << "\n";

    DArrays::MPISession::Finalize();

    return 0;
}