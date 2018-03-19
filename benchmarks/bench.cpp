#include "DArrays.hpp"
#include <iostream>
#include <numeric>
#include <vector>
#include <chrono>
#include <mpi.h>

const int NREP = 10000;
const int N = 30;

int main () {

    DArrays::MPISession::Initialize();
    

    std::array<int, 3> proc_grid_size = {1, 1, 1}; // to be run on 4 processor
    std::array<int, 3> array_size = {N, N, N};
    std::array<int, 3> isperiodic = {false, false, false};
    std::array<int, 3> nhalo_out = {1, 1, 1};
    int nhalo_in = 1;

    DArrays::DArray<double, 3> A(MPI_COMM_WORLD,
                                 proc_grid_size,
                                 array_size,
                                 isperiodic, nhalo_out, nhalo_in);

    // DArrays::DArray<double, 3> B(MPI_COMM_WORLD,
                                //  proc_grid_size,
                                //  array_size,
                                //  isperiodic, nhalo_out, nhalo_in);

    // fill with garbage
    for (auto[i, j, k] : A.indices())
        A(i, j, k) = i*j*k;

    // compute sum of array pieces
    double sum = 0;
    auto t0 = std::chrono::steady_clock::now();
    
    for (int n = 0; n<NREP; n++) {
        for (auto [i, j, k] : A.indices())
            A(i, j, k) = A(i, j, k) - A(i, j, k + 1) - A(i, j, k - 1) -
                                      A(i, j + 1, k) - A(i, j - 1, k) -
                                      A(i + 1, j, k) - A(i - 1, j, k);
    }

    auto t1 = std::chrono::steady_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count() / (double)NREP;
    std::cout << "us\n";

    DArrays::MPISession::Finalize();

    return 0;
}