#include "DArrays.hpp"
#include <array>

int main (int argc, char* argv[]) {

    DArrays::MPISession::Initialize();

    std::array<int, 2> proc_grid_size = {    2,     2};
    std::array<int, 2>     array_size = {   20,    20};
    std::array<int, 2>     isperiodic = {false, false};
    std::array<int, 2>      nhalo_out = {    0,     0};
    int                      nhalo_in =             2;

    DArrays::DArraySpec<2> spec(MPI_COMM_WORLD, 
                                proc_grid_size, 
                                array_size, 
                                isperiodic, nhalo_out, nhalo_in);

    DArrays::DArray<double, 2> u(spec);

    DArrays::MPISession::Finalize();

}