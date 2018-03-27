#pragma once
#include <mpi.h>

namespace DArrays::MPI {

////////////////////////////////////////////////////////
//         INITIALIZE/FINALIZE MPI SESSION            //
////////////////////////////////////////////////////////
inline void Initialize() {
    MPI_Init(nullptr, nullptr);
}

inline void Finalize() {
    MPI_Finalize();
}


////////////////////////////////////////////////////////
//                   MESSAGE TAGS                     //
////////////////////////////////////////////////////////

// I hate seeing 0s and 1s around. This function constructs
// a unique message tag for each halo region.
template<size_t NDIMS>
int message_tag(HaloRegion<NDIMS> region) {
    int tag = 0;
    for ( auto dim : LinRange(NDIMS) )
        tag += std::pow(10, dim+1)*region[dim];
    return tag;
}
