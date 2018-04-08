#pragma once
#include <mpi.h>
#include "subarray.hpp"

// import iterator facilities
using namespace DArrays::Iterators;

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

// I hate seeing 0s and 1s around. Hence, this function 
// constructs a unique message tag for each halo region.
template<size_t NDIMS>
int message_tag(Boundary<NDIMS> boundary) {
    int tag = 0;
    for ( auto dim : LinRange(NDIMS) )
        tag += std::pow(10, dim+1)*static_cast<int>(boundary[dim]);
    return tag;
}

template<size_t NDIMS>
int message_tag(BoundaryTag tag, size_t dim) {
    return 10*dim + static_cast<int>(tag);
}

////////////////////////////////////////////////////////
//                 SEND/RECV HALOREGION               //
////////////////////////////////////////////////////////
template <typename T, size_t NDIMS>
void sendrecv(const HaloRegion<T, NDIMS>& tosend, int dest_rank, 
              const HaloRegion<T, NDIMS>& torecv, int src_rank) {
    MPI_Sendrecv(tosend.parent().data(),
                 1,
                 tosend.type(),
                 dest_rank,
                 0,
                 torecv.parent().data(),
                 1,
                 torecv.type(),
                 src_rank,
                 0,
                 torecv.parent().layout().communicator(), 
                 MPI_STATUS_IGNORE);
}

}