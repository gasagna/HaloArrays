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

////////////////////////////////////////////////////////
//                 SEND/RECV SUBARRAYS                //
////////////////////////////////////////////////////////
template <typename T, size_t NDIMS>
void send(subarray<T, NDIMS> sub, int dest, MESSAGETAG tag) {
    MPI_Datatype type_send;
    MPI_Type_create_subarray(PARENTNDIMS,                // parent number of dimensions
                             sub.parent().size().data(), // parent size
                             sub.size().data(),          // subarray size
                             sub.raw_origin().data(),    // coordinate of start
                             MPI_ORDER_FORTRAN,          // we have column major data
                             T, &type_send);
    MPI_Type_commit(&type_send);
    MPI_Send(sub.parent().data(),
             1,
             type_send,
             dest,
             tag,
             sub.parent().communicator());
    MPI_Type_free(&type_send);
}
}