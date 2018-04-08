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
//                 SEND/RECV SUBARRAYS                //
////////////////////////////////////////////////////////
template <typename T, size_t NDIMS>
void sendrecv(SubArray<T, NDIMS>& tosend, 
              SubArray<T, NDIMS>& torecv, 
              int dest_rank, int src_rank, int msgtag) {

    // create type_send
    MPI_Datatype type_send;
    MPI_Type_create_subarray(NDIMS,                             // parent number of dimensions
                             tosend.parent().raw_size().data(), // parent raw size
                             tosend.size().data(),              // subarray size
                             tosend.raw_origin().data(),        // raw coordinate of start
                             MPI_ORDER_FORTRAN,                 // we have column major data
                             MPI_DOUBLE, &type_send);
    MPI_Type_commit(&type_send);

    // create type_recv
    MPI_Datatype type_recv;
    MPI_Type_create_subarray(NDIMS,                             // parent number of dimensions
                             torecv.parent().raw_size().data(), // parent raw size
                             torecv.size().data(),              // subarray size
                             torecv.raw_origin().data(),        // raw coordinate of start
                             MPI_ORDER_FORTRAN,                 // we have column major data
                             MPI_DOUBLE, &type_recv);
    MPI_Type_commit(&type_recv);

    // actual call
    MPI_Sendrecv(tosend.parent().data(),
                 1,
                 type_send,
                 dest_rank,
                 msgtag,
                 torecv.parent().data(),
                 1,
                 type_recv,
                 src_rank,
                 msgtag,
                 torecv.parent().layout().communicator(), 
                 MPI_STATUS_IGNORE);

    MPI_Type_free(&type_send); MPI_Type_free(&type_recv);
}

}