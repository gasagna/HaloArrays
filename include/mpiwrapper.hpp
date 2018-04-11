#pragma once
// #include "subarray.hpp"

namespace DArrays::MPI {

// ===================================================================== //
// initialize/finalize mpi session            
inline void Initialize() {
    MPI_Init(nullptr, nullptr);
}

inline void Finalize() {
    MPI_Finalize();
}

// ===================================================================== //
// send/recv haloregion               
template <typename T, size_t NDIMS>
void sendrecv(SubArray<T, NDIMS>& tosend, int dest_rank, 
              SubArray<T, NDIMS>& torecv, int src_rank) {
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