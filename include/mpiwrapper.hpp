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
