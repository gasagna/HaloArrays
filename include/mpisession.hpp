#pragma once
#include <mpi.h>

namespace DArrays {
namespace MPISession {

inline void Initialize() {
    MPI_Init(nullptr, nullptr);
}

inline void Finalize() {
    MPI_Finalize();
}

}
}