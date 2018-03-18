namespace DArrays {
namespace MPISession {

void Initialize() {
    MPI_Init(nullptr, nullptr);
}

void Finalize() {
    MPI_Finalize();
}

}
}