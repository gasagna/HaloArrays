#define CATCH_CONFIG_RUNNER
#include <catch.hpp>
#include "DArrays.hpp"


int main(int argc, char* argv[]) {
    DArrays::MPI::Initialize();
    const int result = Catch::Session().run(argc, argv);
    DArrays::MPI::Finalize();
    return result;
}