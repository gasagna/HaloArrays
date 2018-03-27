#pragma once
#include <array>
#include <numeric>
#include <functional>
#include <mpi.h>

namespace DArrays::Topo {



// ===================================================================== //
// DEFINES TOPOLOGY OF THE DISTRIBUTED ARRAY
template <size_t NDIMS>
class DArrayTopology {
private:
    MPI_Comm                       _comm; // communicator
    int                       _comm_size; // size of communicator
    int                       _comm_rank; // my rank within communicator
    std::array<int, NDIMS>    _grid_size; // processor grid size
    std::array<int, NDIMS>  _grid_coords; // processor grid coordinates
    std::array<int, NDIMS>  _is_periodic; // global array size

public:
    DArrayTopology(MPI_Comm comm, 
                   std::array<int, NDIMS> grid_size,
                   std::array<int, NDIMS> is_periodic) 
        : _grid_size   (grid_size )
        , _is_periodic (is_periodic) {

        // get communicator size and my rank
        MPI_Comm_size(comm, &_comm_size);
        MPI_Comm_rank(comm, &_comm_rank);

        // check processor grid size matches with communicator size
        if (std::accumulate(grid_size.begin(),
                            grid_size.end(),
                            1, std::multiplies<int>()) != _comm_size)
            throw std::invalid_argument("incompatible processor count and processor grid size");

        // create communicator with cartesian topology
        MPI_Cart_create(comm, NDIMS, _grid_size.data(),
                        _is_periodic.data(), false, &_comm);

        // get cartesian coordinates of my rank
        MPI_Cart_coords(_comm, _comm_rank, NDIMS, _grid_coords.data());
    }

    // ===================================================================== //
    // GETTERS
    // ~~~ get processor grid size along dimension dim ~~~
    inline int grid_size_along_dim(size_t dim) const {
        return _grid_size[dim];        
    }

    // ~~~ get whether array is periodic along dimension dim ~~~
    inline bool is_periodic_along_dim(size_t dim) const {
        return _is_periodic[dim];
    }

    // ~~~ handle to communicator ~~~
    MPI_Comm& communicator () {
        return _comm;
    }

    // ===================================================================== //
    // QUERIES FOR POSITION ON THE PROCESSOR GRID
    inline bool is_on_boundary(BoundaryTag tag, size_t dim) {
        if (_is_periodic[dim]) return false;
        switch(tag) {
            case BoundaryTag::LEFT:
                return _grid_coords[dim] == 0;
            case BoundaryTag::RIGHT:
                return _grid_coords[dim] == _grid_size[dim] - 1;
        }
    }

    inline bool is_on_boundary() {
        for (auto dim : LinRange(NDIMS))
            if (is_on(Boundary::RIGHT, dim) || is_on(Boundary::LEFT, dim))
                return true;
        return false;
    }

    inline bool has_halo_at(BoundaryRegion<NDIMS> region) {
        for ( auto dim : LinRange(NDIMS) ) {
            if ( is_on_boundary(region[dim], dim) )
                return true;
        }
        return false;
    }
};

}