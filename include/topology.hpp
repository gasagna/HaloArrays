#pragma once
#include <array>
#include <numeric>
#include <functional>
#include <mpi.h>

// import iterator facilities
using namespace DArrays::Iterators;

namespace DArrays::Topology {

// ===================================================================== //
// DEFINES TOPOLOGY OF THE DISTRIBUTED ARRAY
template <size_t NDIMS>
class DArrayTopology {
private:
    // communicator connecting all processor over which the array data is distributed
    MPI_Comm                            _comm;
    // the number of processors in the communicator
    int                            _comm_size;
    // rank of current processor within communicator
public:
    int                            _comm_rank;
private:
    // the size of the processor grid over which data is distributed
    std::array<int, NDIMS>    _proc_grid_size;
    // coordinates of current processor in the grid
    std::array<int, NDIMS>  _proc_grid_coords;
    // whether the processor grid should wrap around
    std::array<int, NDIMS>       _is_periodic;

public:
    DArrayTopology(MPI_Comm comm, 
                   std::array<int, NDIMS> grid_size,
                   std::array<int, NDIMS> is_periodic) 
        : _proc_grid_size   (grid_size )
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
        MPI_Cart_create(comm, NDIMS, _proc_grid_size.data(),
                        _is_periodic.data(), false, &_comm);

        // get cartesian coordinates of my rank
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
        MPI_Cart_coords(_comm, _comm_rank, NDIMS, _proc_grid_coords.data());
    }

    // ~~~ handle to communicator ~~~
    MPI_Comm& communicator () {
        return _comm;
    }

    // ===================================================================== //
    // QUERY WHETHER THE CURRENT PROCESSOR IS ON THE BOUNDARY OF THE PROCESSOR
    // GRID. IF TRUE, THERE IS NO NEIGHBOURING PROCESSOR WE NEED TO SEND/RECV
    // DATA TO/FROM.
    inline bool has_grid_boundary_at(HaloRegion<NDIMS> region) {
        for ( auto dim : LinRange(NDIMS) ) {
            if ( _has_grid_boundary_at(region[dim], dim) )
                return true;
        }
        return false;
    }

    // helper function
    inline bool _has_grid_boundary_at(HaloRegionTag tag, size_t dim) {
        if (_is_periodic[dim]) return false;
        switch(tag) {
            case HaloRegionTag::LEFT:
                return _proc_grid_coords[dim] == 0;
            case HaloRegionTag::RIGHT:
                return _proc_grid_coords[dim] == _proc_grid_size[dim] - 1;
            case HaloRegionTag::CENTER:
                return false;
        }
    }
    
    // ===================================================================== //
    // GET RANK OF NEIGHBOURING PROCESSOR
    // TODO: rather then recomputing this over and over, one could have a 
    // hashmap that stores this information, avoiding recomputations
    int neighbour_proc_rank(HaloRegion<NDIMS> region) {

        // construct the coordinates of the processor we need to talk to
        // based on the halo region we are considering. The CENTER tag
        // corresponds to a zero shift, hence it is not included.         
        std::array<int, NDIMS> target_coords = _proc_grid_coords;
        int target_proc_rank;
        for ( auto dim : LinRange(NDIMS) ) {
            if (region[dim] == HaloRegionTag::LEFT)
                target_coords[dim] -= 1;
            if (region[dim] == HaloRegionTag::RIGHT)
                target_coords[dim] += 1;
        }

        // call to the MPI function
        MPI_Cart_rank(_comm,             
                      target_coords.data(),     
                      &target_proc_rank);

        return target_proc_rank;
    }
};

}