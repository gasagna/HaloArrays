#pragma once
#include <map>
#include <array>
#include <numeric>
#include <iostream>
#include <functional>
#include <mpi.h>

// import iterator facilities
using namespace DArrays::Iterators;

namespace DArrays {

// ===================================================================== //
// DEFINES TOPOLOGY OF THE DISTRIBUTED ARRAY
template <size_t NDIMS>
class DArrayLayout {
// VARIABLES
private:
    std::map<HaloRegion<NDIMS>, int> _rank_of_neighbour_at_map; // ranks of neighbouring processors
    std::array<int, NDIMS>                   _proc_grid_coords; // coordinates of current processor in the grid
    std::array<int, NDIMS>                     _proc_grid_size; // the size of the processor grid over which data is distributed
    std::array<int, NDIMS>                        _is_periodic; // whether the processor grid should wrap around
    int                                             _comm_size; // the number of processors in the communicator
    int                                             _comm_rank; // rank of current processor within communicator
    MPI_Comm                                             _comm; // communicator connecting all processor over which the array data is distributed

// FUNCTIONS
private:
    // ===================================================================== //
    // CONSTRUCT THE COORDINATES OF THE PROCESSOR WE NEED TO TALK TO
    // BASED ON THE HALO REGION WE ARE CONSIDERING. THE CENTER TAG
    // CORRESPONDS TO A ZERO SHIFT, HENCE IT IS NOT INCLUDED. IT IS
    // AN ERROR TO CALL THIS FUNCTION FOR A REGION WHERE THERE IS
    // NO NEIGHBOUR. THIS IS USED AT INITIALISATION ONLY.
    int _rank_of_neighbour_at(HaloRegion<NDIMS> region) {
        std::array<int, NDIMS> target_coords = _proc_grid_coords;
        int target_proc_rank;
        for ( auto dim : LinRange(NDIMS) ) {
            if (region[dim] == HaloRegionTag::LEFT)
                target_coords[dim] -= 1;
            if (region[dim] == HaloRegionTag::RIGHT)
                target_coords[dim] += 1;
        }

        // call to the MPI function
        int ret = MPI_Cart_rank(_comm,             
                      target_coords.data(),     
                      &target_proc_rank);

        return target_proc_rank;
    }

    // ===================================================================== //
    // HELPER FUNCTION
    inline bool _has_neighbour_at(HaloRegionTag tag, size_t dim) {
        if (_is_periodic[dim]) return true;
        switch(tag) {
            case HaloRegionTag::LEFT:
                return _proc_grid_coords[dim] != 0;
            case HaloRegionTag::RIGHT:
                return _proc_grid_coords[dim] != _proc_grid_size[dim] - 1;
            case HaloRegionTag::CENTER:
                return true;
        }
    }

    // ===================================================================== //
    // CHECK WHETHER WE ARE NOT GETTING OUT OF BOUNDS WITH THE DIMENSION
    void _checkboundsdim(size_t dim) {
        if ( dim < 0 or dims >= NDIMS )
            throw std::out_of_range("dimension out of range");
    }

public:
    // ===================================================================== //
    // CONSTRUCTOR
    DArrayLayout(MPI_Comm comm, 
                   std::array<int, NDIMS> proc_grid_size,
                   std::array<int, NDIMS> is_periodic) 
        : _proc_grid_size   (proc_grid_size  )
        , _is_periodic      (is_periodic) {

        // get communicator size and my rank
        MPI_Comm_size(comm, &_comm_size);
        MPI_Comm_rank(comm, &_comm_rank);

        // check processor grid size matches with communicator size
        if (std::accumulate(proc_grid_size.begin(),
                            proc_grid_size.end(),
                            1, std::multiplies<int>()) != _comm_size)
            throw std::invalid_argument("incompatible processor count and processor grid size");

        // create communicator with cartesian topology
        MPI_Cart_create(comm, NDIMS, _proc_grid_size.data(),
                        _is_periodic.data(), false, &_comm);

        // get cartesian coordinates of my rank
        MPI_Cart_coords(_comm, _comm_rank, NDIMS, _proc_grid_coords.data());

        // fill neighbours rank map
        for (auto region : HaloRegions<NDIMS>()) {
            if (has_neighbour_at(region)) {
                _rank_of_neighbour_at_map[region] = _rank_of_neighbour_at(region);
            }
        }
    }

    // ===================================================================== //
    // GET PROCESSOR GRID SIZE ALONG DIMENSION DIM
    inline int grid_size(size_t dim) const {
        #if DARRAY_CONFIG_CHECKBOUNDS
            _checkboundsdim(dim);
        #endif
        return _grid_size[dim];        
    }

    // ===================================================================== //
    // GET WHETHER ARRAY IS PERIODIC ALONG DIMENSION DIM
    inline bool is_periodic(size_t dim) const {
        #if DARRAY_CONFIG_CHECKBOUNDS
            _checkboundsdim(dim);
        #endif
        return _is_periodic[dim];
    }

    // ===================================================================== //
    // GET RANK OF NEIGHBOUR PROCESS SHARING A GIVEN HALO REGION
    inline int rank_of_neighbour_at(HaloRegion<NDIMS> region) const {
        return _rank_of_neighbour_at_map.at(region);
    }

    // ===================================================================== //
    // GET RANK OF CURRENT PROCESSOR
    inline int rank() const {
        return _comm_rank;
    }

    // ===================================================================== //
    // WHETHER THIS PROCESSOR HAS A NEIGHBOUR SHARING A GIVEN HALO REGION
    inline bool has_neighbour_at(HaloRegion<NDIMS> region) {
        for ( auto dim : LinRange(NDIMS) ) {
            if ( !_has_neighbour_at(region[dim], dim) )
                return false;
        }
        return true;
    }
};

}