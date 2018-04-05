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
    std::map<Boundary<NDIMS>, int> _rank_of_neighbour_at_map; // ranks of neighbouring processors
    std::array<int, NDIMS>                      _is_periodic; // whether the processor grid should wrap around
    int                                           _comm_size; // the number of processors in the communicator
    int                                           _comm_rank; // rank of current processor within communicator
    std::array<int, NDIMS>                           _coords; // coordinates of current processor in the grid
    std::array<int, NDIMS>                             _size; // the size of the processor grid over which data is distributed
    MPI_Comm                                           _comm; // communicator connecting all processor over which the array data is distributed

// FUNCTIONS
private:
    // ===================================================================== //
    // CONSTRUCT THE COORDINATES OF THE PROCESSOR WE NEED TO TALK TO
    // BASED ON THE HALO BOUNDARY WE ARE CONSIDERING. THE CENTER TAG
    // CORRESPONDS TO A ZERO SHIFT, HENCE IT IS NOT INCLUDED. IT IS
    // AN ERROR TO CALL THIS FUNCTION FOR A BOUNDARY WHERE THERE IS
    // NO NEIGHBOUR. THIS IS USED AT INITIALISATION ONLY.
    int _rank_of_neighbour_at(Boundary<NDIMS> bnd) {
        // initialise to current coordinates, then modifies as needed
        std::array<int, NDIMS> target_coords = _coords;
        int target_proc_rank;
        for ( auto dim : LinRange(NDIMS) ) {
            if (bnd[dim] == BoundaryTag::LEFT)
                target_coords[dim] -= 1;
            if (bnd[dim] == BoundaryTag::RIGHT)
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
    inline bool _has_neighbour_at(BoundaryTag tag, size_t dim) {
        if (_is_periodic[dim]) return true;
        switch(tag) {
            case BoundaryTag::LEFT:
                return _coords[dim] != 0;
            case BoundaryTag::RIGHT:
                return _coords[dim] != _size[dim] - 1;
            case BoundaryTag::CENTER:
                return true;
        }
    }

    // ===================================================================== //
    // CHECK WHETHER WE ARE NOT GETTING OUT OF BOUNDS WITH THE DIMENSION
    void _checkboundsdim(size_t dim) const {
        if ( dim < 0 or dim >= NDIMS )
            throw std::out_of_range("dimension out of range");
    }

public:
    // ===================================================================== //
    // CONSTRUCTOR
    DArrayLayout(MPI_Comm comm, 
                 std::array<int, NDIMS> size,
                 std::array<int, NDIMS> is_periodic) 
        : _size        (size)
        , _is_periodic (is_periodic) {

        // get communicator size and my rank
        MPI_Comm_size(comm, &_comm_size);
        MPI_Comm_rank(comm, &_comm_rank);

        // check processor grid size matches with communicator size
        if (std::accumulate(size.begin(),
                            size.end(),
                            1, std::multiplies<int>()) != _comm_size)
            throw std::invalid_argument("incompatible processor count and processor grid size");

        // create communicator with cartesian topology
        MPI_Cart_create(comm, NDIMS, _size.data(), _is_periodic.data(), false, &_comm);

        // get cartesian coordinates of my rank
        MPI_Cart_coords(_comm, _comm_rank, NDIMS, _coords.data());

        // fill neighbours rank map
        for (auto bnd : Boundarys<NDIMS>()) {
            if (has_neighbour_at(bnd)) {
                _rank_of_neighbour_at_map[bnd] = _rank_of_neighbour_at(bnd);
            }
        }
    }

    // ===================================================================== //
    // GET PROCESSOR GRID SIZE ALONG DIMENSION DIM
    inline int size(size_t dim) const {
        #if DARRAY_CONFIG_CHECKBOUNDS
            _checkboundsdim(dim);
        #endif
        return _size[dim];        
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
    // GET RANK OF NEIGHBOUR PROCESS SHARING A GIVEN HALO bnd
    inline int rank_of_neighbour_at(Boundary<NDIMS> bnd) const {
        return _rank_of_neighbour_at_map.at(bnd);
    }

    // ===================================================================== //
    // GET RANK OF CURRENT PROCESSOR
    inline int rank() const {
        return _comm_rank;
    }

    // ===================================================================== //
    // WHETHER THIS PROCESSOR HAS A NEIGHBOUR SHARING A GIVEN HALO bnd
    inline bool has_neighbour_at(Boundary<NDIMS> bnd) {
        for ( auto dim : LinRange(NDIMS) ) {
            if ( !_has_neighbour_at(bnd[dim], dim) )
                return false;
        }
        return true;
    }
};

}