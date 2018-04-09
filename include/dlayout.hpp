#pragma once

// ===================================================================== //
// import iterator facilities
using namespace DArrays::Iterators;

namespace DArrays {

// ===================================================================== //
// defines topology of the distributed array
template <size_t NDIMS>
class DArrayLayout {
// VARIABLES
private:
    std::array<int, NDIMS> _is_periodic; // whether the processor grid should wrap around
    int                      _comm_size; // the number of processors in the communicator
    int                      _comm_rank; // rank of current processor within communicator
    std::array<int, NDIMS>      _coords; // coordinates of current processor in the grid
    std::array<int, NDIMS>        _size; // the size of the processor grid over which data is distributed
    MPI_Comm                      _comm; // communicator connecting all processor over which the array data is distributed

    // ===================================================================== //
    // helper
    inline bool _has_neighbour_at(Boundary tag, size_t dim) const {
        #if DARRAY_LAYOUT_CHECKBOUNDS
            _checkdims(dim, NDIMS);
        #endif

        if (_is_periodic[dim])      return true;
        if (tag == Boundary::LEFT)  return _coords[dim] != 0;
        if (tag == Boundary::RIGHT) return _coords[dim] != _size[dim] - 1;

        // CENTER and WILDCARD are neutral
        return true;
    }

public:
    // ===================================================================== //
    // constructor
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
    }

    // ===================================================================== //
    // get communicator
    inline const MPI_Comm& communicator() const {
        return _comm;
    }

    // ===================================================================== //
    // get processor grid size along dimension dim
    inline int size(size_t dim) const {
        #if DARRAY_LAYOUT_CHECKBOUNDS
            _checkdims(dim, NDIMS);
        #endif
        return _size[dim];        
    }

    // ===================================================================== //
    // get whether layout is periodic along dimension dim
    inline bool is_periodic(size_t dim) const {
        #if DARRAY_LAYOUT_CHECKBOUNDS
            _checkdims(dim, NDIMS);
        #endif
        return _is_periodic[dim];
    }

    // ===================================================================== //
    // get rank of neighbour process sharing a given halo region
    inline int rank_of_neighbour_at(const HaloRegionSpec<NDIMS>& halo) const {        
        // return MPI::PROC_NULL if on boundary
        if (!has_neighbour_at(halo))
            return MPI_PROC_NULL;

        // initialise to current coordinates, then modifies as needed
        std::array<int, NDIMS> target_coords = _coords;

        for (auto dim : LinRange(NDIMS)) {
            if (halo[dim] == Boundary::LEFT)  target_coords[dim] -= 1;
            if (halo[dim] == Boundary::RIGHT) target_coords[dim] += 1;
        }

        // call to the MPI function
        int target_proc_rank;
        int ret = MPI_Cart_rank(_comm,             
                    target_coords.data(),     
                    &target_proc_rank);

        return target_proc_rank;
    }
    
    // ===================================================================== //
    // number of processors
    inline int nprocs() const {
        return _comm_size;
    }

    // ===================================================================== //
    // get rank of current processor
    inline int rank() const {
        return _comm_rank;
    }

    // ===================================================================== //
    // whether this processor has a neighbour on given halo
    inline bool has_neighbour_at(const HaloRegionSpec<NDIMS>& halo) const {
        for (auto dim : LinRange(NDIMS))
            if (!_has_neighbour_at(halo[dim], dim))
                return false;
        return true;
    }
};

}