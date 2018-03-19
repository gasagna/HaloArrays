#pragma once
#include <array>
#include <numeric>
#include <functional>
#include <mpi.h>

namespace DArrays {
    
template <size_t NDIMS>
bool _isboundary_l(std::array<int, NDIMS> coords, std::array<int, NDIMS> size, int dim) {
    return coords[dim] == 0;
}

template <size_t NDIMS>
bool _isboundary_r(std::array<int, NDIMS> coords, std::array<int, NDIMS> size, int dim) {
    return coords[dim] == size[dim] - 1;
}

template <size_t NDIMS>
class DArraySpec {
private:
    MPI_Comm                           _comm; // communicator
    int                           _comm_rank; // rank within communicator
    int                           _comm_size; // size of communicator
    std::array<int, NDIMS> _proc_grid_coords; // cartesian coordinates of this process
    std::array<int, NDIMS>   _proc_grid_size; // processor grid size
    std::array<int, NDIMS>  _glob_array_size; // global array size
    std::array<int, NDIMS>       _isperiodic; // global array size
    std::array<int, NDIMS>       _array_size; // local array size
    std::array<int, NDIMS>          _nhalo_l; // number of padding points on 'left'  side (low index)
    std::array<int, NDIMS>          _nhalo_r; // number of padding points on 'right' side (high index)

    // ===================================================================== //
    // CONSTRUCTOR
    DArraySpec(MPI_Comm comm, std::array<int, NDIMS> proc_grid_size,
                              std::array<int, NDIMS> array_size,
                              std::array<int, NDIMS> isperiodic,
                              std::array<int, NDIMS> nhalo_out, int nhalo_in)
        : _proc_grid_size  (proc_grid_size)
        , _glob_array_size (array_size    ) 
        , _isperiodic      (isperiodic    ) {
            // get rank of process and overall size of communicator
            MPI_Comm_size(comm, &_comm_size);
            MPI_Comm_rank(comm, &_comm_rank);

            // check processor grid size matches with communicator size
            if (std::accumulate(proc_grid_size.begin(), 
                                proc_grid_size.end(), 
                                1, std::multiplies<int>()) != _comm_size)
                throw std::invalid_argument("incompatible processor count and processor grid size");
            
            // create communicator with cartesian topology
            MPI_Cart_create(comm, NDIMS, _proc_grid_size.data(), 
                _isperiodic.data(), false, &_comm);

            // get local coordinates within the grid
            MPI_Cart_coords(_comm, _comm_rank, NDIMS, _proc_grid_coords.data());

            // define size and number of halo points of local array, depending on local position
            for (int dim = 0; dim != NDIMS; dim++) {
                _array_size[dim] = _glob_array_size[dim] / _proc_grid_size[dim];
                _nhalo_l[dim] = isperiodic[dim] or _isboundary_l(_proc_grid_coords, _proc_grid_size, dim) ? 
                                    nhalo_in : nhalo_out[dim];
                _nhalo_r[dim] = isperiodic[dim] or _isboundary_l(_proc_grid_coords, _proc_grid_size, dim) ? 
                                    nhalo_in : nhalo_out[dim];
            }
    }

    // ===================================================================== //
    // BOUND CHECKING
    template <typename... INDICES>
    inline void checkbounds(INDICES... indices) {
        return _checkbounds(0, indices...);
    }

    inline void _checkbounds(size_t dim) {}

    template<typename... INDICES>
    inline void _checkbounds(size_t dim, int i, INDICES... indices) {
        if (i < -_nhalo_l[dim] or i > _array_size[dim] + _nhalo_r[dim] - 1) {
            throw std::out_of_range("Out of range");
        }
        _checkbounds(dim+1, indices...);
    }

    // ===================================================================== //
    // INDEXING INTO LINEAR MEMORY BUFFER
    template<typename... INDICES>
    inline size_t tolinearindex(INDICES... indices) {
        return _tolinearindex(0, indices...);
    }

    inline size_t _tolinearindex(size_t dim, int i) {
        return i + _nhalo_l[dim];
    }

    template <typename... INDICES>
    inline size_t _tolinearindex(size_t dim, int i, INDICES... indices)  {
        return _tolinearindex(dim, i) + 
            _array_size[dim]*_tolinearindex(dim+1, indices...);
    }

    // ===================================================================== //
    // UTILITIES 
    // ~~ size of memory buffer, including halo ~~
    inline size_t nelements() const {
        size_t _nelements = 1;
        for (int dim = 0; dim != NDIMS; dim++)
            _nelements *= _nhalo_l[dim] + _array_size[dim] + _nhalo_r[dim];
        return _nelements;
    }

    inline std::array<int, NDIMS> size() {
        return _array_size;
    }
};

}