#pragma once
#include <array>
#include <numeric>
#include <functional>
#include <mpi.h>

#ifndef DARRAY_CONFIG_CHECKBOUNDS
#define DARRAY_CONFIG_CHECKBOUNDS false
#endif

namespace DArrays {

// ===================================================================== //
// HELPER FUNCTIONS
namespace Utils {

template <size_t NDIMS>
bool _is_on_leftboundary(std::array<int, NDIMS> proc_grid_coords, std::array<int, NDIMS> proc_grid_size, int dim) {
    return proc_grid_coords[dim] == 0;
}

template <size_t NDIMS>
bool _is_on_right_boundary(std::array<int, NDIMS> proc_grid_coords, std::array<int, NDIMS> proc_grid_size, int dim) {
    return proc_grid_coords[dim] == proc_grid_size[dim] - 1;
}
}

template <typename T, size_t NDIMS>
class DArray {
private:
    T*                                 _data; // actual data
    MPI_Comm                           _comm; // communicator
    int                           _comm_rank; // rank within communicator
    int                           _comm_size; // size of communicator
    std::array<int, NDIMS>   _proc_grid_size; // processor grid size
    std::array<int, NDIMS>  _glob_array_size; // global array size
    std::array<int, NDIMS>       _isperiodic; // global array size
    std::array<int, NDIMS>       _array_size; // local array size
    std::array<int, NDIMS>          _nhalo_l; // number of halo points on 'left'  side (low index)
    std::array<int, NDIMS>          _nhalo_r; // number of halo points on 'right' side (high index)

    // ===================================================================== //
    // INDEXING INTO LINEAR MEMORY BUFFER
    template<typename... INDICES>
    inline size_t _tolinearindex(INDICES... indices) {
        return __tolinearindex(0, indices...);
    }

    inline size_t __tolinearindex(size_t dim, int i) {
        return i + _nhalo_l[dim];
    }

    template <typename... INDICES>
    inline size_t __tolinearindex(size_t dim, int i, INDICES... indices)  {
        return __tolinearindex(dim, i) + 
            _array_size[dim]*__tolinearindex(dim+1, indices...);
    }

    // ===================================================================== //
    // BOUND CHECKING
    template <typename... INDICES>
    inline void _checkbounds(INDICES... indices) {
        return __checkbounds(0, indices...);
    }

    inline void __checkbounds(size_t dim) {}

    template<typename... INDICES>
    inline void __checkbounds(size_t dim, int i, INDICES... indices) {
        if (i < -_nhalo_l[dim] or i > _array_size[dim] + _nhalo_r[dim] - 1) {
            throw std::out_of_range("Out of range");
        }
        __checkbounds(dim+1, indices...);
    }


public:
    std::array<int, NDIMS> _proc_grid_coords; // cartesian coordinates of this process (ROW MAJOR)

    // ===================================================================== //
    // CONTAINER INTERFACE
    using value_type     = T;
    using DArrayIterator = T*;

    // ===================================================================== //
    // CONSTRUCTOR/DESTRUCTOR
    DArray() = delete;

    DArray(MPI_Comm comm, std::array<int, NDIMS> proc_grid_size,
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
                // depending on location we change the number of halo points
                _nhalo_l[dim] = Utils::_is_on_leftboundary(_proc_grid_coords, _proc_grid_size, dim) ? 
                                        nhalo_out[dim] : nhalo_in;
                _nhalo_r[dim] = Utils::_is_on_right_boundary(_proc_grid_coords, _proc_grid_size, dim) ? 
                                        nhalo_out[dim] : nhalo_in;
                // if periodic we use n_halo_out anyway                                        
                _nhalo_l[dim] = isperiodic[dim] ? nhalo_out[dim] : _nhalo_l[dim];
                _nhalo_r[dim] = isperiodic[dim] ? nhalo_out[dim] : _nhalo_r[dim];
            }

            // allocate memory buffer
            _data = new T[nelements()];
    }

    ~DArray() {
        delete[] _data;
    }

    // ===================================================================== //
    // INDEXING INTO LINEAR MEMORY BUFFER
    // ~~~ linear indexing ~~~
    const inline T& operator [] (size_t i) const { return _data[i]; }
          inline T& operator [] (size_t i)       { return _data[i]; }

    // ~~~ local cartesian indexing - first index is contiguous ~~~
    template <typename... INDICES>
    inline T& operator () (INDICES... indices) {
        static_assert(sizeof...(INDICES) == NDIMS,
                      "Number of indices must match array dimension");
        #if DARRAY_CONFIG_CHECKBOUNDS
            _checkbounds(indices...);
        #endif
        return _data[_tolinearindex(indices...)];
    }

    // ===================================================================== //
    // ITERATION OVER ALL DATA
    DArrayIterator begin() { return _data; }
    DArrayIterator end()   { return _data + nelements(); }

    // ===================================================================== //
    // ITERATOR OVER THE IN-DOMAIN INDICES 
    inline IndexIterator<NDIMS> indices () {
        return IndexIterator<NDIMS>(_array_size);
    }

    // ===================================================================== //
    // UTILITIES 
    // ~~~ size of memory buffer, including halo ~~~
    inline size_t nelements() const {
        size_t _nelements = 1;
        for (int dim = 0; dim != NDIMS; dim++)
            _nelements *= _nhalo_l[dim] + _array_size[dim] + _nhalo_r[dim];
        return _nelements;
    }

    // ~~~ local array size ~~~
    inline std::array<int, NDIMS> size() {
        return _array_size;
    }

    // ===================================================================== //
    // QUERIES FOR POSITION ON THE PROCESSOR GRID
    inline bool is_on_left_boundary(size_t dim) {
        // FIXME: IN CASE WE HAVE A PERIODIC DOMAIN
        return Utils::_is_on_leftboundary(_proc_grid_coords, _proc_grid_size, dim);
    }

    inline bool is_on_right_boundary(size_t dim) {
        // FIXME: IN CASE WE HAVE A PERIODIC DOMAIN
        return Utils::_is_on_right_boundary(_proc_grid_coords, _proc_grid_size, dim);
    }

    inline bool is_on_boundary() {
        // FIXME: IN CASE WE HAVE A PERIODIC DOMAIN
        for (auto dim = 0; dim != NDIMS; dim++)
            if (is_on_left_boundary(dim) || is_on_right_boundary(dim))
                return true;
        return false;
    }

    // ===================================================================== //
    // UPDATE HALO POINTS - specialisations for NDIMS = 1, 2, 3 elsewhere
    void halo_swap();
};

}