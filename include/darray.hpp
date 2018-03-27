#pragma once
#include <array>
#include <mpi.h>

namespace DArrays {

template <typename T, size_t NDIMS>
class DArray {
private:
    T*                                      _data; // actual data
    Topo::DArrayTopology<NDIMS>         _topology; // topologically-aware communicator object
    std::array<int, NDIMS>            _array_size; // global array size
    std::array<int, NDIMS>        _local_arr_size; // local array size
    std::array<int, NDIMS>               _nhalo_left; // number of halo points on 'left'  side (low index)
    std::array<int, NDIMS>               _nhalo_right; // number of halo points on 'right' side (high index)

    // ===================================================================== //
    // INDEXING INTO LINEAR MEMORY BUFFER
    template<typename... INDICES>
    inline size_t _tolinearindex(INDICES... indices) {
        return __tolinearindex(0, indices...);
    }

    inline size_t __tolinearindex(size_t dim, int i) {
        return i + _nhalo_left[dim];
    }

    template <typename... INDICES>
    inline size_t __tolinearindex(size_t dim, int i, INDICES... indices)  {
        return __tolinearindex(dim, i) + 
            _local_arr_size[dim]*__tolinearindex(dim+1, indices...);
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
        if !( is_inbounds(i, dim) )
            throw std::out_of_range("Out of range");
        __checkbounds(dim+1, indices...);
    }

public:

    // ===================================================================== //
    // CONTAINER INTERFACE
    using value_type     = T;
    using DArrayIterator = T*;

    // ===================================================================== //
    // CONSTRUCTOR/DESTRUCTOR
    DArray() = delete;

    DArray(Topo::DArrayTopology<NDIMS> topo, 
           std::array<int, NDIMS> array_size,
           std::array<int, NDIMS> nhalo_out, int nhalo_in)
        : _array_size  (array_size ) 
        , _topology        (topo       ) {
            // define size of local array and number of left/right halo points
            for (auto dim : LinearRange(0, NDIMS))
                _local_arr_size[dim] = _array_size[dim] / topo.grid_size_along_dim(dim);
                _nhalo_left[dim] = topo.is_on_left_boundary(dim) || topo.is_periodic_along_dim(dim) ? 
                                    nhalo_out[dim] : nhalo_in;
                _nhalo_right[dim] = topo.is_on_right_boundary(dim) || topo.is_periodic_along_dim(dim) ? 
                                    nhalo_out[dim] : nhalo_in;
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
        return IndexIterator<NDIMS>(_local_arr_size);
    }

    // ===================================================================== //
    // check whether index i in in bounds along dimension dim
    inline bool is_inbounds(int i, size_t dim) {
        return (i < - _nhalo_left[dim] or i > _local_arr_size[dim] + _nhalo_right[dim] - 1);
    }

    // ===================================================================== //
    // UTILITIES 
    // ~~~ size of memory buffer, including halo ~~~
    inline size_t nelements() const {
        size_t _nelements = 1;
        for (auto dim : LinearRange(0, NDIMS))
            _nelements *= _nhalo_left[dim] + _local_arr_size[dim] + _nhalo_right[dim];
        return _nelements;
    }

    // ~~~ local array size and halo dimension ~~~
    inline std::array<int, NDIMS> size() { 
        return _local_arr_size; 
    }

    inline int size(size_t dim) { 
        return _local_arr_size[dim]; 
    }

    inline int nhalo(HaloTag tag, size_t dim) { 
        switch (tag) {
            case LEFT:   return _nhalo_left[dim];
            case RIGHT:  return _nhalo_right[dim];
            case CENTER: return _local_arr_size[dim];
        }
    }

    // ===================================================================== //
    // UPDATE HALO POINTS - specialisations for NDIMS = 1, 2, 3 elsewhere
    void halo_swap();
};

}