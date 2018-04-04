#pragma once
#include <array>
#include <mpi.h>

namespace DArrays {

template <typename T, size_t NDIMS>
class DArray {
private:
// VARIABLES
    std::array<int, NDIMS> _local_arr_size; // local array size
    std::array<int, NDIMS>   _raw_arr_size; // local array size, including halo points
    std::array<int, NDIMS>    _nhalo_right; // number of halo points on 'right' side (high index)
    std::array<int, NDIMS>     _array_size; // global array size
    std::array<int, NDIMS>     _nhalo_left; // number of halo points on 'left'  side (low index)
    T*                               _data; // actual data

// FUNCTIONS
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
        if !( _is_inbounds(i, dim) )
            throw std::out_of_range("Out of range");
        __checkbounds(dim+1, indices...);
    }

    // ===================================================================== //
    // CHECK WHETHER INDEX I IN IN BOUNDS ALONG DIMENSION DIM, INCLUDING HALO
    inline bool _is_inbounds(int i, size_t dim) {
        return (i < - _nhalo_left[dim] or i > _local_arr_size[dim] + _nhalo_right[dim] - 1);
    }

public:
// VARIABLES
    DArrayLayout<NDIMS>            layout; // topologically-aware communicator object

// FUNCTIONS
    // ===================================================================== //
    // CONTAINER INTERFACE
    using value_type     = T;
    using DArrayIterator = T*;

    // ===================================================================== //
    // CONSTRUCTOR/DESTRUCTOR
    DArray() = delete;

    DArray(DArrayLayout<NDIMS> layout, 
           std::array<int, NDIMS> array_size,
           std::array<int, NDIMS> nhalo_out, int nhalo_in)
        : _array_size (array_size ) 
        , layout      (layout     ) {
            // define size of local array and number of left/right halo points
            for (auto dim : LinRange(NDIMS)) {
                _local_arr_size[dim] = _array_size[dim] / layout.grid_size(dim);
                _nhalo_left[dim]  = layout.has_neighbour_at(HaloRegionTag::LEFT, dim)  || layout.is_periodic(dim) ? 
                                    nhalo_out[dim] : nhalo_in;
                _nhalo_right[dim] = layout.has_neighbour_at(HaloRegionTag::RIGHT, dim) || layout.is_periodic(dim) ? 
                                    nhalo_out[dim] : nhalo_in;

                // full size of the data
                _raw_arr_size[dim] = _local_arr_size[dim] + _nhalo_left[dim] + _nhalo_right[dim];
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
    // UTILITIES 
    // ~~~ size of memory buffer, including halo ~~~
    inline size_t nelements() const {
        return std::reduce(_raw_arr_size.begin(), _raw_arr_size.end(), 
                           1, std::multiplies<>());
    }

    // ~~~ local array size and halo dimension ~~~
    inline std::array<int, NDIMS>& size() const { 
        return _local_arr_size; 
    }

    inline std::array<int, NDIMS>& raw_size() const { 
        return _raw_arr_size; 
    }

    inline int size(size_t dim) { 
        return _local_arr_size[dim]; 
    }

    inline int nhalo(HaloRegionTag tag, size_t dim) { 
        switch (tag) {
            case HaloRegionTag::LEFT   : return _nhalo_left[dim];
            case HaloRegionTag::RIGHT  : return _nhalo_right[dim];
            case HaloRegionTag::CENTER : return _local_arr_size[dim];
        }
    }

};

}