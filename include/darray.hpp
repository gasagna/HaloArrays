#pragma once
#include <array>
#include <mpi.h>

namespace DArrays {

using namespace DArrays::Iterators;

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
        if ( !_is_inbounds(i, dim) )
            throw std::out_of_range("Out of range");
        __checkbounds(dim+1, indices...);
    }

    // ===================================================================== //
    // CHECK WHETHER INDEX I IN IN BOUNDS ALONG DIMENSION DIM, INCLUDING HALO
    inline bool _is_inbounds(int i, size_t dim) {
        return (i >= - _nhalo_left[dim] and i <= _local_arr_size[dim] + _nhalo_right[dim] - 1);
    }

public:
// VARIABLES
    DArrayLayout<NDIMS> layout; // topologically-aware communicator object

// FUNCTIONS
    // ===================================================================== //
    // CONTAINER INTERFACE
    using value_type     = T;

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
                _local_arr_size[dim] = _array_size[dim] / layout.size(dim);
                _nhalo_left[dim]  = layout.has_neighbour_at(BoundaryTag::LEFT, dim)  ? nhalo_in : nhalo_out[dim];
                _nhalo_right[dim] = layout.has_neighbour_at(BoundaryTag::RIGHT, dim) ? nhalo_in : nhalo_out[dim];

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
        #if DARRAY_ARRAY_CHECKBOUNDS
            _checkbounds(indices...);
        #endif
        return _data[_tolinearindex(indices...)];
    }

    // ===================================================================== //
    // ITERATOR OVER THE IN-DOMAIN INDICES 
    inline IndexRange<NDIMS> indices () {
        return IndexRange<NDIMS>(_local_arr_size);
    }

    // ===================================================================== //
    // LOCAL ARRAY SIZE
    inline const std::array<int, NDIMS>& size() const { 
        return _local_arr_size; 
    }

    // ===================================================================== //
    // NUMBER OF HALO POINTS AT A PARTICULAR BOUNDARY
    inline int nhalo(BoundaryTag tag, size_t dim) { 
        switch (tag) {
            case BoundaryTag::LEFT   : return _nhalo_left[dim];
            case BoundaryTag::RIGHT  : return _nhalo_right[dim];
            case BoundaryTag::CENTER : return _local_arr_size[dim];
        }
    }

    // ===================================================================== //
    // LOCAL ARRAY SIZE, INCLUDING HALO ELEMENTS
    inline const std::array<int, NDIMS>& raw_size() const { 
        return _raw_arr_size; 
    }
    
    // ===================================================================== //
    // SIZE OF MEMORY BUFFER, INCLUDING HALO
    inline size_t nelements() const {
        return std::reduce(_raw_arr_size.begin(), _raw_arr_size.end(), 
                           1, std::multiplies<>());
    }

    // ===================================================================== //
    // UPDATE HALO POINTS
    void swap_halo() {
        // loop over the boundary regions and send/recv 
        for ( auto boundary : AllBoundaries<NDIMS>() ) {
            if ( layout.has_neighbour_at(boundary) ) {
                sendrecv(subarray(*this, boundary, BoundaryIntent::SEND), 
                         subarray(*this, boundary, BoundaryIntent::RECV),
                         layout.rank(), layout.rank_of_neighbour_at(boundary),
                         message_tag(boundary));
            }
        }
    }       
};
}