#pragma once
#include <array>
#include "subarray.hpp"
#include "mpiwrapper.hpp"

namespace DArrays {

using namespace DArrays::Iterators;
using namespace DArrays::MPI;

// forward declaration
template <typename T, size_t NDIMS> class HaloRegion;

////////////////////////////////////////////////////////
//                       DArray                       //
////////////////////////////////////////////////////////
template <typename T, size_t NDIMS>
class DArray {
private:
// VARIABLES
    std::array<int, NDIMS>        _local_arr_size; // local array size
    std::array<int, NDIMS>          _raw_arr_size; // local array size, including halo points
    std::array<int, NDIMS>           _nhalo_right; // number of halo points on 'right' side (high index)
    std::array<int, NDIMS>            _array_size; // global array size
    std::array<int, NDIMS>            _nhalo_left; // number of halo points on 'left'  side (low index)
    std::map<int, HaloRegion<T, NDIMS>> _halo_map; // map from integer to halo
    DArrayLayout<NDIMS>                   _layout; // topologically-aware communicator object
    T*                                      _data; // actual data

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

    // ===================================================================== //
    // map from halo specification to an integer
    int _halo_map_hashfun(BoundaryTag tag, HaloIntent intent, size_t dim) {
        #if DARRAY_LAYOUT_CHECKBOUNDS
            _layout.checkbounds(dim);
        #endif
        return (10*dim + static_cast<int>(tag)) * static_cast<int>(intent);
    }

public:
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
        , _layout     (layout     ) {
            // define size of local array and number of left/right halo points
            for (auto dim : LinRange(NDIMS)) {
                _local_arr_size[dim] = _array_size[dim] / _layout.size(dim);
                _nhalo_left[dim]  = _layout.has_neighbour_at(BoundaryTag::LEFT, dim)  ? nhalo_in : nhalo_out[dim];
                _nhalo_right[dim] = _layout.has_neighbour_at(BoundaryTag::RIGHT, dim) ? nhalo_in : nhalo_out[dim];

                // full size of the data
                _raw_arr_size[dim] = _local_arr_size[dim] + _nhalo_left[dim] + _nhalo_right[dim];
            }

            // allocate memory buffer
            _data = new T[nelements()];

            // construct dictionary of the halo regions
            for (auto tag : {BoundaryTag::LEFT, BoundaryTag::RIGHT}) {
                for (auto intent : {HaloIntent::SEND, HaloIntent::RECV}) {
                    for (auto dim : LinRange(NDIMS)) {
                        _halo_map.emplace(_halo_map_hashfun(tag, intent, dim), 
                                   HaloRegion<T, NDIMS>(*this, tag, intent, dim));
                    }
                }
            }
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
    // RAW DATA POINTER
    inline value_type* data() const {
        return _data;
    }

    // ===================================================================== //
    // ITERATE OVER ALL DATA
    inline value_type* begin() { return _data; }
    inline value_type* end()   { return _data + nelements(); }

    // ===================================================================== //
    // LAYOUT
    inline const DArrayLayout<NDIMS>& layout () const {
        return _layout;
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

    inline int size(size_t dim) const { 
        #if DARRAY_LAYOUT_CHECKBOUNDS
            _layout.checkbounds(dim);
        #endif
        return _local_arr_size[dim]; 
    }

    // ===================================================================== //
    // NUMBER OF HALO POINTS AT A PARTICULAR BOUNDARY
    inline int nhalo(BoundaryTag tag, size_t dim) const { 
        #if DARRAY_LAYOUT_CHECKBOUNDS
            _layout.checkbounds(dim);
        #endif
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

    inline HaloRegion<T, NDIMS>& halo(BoundaryTag tag, HaloIntent intent, size_t dim) {
        #if DARRAY_LAYOUT_CHECKBOUNDS
            _layout.checkbounds(dim);
        #endif
        return _halo_map.find(_halo_map_hashfun(tag, intent, dim))->second;
    }

    // ===================================================================== //
    // UPDATE HALO POINTS
    void swap_halo() {
        for ( auto dim : LinRange(NDIMS) ) {
            // +-+---------+-+-+     +-+-+---------+-+-+     +-+-+---------+-+
            // | | (i,j-1) |S| | --> |R| |  (i,j)  |S| | --> |R| | (i,j+1) |S|
            // +-+---------+-+-+     +-+-+---------+-+-+     +-+-+---------+-+
            sendrecv(halo(BoundaryTag::RIGHT, HaloIntent::SEND, dim), 
                     _layout.rank_of_neighbour_at(BoundaryTag::RIGHT, dim),
                     halo(BoundaryTag::LEFT , HaloIntent::RECV, dim), 
                     _layout.rank_of_neighbour_at(BoundaryTag::LEFT,  dim));

            // +-+---------+-+-+     +-+-+---------+-+-+     +-+-+---------+-+
            // | | (i,j-1) | |R| <-- | |S|  (i,j)  | |R| <-- | |S| (i,j+1) | |
            // +-+---------+-+-+     +-+-+---------+-+-+     +-+-+---------+-+
            sendrecv(halo(BoundaryTag::LEFT,  HaloIntent::SEND, dim), 
                     _layout.rank_of_neighbour_at(BoundaryTag::LEFT,  dim),
                     halo(BoundaryTag::RIGHT, HaloIntent::RECV, dim), 
                     _layout.rank_of_neighbour_at(BoundaryTag::RIGHT, dim));
                     
        }
    }    
};
}