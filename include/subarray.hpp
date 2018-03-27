#pragma once

namespace DArray {


////////////////////////////////////////////////////////
//                     SUBARRAY                       //
////////////////////////////////////////////////////////
template <typename T, size_t NDIMS>
class subarray {
private:
    std::array<int, NDIMS>         _size; // size of the subarray
    std::array<int, NDIMS>   _raw_origin; // origin within the raw array
    const DArray<T, NDIMS>&      _parent; // parent array

    // ===================================================================== //
    // FILL SUBARRAY SPECIFICATION DATA
    template <typename... INDICES>
    void _init(INDICES... indices) {
        __init(0, indices...);
    }

    template <typename... INDICES>
    void __init(size_t dim, int from, int to, INDICES... indices) {
        #if DARRAY_CONFIG_CHECKBOUNDS
            _checkparentbounds(dim, from, to);
        #endif
        // size of subarray
        _size[dim]       = to - from + 1;

        // origin of the subarray within the actual underlying data
        _raw_origin[dim] = from + _parent.nhalo_left(dim);

        // process further indices
        __init(dim+1, indices...);
    }

    // base case
    void __init(size_t dim) {}

    // ===================================================================== //
    // BOUND CHECKING FOR SUBARRAY CREATION
    inline void _checkparentbounds(size_t dim, int from, int to) {
        if (from > to)
            throw std::out_of_range("invalid subarray specification");

        if !( _parent.is_inbounds(from, dim) and _parent.is_inbounds(to, dim) ) {
            throw std::out_of_range("invalid subarray specification");
    }

public: 
    // constructor from the indices
    template<typename... INDICES, 
              typename ENABLER = std::enable_if_t< (... && std::is_integral_v<INDICES>) >>>
    subarray(const DArray<T, NDIMS>& parent, INDICES... indices) 
        :_parent (parent) {
            static_assert(sizeof...(INDICES) == 2*NDIMS,
                "Number of indices must match array dimension");
            __init(indices...);
    }

    // constructor from halo region
    subarray(const DArray<T, NDIMS>& parent, HaloRegion<NDIMS> region, HaloIntent intent) {
        // intent == SEND
        for ( auto dim : LinRange(NDIMS) ) {
            _size[dim]       = parent.nhalo(region[dim], dim);
            
            if region[dim] == LEFT
                _raw_origin[dim] = 0;
            if region[dim] == RIGHT
                _raw_origin[dim] = size[dim] - nhalo(RIGHT, dim);            
        }
        
        // intent == RECV
        // for ( auto dim : LinRange(NDIMS) ) {
            // _raw_origin[dim] -= parent.nhalo()
    }

    // return parent darray
    const DArray<T, NDIMS>& parent() {
        return _parent;
    }

    // the size of the subarray
    std::array<int, NDIMS> size() {
        return _size;
    }

    // actual origin of the data
    std::array<int, NDIMS> raw_origin() {
        return _raw_origin;
    }
};

}