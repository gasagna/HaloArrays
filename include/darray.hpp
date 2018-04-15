#pragma once
#include "subarray.hpp"
#include "mpiwrapper.hpp"

namespace DArrays {

using namespace DArrays::Iterators;
using namespace DArrays::MPI;

// forward declaration
template <typename T, size_t NDIMS> class SubArray;

// ===================================================================== //
// DArray
template <typename T, size_t NDIMS>
class DArray {
private:
    std::array<int, NDIMS>            _local_arr_size; // local array size
    std::array<int, NDIMS>              _raw_arr_size; // local array size, including halo points
    std::map<int, SubArray<T, NDIMS>>   _subarray_map; // map from integer to halo
    std::array<int, NDIMS>               _nhalo_right; // number of halo points on 'right' side (high index)
    std::array<int, NDIMS>                _array_size; // global array size
    std::array<int, NDIMS>                _nhalo_left; // number of halo points on 'left'  side (low index)
    DArrayLayout<NDIMS>                       _layout; // topologically-aware communicator object
    T*                                          _data; // actual data

    // ===================================================================== //
    // indexing into linear memory buffer
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
    // bound checking
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
    // check whether index i in in bounds along dimension dim, including halo
    inline bool _is_inbounds(int i, size_t dim) {
        return (i >= - _nhalo_left[dim] and i <= _local_arr_size[dim] + _nhalo_right[dim] - 1);
    }

    // ===================================================================== //
    // index into the dictionary HaloRegionSpec->SubArray
    inline SubArray<T, NDIMS>& _get_subarray(HaloRegionSpec<NDIMS> spec, HaloIntent intent) {
        return _subarray_map.find(spec.hash(intent))->second;
    }

public:
    // ===================================================================== //
    // container interface
    using value_type = T;

    // ===================================================================== //
    // constructor/destructor
    DArray() = delete;

    DArray(DArrayLayout<NDIMS> layout, 
           std::array<int, NDIMS> array_size,
           std::array<int, NDIMS> nhalo_out, 
           std::array<int, NDIMS> nhalo_in)
        : _array_size (array_size ) 
        , _layout     (layout     ) {
            // define size of local array and number of left/right halo points
            for (auto dim : LinRange(NDIMS)) {
                _local_arr_size[dim] = _array_size[dim] / _layout.size(dim);
                _nhalo_left[dim]  = _layout.has_neighbour_at(Boundary::LEFT,  dim) ? nhalo_in[dim] : nhalo_out[dim];
                _nhalo_right[dim] = _layout.has_neighbour_at(Boundary::RIGHT, dim) ? nhalo_in[dim] : nhalo_out[dim];

                // full size of the data
                _raw_arr_size[dim] = _local_arr_size[dim] + _nhalo_left[dim] + _nhalo_right[dim];
            }

            // note that the halo size cannot be larger then the data itself
            for (auto dim : LinRange(NDIMS))
                if (std::max(_nhalo_left[dim], _nhalo_right[dim]) >= _local_arr_size[dim])
                    throw std::invalid_argument("too many halo points for local array size");

            // allocate memory buffer
            _data = new T[nelements()];

            // construct dictionary of the halo regions used for halo swap
            for (auto& spec : std::get<NDIMS>(_halospeclist))
                for (auto intent : {HaloIntent::SEND, HaloIntent::RECV})
                    _subarray_map.emplace(spec.hash(intent), 
                                          SubArray<T, NDIMS>(*this, spec, intent));
    }

    ~DArray() {
        delete[] _data;
    }

    // ===================================================================== //
    // indexing into linear memory buffer
    const inline T& operator [] (size_t i) const { return _data[i]; }
          inline T& operator [] (size_t i)       { return _data[i]; }

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
    // raw data pointer
    inline value_type* data() const {
        return _data;
    }

    // ===================================================================== //
    // iterate over all data
    inline value_type* begin() { return _data; }
    inline value_type* end()   { return _data + nelements(); }

    // ===================================================================== //
    // layout
    inline const DArrayLayout<NDIMS>& layout () const {
        return _layout;
    }

    // ===================================================================== //
    // iterator over the in-domain indices 
    inline IndexRange<NDIMS> indices () {
        return IndexRange<NDIMS>(_local_arr_size);
    }

    // ===================================================================== //
    // local array size
    inline const std::array<int, NDIMS>& size() const { 
        return _local_arr_size; 
    }

    inline int size(size_t dim) const { 
        #if DARRAY_LAYOUT_CHECKBOUNDS
            _checkdims(dim, NDIMS);
        #endif
        return _local_arr_size[dim]; 
    }

    // ===================================================================== //
    // number of halo points at a particular boundary  
    inline int nhalo_points(Boundary bnd, size_t dim) const { 
        #if DARRAY_LAYOUT_CHECKBOUNDS
            _checkdims(dim, NDIMS);
        #endif
        switch (bnd) {
            case Boundary::LEFT     : return _nhalo_left[dim];
            case Boundary::RIGHT    : return _nhalo_right[dim];
            case Boundary::CENTER   : return _local_arr_size[dim];
            case Boundary::WILDCARD : return _raw_arr_size[dim];
        }
    }

    // ===================================================================== //
    // local array size, including halo elements
    inline const std::array<int, NDIMS>& raw_size() const { 
        return _raw_arr_size; 
    }
    
    // ===================================================================== //
    // size of memory buffer, including halo
    inline size_t nelements() const {
        return std::reduce(_raw_arr_size.begin(), _raw_arr_size.end(), 
                           1, std::multiplies<>());
    }
    
    // ===================================================================== //
    // swap halo points with neighbours
    void swap_halo() {
        for (const auto& halo_spec : std::get<NDIMS>(_halospeclist))
            sendrecv(_get_subarray(halo_spec, HaloIntent::SEND),           
                     _layout.rank_of_neighbour_at(halo_spec),
                     _get_subarray(opposite(halo_spec), HaloIntent::RECV), 
                     _layout.rank_of_neighbour_at(opposite(halo_spec)));
    }    
};
}