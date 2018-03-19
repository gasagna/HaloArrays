#pragma once

#ifndef DARRAY_CONFIG_CHECKBOUNDS
#define DARRAY_CONFIG_CHECKBOUNDS false
#endif

namespace DArrays {

template <typename T, size_t NDIMS>
class DArray {
private:
    T*                 _data; // actual data
    DArraySpec<NDIMS>  _spec; // array size object

public:
    // ===================================================================== //
    // CONTAINER INTERFACE
    using value_type     = T;
    using DArrayIterator = T*;

    // ===================================================================== //
    // CONSTRUCTOR/DESTRUCTOR
    DArray() = delete;

    DArray(DArraySpec<NDIMS>& spec)
        : _spec (spec) 
        , _data (new T[spec.nelements()]) {}

    ~DArray() {
        delete[] _data;
    }

    // ===================================================================== //
    // LINEAR INDEXING
    const inline T& operator [] (size_t i) const { return _data[i]; }
          inline T& operator [] (size_t i)       { return _data[i]; }

    // ===================================================================== //
    // LOCAL CARTESIAN INDEXING - FIRST INDEX IS CONTIGUOUS
    template <typename... INDICES>
    inline T& operator () (INDICES... indices) {
        static_assert(sizeof...(INDICES) == NDIMS,
                      "Number of indices must match array dimension");
        #if DARRAY_CONFIG_CHECKBOUNDS
            _spec.checkbounds(indices...);
        #endif
        return _data[_spec.tolinearindex(indices...)];
    }

    // ===================================================================== //
    // ITERATION OVER ALL DATA
    DArrayIterator begin() { return _data; }
    DArrayIterator end()   { return _data + _spec.nelements(); }

    // ===================================================================== //
    // ITERATOR OVER THE IN-DOMAIN INDICES 
    inline IndexIterator<NDIMS> indices () {
        return IndexIterator<NDIMS>(_spec.size())
    }

    // ===================================================================== //
    // UPDATE HALO POINTS - specialisations for NDIMS = 1, 2, 3 elsewhere
    void halo_swap();
};

}