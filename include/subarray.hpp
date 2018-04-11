#pragma once
#include "darray.hpp"

namespace DArrays {

// forward declaration
template <typename T, size_t NDIMS> class DArray;

// ===================================================================== //
// SubArray                     
template <typename T, size_t NDIMS>
class SubArray {
private:
    std::array<int, NDIMS>   _raw_origin; // origin within the raw array
    const DArray<T, NDIMS>&      _parent; // handle to parent array
    std::array<int, NDIMS>         _size; // size of the SubArray
    MPI_Datatype                   _type;

    // ===================================================================== //    
    // init subarray type
    void _init_type(MPI_Datatype type) {
        MPI_Type_create_subarray(NDIMS,
                                 _parent.raw_size().data(),
                                 _size.data(),             
                                 _raw_origin.data(),       
                                 MPI_ORDER_FORTRAN,                
                                 MPI_DOUBLE, &_type);
        MPI_Type_commit(&_type);
    }
    
public:                        
    // ===================================================================== //                
    // constructor from halo region specification and intent
    SubArray(DArray<T, NDIMS>&            parent, 
             const HaloRegionSpec<NDIMS>& spec,  
             HaloIntent                   intent) 
        : _raw_origin ({0})
        , _parent     (parent)
        , _size       ({0}) {
            // construct the size and origin for the case where we want to SEND the data
            for ( auto dim : LinRange(NDIMS) ) {
                _size[dim] = _parent.nhalo_points(spec[dim], dim);
                
                // the raw origin on the underlying array data
                if (spec[dim] == Boundary::LEFT or spec[dim] == Boundary::CENTER)
                    _raw_origin[dim] = _parent.nhalo_points(Boundary::LEFT, dim);

                if (spec[dim] == Boundary::RIGHT) 
                    _raw_origin[dim] = _parent.nhalo_points(Boundary::LEFT, dim) + 
                                       _parent.size(dim) - 
                                       _parent.nhalo_points(Boundary::RIGHT, dim);            
            }
        
            // then shift the origin if we actually need to RECV the data
            if (intent == HaloIntent::RECV) {
                for (auto dim : LinRange(NDIMS)) {
                    switch (spec[dim]) {
                        case Boundary::LEFT :
                            _raw_origin[dim] -= _parent.nhalo_points(spec[dim], dim); break;
                        case Boundary::RIGHT:
                            _raw_origin[dim] += _parent.nhalo_points(spec[dim], dim); break;
                        case Boundary::CENTER:
                            break;
                        case Boundary::WILDCARD:
                            break;
                    }
                }
            }
            
            // create subarray type after everything else
            _init_type(_type); 
    }

    // ===================================================================== //    
    // destructor
    ~SubArray() {
        MPI_Type_free(&_type);
    }

    // ===================================================================== //
    // copy constructor
    SubArray(const SubArray& reg) 
        : _raw_origin (reg.raw_origin())
        , _parent     (reg.parent())
        , _size       (reg.size()) {
            _init_type(_type);
    }

    // ===================================================================== //    
    // return mpi data type for the subarray
    inline const MPI_Datatype& type() const {
        return _type;
    }

    // ===================================================================== //    
    // return parent darray
    inline const DArray<T, NDIMS>& parent() const {
        return _parent;
    }

    // ===================================================================== //
    // the size of the SubArray
    inline const std::array<int, NDIMS>& size() const {
        return _size;
    }

    inline int size(size_t dim) const {
        #if DARRAY_LAYOUT_CHECKBOUNDS
            _checkdims(dim, NDIMS);
        #endif
        return _size[dim];
    }

    // ===================================================================== //
    // actual origin of the data
    inline const std::array<int, NDIMS>& raw_origin() const {
        return _raw_origin;
    }

    inline int raw_origin(size_t dim) const {
        #if DARRAY_LAYOUT_CHECKBOUNDS
            _checkdims(dim, NDIMS);
        #endif
        return _raw_origin[dim];
    }
};

}