#pragma once
#include "darray.hpp"

namespace DArrays {

// forward declaration
template <typename T, size_t NDIMS> class DArray;

////////////////////////////////////////////////////////
//                     HaloRegion                     //
////////////////////////////////////////////////////////
template <typename T, size_t NDIMS>
class HaloRegion {
private:
    std::array<int, NDIMS>   _raw_origin; // origin within the raw array
    const DArray<T, NDIMS>&      _parent; // handle to parent array
    std::array<int, NDIMS>         _size; // size of the HaloRegion
    MPI_Datatype                   _type;

    // ===================================================================== //    
    // CREATE SUBARRAY TYPE
    void _make_type(MPI_Datatype type) {
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
    // CONSTRUCTOR TAG, INTENT AND DIMENSION
    HaloRegion(const DArray<T, NDIMS>& parent, 
               BoundaryTag             tag, 
               HaloIntent              intent,
               size_t                  dim) 
        : _parent     (parent)
        , _size       (parent.raw_size())
        , _raw_origin ({0}) {

            // size along dimension is the halo size
            _size[dim] = _parent.nhalo(tag, dim);

            // origin is initialised to 0, except on the specified dimension
            if (tag == BoundaryTag::LEFT)
                if ( intent == HaloIntent::SEND )
                    _raw_origin[dim] = _parent.nhalo(BoundaryTag::LEFT, dim);

            if (tag == BoundaryTag::RIGHT) {
                if ( intent == HaloIntent::SEND ) {
                    _raw_origin[dim] = _parent.nhalo(BoundaryTag::LEFT, dim) + 
                                       _parent.size(dim) - 
                                       _parent.nhalo(BoundaryTag::RIGHT, dim);     
                } else {
                    _raw_origin[dim] = _parent.nhalo(BoundaryTag::LEFT, dim) + 
                                       _parent.size(dim);     
                }
            }

            // must be called after the other bits have been constructed
            _make_type(_type);
    }
    
    // ===================================================================== //    
    // DESTRUCTOR
    ~HaloRegion() {
        MPI_Type_free(&_type);
    }

    // ===================================================================== //
    // COPY CONSTRUCTOR
    HaloRegion(const HaloRegion& reg) 
        : _raw_origin (reg.raw_origin())
        , _parent     (reg.parent())
        , _size       (reg.size()) {
            _make_type(_type);
    }
                        
    // ===================================================================== //                
    // CONSTRUCTOR FROM BOUNDARY ELEMENT
    HaloRegion(const DArray<T, NDIMS>& parent, 
               const Boundary<NDIMS>&  boundary, 
               HaloIntent&             intent) 
        : _parent (parent) {
            // construct the size and origin for the case where we want to SEND the data
            for ( auto dim : LinRange(NDIMS) ) {
                _size[dim] = _parent.nhalo(boundary[dim], dim);
                
                // note that the raw origin is the coordinate including the halo points
                if (boundary[dim] == BoundaryTag::LEFT)
                    _raw_origin[dim] = _parent.nhalo(BoundaryTag::LEFT, dim);

                if (boundary[dim] == BoundaryTag::RIGHT) 
                    _raw_origin[dim] = _parent.nhalo(BoundaryTag::LEFT, dim) + 
                                       _parent.size(dim) - 
                                       _parent.nhalo(BoundaryTag::RIGHT, dim);            
            }
        
            // then shift the origin it if we actually need to RECV the data
            if ( intent == HaloIntent::RECV ) {
                for ( auto dim : LinRange(NDIMS) ) {
                    switch ( boundary[dim] ) {
                        case BoundaryTag::LEFT :
                            _raw_origin[dim] -= _parent.nhalo(boundary[dim], dim); break;
                        case BoundaryTag::RIGHT:
                            _raw_origin[dim] += _parent.nhalo(boundary[dim], dim); break;
                        case BoundaryTag::CENTER:
                            break;
                    }
                }
            }
            
            // create subarray type
            MPI_Type_create_subarray(NDIMS,
                                     _parent.raw_size().data(),
                                     _size.data(),             
                                     _raw_origin.data(),       
                                     MPI_ORDER_FORTRAN,                
                                     MPI_DOUBLE, &_type);
            MPI_Type_commit(&_type);
    }

    // ===================================================================== //    
    // RETURN MPI DATA TYPE FOR THE SUBARRAY
    inline const MPI_Datatype& type() const {
        return _type;
    }

    // ===================================================================== //    
    // RETURN PARENT DARRAY
    inline const DArray<T, NDIMS>& parent() const {
        return _parent;
    }

    // ===================================================================== //
    // THE SIZE OF THE HALOREGION
    inline const std::array<int, NDIMS>& size() const {
        return _size;
    }

    // ===================================================================== //
    // ACTUAL ORIGIN OF THE DATA
    inline const std::array<int, NDIMS>& raw_origin() const {
        return _raw_origin;
    }
};

}