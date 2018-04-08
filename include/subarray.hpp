#pragma once
#include "darray.hpp"

namespace DArrays {

// forward declaration
template <typename T, size_t NDIMS> class DArray;

////////////////////////////////////////////////////////
//                     SubArray                       //
////////////////////////////////////////////////////////
template <typename T, size_t NDIMS>
class SubArray {
private:
    std::array<int, NDIMS>   _raw_origin; // origin within the raw array
    const DArray<T, NDIMS>&      _parent; // handle to parent array
    std::array<int, NDIMS>         _size; // size of the SubArray

public:
    // constructor from boundary element
    SubArray(const DArray<T, NDIMS>& parent, 
             BoundaryTag tag, 
             size_t dim,
             BoundaryIntent intent) 
        : _parent     (parent)
        , _size       (parent.raw_size())
        , _raw_origin ({0}) {

            // we squeeze the data on this dimension
            _size[dim] = _parent.nhalo(tag, dim);

            // origin is initialised to 0, except on the specified dimension
            if (tag == BoundaryTag::LEFT)
                if ( intent == BoundaryIntent::SEND )
                    _raw_origin[dim] = _parent.nhalo(BoundaryTag::LEFT, dim);

            if (tag == BoundaryTag::RIGHT) {
                if ( intent == BoundaryIntent::SEND ) {
                    _raw_origin[dim] = _parent.nhalo(BoundaryTag::LEFT, dim) + 
                                       _parent.size(dim) - 
                                       _parent.nhalo(BoundaryTag::RIGHT, dim);     
                } else {
                    _raw_origin[dim] = _parent.nhalo(BoundaryTag::LEFT, dim) + 
                                       _parent.size(dim);     
                }
            }
        }
    
    // constructor from boundary element
    SubArray(const DArray<T, NDIMS>& parent, 
             const Boundary<NDIMS>& boundary, 
             const BoundaryIntent& intent) 
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
            if ( intent == BoundaryIntent::RECV ) {
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
    }

    // return parent darray
    inline const DArray<T, NDIMS>& parent() const {
        return _parent;
    }

    // the size of the SubArray
    inline const std::array<int, NDIMS>& size() const {
        return _size;
    }

    // actual origin of the data
    inline const std::array<int, NDIMS>& raw_origin() const {
        return _raw_origin;
    }
};

}