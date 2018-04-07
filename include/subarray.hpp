#pragma once

namespace DArrays {

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