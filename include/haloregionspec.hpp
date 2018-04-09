#pragma once

using namespace DArrays::Iterators;

namespace DArrays {

// ===================================================================== //
// tags for the the left, center and right boundaries
enum class Boundary: int {LEFT = 1, CENTER = 2, RIGHT = 4, WILDCARD = 8};


// ===================================================================== //
// tags for sending or receiving the halo region 
enum class HaloIntent : int {SEND = -1, RECV = 1};


// ===================================================================== //
// HaloRegionSpec                     
template<size_t NDIMS>
class HaloRegionSpec {
private:
    std::array<Boundary, NDIMS> _speclist; // array of boundary tags
    int                            _uhash; // unsigned hash

    void _init_hash() {
        // compute hash
        for (auto dim : LinRange(NDIMS))
            _uhash += std::pow(10, dim)*static_cast<int>(_speclist[dim]);

        // we cannot have as many wildcards as there are dimensions
        if (_uhash == 8 or _uhash == 88 or _uhash == 888)
            throw std::invalid_argument("invalid halo region specification");        
    }

public:
    // constructors
    template <typename... SPECLIST>
    HaloRegionSpec(SPECLIST... speclist) 
        : _speclist ({speclist...}), _uhash (0) { _init_hash(); }
    
    HaloRegionSpec(const std::array<Boundary, NDIMS>& speclist) 
        : _speclist (speclist), _uhash (0) { _init_hash(); }

    // read only array like object
    inline const Boundary operator [] (size_t dim) const {
        #if DARRAY_LAYOUT_CHECKBOUNDS
            _checkdims(dim, NDIMS);
        #endif
        return _speclist[dim];
    }

    // define a unique integer for each halo region specification
    inline int hash(HaloIntent intent) const {
        return _uhash*static_cast<int>(intent);
    }
};

// construct the opposite halo region
template<size_t NDIMS>
inline HaloRegionSpec<NDIMS> opposite(const HaloRegionSpec<NDIMS>& spec) {
    std::array<Boundary, NDIMS> _speclist;
    for (auto dim : LinRange(NDIMS)) {
        // swap the LEFT/RIGHT boundary and keep the rest untouched 
        _speclist[dim] = spec[dim];
        if (spec[dim] == Boundary::LEFT)  _speclist[dim] = Boundary::RIGHT;
        if (spec[dim] == Boundary::RIGHT) _speclist[dim] = Boundary::LEFT;
    }
    return HaloRegionSpec<NDIMS>(_speclist);
}

// ===================================================================== //
// list of all halo regions
// 1D sequence is L, R
static const std::array<HaloRegionSpec<1>, 2> _specs_1d = {{
        HaloRegionSpec<1>(Boundary::LEFT),
        HaloRegionSpec<1>(Boundary::RIGHT)}};


// 2D sequence is: L*, R*, CL, CR
static const std::array<HaloRegionSpec<2>, 4> _specs_2d = {{
        HaloRegionSpec<2>(Boundary::LEFT,   Boundary::WILDCARD), // N
        HaloRegionSpec<2>(Boundary::RIGHT,  Boundary::WILDCARD), // S
        HaloRegionSpec<2>(Boundary::CENTER, Boundary::LEFT),     // W
        HaloRegionSpec<2>(Boundary::CENTER, Boundary::RIGHT)}};  // E

// 3D sequence is: L**, R**, CR*, CL*, CCL, CCR
static const std::array<HaloRegionSpec<3>, 6> _specs_3d = {{
        HaloRegionSpec<3>(Boundary::LEFT,   Boundary::WILDCARD, Boundary::WILDCARD),
        HaloRegionSpec<3>(Boundary::RIGHT,  Boundary::WILDCARD, Boundary::WILDCARD),
        HaloRegionSpec<3>(Boundary::CENTER, Boundary::LEFT,     Boundary::WILDCARD),
        HaloRegionSpec<3>(Boundary::CENTER, Boundary::RIGHT,    Boundary::WILDCARD),
        HaloRegionSpec<3>(Boundary::CENTER, Boundary::CENTER,   Boundary::LEFT),
        HaloRegionSpec<3>(Boundary::CENTER, Boundary::CENTER,   Boundary::RIGHT)}};

// tuple to collect all bits together
static const auto _halospeclist = std::make_tuple(0, _specs_1d, _specs_2d, _specs_3d);

}