#pragma once

namespace DArrays {

// ===================================================================== //
// Allow arbitrary memory orders
template <size_t NDIMS>
struct AbstractMemoryOrder {};

// column major fortran ordering
template <size_t NDIMS>
struct FortranMemoryOrder : AbstractMemoryOrder<NDIMS> {
    inline size_t operator [] (size_t i) { return i; }
};

// row major c ordering
template <size_t NDIMS>
struct CMemoryOrder : AbstractMemoryOrder<NDIMS> {
    inline size_t operator [] (size_t i) { return NDIMS - i - 1; }
};

// arbitrary order
template <size_t NDIMS>
class MemoryOrder : AbstractMemoryOrder<NDIMS> {
    std::array<size_t, NDIMS> _data;

    // check whether we provide a permutation of 0, 1, 2, ... NDIMS
    void _is_permutation() {
        for (auto dim : LinRange(NDIMS))
            if (std::count(_data.begin(), _data.end(), dim) != 1)
                throw std::invalid_argument("argument must be a permutation");
    }

public:
    template <typename... INDICES>
    MemoryOrder(INDICES... indices) {
        static_assert(sizeof...(INDICES) == NDIMS,
                      "Number of indices must match dimension");
        // should check positivity here
        _data = {indices...}; _is_permutation();
    }
    inline size_t operator [] (size_t i) { return _data[i]; }
};

}