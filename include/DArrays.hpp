#pragma once
#include <initializer_list>
#include <functional>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <array>
#include <map>
#include <mpi.h>

namespace DArrays {

// ===================================================================== //
// check whether we are not getting out of bounds with the dimension
inline void _checkdims(size_t dim, size_t NDIMS) {
    // dim < 0 never, for size_t
    if ( dim >= NDIMS )
        throw std::out_of_range("dimension out of range");
}

}

#include "iterators.hpp"
#include "haloregionspec.hpp"
#include "dlayout.hpp"
#include "darray.hpp"
#include "subarray.hpp"
#include "mpiwrapper.hpp"

// DEFAULT CONFIGURATION OPTIONS
#ifndef DARRAY_ARRAY_CHECKBOUNDS
#define DARRAY_ARRAY_CHECKBOUNDS false
#endif

#ifndef DARRAY_LAYOUT_CHECKBOUNDS
#define DARRAY_LAYOUT_CHECKBOUNDS false
#endif