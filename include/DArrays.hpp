#include "iterators.hpp"
#include "boundaries.hpp"
#include "dlayout.hpp"
#include "darray.hpp"
#include "haloregion.hpp"
#include "mpiwrapper.hpp"

// DEFAULT CONFIGURATION OPTIONS
#ifndef DARRAY_ARRAY_CHECKBOUNDS
#define DARRAY_ARRAY_CHECKBOUNDS false
#endif

#ifndef DARRAY_LAYOUT_CHECKBOUNDS
#define DARRAY_LAYOUT_CHECKBOUNDS false
#endif