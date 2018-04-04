#include "iterators.hpp"
#include "haloregion.hpp"
#include "dlayout.hpp"
#include "darray.hpp"
#include "mpiwrapper.hpp"

// DEFAULT CONFIGURATION OPTIONS
#ifndef DARRAY_CONFIG_CHECKBOUNDS
#define DARRAY_CONFIG_CHECKBOUNDS false
#endif