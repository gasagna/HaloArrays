#pragma once
#include<iostream>

namespace DArrays {
namespace Iterators {

template <size_t NDIMS>
class IndexIterator {
private:
    std::array<int, NDIMS> size;

    struct _Iter {
        std::array<int, NDIMS> state;
        std::array<int, NDIMS> size;

        _Iter(std::array<int, NDIMS> size, std::array<int, NDIMS> state)
            : size  (size )
            , state (state) {}

        inline std::array<int, NDIMS>& operator*() {
            return state;
        }

        inline _Iter operator++() {
            state[0]++;
            for (int dim = 0; dim != NDIMS-1; dim++) {
                if (state[dim] == size[dim]) {
                    state[dim] = 0; 
                    state[dim+1]++;
                } else {
                    break;
                }
            }
            return *this;
        }

        inline bool operator != (_Iter &other) const {
            return state != other.state;
        }
    };

public:
    template <typename... NS>
    IndexIterator(NS... ns) {
        static_assert(sizeof...(ns) == NDIMS, "invalid iterator specification");
        size = {ns...};
    }

    _Iter begin() { 
        std::array<int, NDIMS> state = {0};
        return {size, state}; 
    }
    
    _Iter end() {
        std::array<int, NDIMS> state = {0}; state[NDIMS-1] = size[NDIMS-1];
        return {size, state}; 
    }
};

}
}