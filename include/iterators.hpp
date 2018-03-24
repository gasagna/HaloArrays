#pragma once
#include <iostream>
#include <iterator>
#include <numeric>
#include <array>

namespace DArrays::Iterators {

////////////////////////////////////////////////////////
// linear range, from FROM to TO, end points included //
////////////////////////////////////////////////////////
class LinearRange {
private:
    // range iterator
    class RangeIter {
        private:
            int      _state;
            bool _isforward;
        public:
            RangeIter (int state, bool isforward) 
                : _state     (state    )
                , _isforward (isforward) {} 

            inline int operator * () const { 
                return _state; 
            }

            inline RangeIter operator ++ () {
                _isforward ? _state++ : _state--;
                return *this;
            }

            inline bool operator != (RangeIter& other) const {
                return _state != other._state;
            }
        };
    
    // range end points: included
    int _from, _to;
public:
    LinearRange(int from, int to) 
        : _from (from), _to (to) {}
    RangeIter begin() { return {_from, _from < _to}; }
    RangeIter   end() { return {_from < _to ? _to + 1 : _to - 1, _from < _to}; }
};


////////////////////////////////////////////////////////
//                  index range                       //
////////////////////////////////////////////////////////
template <size_t NDIMS>
class IndexRange {
private:
    std::array<int, NDIMS> _size;      // array size

    class _Iter {
    public:
        // ===================================================================== //        
        // ITERATOR TRAITS
        using difference_type = int;
        using value_type = std::array<int, NDIMS>;
        using reference = std::array<int, NDIMS>&;
        using iterator_category = std::random_access_iterator_tag;
        // using pointer = const long*;

    private:
        std::array<int, NDIMS> _state;     // current indices  // e.g. {1, 2, 3}
        std::array<int, NDIMS> _size;      // array sizes      // e.g. {2, 3, 4}
        std::array<int, NDIMS> _size_prod; // product of sizes // e.g. {1, 2, 6}

        // ===================================================================== //
        // ADD LINEARISED INDEX TO THE STATE
        inline void _addlinearindex(difference_type n) {
            div_t divrem;
            for ( auto dim : LinearRange(NDIMS-1, 0) ) {
                divrem = div(n, _size_prod[dim]);
                _state[dim] += divrem.quot;
                n = divrem.rem;
            }
        }

        // ===================================================================== //
        // EXPAND STATE TO LINEARISED INDEX
        inline difference_type _tolinearindex() {
            difference_type n = 0;
            for ( auto dim : LinearRange(0, NDIMS-1) )
                n += _size_prod[dim] * _state[dim];
            return n;
        }

    public:
        // ===================================================================== //
        // CONSTRUCTOR/DESTRUCTOR
        _Iter(std::array<int, NDIMS> size, std::array<int, NDIMS> state)
            : _size       (size )  
            , _state      (state) {
                _size_prod[0] = 1;
                for (auto dim : LinearRange(1, NDIMS-1)) {
                    _size_prod[dim] = _size_prod[dim-1]*_size[dim];
                }
                // std::exclusive_scan(_size.begin(),
                                    // _size.end(),
                                    // _size_prod.begin(), 1, std::multiplies<int>());
            }

        // ===================================================================== //
        // DEREFERENCING
        inline reference operator * () {
            return _state;
        }

        // ===================================================================== //
        // INCREMENT
        inline _Iter operator ++ () {
            _state[0]++;
            // TODO: benchmark this compare to simpler loop. Is the
            // compiler able to unroll this efficiently?
            for ( auto dim : LinearRange(0, NDIMS-2) ) {
                if (_state[dim] == _size[dim]) {
                    _state[dim] = 0; 
                    _state[dim+1]++;
                } else {
                    break;
                }
            }
            return *this;
        }

        // ===================================================================== //
        // ADD/REMOVE LINEAR INDEX
        inline _Iter& operator += (difference_type n) {
            _addlinearindex(n);
            return *this;
        }

        inline _Iter& operator -= (difference_type n) {
            return *this += -n;
        }

        // ===================================================================== //
        // EQUALITY AND COMPARISON
        inline bool operator == (_Iter& other) const {
            return _state == other._state;
        }

        inline bool operator != (_Iter& other) const {
            return !(*this == other);
        }

        inline bool operator < (_Iter& other) const {
            return _tolinearindex() < other._tolinearindex();
        }

        inline bool operator > (_Iter& other) const {
            return _tolinearindex() > other._tolinearindex();
        }
    };

public:
    // ===================================================================== //
    // CONSTRUCTORS
    // from integer arguments
    template <typename... NS,
            typename ENABLER = std::enable_if_t< (... && std::is_integral_v<NS>) >>
    IndexRange(NS... ns) {
        static_assert(sizeof...(ns) == NDIMS, "invalid iterator specification");
        _size = {ns...};
    }

    // from an array of integer sizes
    template<typename T, 
            typename ENABLER = std::enable_if_t< std::is_integral_v<T> >>
    IndexRange(std::array<T, NDIMS> size) 
        : _size (size) {}

    _Iter begin() { 
        std::array<int, NDIMS> _state = {0};
        return {_size, _state}; 
    }
    
    _Iter end() {
        // constuct state for one past the last
        std::array<int, NDIMS> _state = {0}; _state[NDIMS - 1] = _size[NDIMS - 1];
        return {_size, _state};
    }
};

} // namespace Darrays::Iterators