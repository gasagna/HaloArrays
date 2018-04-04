#pragma once
#include <iostream>
#include <iterator>
#include <numeric>
#include <array>
#include <cmath>

namespace DArrays::Iterators {

////////////////////////////////////////////////////////////////
// Integer linear range, from FROM to TO, end point excluded. //
// This has the same semantics of the python range function   //
////////////////////////////////////////////////////////////////
class LinRange {
private:
    // range iterator
    class _LinRangeIter {
        private:
            int _state;
            int  _step;
        public:
            _LinRangeIter (int state, int step) 
                : _state (state)
                , _step  (step ) {} 

            inline int operator * () const { 
                return _state; 
            }

            inline _LinRangeIter operator ++ () {
                _state += _step;
                return *this;
            }

            inline bool operator != (_LinRangeIter& other) const {
                return _state != other._state;
            }
        };
    
    // range end points: included
    int _from, _to, _step;
public:
    LinRange(int from, int to, int step) : _from (from), _to (to), _step (step) {}
    LinRange(int from, int to)           : _from (from), _to (to), _step (1)    {}
    LinRange(int to)                     : _from (0),    _to (to), _step (1)    {}
    _LinRangeIter begin() { return {_from, _step}; }
    _LinRangeIter   end() { return {static_cast<int>(_from + 
                            _step*std::ceil((_to - _from )/(double)_step)), _step}; } 
};


////////////////////////////////////////////////////////
//                  index range                       //
////////////////////////////////////////////////////////
template <size_t NDIMS>
class IndexRange {
private:
    std::array<int, NDIMS> _size;      // array size

    class _IndexRangeIter {
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
        // EXPAND STATE TO/FROM LINEARISED INDEX
        inline difference_type _tolinearindex() const {
            difference_type n = 0;
            for ( auto dim : LinRange(NDIMS) )
                n += _size_prod[dim] * _state[dim];
            return n;
        }

        inline void _fromlinearindex(difference_type n) {
            div_t divrem;            
            for ( auto dim : LinRange(NDIMS-1, -1, -1) ) {
                divrem = div(n, _size_prod[dim]);
                _state[dim] = divrem.quot;
                n = divrem.rem;
            }
        }

    public:
        // ===================================================================== //
        // CONSTRUCTOR/DESTRUCTOR
        _IndexRangeIter(std::array<int, NDIMS> size, std::array<int, NDIMS> state)
            : _size       (size )  
            , _state      (state) {
                // compute product of array sizes
                _size_prod[0] = 1;
                for (auto dim : LinRange(1, NDIMS)) {
                    _size_prod[dim] = _size_prod[dim-1]*_size[dim-1];
                }
            }

        // ===================================================================== //
        // DEREFERENCING
        inline reference operator * () {
            return _state;
        }

        // ===================================================================== //
        // INCREMENT
        inline _IndexRangeIter operator ++ () {
            _state[0]++;
            // TODO: benchmark this compare to simpler loop. Is the
            // compiler able to unroll this efficiently?
            for ( auto dim : LinRange(NDIMS-1) ) {
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
        inline _IndexRangeIter& operator += (difference_type n) {
            _fromlinearindex(_tolinearindex() + n);
            return *this;
        }

        inline _IndexRangeIter& operator -= (difference_type n) {
            _fromlinearindex(_tolinearindex() - n);
            return *this;
        }

        // ===================================================================== //
        // EQUALITY AND COMPARISON
        inline bool operator == (const _IndexRangeIter& other) const {
            return _state == other._state;
        }

        inline bool operator != (const _IndexRangeIter& other) const {
            return !(*this == other);
        }

        inline bool operator < (const _IndexRangeIter& other) const {
            return _tolinearindex() < other._tolinearindex();
        }

        inline bool operator > (const _IndexRangeIter& other) const {
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
        static_assert(sizeof...(ns) == NDIMS, "too many indiced for iterator dimension");
        _size = {ns...};
    }

    // from an array of integer sizes
    template<typename T, 
            typename ENABLER = std::enable_if_t< std::is_integral_v<T> >>
    IndexRange(std::array<T, NDIMS> size) 
        : _size (size) {}

    _IndexRangeIter begin() { 
        std::array<int, NDIMS> _state = {0};
        return {_size, _state}; 
    }
    
    _IndexRangeIter end() {
        // constuct state for one past the last
        std::array<int, NDIMS> _state = {0}; _state[NDIMS - 1] = _size[NDIMS - 1];
        return {_size, _state};
    }
};

} // namespace Darrays::Iterators