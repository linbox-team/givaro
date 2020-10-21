/* ruint/ruint.h - Class definition of ruint<K> from RecInt library

   Copyright Université Grenoble Alpes
Contributors :
Alexis BREUST (alexis.breust@gmail.com 2014)
Christophe CHABOT (christophechabotcc@gmail.com 2011)
Jean-Guillaume Dumas

Time-stamp: <20 Jun 12 10:31:24 Jean-Guillaume.Dumas@imag.fr>

This software is a computer program whose purpose is to provide an fixed precision arithmetic library.

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/

#ifndef RUINT_RUINT_H
#define RUINT_RUINT_H

#include "recdefine.h"

// --------------------------------------------------------------
// ---------------- Declaration of class ruint ------------------

namespace RecInt
{
    /* Basic definition of ruint */
    template <size_t K> class ruint {
    public:
        // this = High|Low
        //ruint<K-1> High, Low;
        ruint<K-1> Low, High;

        // Constructors
        ruint() {}
        ruint(const ruint<K>& r) : Low(r.Low), High(r.High)  {}
        ruint(const ruint<K-1>& rl) : Low(rl) {}
        ruint(const double b) : Low((b < 0)? -b : b) { if (b < 0) *this = -*this; }
        template <typename T, __RECINT_IS_UNSIGNED(T, int) = 0> ruint(const T b) : Low(b) {}
        template <typename T, __RECINT_IS_SIGNED(T, int) = 0>   ruint(const T b) : Low((b < 0)? -b : b)
        { if (b < 0) *this = -*this; }
        template <typename T, __RECINT_IS_NOT_FUNDAMENTAL(T, int) = 0> ruint(const T& b)
        { *this = b.operator ruint<K>(); } // Fix for Givaro::Integer

        ruint(const char* s);

        // type_string
        static const std::string type_string () {
            return "RecInt::ruint<" + std::to_string(K)+ ">";
        }

        // Cast
        // Note: Templated operators and specialization make compilers clang + icpc
        // completely bug (they do not use the template operator somehow)
        // This fix is brutal, but, it works - AB 2015/02/11
        operator bool() const { return (High != 0) || (Low != 0); }
        operator char() const { return char(Low); }
        operator short() const { return short(Low); }
        operator int() const { return int(Low); }
        operator long() const { return long(Low); }
        operator long long() const { return (long long)(Low); }
        operator unsigned char() const { return (unsigned char)(Low); }
        operator unsigned short() const { return (unsigned short)(Low); }
        operator unsigned int() const { return (unsigned int)(Low); }
        operator unsigned long() const { return (unsigned long)(Low); }
        operator unsigned long long() const { return (unsigned long long)(Low); }
        operator float() const { return (float)(Low); }
        operator double() const { return (double)(Low); }
        template <typename T, __RECINT_IS_ARITH(T, int) = 0> operator T() const { return T(Low); }

        // Const reverse iterator
        class cr_iterator {
        public:
            cr_iterator() : _index(-1), _ptr(NULL) {}
            cr_iterator(int index, const ruint<K> *ptr) : _index(index), _ptr(ptr) {}
            cr_iterator& operator++() { --_index; return *this; }
            cr_iterator  operator++(int) { cr_iterator i(*this); --_index; return i; }
            bool operator==(cr_iterator rhs) { return (rhs._index == _index); }
            bool operator!=(cr_iterator rhs) { return (rhs._index != _index); }
            const limb* operator->() { return get_limb_p(*_ptr, (unsigned int)(_index)); }
            const limb& operator*() { return *get_limb_p(*_ptr, (unsigned int)(_index)); }
        private:
            int _index;
            const ruint<K> *_ptr;
        };

        // Const reverse iterator functions
        cr_iterator rbegin() const { return cr_iterator(NBLIMB<K>::value-1, this); }
        cr_iterator rend() const { return cr_iterator(); }
        UDItype size() { return NBLIMB<K>::value; }
        static ruint<K> maxCardinality() ;
        static ruint<K> maxModulus() ;
    };

    // Break points
    template <> class ruint<0> {};
    template <> class ruint<1> {};
    template <> class ruint<2> {};
    template <> class ruint<3> {};
    template <> class ruint<4> {};
    template <> class ruint<5> {}; // 2^32

    /* ruint of size 64 bits */
    template <> class ruint<__RECINT_LIMB_SIZE> {
    public:
        limb Value;

        // Constructors
        ruint() : Value(0u) {}
        ruint(const ruint<__RECINT_LIMB_SIZE>& r) : Value(r.Value) {}
        ruint(const double b) : Value(static_cast<limb>(b)) {}
        template <typename T, __RECINT_IS_UNSIGNED(T, int) = 0> ruint(const T b) : Value(limb(b)) {}
        template <typename T, __RECINT_IS_SIGNED(T, int) = 0> ruint(const T b) : Value(limb(b)) {}
        template <typename T, __RECINT_IS_NOT_FUNDAMENTAL(T, int) = 0> ruint(const T& b)
        { *this = b.operator ruint<__RECINT_LIMB_SIZE>(); } // Fix for Givaro::Integer
        ruint(const char* s);

        // type_string
        static const std::string type_string () {
            return "RecInt::ruint<" + std::to_string(__RECINT_LIMB_SIZE)+ ">";
        }

        // Cast
        // Brutal too, but icc is kind of peaky - AB 2015/02/11
        operator bool() const { return bool(Value); }
        operator char() const { return char(Value); }
        operator short() const { return short(Value); }
        operator int() const { return int(Value); }
        operator long() const { return long(Value); }
        operator long long() const { return (long long)(Value); }
        operator unsigned char() const { return (unsigned char)(Value); }
        operator unsigned short() const { return (unsigned short)(Value); }
        operator unsigned int() const { return (unsigned int)(Value); }
        operator unsigned long() const { return (unsigned long)(Value); }
        operator unsigned long long() const { return (unsigned long long)(Value); }
        operator float() const { return (float)(Value); }
        operator double() const { return (double)(Value); }
        template <typename T, __RECINT_IS_ARITH(T, int) = 0> operator T() const { return T(Value); }

        // Const reverse iterator
        class cr_iterator {
        public:
            cr_iterator() : _ptr(NULL) {}
            cr_iterator(const limb *ptr) : _ptr(ptr) {}
            cr_iterator& operator++() { _ptr = NULL; return *this; }
            cr_iterator  operator++(int) { cr_iterator i(*this); _ptr = NULL; return i; }
            bool operator==(cr_iterator rhs) { return (rhs._ptr == _ptr); }
            bool operator!=(cr_iterator rhs) { return (rhs._ptr != _ptr); }
            const limb* operator->() { return _ptr; }
            const limb& operator*() { return *_ptr; }
        private:
            const limb *_ptr;
        };

        // Const reverse iterator functions
        cr_iterator rbegin() const { return cr_iterator(&Value); }
        cr_iterator rend() const { return cr_iterator(); }
        UDItype size() { return 1; }
        static ruint<__RECINT_LIMB_SIZE> maxCardinality();
        static ruint<__RECINT_LIMB_SIZE> maxModulus();
    };

#if defined(__RECINT_USE_FAST_128)
    /* ruint of size 128 bits */
    template <> class ruint<__RECINT_LIMB_SIZE+1> {
    public:
        // Anonymous structs/unions put definitions in parent scope - A.B. 05-03-2015
        union {
            // FIXME Check BYTE_ORDER, invert Low, High if needed
            struct {
                ruint<__RECINT_LIMB_SIZE> Low;
                ruint<__RECINT_LIMB_SIZE> High;
            };
            __uint128_t Value;
        };

        // Constructors
        ruint() : Value(0u) {}
        ruint(const ruint<__RECINT_LIMB_SIZE+1>& r) : Value(r.Value) {}
        ruint(const double b) : Value(static_cast<__uint128_t>(b)) {}
        template <typename T, __RECINT_IS_UNSIGNED(T, int) = 0> ruint(const T b) : Value(__uint128_t(b)) {}
        template <typename T, __RECINT_IS_SIGNED(T, int) = 0> ruint(const T b) : Value(__uint128_t(b)) {}
        template <typename T, __RECINT_IS_NOT_FUNDAMENTAL(T, int) = 0> ruint(const T& b)
        { *this = b.operator ruint<__RECINT_LIMB_SIZE+1>(); } // Fix for Givaro::Integer
        ruint(const char* s);

        // type_string
        static const std::string type_string () {
            return "RecInt::ruint<" + std::to_string(__RECINT_LIMB_SIZE+1)+ ">";
        }

        // Cast
        operator float() const { return (float)(Value); }
        operator double() const { return (double)(Value); }
        template <typename T, __RECINT_IS_UNSIGNED(T, int) = 0> operator T() const { return T(Value); }
        template <typename T, __RECINT_IS_SIGNED(T, int) = 0> operator T() const
        { T ret = T(Value); if (ret < 0) return T(ret & __RECINT_TYPENOTMAXPOWTWO(T)); else return ret; }

        // Const reverse iterator
        class cr_iterator {
        public:
            cr_iterator() : _index(-1), _ptr(NULL) {}
            cr_iterator(int index, const ruint<__RECINT_LIMB_SIZE+1> *ptr) : _index(index), _ptr(ptr) {}
            cr_iterator& operator++() { --_index; return *this; }
            cr_iterator  operator++(int) { cr_iterator i(*this); --_index; return i; }
            bool operator==(cr_iterator rhs) { return (rhs._index == _index); }
            bool operator!=(cr_iterator rhs) { return (rhs._index != _index); }
            const limb* operator->() { return (_index)? &_ptr->High.Value : &_ptr->Low.Value; }
            const limb& operator*() { return (_index)? _ptr->High.Value : _ptr->Low.Value; }
        private:
            int _index;
            const ruint<__RECINT_LIMB_SIZE+1> *_ptr;
        };

        // Const reverse iterator functions
        cr_iterator rbegin() const { return cr_iterator(NBLIMB<__RECINT_LIMB_SIZE+1>::value-1, this); }
        cr_iterator rend() const { return cr_iterator(); }
        UDItype size() { return NBLIMB<__RECINT_LIMB_SIZE+1>::value; }
        static ruint<__RECINT_LIMB_SIZE+1> maxCardinality();
        static ruint<__RECINT_LIMB_SIZE+1> maxModulus();
    };
#endif

    using ruint64 =  ruint<6>;
    using ruint128 = ruint<7>;
    using ruint256 = ruint<8>;
    using ruint512 = ruint<9>;
}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
