/* ruint/fiddling.h - Bits manipulation for ruint

   Copyright Université Grenoble Alpes
Contributors :
Alexis BREUST (alexis.breust@gmail.com 2014)
Jean-Guillaume DUMAS

Time-stamp: <30 Nov 21 19:05:05 Jean-Guillaume.Dumas@imag.fr>

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


#ifndef RUINT_FIDDLING_H
#define RUINT_FIDDLING_H

#include "ruadd.h" /* operator ++ */

// --------------------------------------------------------------
// ----------------------- DEFINTIONS ---------------------------

namespace RecInt
{
    template <size_t K> ruint<K> operator~(const ruint<K>& c);
    template <size_t K> ruint<K> operator-(const ruint<K>& c);
    template <size_t K> ruint<K>& neg(ruint<K>&, const ruint<K>&);
    template <size_t K> ruint<K>& neg(ruint<K>&);

    template <size_t K> ruint<K>& operator|=(ruint<K>& b, const ruint<K>& c);
    template <size_t K> ruint<K>& operator^=(ruint<K>& b, const ruint<K>& c);
    template <size_t K> ruint<K>& operator&=(ruint<K>& b, const ruint<K>& c);

    template <size_t K, typename T> __RECINT_IS_ARITH(T, ruint<K>&) operator|=(ruint<K>& b, const T& c);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, ruint<K>&) operator^=(ruint<K>& b, const T& c);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, ruint<K>&) operator&=(ruint<K>& b, const T& c);

    template <size_t K> ruint<K> operator|(const ruint<K>& b, const ruint<K>& c);
    template <size_t K> ruint<K> operator^(const ruint<K>& b, const ruint<K>& c);
    template <size_t K> ruint<K> operator&(const ruint<K>& b, const ruint<K>& c);

    template <size_t K, typename T> __RECINT_IS_ARITH(T, ruint<K>) operator|(const ruint<K>& b, const T& c);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, ruint<K>) operator^(const ruint<K>& b, const T& c);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, T) operator&(const ruint<K>& b, const T& c);

    // a = 2^(2^K-1)
    template<size_t K> ruint<K>& max_pow_two(ruint<K>& a);

    // returns a's highest or lowest bit
    template<size_t K> bool highest_bit(const ruint<K>& a);
    template<size_t K> bool lowest_bit(const ruint<K>& a);

    // set a's highest or lowest bit
    template<size_t K> void set_highest_bit(ruint<K>& a);
    template<size_t K> void set_lowest_bit(ruint<K>& a);

    // set a's highest or lowest word
    template<size_t K> void set_highest_word(ruint<K>& a, const limb w);
    template<size_t K> void set_lowest_word(ruint<K>& a, const limb w);
}


// --------------------------------------------------------------
// ----------------------- Operators ----------------------------

namespace RecInt
{
    // Operator ~
    template <size_t K> inline ruint<K> operator~(const ruint<K>& c) {
        ruint<K> b;
        b.High = ~c.High;
        b.Low = ~c.Low;
        return b;
    }
    template <> ruint<__RECINT_LIMB_SIZE> inline operator~(const ruint<__RECINT_LIMB_SIZE>& c) {
        ruint<__RECINT_LIMB_SIZE> b;
        b.Value = ~c.Value;
        return b;
    }

    // Operator neg
    template <size_t K> inline ruint<K>& neg(ruint<K>& r, const ruint<K>& c) {
        r.High = ~c.High;
        r.Low = ~c.Low;
        return ++r;
    }
    template <> inline ruint<__RECINT_LIMB_SIZE>& neg(ruint<__RECINT_LIMB_SIZE>& r, const ruint<__RECINT_LIMB_SIZE>& c) {
        r.Value = ~c.Value;
        return ++r;
    }

    template <size_t K> inline ruint<K>& neg(ruint<K>& r) {
        return neg(r,r);
    }

    // Operator - unary
    template <size_t K> inline ruint<K> operator-(const ruint<K>& c) {
        ruint<K> b; return neg(b,c);
    }

    // Operator |=
    template <size_t K> inline ruint<K>& operator|=(ruint<K>& b, const ruint<K>& c) {
        b.High |= c.High;
        b.Low |= c.Low;
        return b;
    }
    template <> inline ruint<__RECINT_LIMB_SIZE>& operator|=(ruint<__RECINT_LIMB_SIZE>& b, const ruint<__RECINT_LIMB_SIZE>& c) {
        b.Value |= c.Value;
        return b;
    }
    template <size_t K, typename T> inline __RECINT_IS_ARITH(T, ruint<K>&) operator|=(ruint<K>& b, const T& c) {
        b.Low |= c;
        return b;
    }
    template <typename T> inline __RECINT_IS_ARITH(T, ruint<__RECINT_LIMB_SIZE>&) operator|=(ruint<__RECINT_LIMB_SIZE>& b, const T& c) {
        b.Value |= limb(c);
        return b;
    }

    // Operator ^=
    template <size_t K> inline ruint<K>& operator^=(ruint<K>& b, const ruint<K>& c) {
        b.High ^= c.High;
        b.Low ^= c.Low;
        return b;
    }
    template <> inline ruint<__RECINT_LIMB_SIZE>& operator^=(ruint<__RECINT_LIMB_SIZE>& b, const ruint<__RECINT_LIMB_SIZE>& c) {
        b.Value ^= c.Value;
        return b;
    }
    template <size_t K, typename T> inline __RECINT_IS_ARITH(T, ruint<K>&) operator^=(ruint<K>& b, const T& c) {
        b.Low ^= c;
        return b;
    }
    template <typename T> inline __RECINT_IS_ARITH(T, ruint<__RECINT_LIMB_SIZE>&) operator^=(ruint<__RECINT_LIMB_SIZE>& b, const T& c) {
        b.Value ^= limb(c);
        return b;
    }

    // Operator &=
    template <size_t K> inline ruint<K>& operator&=(ruint<K>& b, const ruint<K>& c) {
        b.High &= c.High;
        b.Low &= c.Low;
        return b;
    }
    template <> inline ruint<__RECINT_LIMB_SIZE>& operator&=(ruint<__RECINT_LIMB_SIZE>& b, const ruint<__RECINT_LIMB_SIZE>& c) {
        b.Value &= c.Value;
        return b;
    }
    template <size_t K, typename T> inline __RECINT_IS_ARITH(T, ruint<K>&) operator&=(ruint<K>& b, const T& c) {
        reset(b.High);
        b.Low &= c;
        return b;
    }
    template <typename T> inline __RECINT_IS_ARITH(T, ruint<__RECINT_LIMB_SIZE>&) operator&=(ruint<__RECINT_LIMB_SIZE>& b, const T& c) {
        b.Value &= limb(c);
        return b;
    }

    // Operator |
    template <size_t K> inline ruint<K> operator|(const ruint<K>& b, const ruint<K>& c) {
        ruint<K> a(b);
        return (a |= c);
    }
    template <size_t K, typename T> inline __RECINT_IS_ARITH(T, ruint<K>) operator|(const ruint<K>& b, const T& c) {
        ruint<K> a(b);
        return (a |= c);
    }

    // Operator &
    template <size_t K> inline ruint<K> operator&(const ruint<K>& b, const ruint<K>& c) {
        ruint<K> a(b);
        return (a &= c);
    }
    template <size_t K, typename T> inline __RECINT_IS_ARITH(T, T) operator&(const ruint<K>& b, const T& c) {
        ruint<K> a(b);
        return (a &= c);
    }

    // Operator ^
    template <size_t K> inline ruint<K> operator^(const ruint<K>& b, const ruint<K>& c) {
        ruint<K> a(b);
        return (a ^= c);
    }
    template <size_t K, typename T> inline __RECINT_IS_ARITH(T, ruint<K>) operator^(const ruint<K>& b, const T& c) {
        ruint<K> a(b);
        return (a ^= c);
    }
}


// --------------------------------------------------------------
// -------------------- Extra functions -------------------------

namespace RecInt
{
    // a = 2^(2^K-1)
    template<size_t K> inline ruint<K>& max_pow_two(ruint<K>& a) {
        max_pow_two(a.High);
        reset(a.Low);
        return a;
    }
    template<> inline ruint<__RECINT_LIMB_SIZE>& max_pow_two(ruint<__RECINT_LIMB_SIZE>& a) {
        a.Value = __RECINT_MAXPOWTWO;
        return a;
    }

    // returns a's highest bit
    template<size_t K> inline bool highest_bit(const ruint<K>& a) {
        return highest_bit(a.High);
    }
    template<> inline bool highest_bit(const ruint<__RECINT_LIMB_SIZE>& a) {
        return bool(a.Value & __RECINT_MAXPOWTWO);
    }

    // returns a's lowest bit
    template<size_t K> inline bool lowest_bit(const ruint<K>& a) {
        return lowest_bit(a.Low);
    }
    template<> inline bool lowest_bit(const ruint<__RECINT_LIMB_SIZE>& a) {
        return bool(a.Value & (limb)1);
    }

    // set a's highest bit to 1
    template<size_t K> inline void set_highest_bit(ruint<K>& a) {
        set_highest_bit(a.High);
    }
    template<> inline void set_highest_bit(ruint<__RECINT_LIMB_SIZE>& a) {
        a.Value |= __RECINT_MAXPOWTWO;
    }

    // set a's lowest bit to 1
    template<size_t K> inline void set_lowest_bit(ruint<K>& a) {
        set_lowest_bit(a.Low);
    }
    template<> inline void set_lowest_bit(ruint<__RECINT_LIMB_SIZE>& a) {
        a.Value |= 1;
    }

    // set a's highest word
    template<size_t K> inline void set_highest_word(ruint<K>& a, const limb w) {
        set_highest_word(a.High,w);
    }
    template<> inline void set_highest_word(ruint<__RECINT_LIMB_SIZE>& a, const limb w) {
        a.Value = w;
    }

    // set a's lowest word
    template<size_t K> inline void set_lowest_word(ruint<K>& a, const limb w) {
        set_lowest_word(a.Low,w);
    }
    template<> inline void set_lowest_word(ruint<__RECINT_LIMB_SIZE>& a, const limb w) {
        a.Value = w;
    }


}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
