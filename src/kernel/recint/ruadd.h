/* ruint/arith/add.h - Addition arithmetic functions for ruint

   Copyright Université Grenoble Alpes
Contributors :
Alexis BREUST (alexis.breust@gmail.com 2014)
Christophe CHABOT (christophechabotcc@gmail.com 2011)
Jean-Guillaume DUMAS

Time-stamp: <20 Jun 12 10:28:29 Jean-Guillaume.Dumas@imag.fr>

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


#ifndef RUINT_ARITH_ADD_H
#define RUINT_ARITH_ADD_H

#include "ruruint.h"
#include "rucmp.h"

// --------------------------------------------------------------
// ----------------------- DEFINTIONS ---------------------------

namespace RecInt
{
    template <size_t K> ruint<K>& operator++(ruint<K>&);
    template <size_t K> ruint<K>  operator++(ruint<K>&, int);

    template <size_t K> ruint<K>& operator+=(ruint<K>&, const ruint<K>&);
    template <size_t K, typename T> __RECINT_IS_UNSIGNED(T, ruint<K>&) operator+=(ruint<K>&, const T&);
    template <size_t K, typename T> __RECINT_IS_SIGNED(T, ruint<K>&)   operator+=(ruint<K>&, const T&);

    template <size_t K> ruint<K> operator+(const ruint<K>&, const ruint<K>&);
    template <size_t K, typename T> __RECINT_IS_UNSIGNED(T, ruint<K>) operator+(const ruint<K>&, const T&);
    template <size_t K, typename T> __RECINT_IS_UNSIGNED(T, ruint<K>) operator+(const T&, const ruint<K>&);
    template <size_t K, typename T> __RECINT_IS_SIGNED(T, ruint<K>)   operator+(const ruint<K>&, const T&);
    template <size_t K, typename T> __RECINT_IS_SIGNED(T, ruint<K>)   operator+(const T&, const ruint<K>&);

    // a = b + c    or  a += c  (a, b, c are ruint)
    // r is the carry
    template <size_t K> void add(bool& r, ruint<K>& a, const ruint<K>& b, const ruint<K>& c);
    template <size_t K> void add(bool& r, ruint<K>& a, const ruint<K>& c);
    template <size_t K> void add(ruint<K>& a, const ruint<K>& b, const ruint<K>& c);
    template <size_t K> void add(ruint<K>& a, const ruint<K>& c);

    // a = b + c    or  a += c  (a, b are ruint and c is an integer)
    // r is the carry
    template <size_t K, typename T> __RECINT_IS_ARITH(T, void) add(bool& r, ruint<K>& a, const ruint<K>& b, const T& c);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, void) add(bool& r, ruint<K>& a, const T& c);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, void) add(ruint<K>& a, const ruint<K>& b, const T& c);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, void) add(ruint<K>& a, const T& c);

    // a = b + 1    or  a += 1  (a, b are ruint)
    // r is the carry
    template <size_t K> void add_1(bool& r, ruint<K>& a, const ruint<K>& b);
    template <size_t K> void add_1(bool& r, ruint<K>& a);
    template <size_t K> void add_1(ruint<K>& a, const ruint<K>& b);
    template <size_t K> void add_1(ruint<K>& a);

    // Mostly internal use
    template <size_t K> void add_wc(bool& r, ruint<K>& a, const ruint<K>& b, const ruint<K>& c, const bool& cy);
    template <size_t K> void add_wc(bool& r, ruint<K>& a, const ruint<K>& c, const bool& cy);
    template <size_t K> void add_wc(ruint<K>& a, const ruint<K>& b, const ruint<K>& c, const bool& cy);
    template <size_t K> void add_wc(ruint<K>& a, const ruint<K>& c, const bool& cy);
}


// --------------------------------------------------------------
// ------------------------ Operators ---------------------------

namespace RecInt
{
    // Operator ++
    template <size_t K>
    inline ruint<K>& operator++(ruint<K>& a) {
        add_1(a);
        return a;
    }
    template <size_t K>
    inline ruint<K> operator++(ruint<K>& a, int) {
        ruint<K> temp(a);
        add_1(a);
        return temp;
    }


    // Operator +=
    template <size_t K>
    inline ruint<K>& operator+=(ruint<K>& a, const ruint<K>& b) {
        add(a, b);
        return a;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_UNSIGNED(T, ruint<K>&) operator+=(ruint<K>& a, const T& b) {
        add(a, b);
        return a;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_SIGNED(T, ruint<K>&) operator+=(ruint<K>& a, const T& b) {
        if (b < 0) sub(a, -b);
        else add(a, b);
        return a;
    }


    // Operator +
    template <size_t K>
    inline ruint<K> operator+(const ruint<K>& b, const ruint<K>& c) {
        ruint<K> a;
        add(a, b, c);
        return a;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_UNSIGNED(T, ruint<K>) operator+(const ruint<K>& b, const T& c) {
        ruint<K> a;
        add(a, b, c);
        return a;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_UNSIGNED(T, ruint<K>) operator+(const T& c, const ruint<K>& b) {
        ruint<K> a;
        add(a, b, c);
        return a;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_SIGNED(T, ruint<K>) operator+(const ruint<K>& b, const T& c) {
        ruint<K> a;
        if (c < 0) sub(a, b, -c);
        else add(a, b, c);
        return a;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_SIGNED(T, ruint<K>) operator+(const T& c, const ruint<K>& b) {
        ruint<K> a;
        if (c < 0) sub(a, b, -c);
        else add(a, b, c);
        return a;
    }

}


// --------------------------------------------------------------
// ------------------------- Addition ---------------------------

namespace RecInt
{
    // Add with ruint
    // a = b + c    (r stores the carry)
    template <size_t K>
    inline void add(bool& r, ruint<K>& a, const ruint<K>& b, const ruint<K>& c) {
        bool rl;
        add(rl, a.Low, b.Low, c.Low);
        add_wc(r, a.High, b.High, c.High, rl);
    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline void add(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const ruint<__RECINT_LIMB_SIZE+1>& c) {
        auto bp(b.Value);
        a.Value = b.Value + c.Value;
        r = (a.Value < bp);
    }
#else
    template<>
    inline void add(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const ruint<__RECINT_LIMB_SIZE+1>& c) {
        auto bp(b);
        recint_add_ssaaaa(a.High.Value, a.Low.Value, b.High.Value, b.Low.Value, c.High.Value, c.Low.Value);
        r = (a < bp);
    }
#endif
    template<>
    inline void add(bool& r, ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b, const ruint<__RECINT_LIMB_SIZE>& c) {
        auto bp(b.Value);
        a.Value = b.Value + c.Value;
        r = (a.Value < bp);
    }

    // a += b    (r stores the carry)
    template <size_t K>
    inline void add(bool& r, ruint<K>& a, const ruint<K>& b) {
        bool rl;
        add(rl, a.Low, b.Low);
        add_wc(r, a.High, b.High, rl);
    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline void add(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b) {
        auto bp(b.Value);
        a.Value += b.Value;
        r = (a.Value < bp);
    }
#else
    template<>
    inline void add(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b) {
        auto bp(b);
        recint_add_ssaaaa(a.High.Value, a.Low.Value, a.High.Value, a.Low.Value, b.High.Value, b.Low.Value);
        r = (a < bp);
    }
#endif
    template<>
    inline void add(bool& r, ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b) {
        auto bp(b.Value);
        a.Value += b.Value;
        r = (a.Value < bp);
    }

    // a = b + c    (the carry is lost)
    template <size_t K>
    inline void add(ruint<K>& a, const ruint<K>& b, const ruint<K>& c) {
        bool rl;
        add(rl, a.Low, b.Low, c.Low);
        add_wc(a.High, b.High, c.High, rl);
    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline void add(ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const ruint<__RECINT_LIMB_SIZE+1>& c) {
        a.Value = b.Value + c.Value;
    }
#else
    template<>
    inline void add(ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const ruint<__RECINT_LIMB_SIZE+1>& c) {
        recint_add_ssaaaa(a.High.Value, a.Low.Value, b.High.Value, b.Low.Value, c.High.Value, c.Low.Value);
    }
#endif
    template<>
    inline void add(ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b, const ruint<__RECINT_LIMB_SIZE>& c) {
        a.Value = b.Value + c.Value;
    }

    // a += b    (the carry is lost)
    template <size_t K>
    inline void add(ruint<K>& a, const ruint<K>& b) {
        bool rl;
        add(rl, a.Low, b.Low);
        add_wc(a.High, b.High, rl);
    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline void add(ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b) {
        a.Value += b.Value;
    }
#else
    template<>
    inline void add(ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b) {
        recint_add_ssaaaa(a.High.Value, a.Low.Value, a.High.Value, a.Low.Value, b.High.Value, b.Low.Value);
    }
#endif
    template<>
    inline void add(ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b) {
        a.Value += b.Value;
    }


    // Add with integer
    // a = b + c    (r stores the carry)
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, void) add(bool& r, ruint<K>& a, const ruint<K>& b, const T& c) {
        bool rl;
        add(rl, a.Low, b.Low, c);
        add(r, a.High, b.High, rl);
    }
    // TODO Use __RECINT_USE_FAST_128 here too
    template <typename T>
    inline __RECINT_IS_ARITH(T, void) add(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const T& c) {
        recint_add_ssaaaa(a.High.Value, a.Low.Value, b.High.Value, b.Low.Value, 0, (UWtype)(c));
        r = (a < c);
    }
    template <typename T>
    inline __RECINT_IS_ARITH(T, void) add(bool& r, ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b, const T& c) {
        a.Value = b.Value + limb(c);
        r = (a.Value < limb(c));
    }

    // a += b    (r stores the carry)
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, void) add(bool& r, ruint<K>& a, const T& b) {
        bool rl;
        add(rl, a.Low, b);
        add(r, a.High, rl);
    }
    template <typename T>
    inline __RECINT_IS_ARITH(T, void) add(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a, const T& b) {
        recint_add_ssaaaa(a.High.Value, a.Low.Value, a.High.Value, a.Low.Value, 0, (UWtype)(b));
        r = (a < b);
    }
    template <typename T>
    inline __RECINT_IS_ARITH(T, void) add(bool& r, ruint<__RECINT_LIMB_SIZE>& a, const T& b) {
        a.Value += limb(b);
        r = (a.Value < limb(b));
    }

    // a = b + c    (the carry is lost)
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, void) add(ruint<K>& a, const ruint<K>& b, const T& c) {
        bool r;
        add(r, a, b, c);
    }
    // a += b    (the carry is lost)
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, void) add(ruint<K>& a, const T& b) {
        bool r;
        add(r, a, b);
    }

    // Increment
    // a = b + 1    (r stores the carry)
    template <size_t K>
    inline void add_1(bool& r, ruint<K>& a, const ruint<K>& b) {
        bool rl;
        add_1(rl, a.Low, b.Low);
        add_wc(r, a.High, b.High, rl);
    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline void add_1(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b) {
        a.Value = b.Value + 1;
        r = (a == 0);
    }
#else
    template<>
    inline void add_1(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b) {
        recint_add_ssaaaa(a.High.Value, a.Low.Value, b.High.Value, b.Low.Value, 0, 1);
        r = (a == 0);
    }
#endif
    template<>
    inline void add_1(bool& r, ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b) {
        a.Value = b.Value + 1;
        r = (a.Value == 0);
    }

    // a += 1    (r stores the carry)
    template <size_t K>
    inline void add_1(bool& r, ruint<K>& a) {
        bool rl;
        add_1(rl, a.Low);
        add(r, a.High, rl);
    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline void add_1(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a) {
        ++a.Value;
        r = (a == 0);
    }
#else
    template<>
    inline void add_1(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a) {
        recint_add_ssaaaa(a.High.Value, a.Low.Value, a.High.Value, a.Low.Value, 0, 1);
        r = (a == 0);
    }
#endif
    template<>
    inline void add_1(bool& r, ruint<__RECINT_LIMB_SIZE>& a) {
        a.Value = a.Value + 1;
        r = (a.Value == 0);
    }

    // a = b + 1    (the carry is lost)
    template <size_t K>
    inline void add_1(ruint<K>& a, const ruint<K>& b) {
        bool rl;
        add_1(rl, a.Low, b.Low);
        add_wc(a.High, b.High, rl);
    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline void add_1(ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b) {
        a.Value = b.Value + 1;
    }
#else
    template<>
    inline void add_1(ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b) {
        recint_add_ssaaaa(a.High.Value, a.Low.Value, b.High.Value, b.Low.Value, 0, 1);
    }
#endif
    template<>
    inline void add_1(ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b) {
        a.Value = b.Value + 1;
    }

    // a += 1    (the carry is lost)
    template <size_t K>
    inline void add_1(ruint<K>& a) {
        bool rl;
        add_1(rl, a.Low);
        add(a.High, rl);
    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline void add_1(ruint<__RECINT_LIMB_SIZE+1>& a) {
        ++a.Value;
    }
#else
    template<>
    inline void add_1(ruint<__RECINT_LIMB_SIZE+1>& a) {
        recint_add_ssaaaa(a.High.Value, a.Low.Value, a.High.Value, a.Low.Value, 0, 1);
    }
#endif
    template<>
    inline void add_1(ruint<__RECINT_LIMB_SIZE>& a) {
        a.Value = a.Value + 1;
    }


    // Add with carry
    // a = b + c + cy   (r stores the carry, cy is 0 or 1)
    template <size_t K>
    inline void add_wc(bool& r, ruint<K>& a, const ruint<K>& b, const ruint<K>& c, const bool& cy) {
        bool ret;
        add_wc(ret, a.Low, b.Low, c.Low, cy);
        add_wc(r, a.High, b.High, c.High, ret);
    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline void add_wc(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const ruint<__RECINT_LIMB_SIZE+1>& c, const bool& cy) {
        auto bp(b.Value);
        a.Value = b.Value + c.Value;
        if (cy) {
            ++a.Value;
            r = (a.Value <= bp);
        } else {
            r = (a.Value < bp);
        }
    }
#else
    template<>
    inline void add_wc(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const ruint<__RECINT_LIMB_SIZE+1>& c, const bool& cy) {
        auto bp(b);
        recint_add_ssaaaa(a.High.Value, a.Low.Value, b.High.Value, b.Low.Value, c.High.Value, c.Low.Value);
        if (cy) {
            recint_add_ssaaaa(a.High.Value, a.Low.Value, a.High.Value, a.Low.Value, 0, 1);
            r = (a <= bp);
        } else {
            r = (a < bp);
        }
    }
#endif
    template<>
    inline void add_wc(bool& r, ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b, const ruint<__RECINT_LIMB_SIZE>& c, const bool& cy) {
        auto bp(b.Value);
        a.Value = b.Value + c.Value;
        if (cy) {
            ++a.Value;
            r = (a.Value <= bp);
        } else {
            r = (a.Value < bp);
        }
    }

    // a += b + cy   (r stores the carry, cy is 0 or 1)
    template <size_t K>
    inline void add_wc(bool& r, ruint<K>& a, const ruint<K>& b, const bool& cy) {
        bool ret;
        add_wc(ret, a.Low, b.Low, cy);
        add_wc(r, a.High, b.High, ret);
    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline void add_wc(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const bool& cy) {
        auto bp(b.Value);
        a.Value += b.Value;
        if (cy) {
            ++a.Value;
            r = (a.Value <= bp);
        } else {
            r = (a.Value < bp);
        }
    }
#else
    template<>
    inline void add_wc(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const bool& cy) {
        auto bp(b);
        recint_add_ssaaaa(a.High.Value, a.Low.Value, a.High.Value, a.Low.Value, b.High.Value, b.Low.Value);
        if (cy) {
            recint_add_ssaaaa(a.High.Value, a.Low.Value, a.High.Value, a.Low.Value, 0, 1);
            r = (a <= bp);
        } else {
            r = (a < bp);
        }
    }
#endif
    template<>
    inline void add_wc(bool& r, ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b, const bool& cy) {
        auto bp(b.Value);
        a.Value += b.Value;
        if (cy) {
            ++a.Value;
            r = (a.Value <= bp);
        } else {
            r = (a.Value < bp);
        }
    }

    // a = b + c + cy   (the carry is lost, cy is 0 or 1)
    template <size_t K>
    inline void add_wc(ruint<K>& a, const ruint<K>& b, const ruint<K>& c, const bool& cy) {
        bool ret;
        add_wc(ret, a.Low, b.Low, c.Low, cy);
        add_wc(a.High, b.High, c.High, ret);

    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline void add_wc(ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const ruint<__RECINT_LIMB_SIZE+1>& c, const bool& cy) {
        a.Value = b.Value + c.Value + cy;
    }
#else
    template<>
    inline void add_wc(ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const ruint<__RECINT_LIMB_SIZE+1>& c, const bool& cy) {
        recint_add_ssaaaa(a.High.Value, a.Low.Value, b.High.Value, b.Low.Value, c.High.Value, c.Low.Value);
        if (cy) recint_add_ssaaaa(a.High.Value, a.Low.Value, a.High.Value, a.Low.Value, 0, 1);
    }
#endif
    template<>
    inline void add_wc(ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b, const ruint<__RECINT_LIMB_SIZE>& c, const bool& cy) {
        a.Value = b.Value + c.Value + cy;
    }

    // a += b + cy   (the carry is lost, cy is 0 or 1)
    template <size_t K>
    inline void add_wc(ruint<K>& a, const ruint<K>& b, const bool& cy) {
        bool ret;
        add_wc(ret, a.Low, b.Low, cy);
        add_wc(a.High, b.High, ret);
    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline void add_wc(ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const bool& cy) {
        a.Value += b.Value + cy;
    }
#else
    template<>
    inline void add_wc(ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const bool& cy) {
        recint_add_ssaaaa(a.High.Value, a.Low.Value, a.High.Value, a.Low.Value, b.High.Value, b.Low.Value);
        recint_add_ssaaaa(a.High.Value, a.Low.Value, a.High.Value, a.Low.Value, 0, cy);
    }
#endif
    template<>
    inline void add_wc(ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b, const bool& cy) {
        a.Value += b.Value + cy;
    }
}


#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
