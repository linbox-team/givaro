/* ruint/arith.h - Arithmetic functions for ruint

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


#ifndef RUINT_ARITH_SUB_H
#define RUINT_ARITH_SUB_H

#include "ruruint.h"
#include "rucmp.h"

// --------------------------------------------------------------
// ----------------------- DEFINTIONS ---------------------------

namespace RecInt
{
    template <size_t K> ruint<K>& operator--(ruint<K>&);
    template <size_t K> ruint<K>  operator--(ruint<K>&, int);

    template <size_t K> ruint<K>& operator-=(ruint<K>&, const ruint<K>&);
    template <size_t K, typename T> __RECINT_IS_UNSIGNED(T, ruint<K>&) operator-=(ruint<K>&, const T&);
    template <size_t K, typename T> __RECINT_IS_SIGNED(T, ruint<K>&)   operator-=(ruint<K>&, const T&);

    template <size_t K> ruint<K> operator-(const ruint<K>&, const ruint<K>&);
    template <size_t K, typename T> __RECINT_IS_UNSIGNED(T, ruint<K>) operator-(const ruint<K>&, const T&);
    template <size_t K, typename T> __RECINT_IS_UNSIGNED(T, ruint<K>) operator-(const T&, const ruint<K>&);
    template <size_t K, typename T> __RECINT_IS_SIGNED(T, ruint<K>)   operator-(const ruint<K>&, const T&);
    template <size_t K, typename T> __RECINT_IS_SIGNED(T, ruint<K>)   operator-(const T&, const ruint<K>&);

    // a = b - c    or  a -= c  (a, b, c are ruint)
    // r is the borrow
    template <size_t K> void sub(bool& r, ruint<K>& a, const ruint<K>& b, const ruint<K>& c);
    template <size_t K> void sub(bool& r, ruint<K>& a, const ruint<K>& c);
    template <size_t K> ruint<K>& sub(ruint<K>& a, const ruint<K>& b, const ruint<K>& c);
    template <size_t K> ruint<K>& sub(ruint<K>& a, const ruint<K>& c);

    // a = b - c    or  a -= c  (a, b are ruint and c is integer)
    // r is the borrow
    template <size_t K, typename T> __RECINT_IS_ARITH(T, void) sub(bool& r, ruint<K>& a, const ruint<K>& b, const T& c);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, void) sub(bool& r, ruint<K>& a, const T& c);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, void) sub(ruint<K>& a, const ruint<K>& b, const T& c);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, void) sub(ruint<K>& a, const T& c);

    // a = b - 1    or  a -= 1  (a, b are ruint)
    // r is the borrow
    template <size_t K> void sub_1(bool& r, ruint<K>& a, const ruint<K>& b);
    template <size_t K> void sub_1(bool& r, ruint<K>& a);
    template <size_t K> void sub_1(ruint<K>& a, const ruint<K>& b);
    template <size_t K> void sub_1(ruint<K>& a);

    // Mostly internal use
    template <size_t K> void sub_wc(bool& r, ruint<K>& a, const ruint<K>& b, const ruint<K>& c, const bool& cy);
    template <size_t K> void sub_wc(bool& r, ruint<K>& a, const ruint<K>& c, const bool& cy);
    template <size_t K> void sub_wc(ruint<K>& a, const ruint<K>& b, const ruint<K>& c, const bool& cy);
    template <size_t K> void sub_wc(ruint<K>& a, const ruint<K>& c, const bool& cy);
}


// --------------------------------------------------------------
// ------------------------ Operators ---------------------------

namespace RecInt
{

    // Operator --
    template <size_t K>
    inline ruint<K>& operator--(ruint<K>& a) {
        sub_1(a);
        return a;
    }
    template <size_t K>
    inline ruint<K> operator--(ruint<K>& a, int) {
        ruint<K> temp(a);
        sub_1(a);
        return temp;
    }

    // Operator -=
    template <size_t K>
    inline ruint<K>& operator-=(ruint<K>& a, const ruint<K>& b) {
        sub(a, b);
        return a;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_UNSIGNED(T, ruint<K>&) operator-=(ruint<K>& a, const T& b) {
        sub(a, b);
        return a;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_SIGNED(T, ruint<K>&) operator-=(ruint<K>& a, const T& b) {
        if (b < 0) add(a, -b);
        else sub(a, b);
        return a;
    }

    // Operator -
    template <size_t K>
    inline ruint<K> operator-(const ruint<K>& b, const ruint<K>& c) {
        ruint<K> a;
        sub(a, b, c);
        return a;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_UNSIGNED(T, ruint<K>) operator-(const ruint<K>& b, const T& c) {
        ruint<K> a;
        sub(a, b, c);
        return a;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_UNSIGNED(T, ruint<K>) operator-(const T& c, const ruint<K>& b) {
        ruint<K> a;
        sub(a, b, c);
        return -a;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_SIGNED(T, ruint<K>) operator-(const ruint<K>& b, const T& c) {
        ruint<K> a;
        if (c < 0) add(a, b, -c);
        else sub(a, b, c);
        return a;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_SIGNED(T, ruint<K>) operator-(const T& c, const ruint<K>& b) {
        ruint<K> a;
        if (c < 0) add(a, b, -c);
        else sub(a, b, c);
        return -a;
    }
}


// --------------------------------------------------------------
// ----------------------- Substraction -------------------------

namespace RecInt
{
    // Substract with ruint
    // a = b - c    (r stores the borrow)
    template <size_t K>
    inline void sub(bool& r, ruint<K>& a, const ruint<K>& b, const ruint<K>& c) {
        bool rl;
        sub(rl, a.Low, b.Low, c.Low);
        sub_wc(r, a.High, b.High, c.High, rl);
    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline void sub(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const ruint<__RECINT_LIMB_SIZE+1>& c) {
        r = (b.Value < c.Value);
        a.Value = b.Value - c.Value;
    }
#else
    template<>
    inline void sub(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const ruint<__RECINT_LIMB_SIZE+1>& c) {
        r = (b < c);
        recint_sub_ddmmss(a.High.Value, a.Low.Value, b.High.Value, b.Low.Value, c.High.Value, c.Low.Value);
    }
#endif
    template<>
    inline void sub(bool& r, ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b, const ruint<__RECINT_LIMB_SIZE>& c) {
        r = (b.Value < c.Value);
        a.Value = b.Value - c.Value;
    }

    // a -= b   (r stores the borrow)
    template <size_t K>
    inline void sub(bool& r, ruint<K>& a, const ruint<K>& b) {
        bool rl;
        sub(rl, a.Low, b.Low);
        sub_wc(r, a.High, b.High, rl);
    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline void sub(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b) {
        r = (a < b);
        a.Value -= b.Value;
    }
#else
    template<>
    inline void sub(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b) {
        r = (a < b);
        recint_sub_ddmmss(a.High.Value, a.Low.Value, a.High.Value, a.Low.Value, b.High.Value, b.Low.Value);
    }
#endif
    template<>
    inline void sub(bool& r, ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b) {
        r = (a < b);
        a.Value -= b.Value;
    }

    // a = b - c    (the borrow is lost)
    template <size_t K>
    inline ruint<K>& sub(ruint<K>& a, const ruint<K>& b, const ruint<K>& c) {
        bool rl;
        sub(rl, a.Low, b.Low, c.Low);
        sub_wc(a.High, b.High, c.High, rl);
        return a;
    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline ruint<__RECINT_LIMB_SIZE+1>& sub(ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const ruint<__RECINT_LIMB_SIZE+1>& c) {
        a.Value = b.Value - c.Value;
        return a;
    }
#else
    template<>
    inline ruint<__RECINT_LIMB_SIZE+1>& sub(ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const ruint<__RECINT_LIMB_SIZE+1>& c) {
        recint_sub_ddmmss(a.High.Value, a.Low.Value, b.High.Value, b.Low.Value, c.High.Value, c.Low.Value);
        return a;
    }
#endif
    template<>
    inline ruint<__RECINT_LIMB_SIZE>& sub(ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b, const ruint<__RECINT_LIMB_SIZE>& c) {
        a.Value = b.Value - c.Value;
        return a;
    }

    // a -= b    (the borrow is lost)
    template <size_t K>
    inline ruint<K>& sub(ruint<K>& a, const ruint<K>& b) {
        bool rl;
        sub(rl, a.Low, b.Low);
        sub_wc(a.High, b.High, rl);
        return a;
    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline ruint<__RECINT_LIMB_SIZE+1>& sub(ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b) {
        a.Value -= b.Value;
        return a;
    }
#else
    template<>
    inline ruint<__RECINT_LIMB_SIZE+1>& sub(ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b) {
        recint_sub_ddmmss(a.High.Value, a.Low.Value, a.High.Value, a.Low.Value, b.High.Value, b.Low.Value);
        return a;

    }
#endif
    template<>
    inline ruint<__RECINT_LIMB_SIZE>& sub(ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b) {
        a.Value -= b.Value;
        return a;
    }


    // Substract with integer
    // a = b - c    (r stores the borrow)
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, void) sub(bool& r, ruint<K>& a, const ruint<K>& b, const T& c) {
        bool rl;
        sub(rl, a.Low, b.Low, c);
        sub(r, a.High, b.High, rl);
    }
    // TODO Use __RECINT_USE_FAST_128 optim
    template <typename T>
    inline __RECINT_IS_ARITH(T, void) sub(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const T& c) {
        r = (b < c);
        recint_sub_ddmmss(a.High.Value, a.Low.Value, b.High.Value, b.Low.Value, 0, limb(c));
    }
    template <typename T>
    inline __RECINT_IS_ARITH(T, void) sub(bool& r, ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b, const T& c) {
        r = (b.Value < limb(c));
        a.Value = b.Value - c;
    }

    // a -= b    (r stores the borrow)
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, void) sub(bool& r, ruint<K>& a, const T& b) {
        bool rl;
        sub(rl, a.Low, b);
        sub(r, a.High, rl);
    }
    template <typename T>
    inline __RECINT_IS_ARITH(T, void) sub(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a, const T& b) {
        r = (a < b);
        recint_sub_ddmmss(a.High.Value, a.Low.Value, a.High.Value, a.Low.Value, 0, limb(b));
    }
    template <typename T>
    inline __RECINT_IS_ARITH(T, void) sub(bool& r, ruint<__RECINT_LIMB_SIZE>& a, const T& b) {
        r = (a.Value < limb(b));
        a.Value = a.Value - b;
    }

    // a = b - c    (the borrow is lost)
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, void) sub(ruint<K>& a, const ruint<K>& b, const T& c) {
        bool rl;
        sub(rl, a.Low, b.Low, c);
        sub(a.High, b.High, rl);
    }
    template <typename T>
    inline __RECINT_IS_ARITH(T, void) sub(ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const T& c) {
        recint_sub_ddmmss(a.High.Value, a.Low.Value, b.High.Value, b.Low.Value, 0, limb(c));
    }
    template <typename T>
    inline __RECINT_IS_ARITH(T, void) sub(ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b, const T& c) {
        a.Value = b.Value - limb(c);
    }

    // a -= b    (the borrow is lost)
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, void) sub(ruint<K>& a, const T& b) {
        bool rl;
        sub(rl, a.Low, b);
        sub(a.High, rl);
    }
    template <typename T>
    inline __RECINT_IS_ARITH(T, void) sub(ruint<__RECINT_LIMB_SIZE+1>& a, const T& b) {
        recint_sub_ddmmss(a.High.Value, a.Low.Value, a.High.Value, a.Low.Value, 0, limb(b));
    }
    template <typename T>
    inline __RECINT_IS_ARITH(T, void) sub(ruint<__RECINT_LIMB_SIZE>& a, const T& b) {
        a.Value = a.Value - b;
    }


    // Decrement
    // a = b - 1    (r stores the borrow)
    template <size_t K>
    inline void sub_1(bool& r, ruint<K>& a, const ruint<K>& b) {
        bool rl;
        sub_1(rl, a.Low, b.Low);
        sub(r, a.High, b.High, rl);
    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline void sub_1(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b) {
        r = (b.Value == 0);
        a.Value = b.Value - 1;
    }
#else
    template<>
    inline void sub_1(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b) {
        r = (b == 0);
        recint_sub_ddmmss(a.High.Value, a.Low.Value, b.High.Value, b.Low.Value, 0, 1);
    }
#endif
    template<>
    inline void sub_1(bool& r, ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b) {
        r = (b.Value == 0);
        a.Value = b.Value - 1;
    }

    // a -= 1    (r stores the borrow)
    template <size_t K>
    inline void sub_1(bool& r, ruint<K>& a) {
        bool rl;
        sub_1(rl, a.Low);
        sub(r, a.High, rl);
    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline void sub_1(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a) {
        r = (a.Value == 0);
        a.Value = a.Value - 1;
    }
#else
    template<>
    inline void sub_1(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a) {
        r = (a == 0);
        recint_sub_ddmmss(a.High.Value, a.Low.Value, a.High.Value, a.Low.Value, 0, 1);
    }
#endif
    template<>
    inline void sub_1(bool& r, ruint<__RECINT_LIMB_SIZE>& a) {
        r = (a.Value == 0);
        a.Value = a.Value - 1;
    }

    // a = b - 1    (the borrow is lost)
    template <size_t K>
    inline void sub_1(ruint<K>& a, const ruint<K>& b) {
        bool rl;
        sub_1(rl, a.Low, b.Low);
        sub(a.High, b.High, rl);
    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline void sub_1(ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b) {
        a.Value = b.Value - 1;
    }
#else
    template<>
    inline void sub_1(ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b) {
        recint_sub_ddmmss(a.High.Value, a.Low.Value, b.High.Value, b.Low.Value, 0, 1);
    }
#endif
    template<>
    inline void sub_1(ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b) {
        a.Value = b.Value - 1;
    }

    // a -= 1    (the borrow is lost)
    template <size_t K>
    inline void sub_1(ruint<K>& a) {
        bool rl;
        sub_1(rl, a.Low);
        sub(a.High, rl);
    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline void sub_1(ruint<__RECINT_LIMB_SIZE+1>& a) {
        --a.Value;
    }
#else
    template<>
    inline void sub_1(ruint<__RECINT_LIMB_SIZE+1>& a) {
        recint_sub_ddmmss(a.High.Value, a.Low.Value, a.High.Value, a.Low.Value, 0, 1);
    }
#endif
    template<>
    inline void sub_1(ruint<__RECINT_LIMB_SIZE>& a) {
        --a.Value;
    }


    // Sub with carry
    // a = b - c - cy   (r stores the borrow, cy is 0 or 1)
    template <size_t K>
    inline void sub_wc(bool& r, ruint<K>& a, const ruint<K>& b, const ruint<K>& c, const bool& cy) {
        bool ret;
        sub_wc(ret, a.Low, b.Low, c.Low, cy);
        sub_wc(r, a.High, b.High, c.High, ret);
    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline void sub_wc(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const ruint<__RECINT_LIMB_SIZE+1>& c, const bool& cy) {
        if (cy) r = (b.Value <= c.Value);
        else r = (b.Value < c.Value);
        a.Value = b.Value - c.Value - cy;
    }
#else
    template<>
    inline void sub_wc(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const ruint<__RECINT_LIMB_SIZE+1>& c, const bool& cy) {
        if (cy) r = (b <= c);
        else r = (b < c);
        recint_sub_ddmmss(a.High.Value, a.Low.Value, b.High.Value, b.Low.Value, c.High.Value, c.Low.Value);
        if (cy) recint_sub_ddmmss(a.High.Value, a.Low.Value, a.High.Value, a.Low.Value, 0, 1);
    }
#endif
    template<>
    inline void sub_wc(bool& r, ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b, const ruint<__RECINT_LIMB_SIZE>& c, const bool& cy) {
        if (cy) r = (b.Value <= c.Value);
        else r = (b.Value < c.Value);
        a.Value = b.Value - c.Value - cy;
    }

    // a -= b + cy   (r stores the borrow, cy is 0 or 1)
    template <size_t K>
    inline void sub_wc(bool& r, ruint<K>& a, const ruint<K>& b, const bool& cy) {
        bool ret;
        sub_wc(ret, a.Low, b.Low, cy);
        sub_wc(r, a.High, b.High, ret);
    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline void sub_wc(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const bool& cy) {
        if (cy) r = (a.Value <= b.Value);
        else r = (a.Value < b.Value);
        a.Value = a.Value - b.Value - cy;
    }
#else
    template<>
    inline void sub_wc(bool& r, ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const bool& cy) {
        if (cy) r = (a <= b);
        else r = (a < b);
        recint_sub_ddmmss(a.High.Value, a.Low.Value, a.High.Value, a.Low.Value, b.High.Value, b.Low.Value);
        if (cy) recint_sub_ddmmss(a.High.Value, a.Low.Value, a.High.Value, a.Low.Value, 0, 1);
    }
#endif
    template<>
    inline void sub_wc(bool& r, ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b, const bool& cy) {
        if (cy) r = (a.Value <= b.Value);
        else r = (a.Value < b.Value);
        a.Value = a.Value - b.Value - cy;
    }

    // a = b - c - cy   (the borrow is lost, cy is 0 or 1)
    template <size_t K>
    inline void sub_wc(ruint<K>& a, const ruint<K>& b, const ruint<K>& c, const bool& cy) {
        bool ret;
        sub_wc(ret, a.Low, b.Low, c.Low, cy);
        sub_wc(a.High, b.High, c.High, ret);
    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline void sub_wc(ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const ruint<__RECINT_LIMB_SIZE+1>& c, const bool& cy) {
        a.Value = b.Value - c.Value - cy;
    }
#else
    template<>
    inline void sub_wc(ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const ruint<__RECINT_LIMB_SIZE+1>& c, const bool& cy) {
        recint_sub_ddmmss(a.High.Value, a.Low.Value, b.High.Value, b.Low.Value, c.High.Value, c.Low.Value);
        recint_sub_ddmmss(a.High.Value, a.Low.Value, a.High.Value, a.Low.Value, 0, cy);
    }
#endif
    template<>
    inline void sub_wc(ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b, const ruint<__RECINT_LIMB_SIZE>& c, const bool& cy) {
        a.Value = b.Value - c.Value - cy;
    }

    // a -= b + cy   (the borrow is lost, cy is 0 or 1)
    template <size_t K>
    inline void sub_wc(ruint<K>& a, const ruint<K>& b, const bool& cy) {
        bool ret;
        sub_wc(ret, a.Low, b.Low, cy);
        sub_wc(a.High, b.High, ret);
    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline void sub_wc(ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const bool& cy) {
        a.Value = a.Value - b.Value - cy;
    }
#else
    template<>
    inline void sub_wc(ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const bool& cy) {
        recint_sub_ddmmss(a.High.Value, a.Low.Value, a.High.Value, a.Low.Value, 0, cy);
        recint_sub_ddmmss(a.High.Value, a.Low.Value, a.High.Value, a.Low.Value, b.High.Value, b.Low.Value);
    }
#endif
    template<>
    inline void sub_wc(ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b, const bool& cy) {
        a.Value = a.Value - b.Value - cy;
    }
}
#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
