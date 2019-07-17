/* rint/arith.h - Arithmetic functions for rint

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


#ifndef RINT_ARITH_SUB_H
#define RINT_ARITH_SUB_H

#include "rrint.h"
#include "rusub.h"

// --------------------------------------------------------------
// ----------------------- DEFINTIONS ---------------------------

namespace RecInt
{
    template <size_t K> rint<K>& operator--(rint<K>&);
    template <size_t K> rint<K>  operator--(rint<K>&, int);

    template <size_t K> rint<K>& operator-=(rint<K>&, const rint<K>&);
    template <size_t K, typename T> __RECINT_IS_UNSIGNED(T, rint<K>&) operator-=(rint<K>&, const T&);
    template <size_t K, typename T> __RECINT_IS_SIGNED(T, rint<K>&)   operator-=(rint<K>&, const T&);

    template <size_t K> rint<K> operator-(const rint<K>&, const rint<K>&);
    template <size_t K, typename T> __RECINT_IS_UNSIGNED(T, rint<K>) operator-(const rint<K>&, const T&);
    template <size_t K, typename T> __RECINT_IS_UNSIGNED(T, rint<K>) operator-(const T&, const rint<K>&);
    template <size_t K, typename T> __RECINT_IS_SIGNED(T, rint<K>)   operator-(const rint<K>&, const T&);
    template <size_t K, typename T> __RECINT_IS_SIGNED(T, rint<K>)   operator-(const T&, const rint<K>&);

    // a = b - c    or  a -= c  (a, b, c are rint)
    // r is the borrow
    template <size_t K> void sub(bool& r, rint<K>& a, const rint<K>& b, const rint<K>& c);
    template <size_t K> void sub(bool& r, rint<K>& a, const rint<K>& c);
    template <size_t K> rint<K>& sub(rint<K>& a, const rint<K>& b, const rint<K>& c);
    template <size_t K> rint<K>& sub(rint<K>& a, const rint<K>& c);

    // a = b - c    or  a -= c  (a, b are rint and c is integer)
    // r is the borrow
    template <size_t K, typename T> __RECINT_IS_ARITH(T, void) sub(bool& r, rint<K>& a, const rint<K>& b, const T& c);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, void) sub(bool& r, rint<K>& a, const T& c);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, void) sub(rint<K>& a, const rint<K>& b, const T& c);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, void) sub(rint<K>& a, const T& c);

    // a = b - 1    or  a -= 1  (a, b are rint)
    // r is the borrow
    template <size_t K> void sub_1(bool& r, rint<K>& a, const rint<K>& b);
    template <size_t K> void sub_1(bool& r, rint<K>& a);
    template <size_t K> void sub_1(rint<K>& a, const rint<K>& b);
    template <size_t K> void sub_1(rint<K>& a);
}


// --------------------------------------------------------------
// ------------------------ Operators ---------------------------

namespace RecInt
{

    // Operator --
    template <size_t K>
    inline rint<K>& operator--(rint<K>& a) {
        sub_1(a);
        return a;
    }
    template <size_t K>
    inline rint<K> operator--(rint<K>& a, int) {
        rint<K> temp(a);
        sub_1(a);
        return temp;
    }

    // Operator -=
    template <size_t K>
    inline rint<K>& operator-=(rint<K>& a, const rint<K>& b) {
        sub(a, b);
        return a;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_UNSIGNED(T, rint<K>&) operator-=(rint<K>& a, const T& b) {
        sub(a, b);
        return a;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_SIGNED(T, rint<K>&) operator-=(rint<K>& a, const T& b) {
        sub(a, b);
        return a;
    }

    // Operator -
    template <size_t K>
    inline rint<K> operator-(const rint<K>& b, const rint<K>& c) {
        rint<K> a;
        sub(a, b, c);
        return a;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_UNSIGNED(T, rint<K>) operator-(const rint<K>& b, const T& c) {
        rint<K> a;
        sub(a, b, c);
        return a;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_UNSIGNED(T, rint<K>) operator-(const T& c, const rint<K>& b) {
        rint<K> a;
        sub(a, b, c);
        return -a;
    }
}


// --------------------------------------------------------------
// ----------------------- Substraction -------------------------

namespace RecInt
{
    // Substract with rint
    // a = b - c    (r stores the borrow)
    template <size_t K>
    inline void sub(bool& r, rint<K>& a, const rint<K>& b, const rint<K>& c) {
        sub(r, a.Value, b.Value, c.Value);
    }

    // a -= b   (r stores the borrow)
    template <size_t K>
    inline void sub(bool& r, rint<K>& a, const rint<K>& b) {
        sub(r, a.Value, b.Value);
    }

    // a = b - c    (the borrow is lost)
    template <size_t K>
    inline rint<K>& sub(rint<K>& a, const rint<K>& b, const rint<K>& c) {
        sub(a.Value, b.Value, c.Value);
        return a;
    }

    // a -= b    (the borrow is lost)
    template <size_t K>
    inline rint<K>& sub(rint<K>& a, const rint<K>& b) {
        sub(a.Value, b.Value);
        return a;
    }


    // Substract with integer
    // a = b - c    (r stores the borrow)
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, void) sub(bool& r, rint<K>& a, const rint<K>& b, const T& c) {
        sub(r, a.Value, b.Value, c);
    }

    // a -= b    (r stores the borrow)
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, void) sub(bool& r, rint<K>& a, const T& b) {
        sub(r, a.Value, b);
    }

    // a = b - c    (the borrow is lost)
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, void) sub(rint<K>& a, const rint<K>& b, const T& c) {
        sub(a.Value, b.Value, c);
    }

    // a -= b    (the borrow is lost)
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, void) sub(rint<K>& a, const T& b) {
        sub(a.Value, b);
    }


    // Decrement
    // a = b - 1    (r stores the borrow)
    template <size_t K>
    inline void sub_1(bool& r, rint<K>& a, const rint<K>& b) {
        sub_1(r, a.Value, b.Value);
    }

    // a -= 1    (r stores the borrow)
    template <size_t K>
    inline void sub_1(bool& r, rint<K>& a) {
        sub_1(r, a.Value);
    }

    // a = b - 1    (the borrow is lost)
    template <size_t K>
    inline void sub_1(rint<K>& a, const rint<K>& b) {
        sub_1(a.Value, b.Value);
    }

    // a -= 1    (the borrow is lost)
    template <size_t K>
    inline void sub_1(rint<K>& a) {
        sub_1(a.Value);
    }
}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
