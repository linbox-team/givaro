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


#ifndef RINT_ARITH_ADD_H
#define RINT_ARITH_ADD_H

#include "rrint.h"
#include "ruadd.h"

// --------------------------------------------------------------
// ----------------------- DEFINTIONS ---------------------------

namespace RecInt
{
    template <size_t K> rint<K>& operator++(rint<K>&);
    template <size_t K> rint<K>  operator++(rint<K>&, int);

    template <size_t K> rint<K>& operator+=(rint<K>&, const rint<K>&);
    template <size_t K, typename T> __RECINT_IS_UNSIGNED(T, rint<K>&) operator+=(rint<K>&, const T&);
    template <size_t K, typename T> __RECINT_IS_SIGNED(T, rint<K>&)   operator+=(rint<K>&, const T&);

    template <size_t K> rint<K> operator+(const rint<K>&, const rint<K>&);
    template <size_t K, typename T> __RECINT_IS_UNSIGNED(T, rint<K>) operator+(const rint<K>&, const T&);
    template <size_t K, typename T> __RECINT_IS_UNSIGNED(T, rint<K>) operator+(const T&, const rint<K>&);
    template <size_t K, typename T> __RECINT_IS_SIGNED(T, rint<K>)   operator+(const rint<K>&, const T&);
    template <size_t K, typename T> __RECINT_IS_SIGNED(T, rint<K>)   operator+(const T&, const rint<K>&);

    // a = b + c    or  a += c  (a, b, c are ruint)
    // r is the carry
    template <size_t K> void add(bool& r, rint<K>& a, const rint<K>& b, const rint<K>& c);
    template <size_t K> void add(bool& r, rint<K>& a, const rint<K>& c);
    template <size_t K> void add(rint<K>& a, const rint<K>& b, const rint<K>& c);
    template <size_t K> void add(rint<K>& a, const rint<K>& c);

    // a = b + c    or  a += c  (a, b are ruint and c is an integer)
    // r is the carry
    template <size_t K, typename T> __RECINT_IS_ARITH(T, void) add(bool& r, rint<K>& a, const rint<K>& b, const T& c);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, void) add(bool& r, rint<K>& a, const T& c);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, void) add(rint<K>& a, const rint<K>& b, const T& c);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, void) add(rint<K>& a, const T& c);

    // a = b + 1    or  a += 1  (a, b are ruint)
    // r is the carry
    template <size_t K> void add_1(bool& r, rint<K>& a, const rint<K>& b);
    template <size_t K> void add_1(bool& r, rint<K>& a);
    template <size_t K> void add_1(rint<K>& a, const rint<K>& b);
    template <size_t K> void add_1(rint<K>& a);
}


// --------------------------------------------------------------
// ------------------------ Operators ---------------------------

namespace RecInt
{
    // Operator ++
    template <size_t K>
    inline rint<K>& operator++(rint<K>& a) {
        add_1(a);
        return a;
    }
    template <size_t K>
    inline rint<K> operator++(rint<K>& a, int) {
        rint<K> temp(a);
        add_1(a);
        return temp;
    }


    // Operator +=
    template <size_t K>
    inline rint<K>& operator+=(rint<K>& a, const rint<K>& b) {
        add(a, b);
        return a;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_UNSIGNED(T, rint<K>&) operator+=(rint<K>& a, const T& b) {
        add(a, b);
        return a;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_SIGNED(T, rint<K>&) operator+=(rint<K>& a, const T& b) {
        add(a, b);
        return a;
    }


    // Operator +
    template <size_t K>
    inline rint<K> operator+(const rint<K>& b, const rint<K>& c) {
        rint<K> a;
        add(a, b, c);
        return a;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, rint<K>) operator+(const rint<K>& b, const T& c) {
        rint<K> a;
        add(a, b, c);
        return a;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, rint<K>) operator+(const T& c, const rint<K>& b) {
        rint<K> a;
        add(a, b, c);
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
    inline void add(bool& r, rint<K>& a, const rint<K>& b, const rint<K>& c) {
        add(r, a.Value, b.Value, c.Value);
    }

    // a += b    (r stores the carry)
    template <size_t K>
    inline void add(bool& r, rint<K>& a, const rint<K>& b) {
        add(r, a.Value, b.Value);
    }

    // a = b + c    (the carry is lost)
    template <size_t K>
    inline void add(rint<K>& a, const rint<K>& b, const rint<K>& c) {
        add(a.Value, b.Value, c.Value);
    }

    // a += b    (the carry is lost)
    template <size_t K>
    inline void add(rint<K>& a, const rint<K>& b) {
        add(a.Value, b.Value);
    }


    // Add with integer
    // a = b + c    (r stores the carry)
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, void) add(bool& r, rint<K>& a, const rint<K>& b, const T& c) {
        add(r, a.Value, b.Value, c);
    }

    // a += b    (r stores the carry)
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, void) add(bool& r, rint<K>& a, const T& b) {
        add(r, a.Value, b);
    }

    // a = b + c    (the carry is lost)
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, void) add(rint<K>& a, const rint<K>& b, const T& c) {
        add(a.Value, b.Value, c);
    }

    // a += b    (the carry is lost)
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, void) add(rint<K>& a, const T& b) {
        add(a.Value, b);
    }


    // Increment
    // a = b + 1    (r stores the carry)
    template <size_t K>
    inline void add_1(bool& r, rint<K>& a, const rint<K>& b) {
        add_1(r, a.Value, b.Value);
    }

    // a += 1    (r stores the carry)
    template <size_t K>
    inline void add_1(bool& r, rint<K>& a) {
        add_1(r, a.Value);
    }

    // a = b + 1    (the carry is lost)
    template <size_t K>
    inline void add_1(rint<K>& a, const rint<K>& b) {
        add_1(a.Value, b.Value);
    }

    // a += 1    (the carry is lost)
    template <size_t K>
    inline void add_1(rint<K>& a) {
        add_1(a.Value);
    }
}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
