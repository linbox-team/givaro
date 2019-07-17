/* rmint/arith.h - Arithmetic functions for rmint

   Copyright Université Grenoble Alpes
Contributors :
Alexis BREUST (alexis.breust@gmail.com 2014)
Christophe CHABOT (christophechabotcc@gmail.com 2011)


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


#ifndef RMINT_COMMON_ARITH_DIV_H
#define RMINT_COMMON_ARITH_DIV_H

/** NOTE : For this common file, either basic/reduc.h or mg/reduc.h
  has to be pre-included. **/

#include "rudiv.h"

// --------------------------------------------------------------
// ----------------------- DEFINTIONS ---------------------------

namespace RecInt
{
    template <size_t K, size_t MG> rmint<K, MG>& operator%=(rmint<K, MG>&, const rmint<K, MG>&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, rmint<K, MG>&) operator%=(rmint<K, MG>&, const T&);

    template <size_t K, size_t MG> rmint<K, MG>& operator/=(rmint<K, MG>&, const rmint<K, MG>&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, rmint<K, MG>&) operator/=(rmint<K, MG>&, const T&);

    template <size_t K, size_t MG> rmint<K, MG> operator/(const rmint<K, MG>&, const rmint<K, MG>&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, rmint<K, MG>) operator/(const rmint<K, MG>&, const T&);

    template <size_t K, size_t MG> rmint<K, MG> operator%(const rmint<K, MG>&, const rmint<K, MG>&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, rmint<K, MG>) operator%(const rmint<K, MG>&, const T&);

    template <size_t K, size_t MG> void div(rmint<K, MG>&, const rmint<K, MG>&, const rmint<K, MG>&);
    template <size_t K, size_t MG> void div(rmint<K, MG>&, const rmint<K, MG>&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, void) div(rmint<K, MG>&, const rmint<K, MG>&, const T&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, void) div(rmint<K, MG>&, const T&);

    template <size_t K, size_t MG> void mod(rmint<K, MG>&, const rmint<K, MG>&, const rmint<K, MG>&);
    template <size_t K, size_t MG> void mod(rmint<K, MG>&, const rmint<K, MG>&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, void) mod(rmint<K, MG>&, const rmint<K, MG>&, const T&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, void) mod(rmint<K, MG>&, const T&);
}


// --------------------------------------------------------------
// ------------------------ Operators ---------------------------

namespace RecInt
{
    // Operator %=
    template <size_t K, size_t MG>
    inline rmint<K, MG>& operator%=(rmint<K, MG>& a, const rmint<K, MG>& b) {
        mod(a, b);
        return a;
    }

    template <size_t K, size_t MG, typename T>
    inline __RECINT_IS_ARITH(T, rmint<K, MG>&) operator%=(rmint<K, MG>& a, const T& b) {
        mod(a, b);
        return a;
    }

    // Operator /=
    template <size_t K, size_t MG>
    inline rmint<K, MG>& operator/=(rmint<K, MG>& a, const rmint<K, MG>& b) {
        div(a, b);
        return a;
    }

    template <size_t K, size_t MG, typename T>
    inline __RECINT_IS_ARITH(T, rmint<K, MG>&) operator/=(rmint<K, MG>& a, const T& b) {
        div(a, b);
        return a;
    }

    // Operator /
    template <size_t K, size_t MG>
    inline rmint<K, MG> operator/(const rmint<K, MG>& b, const rmint<K, MG>& c) {
        rmint<K, MG> a;
        div(a, b, c);
        return a;
    }

    template <size_t K, size_t MG, typename T>
    inline __RECINT_IS_ARITH(T, rmint<K, MG>) operator/(const rmint<K, MG>& b, const T& c) {
        rmint<K, MG> a;
        div(a, b, c);
        return a;
    }

    // Operator %
    template <size_t K, size_t MG>
    inline rmint<K, MG> operator%(const rmint<K, MG>& b, const rmint<K, MG>& c) {
        rmint<K, MG> a;
        mod(a, b, c);
        return a;
    }

    template <size_t K, size_t MG, typename T>
    inline __RECINT_IS_ARITH(T, rmint<K, MG>) operator%(const rmint<K, MG>& b, const T& c) {
        rmint<K, MG> a;
        mod(a, b, c);
        return a;
    }
}


// --------------------------------------------------------------
// ------------------------- Division ---------------------------

namespace RecInt
{
    // a = b * c^(-1) mod a.p
    template <size_t K, size_t MG>
    inline void div(rmint<K, MG>& a, const rmint<K, MG>& b, const rmint<K, MG>& c) {
        rmint<K, MG> ci;
        inv(ci, c);
        if (ci.Value == 0) reset(a.Value);
        else mul(a, b, ci);
    }

    // a *= b^(-1) mod a.p
    template <size_t K, size_t MG>
    inline void div(rmint<K, MG>& a, const rmint<K, MG>& b) {
        div(a, a, b);
    }

    // a = b * c^(-1) mod a.p
    template <size_t K, size_t MG, typename T>
    inline __RECINT_IS_ARITH(T, void) div(rmint<K, MG>& a, const rmint<K, MG>& b, const T& c) {
        rmint<K, MG> cr(c);
        div(a, b, cr);
    }

    // a *= b^(-1) mod a.p
    template <size_t K, size_t MG, typename T>
    inline __RECINT_IS_ARITH(T, void) div(rmint<K, MG>& a, const T& b) {
        rmint<K, MG> br(b);
        div(a, a, br);
    }
}


// --------------------------------------------------------------
// ------------------------- Module ---------------------------

namespace RecInt
{
    // a = b mod c
    template <size_t K, size_t MG>
    inline void mod(rmint<K, MG>& a, const rmint<K, MG>& b, const rmint<K, MG>& c) {
        div_r(a.Value, get_ruint(b), get_ruint(c));
        get_ready(a);
    }

    // a = a mod b
    template <size_t K, size_t MG>
    inline void mod(rmint<K, MG>& a, const rmint<K, MG>& b) {
        div_r(a.Value, get_ruint(a), get_ruint(b));
        get_ready(a);
    }

    // a = b mod c
    template <size_t K, size_t MG, typename T>
    inline __RECINT_IS_ARITH(T, void) mod(rmint<K, MG>& a, const rmint<K, MG>& b, const T& c) {
        rmint<K, MG> cr(c);
        mod(a, b, cr);
    }

    // a = a mod b
    template <size_t K, size_t MG, typename T>
    inline __RECINT_IS_ARITH(T, void) mod(rmint<K, MG>& a, const T& b) {
        rmint<K, MG> br(b);
        mod(a, a, br);
    }
}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
