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


#ifndef RMINT_COMMON_ARITH_ADD_H
#define RMINT_COMMON_ARITH_ADD_H

/** NOTE : For this common file, either basic/reduc.h or mg/reduc.h
  has to be pre-included. **/

// --------------------------------------------------------------
// ----------------------- DEFINTIONS ---------------------------

namespace RecInt
{
    template <size_t K, size_t MG> rmint<K, MG>& operator++(rmint<K, MG>&);
    template <size_t K, size_t MG> rmint<K, MG>  operator++(rmint<K, MG>&, int);

    template <size_t K, size_t MG> rmint<K, MG>& operator+=(rmint<K, MG>&, const rmint<K, MG>&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, rmint<K, MG>&) operator+=(rmint<K, MG>&, const T&);

    template <size_t K, size_t MG> rmint<K, MG> operator+(const rmint<K, MG>&, const rmint<K, MG>&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, rmint<K, MG>) operator+(const rmint<K, MG>&, const T&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, rmint<K, MG>) operator+(const T&, const rmint<K, MG>&);

    template <size_t K, size_t MG> rmint<K, MG>& add(rmint<K, MG>&, const rmint<K, MG>&, const rmint<K, MG>&);
    template <size_t K, size_t MG> rmint<K, MG>& add(rmint<K, MG>&, const rmint<K, MG>&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, rmint<K, MG>&) add(rmint<K, MG>&, const rmint<K, MG>&, const T&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, rmint<K, MG>&) add(rmint<K, MG>&, const T&);
}


// --------------------------------------------------------------
// ------------------------ Operators ---------------------------

namespace RecInt
{
    // Operator ++
    template <size_t K, size_t MG>
    inline rmint<K, MG>& operator++(rmint<K, MG>& a) {
        add(a, 1);
        return a;
    }

    template <size_t K, size_t MG>
    inline rmint<K, MG> operator++(rmint<K, MG>& a, int) {
        rmint<K, MG> temp(a);
        add(a, 1);
        return temp;
    }

    // Operator +=
    template <size_t K, size_t MG>
    inline rmint<K, MG>& operator+=(rmint<K, MG>& a, const rmint<K, MG>& b) {
        return add(a, b);
    }

    template <size_t K, size_t MG, typename T>
    inline __RECINT_IS_ARITH(T, rmint<K, MG>&) operator+=(rmint<K, MG>& a, const T& b) {
        return add(a, b);
    }

    // Operator +
    template <size_t K, size_t MG>
    inline rmint<K, MG> operator+(const rmint<K, MG>& b, const rmint<K, MG>& c) {
        rmint<K, MG> a;
        return add(a, b, c);
    }

    template <size_t K, size_t MG, typename T>
    inline __RECINT_IS_ARITH(T, rmint<K, MG>) operator+(const rmint<K, MG>& b, const T& c) {
        rmint<K, MG> a;
        return add(a, b, c);
    }

    template <size_t K, size_t MG, typename T>
    inline __RECINT_IS_ARITH(T, rmint<K, MG>) operator+(const T& c, const rmint<K, MG>& b) {
        rmint<K, MG> a;
        return add(a, b, c);
    }
}


// --------------------------------------------------------------
// ------------------------- Addition ---------------------------

namespace RecInt
{
    // a = b + c mod a.p
    template <size_t K, size_t MG>
    inline rmint<K, MG>& add(rmint<K, MG>& a, const rmint<K, MG>& b, const rmint<K, MG>& c) {
        bool r;
        add(r, a.Value, b.Value, c.Value);
        if (r || a.Value >= a.p) sub(a.Value, a.p);
        return a;
    }

    // a += b mod a.p
    template <size_t K, size_t MG>
    inline rmint<K, MG>& add(rmint<K, MG>& a, const rmint<K, MG>& b) {
        bool r;
        add(r, a.Value, b.Value);
        if (r || a.Value >= a.p) sub(a.Value, a.p);
        return a;
    }

    // a = b + c mod a.p
    template <size_t K, size_t MG, typename T>
    inline __RECINT_IS_ARITH(T, rmint<K, MG>&) add(rmint<K, MG>& a, const rmint<K, MG>& b, const T& c) {
        bool r;
        rmint<K, MG> cp(c);
        add(r, a.Value, b.Value, cp.Value);
        if (r || a.Value >= a.p) sub(a.Value, a.p);
        return a;
    }

    // a += b mod a.p
    template <size_t K, size_t MG, typename T>
    inline __RECINT_IS_ARITH(T, rmint<K, MG>&) add(rmint<K, MG>& a, const T& b) {
        bool r;
        rmint<K, MG> bp(b);
        add(r, a.Value, bp.Value);
        if (r || a.Value >= a.p) sub(a.Value, a.p);
        return a;
    }
}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
