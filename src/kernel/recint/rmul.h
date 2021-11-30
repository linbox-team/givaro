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


#ifndef RINT_ARITH_MUL_H
#define RINT_ARITH_MUL_H

#include "rrint.h"
#include "rumul.h"

// --------------------------------------------------------------
// ----------------------- DEFINTIONS ---------------------------

namespace RecInt
{
    template <size_t K> rint<K>& operator*=(rint<K>&, const rint<K>&);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, rint<K>&) operator*=(rint<K>&, const T&);

    template <size_t K> rint<K> operator*(const rint<K>&, const rint<K>&);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, rint<K>) operator*(const rint<K>&, const T&);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, rint<K>) operator*(const T&, const rint<K>&);

    // a = (b*c).Low    or a = (a*c).Low
    // The higher part is lost
    template <size_t K> rint<K>& mul(rint<K>& a, const rint<K>& b, const rint<K>& c);
    template <size_t K> rint<K>& mul(rint<K>& a, const rint<K>& c);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, rint<K>&) mul(rint<K>& a, const rint<K>& b, const T& c);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, rint<K>&) mul(rint<K>& a, const T& c);

    template <size_t K> void lmul(rint<K+1>&, const rint<K>&, const rint<K>&);

        // a += b*c
    template <size_t K> void addmul(rint<K>&, const rint<K>&, const rint<K>&);



    // a = b*b
    template <size_t K> void lsquare(rint<K+1>& a, const rint<K>& b);
}


// --------------------------------------------------------------
// ------------------------ Operators ---------------------------

namespace RecInt
{
    // Operator *=
    template <size_t K>
    inline rint<K>& operator*=(rint<K>& a, const rint<K>& b) {
        return mul(a, b);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, rint<K>&) operator*=(rint<K>& a, const T& b) {
        return mul(a, b);
    }

    // Operator *
    template <size_t K>
    inline rint<K> operator*(const rint<K>& b, const rint<K>& c) {
        rint<K> a;
        return mul(a, b, c);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, rint<K>) operator*(const rint<K>& b, const T& c) {
        rint<K> a;
        return mul(a, b, c);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, rint<K>) operator*(const T& c, const rint<K>& b) {
        rint<K> a;
        return mul(a, b, c);
    }
}


// --------------------------------------------------------------
// --------------------- Multiplication -------------------------

namespace RecInt
{
    // al = (b*c).Low
    // The higher part is lost
    // Note: this function is safe, al is correctly computed
    // even if b, c are really al
    template <size_t K>
    inline rint<K>& mul(rint<K>& al, const rint<K>& b, const rint<K>& c) {
        mul(al.Value, b.Value, c.Value);
        return al;
    }

    // al = (al*c).Low
    // The higher part is lost
    // Note: this function is safe, al is correctly computed
    // even if c is really al
    template <size_t K>
    inline rint<K>& mul(rint<K>& al, const rint<K>& c) {
        mul(al.Value, c.Value);
        return al;
    }

    // a = (b*c).Low
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, rint<K>&) mul(rint<K>& a, const rint<K>& b, const T& c) {
        mul(a.Value, b.Value, c);
        return a;
    }

    // a = (a*b).Low
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, rint<K>&) mul(rint<K>& a, const T& b) {
        mul(a.Value, b);
        return a;
    }


        // a += b*c
    template <size_t K> inline void addmul(rint<K>& a, const rint<K>& b, const rint<K>& c) {
        rint<K> tmp; mul(tmp, b, c);
        add(a,tmp);
    }

    template <size_t K> inline void lmul(rint<K+1>& a, const rint<K>& b, const rint<K>& c) {
        if (b.isPositive()) {
            if (c.isPositive()) {
                lmul(a.Value, b.Value, c.Value);
            } else {
                lmul(a.Value, b.Value, (-c).Value);
                neg(a);
            }
        } else {
            if (c.isPositive()) {
                lmul(a.Value, (-b).Value, c.Value);
                neg(a);
            } else {
                lmul(a.Value, (-b).Value, (-c).Value);
            }
        }
    }


}


// --------------------------------------------------------------
// -------------------------- Square -----------------------------

namespace RecInt
{
    // a = b*b
    template <size_t K>
    inline void lsquare(rint<K+1>& a, const rint<K>& b) {
        lsquare(a.Value, b.Value);
    }
}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
