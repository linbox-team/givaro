/* rmint/tools.h - Modular arithmetic functions for ruint

   Copyright Université Grenoble Alpes
Contributors :
Alexis BREUST (alexis.breust@gmail.com, 2014)
Christophe CHABOT (christophechabotcc@gmail.com, 2011)

This software is a cmputer program whose purpose is to provide an fixed precision arithmetic library.

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
that may mean  that it is cmplicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth cmputer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/


#ifndef RUINT_ARITH_INV_MOD_H
#define RUINT_ARITH_INV_MOD_H

#include "rutools.h" /* mod_n() */

#ifdef __GIVARO_DEBUG
#include "givaro/giverror.h"
#endif
// --------------------------------------------------------------
// ----------------------- DEFINITIONS --------------------------

namespace RecInt
{
    // If c and d are relatively primes:
    // x * c = 1 mod d
    // y * d = 1 mod c
    template <size_t K> void bezout_mod(ruint<K>& x, ruint<K>& y, const ruint<K>& c, const ruint<K>& d);

    // a = b^{-1} mod c (if b is not invertible, a = 0)
    template <size_t K> ruint<K>& inv_mod(ruint<K>& a, const ruint<K>& b, const ruint<K>& c);
}


// --------------------------------------------------------------
// --------------------- Implementation -------------------------

namespace RecInt
{
    // If c and d are relatively prime
    // Compute lastx and lasty such that
    // lastx * c = 1 mod d
    // lasty * d = 1 mod c
    template <size_t K>
    inline void bezout_mod(ruint<K>& lastx, ruint<K>& lasty, const ruint<K>& c, const ruint<K>& d) {
        ruint<K+1> resmul;
        ruint<K> x(0), y(1), a, b, q, r, temp;
        bool ret;

        copy(a, c);
        copy(b, d);

        lastx = 1;
        reset(lasty);

        while (b != 0) {
            // Classical Euclidean algorithm
            div(q, r, a, b);
            copy(a, b);
            copy(b, r);

            // x(i+1) = x(i-1) - q * x(i) mod d
            // temp = q * x mod d
            lmul(resmul, q, x);
            mod_n(temp, resmul, d);

            // temp = -temp mod d
            if (temp != 0) sub(temp, c, temp);

            // temp += lastx mod d
            add(ret, temp, lastx);
            if (ret || cmp(temp, d) >= 0) sub(temp, d);

            // lastx = x and x = temp
            copy(lastx, x);
            copy(x, temp);

            // y(i+1) = y(i-1) - q * y(i) mod c
            // temp = q * y mod c
            lmul(resmul, q, y);
            mod_n(temp, resmul, c);

            // temp = -temp mod c
            if (temp != 0) sub(temp, c, temp);

            // temp += lasty mod c
            add(ret, temp, lasty);
            if (ret || cmp(temp, c) >= 0) sub(temp, c);

            // lasty = y and y = temp
            copy(lasty, y);
            copy(y, temp);
        }
    }

    // a = b^(-1) mod c if b invertible
    template <size_t K>
    inline ruint<K>& inv_mod(ruint<K>& a, const ruint<K>& b, const ruint<K>& c) {
        ruint<K+1> resmul;
        ruint<K> x(0), a2, b2, q, r, temp;
        bool ret;

        copy(a2, b);
        copy(b2, c);
        a = 1;

        while (b2 != 0) {
            div(q, r, a2, b2);
            copy(a2, b2);
            copy(b2, r);

            // temp = q * x mod c
            lmul(resmul, q, x);
            mod_n(temp, resmul, c);
            // temp = -temp mod c
            if (temp != 0) sub(temp, c, temp);
            // temp += a mod c
            add(ret, temp, a);
            if (ret || temp >= c) sub(temp, c);

            copy(a, x);
            copy(x, temp);
        }
#ifdef __GIVARO_DEBUG
        if ( a2 != 1) {
            throw Givaro::GivMathDivZero("*** Error: division by zero, in operator RecInt::inv_mod in ruinvmod.h") ;
        }
#endif
        return a;
    }
}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
