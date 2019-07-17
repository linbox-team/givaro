/* ruint/misc.h - Miscellaneous functions for ruint

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


#ifndef RUINT_TOOLS_H
#define RUINT_TOOLS_H

#include "ruruint.h"

// --------------------------------------------------------------
// ----------------------- DEFINTIONS ---------------------------

namespace RecInt
{
    // d is the number of empty bits behind b
    template <size_t K> void normalization(UDItype& d, const ruint<K>& b);

    //-------- Modular algorithms ----------

    // a = b mod n or a = a mod n
    template <size_t K> void mod_n(ruint<K>& a, const ruint<K+1>& b, const ruint<K>& n);
    template <size_t K> void mod_n(ruint<K>& a, const ruint<K>& b,   const ruint<K>& n);
    template <size_t K> void mod_n(ruint<K>& a, const ruint<K>& n);
}


// --------------------------------------------------------------
// ---------------------- Normalization -------------------------

namespace RecInt
{
    // d is the number of empty bits behind b
    template <size_t K> inline void normalization(UDItype& d, const ruint<K>& b) {
        typename ruint<K>::cr_iterator it(b.rbegin());
        limb mask;

        for (d = 0; it != b.rend(); it++) {
            if (*it == 0) d += __RECINT_LIMB_BITS;
            else {
                for (mask = __RECINT_MAXPOWTWO; mask != 0; mask >>= 1) {
                    if ((*it) & mask) return;
                    else ++d;
                }
            }
        }
    }
}


// --------------------------------------------------------------
// ------------------- Modular algorithms -----------------------

namespace RecInt
{
    // a = b mod n
    template <size_t K>
    inline void mod_n(ruint<K>& a, const ruint<K>& b, const ruint<K>& n) {
        a = b % n;
    }

    // a = b mod n
    template <size_t K>
    inline void mod_n(ruint<K>& a, const ruint<K+1>& b, const ruint<K>& n) {
        UDItype d;
        ruint<K+2> bb;
        ruint<K> nn, q, r;

        // find number of zeros behind n
        normalization(d, n);
        // bb = b << d and nn = n << d
        left_shift(bb, b, d);
        left_shift(nn, n, d);
        // (bb.High.Low|bb.Low.High) =  q * nn + r <=> r = (bb.High.Low|bb.Low.High) mod nn
        div_2_1(q, r, bb.High.Low, bb.Low.High, nn);
        // (r|bb.Low.Low) =  q * nn + a <=> a = (r|bb.Low.Low) mod nn
        div_2_1(q, a, r, bb.Low.Low, nn);
        // a >>= d
        right_shift(a, a, d);
    }

    // a = a mod n
    template <size_t K>
    inline void mod_n(ruint<K>& a, const ruint<K>& n) {
        a %= n;
    }
}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
