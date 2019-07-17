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


#ifndef RMINT_MG_ARITH_EXP_H
#define RMINT_MG_ARITH_EXP_H

#include "rmgrmint.h"
#include "rmgreduc.h"

// --------------------------------------------------------------
// ----------------------- DEFINTIONS ---------------------------

namespace RecInt
{
    template <size_t K> void exp(rmint<K, MGA>&, const rmint<K, MGA>&, const UDItype&);
    template <size_t K> void exp(rmint<K, MGA>&, const rmint<K, MGA>&, const ruint<K>&);
}


// --------------------------------------------------------------
// --------------------- Implementation -------------------------

namespace RecInt
{
    // a = b^c mod a.p
    template <size_t K>
    inline void exp(rmint<K, MGA>& a, const rmint<K, MGA>& b, const ruint<K>& c) {
        limb **tab, **originalTab;
        originalTab = (limb**)malloc(NBLIMB<K>::value * sizeof(limb*)); // TODO Cleaner avec iterateurs
        pointers_list(originalTab, c);

        rmint<K, MGA> *g;
        g = (rmint<K, MGA>*)malloc(16*sizeof(rmint<K, MGA>));

        int i, j, NB_WIN;

        UDItype exp, mask=0xf;
        NB_WIN=__RECINT_LIMB_BITS/4;

        copy(g[0].Value, rmint<K, MG_ACTIVE>::r);
        for (i=1;i<16;i++) mul(g[i], g[i-1], b);

        copy(a.Value, rmint<K, MG_ACTIVE>::r);
        tab=&originalTab[NBLIMB<K>::value - 1];

        for (i=NBLIMB<K>::value - 1;i>0;i--) {
            exp=**tab;
            for (j=NB_WIN-1;j>=0;j--) {
                mul(a, a, g[(exp >> (j<<2)) & mask]);
                square(a, a);
                square(a, a);
                square(a, a);
                square(a, a);
            }
            tab--;
        }
        exp=**tab;
        for (j=NB_WIN-1;j>0;j--) {
            mul(a, a, g[(exp >> (j<<2)) & mask]);
            square(a, a);
            square(a, a);
            square(a, a);
            square(a, a);
        }
        mul(a, a, g[exp & mask]);

        free(originalTab);
        free(g);
    }

    // a = b^c mod a.p
    template <size_t K>
    inline void exp(rmint<K, MGA>& a, const rmint<K, MGA>& b, const UDItype& c) {
        UDItype exp = c;
        rmint<K, MGA> x;

        // x = b
        copy(x, b);
        // a = 1
        copy(a.Value, rmint<K, MG_ACTIVE>::r);

        while (exp != 0) {
            if (exp & 1) mul(a, a, x);
            mul(x, x, x);
            exp = exp >> 1;
        }
    }
}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
