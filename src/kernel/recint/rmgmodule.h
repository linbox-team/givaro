/* rmint/mg/module.h - Module functions for rmint<K,MGA>

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


#ifndef RMINT_MG_MODULE_H
#define RMINT_MG_MODULE_H

#include "rmgrmint.h"
#include "ruruint.h"
#include "rumanip.h" /* copy() */
#include "rufiddling.h" /* bits manip */
#include "rumul.h" /* lmul() */
#include "ruadd.h" /* add() */


// --------------------------------------------------------------
// ----------------------- DEFINITIONS --------------------------

namespace RecInt
{
    //--------- Arazi & Qi ----------

    // If a is odd, u = inv(a) mod 2^(2^K)
    template <size_t K> ruint<K>& arazi_qi(ruint<K>& u, const ruint<K>& a);
}


// ---------------------------------------------------------------------
// ----------------- Arazi & Qi with fast basecase ---------------------
// -- See [JG Dumas, On Newton-Raphson iteration for                  --
// --     multiplicative inverses modulo prime powers.                --
// --     IEEE Transactions on Computers, 63(8), pp 2106-2109, 2014]. --

namespace RecInt
{
    // If a is odd, ia = inv(a) mod 2^(2^K)
    template <size_t K>
    inline ruint<K>& arazi_qi(ruint<K>& u, const ruint<K>& a) {
        ruint<K-1> t1, t2;

        // Get previous u - we're doing step i = m/2
        arazi_qi(u.Low, a.Low);

        // b = a & (2^i - 1)
        // t1 = u * b
        // t1 >>= i (NOTE: t2 is junk)
        lmul(t1, t2, u.Low, a.Low);

        // c = (a >> i) & (2^i - 1)
        // t2 = (u * c) & (2^i - 1)
        mul(t2, u.Low, a.High);

        // t1 += t2
        add(t1, t2);

        // t1 *= u
        // t1 &= (2^i - 1)
        mul(t1, u.Low);

        // t1 = 2^i - t1
        // t1 <<= i
        // u |= t1
        copy(u.High, -t1);

        return u;
    }

    // If a is odd, u = inv(a) mod 2^64
    template <>
    inline ruint<__RECINT_LIMB_SIZE>& arazi_qi(ruint<__RECINT_LIMB_SIZE>& u, const ruint<__RECINT_LIMB_SIZE>& a) {
        if (a.Value == 1) {
            u.Value = 1;
            return u;
        }

        UDItype amone(a.Value-1);
        u.Value = 1;

        for (size_t i = 2; i < __RECINT_LIMB_BITS; i <<= 1) {
            amone *= amone;
            u.Value *= ++amone;
            --amone;
        }

        u.Value *= (2 - a.Value); // 2-a
        return u;
    }
}


// --------------------------------------------------------------
// ------------------------- Module -----------------------------

namespace RecInt
{
    // Initialize the module of rmint to p (must be odd !)
    template <size_t K>
    inline void rmint<K, MGA>::init_module(const ruint<K>& _p) {
        // p is the new module
        copy(rmint<K, MGA>::p, _p);

        // p1 = -inv(p) mod 2^(2^K)
        arazi_qi(rmint<K, MGA>::p1, -_p);

        // r = 2^(2^K) mod p
        div_r(rmint<K, MGA>::r, -_p, _p);
    }

    // Get the module of rmint to p
    template <size_t K>
    inline void rmint<K, MGA>::get_module(ruint<K>& _p) {
        copy(_p, rmint<K, MGA>::p);
    }
}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
