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


#ifndef RUINT_ARITH_ADDMUL_H
#define RUINT_ARITH_ADDMUL_H

#include "ruruint.h"

#include "ruadd.h"
#include "rumul.h"

// --------------------------------------------------------------
// ----------------------- DEFINTIONS ---------------------------

namespace RecInt
{
    // a = b*c + d  or  (ahal) = b*c + d    (ah, al, b, c, d are ruint<K> and a is ruint<K+1>)
    // r is the carry
    template <size_t K> void laddmul(bool& r, ruint<K+1>& a, const ruint<K>& b, const ruint<K>& c, const ruint<K>& d);
    template <size_t K> void laddmul(bool& r, ruint<K>& ah, ruint<K>& al, const ruint<K>& b, const ruint<K>& c, const ruint<K>& d);
    template <size_t K> void laddmul(ruint<K+1>& a, const ruint<K>& b, const ruint<K>& c, const ruint<K>& d);
    template <size_t K> void laddmul(ruint<K>& ah, ruint<K>& al, const ruint<K>& b, const ruint<K>& c, const ruint<K>& d);

    // a = b*c + d  or  (ahal) = b*c + d    (ah, al, b, c are ruint<K> and a, d are ruint<K+1>)
    // r is the carry
    template <size_t K> void laddmul(bool& r, ruint<K+1>& a, const ruint<K>& b, const ruint<K>& c, const ruint<K+1>& d);
    template <size_t K> void laddmul(bool& r, ruint<K>& ah, ruint<K>& al, const ruint<K>& b, const ruint<K>& c, const ruint<K+1>& d);

    // a += b*c  (a, b are ruint and is ruint or UDItype)
    // The higher part is lost
    template <size_t K> void addmul(ruint<K>& a, const ruint<K>& b, const ruint<K>& c);
    template <size_t K> void addmul(ruint<K>& a, const ruint<K>& b, const UDItype& c);
}


// --------------------------------------------------------------
// -------------------------- Addmul ----------------------------

namespace RecInt
{
    // a = b*c + d (r stores the carry)
    // Note: this function is safe, ah|al is correctly computed
    // even if b, c, d are really ah or al
    template <size_t K>
    inline void laddmul(bool& r, ruint<K+1>& a, const ruint<K>& b, const ruint<K>& c, const ruint<K>& d) {
        laddmul(r, a.High, a.Low, b, c, d);
    }

    // (ah|al) = b*c + d (r stores the carry)
    // Note: this function is safe, ah|al is correctly computed
    // even if b, c, d are really ah or al
    template <size_t K>
    inline void laddmul(bool& r, ruint<K>& ah, ruint<K>& al, const ruint<K>& b, const ruint<K>& c, const ruint<K>& d) {
        bool rlow, rlow2, rmid, rhigh;
        ruint<K> bcmid, blcld;

        // Low part
        laddmul(rlow, blcld, b.Low, c.Low, d);

        // Middle part
        lmul(bcmid, b.High, c.Low);
        laddmul(rmid, bcmid, b.Low, c.High, bcmid);

        // High part
        laddmul(rhigh, ah, b.High, c.High, bcmid.High);

        // Below, we do not need b, c, d anymore, go fill ah|al no problem
        copy(al.Low, blcld.Low);
        add(rlow2, al.High, blcld.High, bcmid.Low);

        if (rlow)  add_1(rlow, ah);
        if (rlow2) add_1(rlow2, ah);
        if (rmid)  add_1(rmid, ah.High);

        r = (rlow || rlow2 || rmid || rhigh);
    }

    template <>
    inline void laddmul(bool& r, ruint<__RECINT_LIMB_SIZE>& ah, ruint<__RECINT_LIMB_SIZE>& al, const ruint<__RECINT_LIMB_SIZE>& b, const ruint<__RECINT_LIMB_SIZE>& c, const ruint<__RECINT_LIMB_SIZE>& d) {
        auto dp(d.Value);
        recint_umul_ppmm(ah.Value, al.Value, b.Value, c.Value);
        recint_add_ssaaaa(ah.Value, al.Value, ah.Value, al.Value, 0, dp);
        r = ((ah.Value == 0) && (al.Value < dp));
    }

    // a = b*c + d (the carry is lost)
    template <size_t K>
    inline void laddmul(ruint<K+1>& a, const ruint<K>& b, const ruint<K>& c, const ruint<K>& d) {
        laddmul(a.High, a.Low, b, c, d);
    }

    // (ah|al) = b*c + d (the carry is lost)
    // Note: this function is safe, ah|al is correctly computed
    // even if b, c, d are really ah or al
    template <size_t K>
    inline void laddmul(ruint<K>& ah, ruint<K>& al, const ruint<K>& b, const ruint<K>& c, const ruint<K>& d) {
        bool rlow, rlow2, rmid;
        ruint<K> bcmid, blcld;

        // Low part
        laddmul(rlow, blcld, b.Low, c.Low, d);

        // Middle part
        lmul(bcmid, b.High, c.Low);
        laddmul(rmid, bcmid, b.Low, c.High, bcmid);

        // High part
        laddmul(ah, b.High, c.High, bcmid.High);

        // Below, we do not need b, c, d anymore, go fill ah|al no problem
        copy(al.Low, blcld.Low);
        add(rlow2, al.High, blcld.High, bcmid.Low);

        if (rlow)  add_1(ah);
        if (rlow2) add_1(ah);
        if (rmid)  add_1(ah.High);
    }

    template <>
    inline void laddmul(ruint<__RECINT_LIMB_SIZE>& ah, ruint<__RECINT_LIMB_SIZE>& al, const ruint<__RECINT_LIMB_SIZE>& b, const ruint<__RECINT_LIMB_SIZE>& c, const ruint<__RECINT_LIMB_SIZE>& d) {
        auto dp(d.Value);
        recint_umul_ppmm(ah.Value, al.Value, b.Value, c.Value);
        recint_add_ssaaaa(ah.Value, al.Value, ah.Value, al.Value, 0, dp);
    }

    // a = b*c + d (r stores the carry)
    template <size_t K>
    inline void laddmul(bool& r, ruint<K+1>& a, const ruint<K>& b, const ruint<K>& c, const ruint<K+1>& d) {
        laddmul(r, a.High, a.Low, b, c, d);
    }

    // (ah|al) = b*c + d (r stores the carry)
    // Note: this function is safe, ah|al is correctly computed
    // even if b, c, d.High, d.Low are really ah or al
    template <size_t K>
    inline void laddmul(bool& r, ruint<K>& ah, ruint<K>& al, const ruint<K>& b, const ruint<K>& c, const ruint<K+1>& d) {
        bool rlow, rlow2, rmid, rmid2, rhigh;
        ruint<K> bcmid, blcldl;

        // Low part
        laddmul(rlow, blcldl, b.Low, c.Low, d.Low);

        // Middle part
        lmul(bcmid, b.High, c.Low);
        laddmul(rmid, bcmid, b.Low, c.High, bcmid);

        // High part
        laddmul(rhigh, ah, b.High, c.High, d.High);

        // Below, we do not need b, c, d anymore, go fill ah|al no problem
        copy(al.Low, blcldl.Low);
        add(rlow2, al.High, blcldl.High, bcmid.Low);
        add(rmid2, ah.Low, bcmid.High);

        if (rlow)  add_1(rlow,  ah);
        if (rlow2) add_1(rlow2, ah);
        if (rmid)  add_1(rmid,  ah.High);
        if (rmid2) add_1(rmid2, ah.High);

        r = (rlow || rlow2 || rmid || rmid2 || rhigh);
    }

    template <> inline void
    laddmul(bool& r, ruint<__RECINT_LIMB_SIZE>& ah, ruint<__RECINT_LIMB_SIZE>& al, const ruint<__RECINT_LIMB_SIZE>& b, const ruint<__RECINT_LIMB_SIZE>& c, const ruint<__RECINT_LIMB_SIZE+1>& d) {
        auto dph(d.High.Value);
        auto dpl(d.Low.Value);

        recint_umul_ppmm(ah.Value, al.Value, b.Value, c.Value);
        recint_add_ssaaaa(ah.Value, al.Value, ah.Value, al.Value, dph, dpl);
        r = ((ah.Value < dph) || ((ah.Value == dph) && (al.Value < dpl)));
    }

    // a += b*c
    // The higher part is lost
    template <size_t K>
    inline void addmul(ruint<K>& a, const ruint<K>& b, const ruint<K>& c) {
        ruint<K> bc;
        mul(bc, b, c);
        add(a, bc);
    }
#if defined(__RECINT_USE_FAST_128)
    template <>
    inline void addmul(ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE+1>& b, const ruint<__RECINT_LIMB_SIZE+1>& c) {
        a.Value += b.Value * c.Value;
    }
#endif
    template <>
    inline void addmul(ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b, const ruint<__RECINT_LIMB_SIZE>& c) {
        a.Value += b.Value * c.Value;
    }

    // a += b*c with c an integer
    // The higher part is lost
    template <size_t K>
    inline void addmul(ruint<K>& a, const ruint<K>& b, const UDItype& c) {
        ruint<K> bc;
        mul(bc, b, c);
        add(a, bc);
    }
}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
