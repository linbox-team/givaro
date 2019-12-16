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


#ifndef RUINT_ARITH_MUL_H
#define RUINT_ARITH_MUL_H

#include "ruruint.h"
#include "rucmp.h"

// --------------------------------------------------------------
// ----------------------- DEFINTIONS ---------------------------

namespace RecInt
{
    template <size_t K> ruint<K>& operator*=(ruint<K>&, const ruint<K>&);
    template <size_t K, typename T> __RECINT_IS_UNSIGNED(T, ruint<K>&) operator*=(ruint<K>&, const T&);
    template <size_t K, typename T> __RECINT_IS_SIGNED(T, ruint<K>&)   operator*=(ruint<K>&, const T&);

    template <size_t K> ruint<K> operator*(const ruint<K>&, const ruint<K>&);
    template <size_t K, typename T> __RECINT_IS_UNSIGNED(T, ruint<K>) operator*(const ruint<K>&, const T&);
    template <size_t K, typename T> __RECINT_IS_UNSIGNED(T, ruint<K>) operator*(const T&, const ruint<K>&);
    template <size_t K, typename T> __RECINT_IS_SIGNED(T, ruint<K>)   operator*(const ruint<K>&, const T&);
    template <size_t K, typename T> __RECINT_IS_SIGNED(T, ruint<K>)   operator*(const T&, const ruint<K>&);

    // a = (ahal) = b*c with naive method
    template <size_t K> void lmul_naive(ruint<K+1>& a, const ruint<K>& b, const ruint<K>& c);
    template <size_t K> void lmul_naive(ruint<K>& ah, ruint<K>& al, const ruint<K>& b, const ruint<K>& c);

    // a = (ahal) = b*c with Karatsuba method
    template <size_t K> void lmul_kara(ruint<K+1>& a, const ruint<K>& b, const ruint<K>& c);
    template <size_t K> void lmul_kara(ruint<K>& ah, ruint<K>& al, const ruint<K>& b, const ruint<K>& c);

    // Choose between naive and Karatsuba methods according to __RECINT_THRESHOLD_KARA constant
    // a = (ahal) = b*c
    template <size_t K> void lmul(ruint<K+1>& a, const ruint<K>& b, const ruint<K>& c);
    template <size_t K> void lmul(ruint<K>& ah, ruint<K>& al, const ruint<K>& b, const ruint<K>& c);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, void) lmul(ruint<K+1>& a, const ruint<K>& b, const T& c);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, void) lmul(limb& ah, ruint<K>& al, const ruint<K>& b, const T& c);

    // a = (b*c).Low    or a = (a*c).Low
    // The higher part is lost
    template <size_t K> ruint<K>& mul(ruint<K>& a, const ruint<K>& b, const ruint<K>& c);
    template <size_t K> ruint<K>& mul(ruint<K>& a, const ruint<K>& c);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, ruint<K>&) mul(ruint<K>& a, const ruint<K>& b, const T& c);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, ruint<K>&) mul(ruint<K>& a, const T& c);

    // a = b*b
    template <size_t K> ruint<K+1>& lsquare(ruint<K+1>& a, const ruint<K>& b);
    template <size_t K> ruint<K>& square(ruint<K>& al, const ruint<K>& b);
}


// --------------------------------------------------------------
// ------------------------ Operators ---------------------------

namespace RecInt
{
    // Operator *=
    template <size_t K>
    inline ruint<K>& operator*=(ruint<K>& a, const ruint<K>& b) {
        return mul(a, b);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_UNSIGNED(T, ruint<K>&) operator*=(ruint<K>& a, const T& b) {
        return mul(a, b);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_SIGNED(T, ruint<K>&) operator*=(ruint<K>& a, const T& b) {
        if (b < 0) {
            mul(a, -b);
            return (a = -a);
        } else return mul(a, b);
    }

    // Operator *
    template <size_t K>
    inline ruint<K> operator*(const ruint<K>& b, const ruint<K>& c) {
        ruint<K> a;
        return mul(a, b, c);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_UNSIGNED(T, ruint<K>) operator*(const ruint<K>& b, const T& c) {
        ruint<K> a;
        return mul(a, b, c);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_UNSIGNED(T, ruint<K>) operator*(const T& c, const ruint<K>& b) {
        ruint<K> a;
        return mul(a, b, c);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_SIGNED(T, ruint<K>) operator*(const ruint<K>& b, const T& c) {
        ruint<K> a;
        if (c < 0) return -mul(a, b, -c);
        else return mul(a, b, c);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_SIGNED(T, ruint<K>) operator*(const T& c, const ruint<K>& b) {
        ruint<K> a;
        if (c < 0) return -mul(a, b, -c);
        else return mul(a, b, c);
    }
}


// --------------------------------------------------------------
// --------------------- Multiplication -------------------------

namespace RecInt
{
    // a = ahal = b*c   with naive method
    // Note: this function is safe, ah|al is correctly computed
    // even if b, c are really ah or al

    template <size_t K>
    inline void lmul_naive(ruint<K>& ah, ruint<K>& al, const ruint<K>& b, const ruint<K>& c) {
        bool rmid, rlow;
        ruint<K> bcmid, blcl;

        // Low part
        lmul_naive(blcl, b.Low, c.Low);

        // Middle part
        lmul_naive(bcmid, b.High, c.Low);
        laddmul(rmid, bcmid, b.Low, c.High, bcmid);

        // High part
        laddmul(ah, b.High, c.High, bcmid.High);

        // Below, we do not need b, c, d anymore, go fill ah|al no problem
        copy(al.Low, blcl.Low);
        add(rlow, al.High, blcl.High, bcmid.Low);

        if (rlow)  add_1(ah);
        if (rmid)  add_1(rmid, ah.High);
    }

    template<>
    inline void lmul_naive(ruint<__RECINT_LIMB_SIZE>& ah, ruint<__RECINT_LIMB_SIZE>& al, const ruint<__RECINT_LIMB_SIZE>& b, const ruint<__RECINT_LIMB_SIZE>& c) {
        recint_umul_ppmm(ah.Value, al.Value, b.Value, c.Value);
    }
    template <size_t K>
    inline void lmul_naive(ruint<K+1>& a, const ruint<K>& b, const ruint<K>& c) {
        lmul_naive(a.High, a.Low, b, c);
    }

    // a = ahal = b*c   with Karatsuba method
    // Note: FIXME NOT safe - if al or ah is the variable than b or c, then pb
    template <size_t K>
    inline void lmul_kara(ruint<K>& ah, ruint<K>& al, const ruint<K>& b, const ruint<K>& c) {
        ruint<K-1> bb, cc;
        ruint<K> bc;
        bool rb, rc, r, rt1 = false, rt2 = false, rt3, rt4, rt5, rt6;

        add(rb, bb, b.High, b.Low);
        add(rc, cc, c.High, c.Low);
        lmul(ah, b.High, c.High);
        lmul(al, b.Low, c.Low);
        lmul(bc, bb, cc);

        if (rb) add(rt1, bc.High, cc);
        if (rc) add(rt2, bc.High, bb);

        sub(rt3, bc, ah);
        sub(rt4, bc, al);
        r = (rb&rc)+rt1+rt2-rt3-rt4;
        add(rt5, al.High, bc.Low);
        if (rt5) add_1(ah);
        add(rt6, ah.Low, bc.High);
        if (rt6 || r) add(ah.High, rt6 + r);
    }
    template<>
    inline void lmul_kara(ruint<__RECINT_LIMB_SIZE>& ah, ruint<__RECINT_LIMB_SIZE>& al, const ruint<__RECINT_LIMB_SIZE>& b, const ruint<__RECINT_LIMB_SIZE>& c) {
        lmul_naive(ah, al, b, c);
    }
    template <size_t K>
    inline void lmul_kara(ruint<K+1>& a, const ruint<K>& b, const ruint<K>& c) {
        lmul_kara(a.High, a.Low, b, c);
    }

    // Choose between naive and Karatsuba methods according to __RECINT_THRESHOLD_KARA constant
    // a = ahal = b*c
    template <size_t K>
    inline void lmul(ruint<K>& ah, ruint<K>& al, const ruint<K>& b, const ruint<K>& c) {
        if (K < __RECINT_THRESHOLD_KARA) lmul_naive(ah, al, b, c);
        else lmul_kara(ah, al, b, c);
    }
    template<>
    inline void lmul(ruint<__RECINT_LIMB_SIZE>& ah, ruint<__RECINT_LIMB_SIZE>& al, const ruint<__RECINT_LIMB_SIZE>& b, const ruint<__RECINT_LIMB_SIZE>& c) {
        lmul_naive(ah, al, b, c);
    }
    template <size_t K>
    inline void lmul(ruint<K+1>& a, const ruint<K>& b, const ruint<K>& c) {
        lmul(a.High, a.Low, b, c);
    }

    // ret|a = b*c
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, void) lmul(limb& ret, ruint<K>& a, const ruint<K>& b, const T& c) {
        limb retl;
        bool ret_temp;
        lmul(retl, a.Low, b.Low, c);
        lmul(ret, a.High, b.High, c);
        add(ret_temp, a.High, retl);
        ret += ret_temp;
    }
    template<typename T>
    inline __RECINT_IS_ARITH(T, void) lmul(limb& ret, ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b, const T& c) {
        recint_umul_ppmm(ret, a.Value, b.Value, limb(c));
    }

    // a = b*c
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, void) lmul(ruint<K+1>& a, const ruint<K>& b, const T& c) {
        limb ret;
        lmul(ret, a.Low, b, c);
        a.High = ret;
    }

    // al = (b*c).Low
    // The higher part is lost
    // Note: this function is safe, al is correctly computed
    // even if b, c are really al
    template <size_t K>
    inline ruint<K>& mul(ruint<K>& al, const ruint<K>& b, const ruint<K>& c) {
        ruint<K-1> bcmid;
        mul(bcmid, b.Low, c.High);
        addmul(bcmid, b.High, c.Low);
        lmul(al, b.Low, c.Low);
        add(al.High, bcmid);
        return al;
    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline ruint<__RECINT_LIMB_SIZE+1>& mul(ruint<__RECINT_LIMB_SIZE+1>& al, const ruint<__RECINT_LIMB_SIZE+1>& b, const ruint<__RECINT_LIMB_SIZE+1>& c) {
        al.Value = b.Value * c.Value;
        return al;
    }
#endif
    template<>
    inline ruint<__RECINT_LIMB_SIZE>& mul(ruint<__RECINT_LIMB_SIZE>& al, const ruint<__RECINT_LIMB_SIZE>& b, const ruint<__RECINT_LIMB_SIZE>& c) {
        al.Value = b.Value * c.Value;
        return al;
    }

    // al = (al*c).Low
    // The higher part is lost
    // Note: this function is safe, al is correctly computed
    // even if c is really al
    template <size_t K>
    inline ruint<K>& mul(ruint<K>& al, const ruint<K>& c) {
        ruint<K-1> acmid;
        mul(acmid, al.Low, c.High);
        addmul(acmid, al.High, c.Low);
        lmul(al, al.Low, c.Low);
        add(al.High, acmid);
        return al;
    }
#if defined(__RECINT_USE_FAST_128)
    template<>
    inline ruint<__RECINT_LIMB_SIZE+1>& mul(ruint<__RECINT_LIMB_SIZE+1>& al, const ruint<__RECINT_LIMB_SIZE+1>& c) {
        al.Value *= c.Value;
        return al;
    }
#endif
    template<>
    inline ruint<__RECINT_LIMB_SIZE>& mul(ruint<__RECINT_LIMB_SIZE>& al, const ruint<__RECINT_LIMB_SIZE>& c) {
        al.Value *= c.Value;
        return al;
    }

    // a = (b*c).Low
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, ruint<K>&) mul(ruint<K>& a, const ruint<K>& b, const T& c) {
        limb ret;
        lmul(ret, a, b, c);
        return a;
    }

    // a = (a*b).Low
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, ruint<K>&) mul(ruint<K>& a, const T& b) {
        limb ret;
        lmul(ret, a, a, b);
        return a;
    }
}


// --------------------------------------------------------------
// -------------------------- Square -----------------------------

namespace RecInt
{
    // a = b*b
    template <size_t K>
    inline ruint<K+1>& lsquare(ruint<K+1>& a, const ruint<K>& b) {
        bool rbb, ralb, rbah;
        ruint<K> bhbl;

        // a = bh*bh*B*B + 2*B*bh*bl + bl*bl
        lmul(bhbl, b.High, b.Low);
        lsquare(a.High, b.High);
        lsquare(a.Low, b.Low);

        rbb = highest_bit(bhbl);
        left_shift_1(bhbl, bhbl);

        add(ralb, a.Low.High, bhbl.Low);
        add(rbah, a.High.Low, bhbl.High);

        if (ralb) add_1(a.High);
        if (rbah || rbb) add(a.High.High, rbah + rbb);
        return a;
    }
    template<>
    inline ruint<__RECINT_LIMB_SIZE+1>& lsquare(ruint<__RECINT_LIMB_SIZE+1>& a, const ruint<__RECINT_LIMB_SIZE>& b) {
        recint_umul_ppmm(a.High.Value, a.Low.Value, b.Value, b.Value);
        return a;
    }

    // <>|al = b*b
    template <size_t K>
    inline ruint<K>& square(ruint<K>& al, const ruint<K>& b) {
        ruint<K-1> bhbll;

        // al = (2*bh*bl).Low * B + bl*bl
        mul(bhbll, b.High, b.Low);
        lsquare(al, b.Low);

        left_shift_1(bhbll, bhbll);
        add(al.High, bhbll);

        return al;
    }
    template<>
    inline ruint<__RECINT_LIMB_SIZE>& square(ruint<__RECINT_LIMB_SIZE>& al, const ruint<__RECINT_LIMB_SIZE>& b) {
        al.Value = b.Value * b.Value;
        return al;
    }
}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
