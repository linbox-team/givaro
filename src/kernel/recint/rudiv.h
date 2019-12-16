/* ruint/arith.h - Arithmetic functions for ruint

   Copyright Université Grenoble Alpes
Contributors :
Alexis BREUST (alexis.breust@gmail.com 2014)
Christophe CHABOT (christophechabotcc@gmail.com 2011)
Jean-Guillaume DUMAS

Time-stamp: <20 Jun 12 10:28:29 Jean-Guillaume.Dumas@imag.fr>

This software is a computer program whose purpose is to provide an
fixed precision arithmetic library.

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


#ifndef RUINT_ARITH_DIV_H
#define RUINT_ARITH_DIV_H

#include "ruruint.h"
#include "rucmp.h"

#include "rushift.h" /* right_shift() */

// --------------------------------------------------------------
// ----------------------- DEFINTIONS ---------------------------

namespace RecInt
{
    template <size_t K> ruint<K>& operator%=(ruint<K>&, const ruint<K>&);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, ruint<K>&) operator%=(ruint<K>&, const T&);

    template <size_t K> ruint<K>& operator/=(ruint<K>&, const ruint<K>&);
    template <size_t K, typename T> __RECINT_IS_UNSIGNED(T, ruint<K>&) operator/=(ruint<K>&, const T&);
    template <size_t K, typename T> __RECINT_IS_SIGNED(T, ruint<K>&)   operator/=(ruint<K>&, const T&);

    template <size_t K> ruint<K> operator%(const ruint<K>&, const ruint<K>&);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, ruint<K>) operator%(const ruint<K>&, const T&);

    template <size_t K> ruint<K> operator/(const ruint<K>&, const ruint<K>&);
    template <size_t K, typename T> __RECINT_IS_UNSIGNED(T, ruint<K>) operator/(const ruint<K>&, const T&);
    template <size_t K, typename T> __RECINT_IS_SIGNED(T, ruint<K>)   operator/(const ruint<K>&, const T&);

    // Euclidean division of the 3-ruint integer (a2|a1|a0) by the 2-ruint integer (b1|b0)
    // the 1-ruint quotient is stored in q
    // the 2-ruint remainder is stored in (r1|r0)
    template <size_t K> void div_3_2(ruint<K>& q, ruint<K>& r1, ruint<K>& r0,
                                     const ruint<K>& a2, const ruint<K>& a1, const ruint<K>& a0,
                                     const ruint<K>& b1, const ruint<K>& b0);

    // Euclidean division of the 2-ruint integer (a1|a0) by b
    // q stores the quotient and r the remainder
    template <size_t K> void div_2_1(ruint<K>& q, ruint<K>& r,
                                     const ruint<K>& a1, const ruint<K>& a0,
                                     const ruint<K>& b);

    // computes (q, r) such that a = q*b + r (0 <= r < b)
    template <size_t K> void div(ruint<K>& q, ruint<K>& r, const ruint<K>& a, const ruint<K>& b);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, void) div(ruint<K>& q, T& r, const ruint<K>& a, const T& b);

    // q = floor(a/b)
    template <size_t K> ruint<K>& div_q(ruint<K>& q, const ruint<K>& a, const ruint<K>& b);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, ruint<K>&) div_q(ruint<K>& q, const ruint<K>& a, const T& b);

    // r = a mod b
    template <size_t K> ruint<K>& div_r(ruint<K>& r, const ruint<K>& a, const ruint<K>& b);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, T&) div_r(T& r, const ruint<K>& a, const T& b);
}


// --------------------------------------------------------------
// ------------------------ Operators ---------------------------

namespace RecInt
{
    // Operator %=
    template <size_t K>
    inline ruint<K>& operator%=(ruint<K>& a, const ruint<K>& b) {
        return div_r(a, a, b);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, ruint<K>&) operator%=(ruint<K>& a, const T& b) {
        T aa;
        div_r(aa, a, b);
        return (a = aa);
    }

    // Operator /=
    template <size_t K>
    inline ruint<K>& operator/=(ruint<K>& a, const ruint<K>& b) {
        return div_q(a, a, b);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_UNSIGNED(T, ruint<K>&) operator/=(ruint<K>& a, const T& b) {
        return div_q(a, a, b);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_SIGNED(T, ruint<K>&) operator/=(ruint<K>& a, const T& b) {
        if (b < 0) {
            div_q(a, a, -b);
            return (a = -a);
        } else return div_q(a, a, b);
    }

    // Operator %
    template <size_t K>
    inline ruint<K> operator%(const ruint<K>& b, const ruint<K>& c) {
        ruint<K> a;
        div_r(a, b, c);
        return a;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, ruint<K>) operator%(const ruint<K>& b, const T& c) {
        ruint<K> a;
        T aa;
        div_r(aa, b, c);
        return (a = aa);
    }

    // Operator /
    template <size_t K>
    inline ruint<K> operator/(const ruint<K>& b, const ruint<K>& c) {
        ruint<K> a;
        div_q(a, b, c);
        return a;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_UNSIGNED(T, ruint<K>) operator/(const ruint<K>& b, const T& c) {
        ruint<K> a;
        return div_q(a, b, c);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_SIGNED(T, ruint<K>) operator/(const ruint<K>& b, const T& c) {
        ruint<K> a;
        if (c < 0) {
            div_q(a, b, -c);
            return (a = -a);
        } else return div_q(a, b, c);
    }
}


// --------------------------------------------------------------
// ------------------------ Division ---------------------------

namespace RecInt
{
    // Euclidean division of the 3-ruint integer (a2|a1|a0) by the 2-ruint integer (b1|b0)
    // the 1-ruint quotient is stored in q
    // the 2-ruint remainder is stored in (r1|r0)
    template <size_t K>  inline void div_3_2(ruint<K>& q, ruint<K>& r1, ruint<K>& r0,
                                             const ruint<K>& a2, const ruint<K>& a1, const ruint<K>& a0,
                                             const ruint<K>& b1, const ruint<K>& b0) {
        ruint<K> c, d1, d0;
        bool ret_sub, ret1 = false;

        if (a2 < b1) {
            div_2_1(q, c, a2, a1, b1);
        } else {
            fill_with_1(q);
            add(ret1, c, a1, b1);
        }

        lmul(d1, d0, q, b0);
        sub(ret_sub, r0, a0, d0);
        sub_wc(r1, c, d1, ret_sub);

        if ((ret1 == 0) && (d1 > c || (d1 == c && d0 > a0))) {
            bool ret;
            sub_1(q);
            add(ret_sub, r0, b0);
            add_wc(ret, r1, b1, ret_sub);

            if (!ret) {
                sub_1(q);
                add(ret_sub, r0, b0);
                add_wc(r1, b1, ret_sub);
            }
        }
    }
    template <>
    inline void div_3_2(ruint<__RECINT_LIMB_SIZE>& q, ruint<__RECINT_LIMB_SIZE>& r1, ruint<__RECINT_LIMB_SIZE>& r0,
                        const ruint<__RECINT_LIMB_SIZE>& a2, const ruint<__RECINT_LIMB_SIZE>& a1, const ruint<__RECINT_LIMB_SIZE>& a0,
                        const ruint<__RECINT_LIMB_SIZE>& b1, const ruint<__RECINT_LIMB_SIZE>& b0) {
        limb c, d1, d0;
        bool ret = false;

        if (a2.Value < b1.Value) {
            recint_udiv_qrnnd(q.Value, c, a2.Value, a1.Value, b1.Value);
        } else {
            q.Value = __RECINT_MINUSONE;
            c = a1.Value + b1.Value;
            if (c < a1.Value)
                ret = true;
        }

        recint_umul_ppmm(d1, d0, q.Value, b0.Value);
        recint_sub_ddmmss(r1.Value, r0.Value, c, a0.Value, d1, d0);

        if (!ret && ((d1 > c) || ((d1 == c) && (d0 > a0.Value)))) {
            q.Value--;
            r0.Value += b0.Value;
            r1.Value += b1.Value;
            if (r0.Value < b0.Value)
                r1.Value++;

            if ((r1.Value > b1.Value) || ((r1.Value == b1.Value) && (r0.Value >= b0.Value))) {
                q.Value--;
                r0.Value += b0.Value;
                r1.Value += b1.Value;
                if (r0.Value<b0.Value)
                    r1.Value++;
            }
        }
    }

    // Euclidean division of the 2-ruint integer (a1|a0) by b
    // q stores the quotient and r the remainder
    template <size_t K>
    inline void div_2_1(ruint<K>& q, ruint<K>& r,
                        const ruint<K>& ah, const ruint<K>& al,
                        const ruint<K>& b) {
        ruint<K> s;
        div_3_2(q.High, s.High, s.Low, ah.High, ah.Low, al.High, b.High, b.Low);
        div_3_2(q.Low, r.High, r.Low, s.High, s.Low, al.Low, b.High, b.Low);
    }
    template <>
    inline void div_2_1(ruint<__RECINT_LIMB_SIZE>& q, ruint<__RECINT_LIMB_SIZE>& r,
                        const ruint<__RECINT_LIMB_SIZE>& ah, const ruint<__RECINT_LIMB_SIZE>& al,
                        const ruint<__RECINT_LIMB_SIZE>& b) {
        recint_udiv_qrnnd(q.Value, r.Value, ah.Value, al.Value, b.Value);
    }

    // computes (q, r) such that a = q*b + r (0 <= r < b)
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, void) div(ruint<K>& q, T& r, const ruint<K>& a, const T& b) {
        if (b == 2) {
            bool z;
            right_shift_1(z, q, a);
            r = (z)? 1: 0;
        } else {
            ruint<K> bb(b), rr;
            div(q, rr, a, bb);
            r = static_cast<T>(rr);
        }
    }

    // computes (q, r) such that a = q*b + r (0 <= r < b)
    template <size_t K>
    inline void div(ruint<K>& q, ruint<K>& r, const ruint<K>& a, const ruint<K>& b) {
        UDItype d;
        ruint<K+1> aa;
        ruint<K> bb;

        normalization(d, b);
        left_shift(aa, a, d);
        left_shift(bb, b, d);
        div_2_1(q, r, aa.High, aa.Low, bb);
        right_shift(r, r, d);
    }

    // computes (q, r) such that a = q*b + r (0 <= r < b)
    // Note: with the correct option (-O), a good compiler should compute a/b and a%b with one call
    struct sdiv { limb quot; limb rem; };
    inline void udiv_qrnd(limb& q, limb& r, const limb& a, const limb& b) {
        sdiv x{a/b, a%b};
        q = x.quot; r = x.rem;
    }
    template <>
    inline void div(ruint<__RECINT_LIMB_SIZE>& q, ruint<__RECINT_LIMB_SIZE>& r, const ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b) {
        udiv_qrnd(q.Value, r.Value, a.Value, b.Value);
    }

    // q = floor(a/b)
    template <size_t K>
    inline ruint<K>& div_q(ruint<K>& q, const ruint<K>& a, const ruint<K>& b) {
        ruint<K> r;
        div(q, r, a, b);
        return q;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, ruint<K>&) div_q(ruint<K>& q, const ruint<K>& a, const T& b) {
        ruint<K> r, bb(b);
        div(q, r, a, bb);
        return q;
    }

    // r = a mod b
    template <size_t K>
    inline ruint<K>& div_r(ruint<K>& r, const ruint<K>& a, const ruint<K>& b) {
        ruint<K> q;
        div(q, r, a, b);
        return r;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, T&) div_r(T& r, const ruint<K>& a, const T& b) {
        ruint<K> q;
        div(q, r, a, b);
        return r;
    }
}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
