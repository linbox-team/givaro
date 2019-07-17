/* rmint/montgomery.h - Montgomery functions for rmint

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


#ifndef RMINT_MG_REDUC_H
#define RMINT_MG_REDUC_H

#include "rutools.h" /* mod_n() */

// --------------------------------------------------------------
// ----------------------- DEFINITIONS --------------------------

namespace RecInt
{
    //-------- Reduction algorithms ----------

    // a is set to Montgomery representation of b
    template <size_t K> rmint<K, MGA>& to_mg(rmint<K, MGA>& a, const rmint<K, MGA>& b);
    template <size_t K> rmint<K, MGA>& to_mg(rmint<K, MGA>& a);

    // t = c * r1 mod p
    // where r1 = r^(-1) mod p, r = 2^(2^k) and p is the module of rmint
    template <size_t K> rmint<K, MGA>& reduction(rmint<K, MGA>& t, const ruint<K+1>& c);
    template <size_t K> rmint<K, MGA>& reduction(rmint<K, MGA>& t, const ruint<K>& c);
    template <size_t K> rmint<K, MGA>& reduction(rmint<K, MGA>& t, const rmint<K, MGA>& c);
    template <size_t K> rmint<K, MGA>& reduction(rmint<K, MGA>& t);

    //---------- Compatibility functions ----------

    // Returns a.Value, demontgomerized if necessary
    template <size_t K> ruint<K> get_ruint(const rmint<K, MGA>& a);

    // Reduction or montgomerizing
    template <size_t K> rmint<K, MGA>& get_ready(rmint<K, MGA>& a);
}


// --------------------------------------------------------------
// --------------------- Montgomerizing -------------------------

namespace RecInt
{
    // a is Montgomery representation of b
    template <size_t K>
    inline rmint<K, MGA>& to_mg(rmint<K, MGA>& a, const rmint<K, MGA>& b) {
        ruint<K+1> res;
        copy(res.High, b.Value);
        mod_n(a.Value, res, rmint<K, MGA>::p);
        return a;
    }

    template <size_t K>
    inline rmint<K, MGA>& to_mg(rmint<K, MGA>& a) {
        return to_mg(a, a);
    }
}


// --------------------------------------------------------------
// -------------------- Demontgomerizing ------------------------

namespace RecInt
{
    // t = a * r1 mod p
    // where r1 = r^(-1) mod p, r = 2^(2^k) and p is the module of rmint
    template <size_t K>
    inline rmint<K, MGA>& reduction(rmint<K, MGA>& t, const ruint<K+1>& a) {
        bool r;
        ruint<K> a0;

        // m = a * p1 mod r
        mul(a0, a.Low, rmint<K, MGA>::p1);
        // t|a0 = a * (p1 * p + 1)
        laddmul(r, t.Value, a0, a0, rmint<K, MGA>::p, a);

        if (r || t.Value >= rmint<K, MGA>::p) sub(t.Value, rmint<K, MGA>::p);
        return t;
    }

    // t = a * r1 mod p
    // where r1 = r^(-1) mod p, r = 2^(2^k)
    // Note: this function is safe, the result is correct
    // even if &a == &t.Value
    template <size_t K>
    inline rmint<K, MGA>& reduction(rmint<K, MGA>& t, const ruint<K>& a) {
        bool r;
        ruint<K> a0;

        // m = a * p1 mod r
        mul(a0, a, rmint<K, MGA>::p1);
        // t|a0 = a * (p1 * p + 1)
        laddmul(r, t.Value, a0, a0, rmint<K, MGA>::p, a);

        if (r || t.Value >= rmint<K, MGA>::p) sub(t.Value, rmint<K, MGA>::p);
        return t;
    }

    // t = a * r1 mod a.p
    // where r1 = r^(-1) mod p, r = 2^(2^k)
    template <size_t K>
    inline rmint<K, MGA>& reduction(rmint<K, MGA>& t, const rmint<K, MGA>& a) {
        return reduction(t, a.Value);
    }

    // t = t * r1 mod a.p
    // where r1 = r^(-1) mod p, r = 2^(2^k)
    template <size_t K>
    inline rmint<K, MGA>& reduction(rmint<K, MGA>& t) {
        return reduction(t, t.Value);
    }
}


// --------------------------------------------------------------
// ------------------ Compatibily functions ---------------------

namespace RecInt
{
    // Returns a.Value, demontgomerized if necessary
    template <size_t K>
    inline ruint<K> get_ruint(const rmint<K, MGA>& a) {
        rmint<K, MGA> ap(a);
        return reduction(ap).Value;
    }

    // Reduction or montgomerizing
    template <size_t K>
    inline rmint<K, MGA>& get_ready(rmint<K, MGA>& a) {
        return to_mg(a);
    }
}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
