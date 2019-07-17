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


#ifndef RMINT_BASIC_REDUC_H
#define RMINT_BASIC_REDUC_H

#include "rutools.h" /* mod_n() */

// --------------------------------------------------------------
// ----------------------- DEFINITIONS --------------------------

namespace RecInt
{
    //-------- Reduction algorithms ----------

    // t = c mod p
    template <size_t K> rmint<K, MGI>& reduction(rmint<K, MGI>& t, const ruint<K+1>& c);
    template <size_t K> rmint<K, MGI>& reduction(rmint<K, MGI>& t, const ruint<K>& c);
    template <size_t K> rmint<K, MGI>& reduction(rmint<K, MGI>& t, const rmint<K, MGI>& c);
    template <size_t K> rmint<K, MGI>& reduction(rmint<K, MGI>& t);

    //---------- Compatibility functions ----------

    // Returns a.Value, demontgomerized if necessary
    template <size_t K> ruint<K> get_ruint(const rmint<K, MGI>& a);

    // Reduction or montgomerizing
    template <size_t K> rmint<K, MGI>& get_ready(rmint<K, MGI>& a);
}


// --------------------------------------------------------------
// ------------------ Reduction algorithm -----------------------

namespace RecInt
{
    // t = c mod p
    template <size_t K>
    inline rmint<K, MGI>& reduction(rmint<K, MGI>& t, const ruint<K+1>& a) {
        mod_n(t.Value, a, t.p);
        return t;
    }

    // t = c mod p
    template <size_t K>
    inline rmint<K, MGI>& reduction(rmint<K, MGI>& t, const ruint<K>& a) {
        mod_n(t.Value, a, t.p);
        return t;
    }

    // t = c mod p
    template <size_t K>
    inline rmint<K, MGI>& reduction(rmint<K, MGI>& t, const rmint<K, MGI>& c) {
        mod_n(t.Value, c.Value, t.p);
        return t;
    }

    // t = t mod p
    template <size_t K>
    inline rmint<K, MGI>& reduction(rmint<K, MGI>& t) {
        mod_n(t.Value, t.p);
        return t;
    }
}


// --------------------------------------------------------------
// ------------------ Compatibily functions ---------------------

namespace RecInt
{
    // Returns a.Value, demontgomerized if necessary
    template <size_t K>
    inline ruint<K> get_ruint(const rmint<K, MGI>& a) {
        return a.Value;
    }

    // Reduction or montgomerizing
    template <size_t K>
    inline rmint<K, MGI>& get_ready(rmint<K, MGI>& a) {
        mod_n(a.Value, a.p);
        return a;
    }
}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
