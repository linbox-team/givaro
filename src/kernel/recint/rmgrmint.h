/* rmint/rmint_mg.h - Class definition of rmint<K,MG_ACTIVE> from RecInt library

   Copyright Université Grenoble Alpes
Contributors :
Alexis BREUST (alexis.breust@gmail.com 2014)
Christophe CHABOT (christophechabotcc@gmail.com 2011)
Jean-Guillaume Dumas

Time-stamp: <20 Jun 12 10:31:24 Jean-Guillaume.Dumas@imag.fr>

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

#ifndef RMINT_MG_RMINT_H
#define RMINT_MG_RMINT_H

/* If not previously defined, MG_DEFAULT is set
   for this file to run well. */
#if not defined(MG_DEFAULT)
#define MG_DEFAULT MG_ACTIVE
#endif

#include "rint.h"
#include "ruruint.h"
#include "rutools.h" /* mod_n() */
#include "rmdefine.h"
#include "rmgreduc.h"

// --------------------------------------------------------------
// ------- Declaration of class rmint (with Montgomery) ---------

namespace RecInt
{
    /* For modular calculus */
    template <size_t K> class rmint<K, MGA> {
    public:
        // p is the module (must be odd and > 1)
        static ruint<K> p;
        // p1 = -inv(p) mod 2^(2^K)
        static ruint<K> p1;
        // r = 2^(2^K) mod p
        static ruint<K> r;
        // Current value (always < p)
        ruint<K> Value;

        // Constructors
        rmint() : Value(0) {}
        rmint(const ruint<K>& c) : Value(c) { to_mg(*this); }
        rmint(const rint<K>& c) : Value( c.isNegative() ? (-c).Value : c.Value) {
            to_mg(*this); if (c.isNegative()) neg(*this); }
        rmint(const rmint<K, MGI>& c) : Value(c.Value) { to_mg(*this); }
        rmint(const rmint<K, MGA>& c) : Value(c.Value) {}
        template <typename T, __RECINT_IS_UNSIGNED(T, int) = 0> rmint(const T b) : Value(b) { to_mg(*this); }
        template <typename T, __RECINT_IS_SIGNED(T, int) = 0>   rmint(const T b) : Value((b < 0)? -b : b)
        { mod_n(Value, p); if (b < 0) sub(Value, p, Value); to_mg(*this); }

        rmint<K, MGA>& random();

        // Cast
        template <typename T, __RECINT_IS_ARITH(T, int) = 0> operator T() const { return T(get_ruint(*this)); }

        // Module functions
        static void init_module(const ruint<K>& p);
        static void get_module(ruint<K>& p);
    };

    /* Declarations for rmint<K, MG_ACTIVE> module */
    template <size_t K> ruint<K>   rmint<K, MGA>::p    = 0;
    template <size_t K> ruint<K>   rmint<K, MGA>::p1   = 0;
    template <size_t K> ruint<K>   rmint<K, MGA>::r    = 0;
}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
