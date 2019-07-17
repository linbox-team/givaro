/* misc/convert_gmp.h - Conversion functions between r(u/m)int and mpz_class from GMP

   Copyright Université Grenoble Alpes
Contributors :
Alexis BREUST (alexis.breust@gmail.com 2014)
Christophe CHABOT (christophechabotcc@gmail.com 2011)


This software is a computer program whose purpose is to provide an fixed precision arithmetic library.

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software.	You can	use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and	rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty	and the software's author,	the holder of the
economic rights,	and the successive licensors	have only	limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,	using,	modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean	that it is complicated to manipulate,	and	that	also
therefore means	that it is reserved for developers	and	experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,	more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/


#ifndef RINT_CONVERT_H
#define RINT_CONVERT_H

#include <cstddef> // required by gmp versions <= 5.1.3
#include <gmpxx.h>

#include "rrint.h"
#include "ruconvert.h"
#include "rufiddling.h" // Unary operator

// --------------------------------------------------------------
// ----------------------- DEFINTIONS ---------------------------

namespace RecInt
{
    // Converts a mpz_class to a ruint<K> and vice-versa
    template <size_t K> rint<K>& mpz_to_rint(ruint<K>&, const mpz_class&);
    template <size_t K> mpz_class& rint_to_mpz(mpz_class&, const rint<K>&);
}


// --------------------------------------------------------------
// --------------------- Convert ruint --------------------------

namespace RecInt
{
    // Convert a GMP integer into a ruint
    template <size_t K>
    inline rint<K>& mpz_to_rint(rint<K>& a, const mpz_class& b) {
        // If negative, get positive one and reinverse
        if (b < 0) {
            mpz_to_ruint(a.Value, -b);
            a.Value = -a.Value;
        }
        // If positive, no problem
        else {
            mpz_to_ruint(a.Value, b);
        }

        return a;
    }

    // Convert a ruint into a GMP integer
    template <size_t K>
    inline mpz_class& rint_to_mpz(mpz_class& a, const rint<K>& b) {
        // If negative, get positive one and reinverse
        if (b.isNegative()) {
            ruint_to_mpz(a, -b.Value);
            a = -a;
        }
        // If positive, no problem
        else {
            ruint_to_mpz(a, b.Value);
        }

        return a;
    }

    // Convert a rint into a GMP integer
    template <size_t K>
    inline mpz_ptr rint_to_mpz_t(mpz_ptr a, const rint<K>& b) {
        // TODO Optimize...
        mpz_class r;
        RecInt::rint_to_mpz(r, b);
        mpz_init_set(a, r.get_mpz_t()) ;
        return a;
    }

    // Convert a rint into a GMP integer
    template <size_t K>
    inline rint<K>& mpz_t_to_rint(rint<K>& a, mpz_srcptr b) {
        // TODO Optimize...
        mpz_class r(b);
        return RecInt::mpz_to_rint(a, r);
    }
}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
