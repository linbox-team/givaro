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


#ifndef RMINT_CONVERT_GMP_H
#define RMINT_CONVERT_GMP_H

#include <cstddef> // required by gmp versions <= 5.1.3
#include <gmpxx.h>

#include "ruconvert.h"

/** NOTE : For this common file, either basic/reduc.h or mg/reduc.h
  has to be pre-included. **/

// --------------------------------------------------------------
// ----------------------- DEFINTIONS ---------------------------

namespace RecInt
{
    // Converts a mpz_class to a rmint<K,MG> and vice-versa
    template <size_t K, size_t MG> rmint<K, MG>& mpz_to_rmint(rmint<K, MG>&, const mpz_class&);
    template <size_t K, size_t MG> mpz_class& rmint_to_mpz(mpz_class&, const rmint<K, MG>&);
}


// --------------------------------------------------------------
// --------------------- Convert rmint --------------------------

namespace RecInt
{
    // Convert a GMP integer into a rmint
    template <size_t K, size_t MG>
    inline rmint<K, MG>& mpz_to_rmint(rmint<K, MG>& a, const mpz_class& b) {
        mpz_to_ruint(a.Value, b);
        return get_ready(a);
    }

    // Convert a rmint into a GMP integer
    template <size_t K, size_t MG>
    inline mpz_class& rmint_to_mpz(mpz_class& a, const rmint<K, MG>& b) {
        return ruint_to_mpz(a, get_ruint(b));
    }
}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
