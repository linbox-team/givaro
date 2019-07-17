/* ruint/misc.h - Miscellaneous functions for ruint

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


#ifndef RUINT_RANDOM_H
#define RUINT_RANDOM_H

#include <random>

#include "ruruint.h"

// --------------------------------------------------------------
// ----------------------- DEFINTIONS ---------------------------

namespace RecInt
{
    /* Define the global variable for random */
    static std::mt19937_64 rand_gen;

    // a is set to a random value (0 <= a < 2^(2^K))
    template <size_t K> ruint<K>& rand(ruint<K>& a);

    // Set the random seed for RecInt library
    void srand(const limb s);
}


// --------------------------------------------------------------
// ------------------------- Random -----------------------------

namespace RecInt
{

    // a is set to a random value
    template <size_t K> inline ruint<K>& rand(ruint<K>& a) {
        rand(a.High);
        rand(a.Low);
        return a;
    }
    template <> inline ruint<__RECINT_LIMB_SIZE>& rand(ruint<__RECINT_LIMB_SIZE>& a) {
        a.Value = rand_gen();
        return a;
    }

    // Set the random seed for RecInt library
    inline void srand(const limb s) {
        std::srand((unsigned int)(s));
        rand_gen.seed(s);
    }
}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
