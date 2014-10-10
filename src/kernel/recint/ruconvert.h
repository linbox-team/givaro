/* misc/convert_gmp.h - Conversion functions between r(u/m)int and mpz_class from GMP

Copyright Université Joseph Fourier - Grenoble
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


#ifndef RUINT_CONVERT_H
#define RUINT_CONVERT_H

#include <gmpxx.h>

#include "rumanip.h" /* reset() */

// --------------------------------------------------------------
// ----------------------- DEFINTIONS ---------------------------

namespace RecInt
{
    // Converts a mpz_class to a ruint<K> and vice-versa
    template <size_t K> ruint<K>& mpz_to_ruint(ruint<K>&, const mpz_class&);
    template <size_t K> mpz_class& ruint_to_mpz(mpz_class&, const ruint<K>&);

}


// --------------------------------------------------------------
// --------------------- Convert ruint --------------------------

namespace RecInt
{
    // Convert a GMP integer into a ruint
    template <size_t K>
    inline ruint<K>& mpz_to_ruint(ruint<K>& a, const mpz_class& b) {
	    unsigned int i;
	    mpz_class c(b);
        limb l;

        reset(a);
	    for (i = 0; i < NBLIMB<K>::value; i++) {
		    l = c.get_ui(); c >>= 32;
		    l |= (limb(c.get_ui()) << 32);
		    set_limb(a, l, i); c >>= 32;
	    }

        return a;
    }

    template <>
    inline ruint<LIMB_SIZE>& mpz_to_ruint(ruint<LIMB_SIZE>& a, const mpz_class& b) {
	    mpz_class c(b);
        limb l;
        
	    l = c.get_ui(); c >>= 32;
	    l |= (limb(c.get_ui()) << 32);
	    a.Value = l;
        
        return a;
    }

    // Convert a ruint into a GMP integer
    template <size_t K>
    inline mpz_class& ruint_to_mpz(mpz_class& a, const ruint<K>& b) {
        a = 0;
        for( auto it(b.rbegin()); it != b.rend(); ++it) {
		    // GMP does not handle uint64_t, need to break it
            a <<= 32;
            a += USItype((*it) >> 32);
            a <<= 32;
		    a += USItype(*it);
        }
        
        return a;
    }
}

#endif

