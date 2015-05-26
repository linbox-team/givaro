/* ruint/misc.h - Miscellaneous functions for ruint

Copyright Université Joseph Fourier - Grenoble
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


#ifndef RUINT_MANIP_H
#define RUINT_MANIP_H

#include "ruruint.h"

// --------------------------------------------------------------
// ----------------------- DEFINTIONS ---------------------------

namespace RecInt
{
    //---------- Basic manipulation -----------

    // a = 0
    template <size_t K> void reset(ruint<K>& a);

    // a = -1 (a is filled with 1)
    template <size_t K> void fill_with_1(ruint<K>& a);

    // a = b
    template <size_t K> void copy(ruint<K>& a, const ruint<K>& b);

    //---------- Limb manipulation -----------

    // p is the ordered list of all limbs of a
    template <size_t K> void pointers_list(limb** p, const ruint<K>& a);

    // Set the index-th limb of a to b
    template <size_t K> void set_limb(ruint<K>& a, const limb& b, unsigned int index);
    
    // Get the index-th limb of a
    template <size_t K> limb ms_limb(const ruint<K>& a); // Most significant limb
    template <size_t K> limb get_limb(const ruint<K>& a, unsigned int index);
    template <size_t K> const limb* get_limb_p(const ruint<K>& a, unsigned int index);
}


// --------------------------------------------------------------
// ------------------- Basic manipulation -----------------------

namespace RecInt
{
    // a = 0
    template <size_t K> inline void reset(ruint<K>& a) {
        reset(a.High);
        reset(a.Low);
    }
    template <> inline void reset(ruint<__RECINT_LIMB_SIZE>& a) {
        a.Value = 0;
    }

    // a = -1 (a is filled with 1)
    template <size_t K> inline void fill_with_1(ruint<K>& a) {
        fill_with_1(a.High);
        fill_with_1(a.Low);
    }
    template <>  inline void fill_with_1(ruint<__RECINT_LIMB_SIZE>& a) {
        a.Value = __RECINT_MINUSONE;
    }

    // a = b
    template <size_t K> inline void copy(ruint<K>& a, const ruint<K>& b) {
        if (&a == &b) return;
        copy(a.High, b.High);
        copy(a.Low, b.Low);  
    }
    template <> inline void copy(ruint<__RECINT_LIMB_SIZE>& a, const ruint<__RECINT_LIMB_SIZE>& b) {
        a.Value = b.Value;
    }
}


// --------------------------------------------------------------
// -------------------- Limb manipulation -----------------------

namespace RecInt
{
    // tab is the ordered list of all limbs of a
    template <size_t K> inline void pointers_list(limb **tab, const ruint<K>& a) {
        pointers_list(tab, a.Low);
        pointers_list(&tab[NBLIMB<K-1>::value], a.High);
    }
    template <> inline void pointers_list(limb **tab, const ruint<__RECINT_LIMB_SIZE>& a) {
        const limb* toto = static_cast<const limb*>( &(a.Value) );
        *tab = const_cast<limb*>(toto);
    }


    // Set the index-th limb of a to b
    template <size_t K> inline void set_limb(ruint<K>& a, const limb& b, unsigned int index) {
        if (index < NBLIMB<K-1>::value) set_limb(a.Low, b, index);
        else set_limb(a.High, b, index-NBLIMB<K-1>::value);    
    }
    template <> inline void set_limb(ruint<__RECINT_LIMB_SIZE>& a, const limb& b, unsigned int index) {
        if (index == 0) a.Value = b;
    }


    // Get the index-th limb of a
    template <size_t K> inline limb ms_limb(const ruint<K>& a) {
        return ms_limb(a.High);  
    }
    template <> inline limb ms_limb(const ruint<__RECINT_LIMB_SIZE>& a) {
        return a.Value;
    }
    
    template <size_t K> inline limb get_limb(const ruint<K>& a, unsigned int index) {
        if (index < NBLIMB<K-1>::value) return get_limb(a.Low, index);
        else return get_limb(a.High, index - NBLIMB<K-1>::value);    
    }
    template <> inline limb get_limb(const ruint<__RECINT_LIMB_SIZE>& a, unsigned int) {
        return a.Value;
    }
    
    template <size_t K> inline const limb* get_limb_p(const ruint<K>& a, unsigned int index) {
        if (index < NBLIMB<K-1>::value) return get_limb_p(a.Low, index);
        else return get_limb_p(a.High, index - NBLIMB<K-1>::value);    
    }
    template <> inline const limb* get_limb_p(const ruint<__RECINT_LIMB_SIZE>& a, unsigned int) {
        return &(a.Value);
    }
}

#endif

