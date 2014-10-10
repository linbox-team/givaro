/* ruint/arith.h - Arithmetic functions for ruint

Copyright Université Joseph Fourier - Grenoble
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


#ifndef RUINT_ARITH_ADDMUL_H
#define RUINT_ARITH_ADDMUL_H

#include "ruruint.h"

#include "ruadd.h"
#include "rumul.h"

// --------------------------------------------------------------
// ----------------------- DEFINTIONS ---------------------------

namespace RecInt
{
    // a = b*c + d  or  (ahal) = b*c + d    (ah, al, b, c, d are ruint<K> and a is ruint<K+1>)
    // r is the carry
    template <size_t K> void laddmul(bool& r, ruint<K+1>& a, const ruint<K>& b, const ruint<K>& c, const ruint<K>& d);
    template <size_t K> void laddmul(bool& r, ruint<K>& ah, ruint<K>& al, const ruint<K>& b, const ruint<K>& c, const ruint<K>& d);
    template <size_t K> void laddmul(ruint<K+1>& a, const ruint<K>& b, const ruint<K>& c, const ruint<K>& d);
    template <size_t K> void laddmul(ruint<K>& ah, ruint<K>& al, const ruint<K>& b, const ruint<K>& c, const ruint<K>& d);

    // a = b*c + d  or  (ahal) = b*c + d    (ah, al, b, c are ruint<K> and a, d are ruint<K+1>)
    // r is the carry
    template <size_t K> void laddmul(bool& r, ruint<K+1>& a, const ruint<K>& b, const ruint<K>& c, const ruint<K+1>& d);
    template <size_t K> void laddmul(bool& r, ruint<K>& ah, ruint<K>& al, const ruint<K>& b, const ruint<K>& c, const ruint<K+1>& d);

    // a += b*c  (a, b, c are ruint)
    // r is the carry
    // The higher part is lost
    template <size_t K> void addmul(bool& r, ruint<K>& a, const ruint<K>& b, const ruint<K>& c);
    template <size_t K> void addmul(ruint<K>& a, const ruint<K>& b, const ruint<K>& c);

    // a += b*c  (a, b are ruint and c is UDItype)
    // r is the carry
    // The higher part is lost
    template <size_t K> void addmul(bool& r, ruint<K>& a, const ruint<K>& b, const UDItype& c);
    template <size_t K> void addmul(ruint<K>& a, const ruint<K>& b, const UDItype& c);
}


// --------------------------------------------------------------
// -------------------------- Addmul ----------------------------

namespace RecInt
{
    // a = b*c + d (r stores the carry)
    template <size_t K>
    inline void laddmul(bool& r, ruint<K+1>& a, const ruint<K>& b, const ruint<K>& c, const ruint<K>& d) {
        laddmul(r, a.High, a.Low, b, c, d);
    }

    // (ah|al) = b*c + d (r stores the carry)
    template <size_t K>
    inline void laddmul(bool& r, ruint<K>& ah, ruint<K>& al, const ruint<K>& b, const ruint<K>& c, const ruint<K>& d) {
        unsigned int ret1, ret2, ret3, ret4;
        ruint<K> bhcl, blch;

        laddmul(ret1, al, b.Low, c.Low, d);
        lmul(bhcl, b.High, c.Low);
        laddmul(ret2, blch, b.Low, c.High, bhcl);
        laddmul(r, ah, b.High, c.High, blch.High);
        add(ret3, al.High, blch.Low);
        add(ret4, ah, ret1+ret3);
        r+=ret4;
        add(ret4, ah.High, ret2);
        r+=ret4;
    }

    template <size_t K>
    inline void laddmul(bool& r, ruint<LIMB_SIZE>& ah, ruint<LIMB_SIZE>& al, const ruint<LIMB_SIZE>& b, const ruint<LIMB_SIZE>& c, const ruint<LIMB_SIZE>& d) {
        umul_ppmm(ah.Value, al.Value, b.Value, c.Value);
        add_ssaaaa(ah.Value, al.Value, ah.Value, al.Value, 0, d.Value);
        r = ((ah.Value == 0) && (al.Value < d.Value));
    }

    // a = b*c + d (the carry is lost)
    template <size_t K>
    inline void laddmul(ruint<K+1>& a, const ruint<K>& b, const ruint<K>& c, const ruint<K>& d) {
        laddmul(a.High, a.Low, b, c, d);
    }

    // (ah|al) = b*c + d (the carry is lost)
    template <size_t K>
    inline void laddmul(ruint<K>& ah, ruint<K>& al, const ruint<K>& b, const ruint<K>& c, const ruint<K>& d) {
        bool ret1, ret2, ret3;
        ruint<K> bhcl, blch;

        laddmul(ret1, al, b.Low, c.Low, d);
        lmul(bhcl, b.High, c.Low);
        laddmul(ret2, blch, b.Low, c.High, bhcl);
        laddmul(ah, b.High, c.High, blch.High);
        add(ret3, al.High, blch.Low);
        if (ret1 || ret3) add(ah, ret1 + ret3);
        if (ret2) add_1(ah.High);
    }

    template <> 
    inline void laddmul(ruint<LIMB_SIZE>& ah, ruint<LIMB_SIZE>& al, const ruint<LIMB_SIZE>& b, const ruint<LIMB_SIZE>& c, const ruint<LIMB_SIZE>& d) {
        umul_ppmm(ah.Value, al.Value, b.Value, c.Value);
        add_ssaaaa(ah.Value, al.Value, ah.Value, al.Value, 0, d.Value);
    }

    // a = b*c + d (r stores the carry)
    template <size_t K>
    inline void laddmul(bool& r, ruint<K+1>& a, const ruint<K>& b, const ruint<K>& c, const ruint<K+1>& d) {
        laddmul(r, a.High, a.Low, b, c, d);
    }

    // (ah|al) = b*c + d (r stores the carry)
    template <size_t K>
    inline void laddmul(bool& r, ruint<K>& ah, ruint<K>& al, const ruint<K>& b, const ruint<K>& c, const ruint<K+1>& d) {
        bool ret1, ret2, ret3, ret4;
        ruint<K> bhcl, blch;

        laddmul(ret1, al, b.Low, c.Low, d.Low);
        lmul(bhcl, b.High, c.Low);
        laddmul(ret2, blch, b.Low, c.High, bhcl);
        laddmul(r, ah, b.High, c.High, d.High);
        add(ret3, al.High, blch.Low);
        add(ret4, ah, ret1+ret3);
        r += ret4;
        add(ret3, ah.Low, blch.High);
        if (ret2 || ret3) {
            add(ret4, ah.High, ret2 + ret3);
            r += ret4;
        }
    }

    template <> inline void
    laddmul(bool& r, ruint<LIMB_SIZE>& ah, ruint<LIMB_SIZE>& al, const ruint<LIMB_SIZE>& b, const ruint<LIMB_SIZE>& c, const ruint<LIMB_SIZE+1>& d) {
        umul_ppmm(ah.Value, al.Value, b.Value, c.Value);
        add_ssaaaa(ah.Value, al.Value, ah.Value, al.Value, d.High.Value, d.Low.Value);
        r = ((ah.Value < d.High.Value) || ((ah.Value == d.High.Value) && (al.Value < d.Low.Value)));
    }

    // a += b*c (r stores the carry)
    // The higher part is lost
    template <size_t K>
    inline void addmul(bool& r, ruint<K>& a, const ruint<K>& b, const ruint<K>& c) {
        ruint<K> bc;
        mul(bc, b, c);
        add(r, a, bc);
    }

    // a += b*c  (the carry is lost)
    // The higher part is lost
    template <size_t K>
    inline void addmul(ruint<K>& a, const ruint<K>& b, const ruint<K>& c) {
        ruint<K> bc;
        mul(bc, b, c);
        add(a, bc);
    }

    // a += b*c  (r stores the carry)
    // The higher part is lost
    template <size_t K>
    inline void addmul(bool& r, ruint<K>& a, const ruint<K>& b, const UDItype& c) {
        ruint<K> bc;
        mul(bc, b, c);
        add(r, a, bc);
    }

    // a += b*c  (the carry is lost)
    // The higher part is lost
    template <size_t K>
    inline void addmul(ruint<K>& a, const ruint<K>& b, const UDItype& c) {
        ruint<K> bc;
        mul(bc, b, c);
        add(a, bc);
    }
}

#endif

