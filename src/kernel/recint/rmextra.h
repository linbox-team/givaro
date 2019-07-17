/* rmint/arith.h - Arithmetic functions for rmint

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


#ifndef RMINT_COMMON_EXTRA_H
#define RMINT_COMMON_EXTRA_H

// --------------------------------------------------------------
// ----------------------- DEFINITIONS --------------------------

namespace RecInt
{
    // returns true iff b exists such that b^2 = a mod a.p
    template <size_t K, size_t MG> bool is_quadratic_residue(const rmint<K, MG>& a);

    // computes r such that r*r = a mod a.p (if a is not a quadratic residue, r = 0)
    template <size_t K, size_t MG> void square_root(rmint<K, MG>& r, const rmint<K, MG>& a);
}


// --------------------------------------------------------------
// --------------------- Implementation -------------------------

namespace RecInt
{
    // Let a.p be an odd prime
    // returns true iff b exists such that b^2 = a mod a.p
    template <size_t K, size_t MG>
    inline bool is_quadratic_residue(const rmint<K, MG>& a) {
        rmint<K, MG> aa;
        ruint<K> pp(a.p - 1);

        right_shift_1(pp, pp);
        exp(aa, a, pp);

        return (aa == 1);
    }

    // computes r such that r*r = a mod a.p (if a is not a quadratic residue, r = 0)
    template <size_t K, size_t MG>
    inline void square_root(rmint<K, MG>& r, const rmint<K, MG>& a) {
        if (!is_quadratic_residue(a)) {
            // a is not a square
            reset(r);
            return;
        }

        ruint<K> temp, pp, t;
        rmint<K, MG> tempmod, ppmod;
        UDItype s(0), tempUDItype(0);

        // Compute t and s such that p - 1 = 2^s * t , where t is odd
        pp = a.p - 1;
        div(pp, tempUDItype, pp, UDItype(2));

        while(tempUDItype == 0) {
            s += 1;
            copy(t, pp);
            div(pp, tempUDItype, pp, UDItype(2));
        }

        if (s == 1) {
            exp(r, a, (a.p + 1)/4);
        } else if (s == 2) {
            exp(tempmod, a, (a.p-1)/4);
            if (tempmod == 1) {
                exp(r, a, (a.p + 3)/8);
            } else {
                exp(ppmod, 4*a, (a.p - 5)/8);
                r = 2*a*ppmod;
            }
        } else if (s==3) {
            rmint<K, MG> S, d, h, z, i;
            exp(S, 2*a, (a.p - 1)/4);
            if (S == 1) {
                do { rand(d); }
                while (is_quadratic_residue(d));
            } else {
                do { rand(d); }
                while (!is_quadratic_residue(d));
            }
            mul(h, d, d);
            h *= 2*a;
            exp(z, h, (a.p - 9)/16);
            mul(i, z, z);
            i *= h;
            r = (i-1)*z*d*a;
        } else {
            rmint<K, MG> d, z, b, y;
            UDItype s1;

            do { rand(d); }
            while (is_quadratic_residue(d));
            exp(z, d, t);
            exp(r, a, (t+1)/2);
            exp(b, a, t);

            do {
                s1 = 0;
                copy(tempmod, b);
                while (tempmod != 1) {
                    mul(tempmod, tempmod, tempmod);
                    s1 += 1;
                }

                if (s1 != 0) {
                    exp(y, z, (UDItype)(1 << (s-s1-1)));
                    r *= y;
                    mul(z, y, y);
                    b *= z;
                    s = s1;
                }
            } while (s1 != 0);
        }
    }
}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
