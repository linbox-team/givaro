/* ruint/display.h - Display functions for ruint

   Copyright Université Grenoble Alpes
Contributors :
Alexis BREUST (alexis.breust@gmail.com 2014)
Christophe CHABOT (christophechabotcc@gmail.com 2011)
Jean-Guillaume DUMAS

Time-stamp: <19 Jun 12 18:22:16 Jean-Guillaume.Dumas@imag.fr>


This software is a computer program whose purpose is to provide a
fixed precision arithmetic library.

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


#ifndef RUINT_DISPLAY_H
#define RUINT_DISPLAY_H

#include <iostream> /* For streams */
#include <iomanip> /* For setw and so on */

#include "ruruint.h"
#include "ruconvert.h"
#include "rucmp.h"

// --------------------------------------------------------------
// ----------------------- DEFINTIONS ---------------------------

namespace RecInt
{
    // Prints a ruint
    template <size_t K> inline std::ostream& operator<<(std::ostream&, const ruint<K>&);

    // Reads a ruint
    template <size_t K> inline std::istream& operator>>(std::istream&, ruint<K>&);
}


// --------------------------------------------------------------
// ----------------------- ruint print --------------------------

namespace RecInt
{
    // Out stream
    template <size_t K>
    inline std::ostream& operator<<(std::ostream& out, const ruint<K>& a) {
        std::ios_base::fmtflags ff(out.flags());
        if (ff & out.hex) return display_hex(out, a);
        else return display_dec(out, a);
        out.flags(ff);
    }

    template <>
    inline std::ostream& operator<<(std::ostream& out, const ruint<__RECINT_LIMB_SIZE>& a) {
        return out << a.Value;
    }

    // Display the ruint in decimal mode
    template <size_t K>
    inline std::ostream& display_dec(std::ostream& out, const ruint<K>& a) {
        ruint<K> b(a);
        char result[1024];
        limb m(0), ten(10);
        int i;

        if (b == 0) out << '0';

        for (i = 0; b != 0 && i < 1024; i++) {
            div(b, m, b, ten);
            result[i] = char('0' + m);
        }

        for (i--; i >= 0; i--) out << result[i];
        return out;
    }

    // Display the ruint a in hexadecimal mode
    template <size_t K>
    inline std::ostream& display_hex(std::ostream& out, const ruint<K>& a) {
        return display_hex(display_hex(out, a.High), a.Low);
    }

    template <>
    inline std::ostream& display_hex(std::ostream& out, const ruint<__RECINT_LIMB_SIZE>& a) {
        return out << std::setw(__RECINT_LIMB_BITS/4) << std::setfill('0') << a.Value;
    }

    // Display the ruint as raw (usually used for binary streams)
    template <size_t K>
    inline std::ostream& write_raw(std::ostream& out, const ruint<K>& a) {
        return out.write(reinterpret_cast<const char*>(&a), sizeof(a));
    }

    // Reads a ruint
    template <size_t K> inline std::istream& operator>>(std::istream& is, ruint<K>& a) {
        mpz_class g;
        is >> g;
        mpz_to_ruint(a, g);
        return is;
    }

    template <size_t K>
    inline std::istream& read_raw(std::istream& is, ruint<K>& a) {
        return is.read(reinterpret_cast<char*>(&a), sizeof(a));
    }
}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
