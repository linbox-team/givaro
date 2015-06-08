/* rint/display.h - Display functions for rint

Copyright Université Joseph Fourier - Grenoble
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


#ifndef RINT_DISPLAY_H
#define RINT_DISPLAY_H

#include "rrint.h"
#include "rconvert.h"
#include "rudisplay.h"

// --------------------------------------------------------------
// ----------------------- DEFINTIONS ---------------------------

namespace RecInt
{
    // Prints a rint
    template <size_t K> inline std::ostream& operator<<(std::ostream&, const rint<K>&);
}


// --------------------------------------------------------------
// ----------------------- rint print --------------------------

namespace RecInt
{
    // Out stream
    template <size_t K>
    inline std::ostream& operator<<(std::ostream& out, const rint<K>& a) {
        std::ios_base::fmtflags ff(out.flags());
        if (ff & out.hex) return display_hex(out, a.Value);
        else return display_dec(out, a);
        out.flags(ff);
    }

    // Display the rint in decimal mode
    template <size_t K>
    inline std::ostream& display_dec(std::ostream& out, const rint<K>& a) {
        if (a.isNegative()) {
            out << '-';
            display_dec(out, (-a).Value);
        }
        else {
            display_dec(out, a.Value);
        }
        return out;
    }

    // Reads a rint
    template <size_t K> inline std::istream& operator>>(std::istream& is, rint<K>& a) {
    	mpz_class g;
    	is >> g;
    	mpz_to_rint(a, g);
    	return is;
    }
}

#endif

