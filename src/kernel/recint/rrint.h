/* rint/rint.h - Class definition of rint<K> from RecInt library

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

#ifndef RINT_RINT_H
#define RINT_RINT_H

#include "recdefine.h"
#include "ruruint.h"
#include "givaro/givtypestring.h"
#include "rumanip.h" // ms_limb
#include "rutools.h" // mod_n

// --------------------------------------------------------------
// ---------------- Declaration of class rint ------------------

namespace RecInt
{
    /* Basic definition of rint */
    template <size_t K> class rint {
    public:
        // A rint is stored as a ruint
        ruint<K> Value;

        // Constructors
        rint() {}
        rint(const rint<K>& r) : Value(r.Value) {}
        rint(const rint<K-1>& rl) : Value(rl.Value)
        { if (rl < 0) { Value.Low = -Value.Low; Value = -Value; } }
        rint(const ruint<K>& r) : Value(r) {}
        template <typename T> rint(const T& b) : Value(b) {}

        // type_string
        static const std::string type_string () {
            return "RecInt::rint<" + std::to_string(K)+ ">";
        }

        // Cast
        template <typename T> operator T() const { return static_cast<T>(Value); }

        // Quick sign evaluation
        // *this < 0
        inline bool isNegative() const { return ms_limb(Value) & __RECINT_MAXPOWTWO; }
        // *this >= 0
        inline bool isPositive() const { return !(isNegative()); }

		// max Cardinality : 2^( (2^K-1) / 2 )
        static rint<K> maxCardinality() {
            rint<K> max( ruint<K>::maxFFLAS() );
           return max;
        }

		// max Elements
        static rint<K> maxElement() {
            return ruint<K>::maxElement()/2;
        }

		// max Cardinality for fflas-ffpack : supports (a*b+c*d)
        // 2^(2^(K-1)-1)
        static rint<K> maxFFLAS();
    };

    using rint64 =  rint<6>;
    using rint128 = rint<7>;
    using rint256 = rint<8>;
    using rint512 = rint<9>;


}

namespace RecInt
{
    // a = 0
    template <size_t K> inline void reset(rint<K>& a) {
        reset(a.Value);
    }

        // a = b
    template <size_t K> inline void copy(rint<K>& a, const rint<K>& b) {
        copy(a.Value, b.Value);
    }

}

namespace std
{
    template <size_t K> struct make_signed<RecInt::rint<K>> {
        typedef RecInt::rint<K> type;
    };
    template <size_t K> struct make_signed<RecInt::ruint<K>> {
        typedef RecInt::rint<K> type;
    };

}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
