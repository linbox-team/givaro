/* ruint/shift.h - Shift functions for ruint

   Copyright Université Grenoble Alpes
Contributors :
Alexis BREUST (alexis.breust@gmail.com 2014)
Jean-Guillaume DUMAS

Time-stamp: <20 Jun 12 10:28:28 Jean-Guillaume.Dumas@imag.fr>

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


#ifndef RUINT_SHIFT_H
#define RUINT_SHIFT_H

#include "ruruint.h"

// --------------------------------------------------------------
// ----------------------- DEFINTIONS ---------------------------

namespace RecInt
{
    template <size_t K, typename T> ruint<K> operator<<(const ruint<K>&, const T&);
    template <size_t K, typename T> ruint<K> operator>>(const ruint<K>&, const T&);

    template <size_t K, typename T> ruint<K>& operator<<=(ruint<K>&, const T&);
    template <size_t K, typename T> ruint<K>& operator>>=(ruint<K>&, const T&);

    // a = b << c
    template <size_t K, typename T> ruint<K>& left_shift(ruint<K>& a, const ruint<K>& b, const T& c);

    // a = b << 1   (z is the lost bit)
    template <size_t K> ruint<K>& left_shift_1(bool& z, ruint<K>& a, const ruint<K>& b);
    template <size_t K> ruint<K>& left_shift_1(ruint<K>& a, const ruint<K>& b);

    // a = b >> c
    template <size_t K, typename T> ruint<K>& right_shift(ruint<K>& a, const ruint<K>& b, const T& c);

    // a = b >> 1   (z is the lost bit)
    template <size_t K> ruint<K>& right_shift_1(bool& z, ruint<K>& a, const ruint<K>& b);
    template <size_t K> ruint<K>& right_shift_1(ruint<K>& a, const ruint<K>& b);

    // Internal use
    template <size_t K, typename T> void left_shift(ruint<K+1>& b, const ruint<K>& a, const T& m);
}


// --------------------------------------------------------------
// ----------------------- Operators ----------------------------

namespace RecInt
{
    // Left shift operator
    template <size_t K, typename T>
    inline ruint<K> operator<<(const ruint<K>& a, const T& d) {
        ruint<K> c;
        return left_shift(c, a, d);
    }

    template <size_t K, typename T>
    inline ruint<K>& operator<<=(ruint<K>& b, const T& c) {
        ruint<K> bp(b);
        return left_shift(b, bp, c);
    }

    // Right shift operator
    template <size_t K, typename T>
    inline ruint<K> operator>>(const ruint<K>& a, const T& d) {
        ruint<K> c;
        return right_shift(c, a, d);
    }

    template <size_t K, typename T>
    inline ruint<K>& operator>>=(ruint<K>& b, const T& c) {
        ruint<K> bp(b);
        return right_shift(b, bp, c);
    }
}


// --------------------------------------------------------------
// ----------------------- Left shift ---------------------------

namespace RecInt
{
    // b = a << d
    template <size_t K, typename T>
    inline ruint<K>& left_shift(ruint<K>& b, const ruint<K>& a, const T& d) {
        const DItype defect((DItype)NBBITS<K-1>::value - (DItype)d);

        if (d == 0) {
            copy(b, a);
        } else if (d == 1) {
            left_shift_1(b, a);
        } else if (d > T(NBBITS<K>::value)) {
            reset(b);
        } else if (defect > 0) {
            // b = a << d <=> b.High = (a.High << d) + (a.Low >> (NBBITS<K-1> - d))
            // and b.Low = a.Low << d
            ruint<K-1> ahd;
            left_shift(ahd, a.High, d);
            right_shift(b.High, a.Low, defect);
            b.High |= ahd;
            left_shift(b.Low, a.Low, d);
        } else if (defect < 0) {
            // b = a << d <=> b.High = a.Low << (d - NBBITS<K-1>)
            // and b.Low = 0
            left_shift(b.High, a.Low, -defect);
            reset(b.Low);
        } else {
            // b = a << d <=> b.High = a.Low
            // and b.Low = 0
            copy(b.High, a.Low);
            reset(b.Low);
        }

        return b;
    }

    template <typename T>
    inline ruint<__RECINT_LIMB_SIZE>& left_shift(ruint<__RECINT_LIMB_SIZE>& b, const ruint<__RECINT_LIMB_SIZE>& a, const T& d) {
        if (d == 0) b.Value = a.Value;
        else if (d < __RECINT_LIMB_BITS)  b.Value = a.Value << d;
        else b.Value = 0;
        return b;
    }

    // b = a << 1   (z is the lost bit)
    template <size_t K>
    inline ruint<K>& left_shift_1(bool& z, ruint<K>& b, const ruint<K>& a) {
        bool zl;

        left_shift_1(z, b.High, a.High);
        left_shift_1(zl, b.Low, a.Low);
        if (zl) set_lowest_bit(b.High);

        return b;
    }

    template <>
    inline ruint<__RECINT_LIMB_SIZE>& left_shift_1(bool& z, ruint<__RECINT_LIMB_SIZE>& b, const ruint<__RECINT_LIMB_SIZE>& a) {
        z = (a.Value & __RECINT_MAXPOWTWO);
        b.Value = a.Value << 1;
        return b;
    }

    template <size_t K>
    inline ruint<K>& left_shift_1(ruint<K>& b, const ruint<K>& a) {
        bool z;
        return left_shift_1(z, b, a);
    }
}


// --------------------------------------------------------------
// ---------------------- Right shift ---------------------------

namespace RecInt
{
    // b = a >> d
    template <size_t K, typename T>
    inline ruint<K>& right_shift(ruint<K>& b, const ruint<K>& a, const T& d) {
        const DItype defect((DItype)NBBITS<K-1>::value - (DItype)d);

        if (d == 0) {
            copy(b, a);
        } else if (d == 1) {
            right_shift_1(b, a);
        } else if (d > T(NBBITS<K>::value)) {
            reset(b);
        } else if (defect > 0) {
            // b = a >> d <=> b.High = a.High >> d
            // and b.Low = (a.Low >> d) + (a.High << (NBBITS<K-1> - d))
            ruint<K-1> ald;
            right_shift(ald, a.Low, d);
            left_shift(b.Low, a.High, defect);
            right_shift(b.High, a.High, d);
            b.Low |= ald;
        } else if (defect < 0) {
            // b = a >> d <=> b.High = 0
            // and b.Low = a.High >> (d - NBBITS<K-1>)
            right_shift(b.Low, a.High, -defect);
            reset(b.High);
        } else {
            // b = a >> d <=> b.High = 0
            // and b.Low = a.High
            copy(b.Low, a.High);
            reset(b.High);
        }

        return b;
    }

    template <typename T>
    inline ruint<__RECINT_LIMB_SIZE>& right_shift(ruint<__RECINT_LIMB_SIZE>& b, const ruint<__RECINT_LIMB_SIZE>& a, const T& d) {
        if (d == 0) b.Value = a.Value;
        else if (d < __RECINT_LIMB_BITS)  b.Value = a.Value >> d;
        else b.Value = 0;
        return b;
    }

    // b = a >> 1 (z is the lost bit)
    template <size_t K>
    inline ruint<K>& right_shift_1(bool& z, ruint<K>& b, const ruint<K>& a) {
        bool zh;

        right_shift_1(zh, b.High, a.High);
        right_shift_1(z, b.Low, a.Low);
        if (zh) set_highest_bit(b.Low);

        return b;
    }

    template <>
    inline ruint<__RECINT_LIMB_SIZE>& right_shift_1(bool& z, ruint<__RECINT_LIMB_SIZE>& b, const ruint<__RECINT_LIMB_SIZE>& a) {
        z = (a.Value & 1);
        b.Value = a.Value >> 1;
        return b;
    }

    template <size_t K>
    inline ruint<K>& right_shift_1(ruint<K>& b, const ruint<K>& a) {
        bool z;
        return right_shift_1(z, b, a);
    }
}


// --------------------------------------------------------------
// ------------------------- Extra ------------------------------

namespace RecInt
{
    // b = a << d
    template <size_t K, typename T>
    inline void left_shift(ruint<K+1>& b, const ruint<K>& a, const T& d) {
        const DItype defect((DItype)NBBITS<K>::value - (DItype)d);

        if (d == 0) {
            copy(b.Low, a);
            reset(b.High);
        } else if (d > NBBITS<K+1>::value) {
            reset(b);
        }  else if (defect > 0) {
            // b = a << d <=> b.High = (a >> (NBBITS<K> - d))
            // and b.Low = a << d
            right_shift(b.High, a, defect);
            left_shift(b.Low, a, d);
        } else if (defect < 0) {
            // b = a << d <=> b.High = a << (d - NBBITS<K>)
            // and b.Low = 0
            left_shift(b.High, a, -defect);
            reset(b.Low);
        } else {
            // b = a << d <=> b.High = a
            // and b.Low = 0
            copy(b.High, a);
            reset(b.Low);
        }
    }
}

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
