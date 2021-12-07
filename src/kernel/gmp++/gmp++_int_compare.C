// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_compare.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama, T. Gautier
// $Id: gmp++_int_compare.C,v 1.6 2009-09-17 14:28:22 jgdumas Exp $
// ==========================================================================

/*! @file gmp++/gmp++_int_compare.C
 * @brief routines to compare integers.
 */

#ifndef __GIVARO_gmpxx_gmpxx_int_compare_C
#define __GIVARO_gmpxx_gmpxx_int_compare_C

#ifndef __GIVARO_INLINE_ALL
#include "gmp++/gmp++.h"
#endif
#include <cstdlib>

namespace Givaro {

    int32_t compare(const Integer &a, const Integer& b)
    {
        return mpz_cmp ( (mpz_srcptr)&a.gmp_rep, (mpz_srcptr)&b.gmp_rep );
    }

    // absCompare
    int32_t absCompare(const Integer &a, const Integer &b)
    {
        return mpz_cmpabs( (mpz_srcptr)&(a.gmp_rep), (mpz_srcptr)&(b.gmp_rep));
    }

    int32_t absCompare(const Integer &a, const double b)
    {
        return mpz_cmpabs_d( (mpz_srcptr)&(a.gmp_rep), b);
    }

    int32_t absCompare(const Integer &a, const float b)
    {
        return mpz_cmpabs_d( (mpz_srcptr)&(a.gmp_rep), (double)b);
    }

    int32_t absCompare(const Integer &a, const uint64_t b)
    {
#if __GIVARO_SIZEOF_LONG < 8
        return absCompare( a, Integer(b));
#else
        return mpz_cmpabs_ui( (mpz_srcptr)&(a.gmp_rep), b);
#endif
    }

    int32_t absCompare(const Integer &a, const uint32_t b)
    {
#if __GIVARO_SIZEOF_LONG == 4
        return mpz_cmpabs_ui( (mpz_srcptr)&(a.gmp_rep), b);
#else
        return mpz_cmpabs_ui( (mpz_srcptr)&(a.gmp_rep), (uint64_t)b);
#endif
    }

    int32_t absCompare(const Integer &a, const int64_t b)
    {
#if __GIVARO_SIZEOF_LONG < 8
        return absCompare( a, Integer(b));
#else
        return mpz_cmpabs_ui( (mpz_srcptr)&(a.gmp_rep), (uint64_t) std::abs(b));
#endif
    }

    int32_t absCompare(const Integer &a, const int32_t b)
    {
#if __GIVARO_SIZEOF_LONG == 4
        return mpz_cmpabs_ui( (mpz_srcptr)&(a.gmp_rep), std::abs(b));
#else
        return mpz_cmpabs_ui( (mpz_srcptr)&(a.gmp_rep), (uint64_t) std::abs(b));
#endif
    }

    // Operator !=
    int32_t Integer::operator != (const Integer & l) const
    {
        return mpz_cmp((mpz_srcptr)&gmp_rep,  (mpz_srcptr)l.get_mpz_const()) != 0;
    }

    int32_t Integer::operator != (const double l) const
    {
        return mpz_cmp_d((mpz_srcptr)&gmp_rep,  l) != 0;
    }

    int32_t Integer::operator != (const float l) const
    {
        return mpz_cmp_d((mpz_srcptr)&gmp_rep, (float)  l) != 0;
    }

    int32_t Integer::operator != (const int32_t l) const
    {
        return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) != 0;
    }

    int32_t Integer::operator != (const uint32_t l) const
    {
#if __GIVARO_SIZEOF_LONG == 4
        return mpz_cmp_ui((mpz_srcptr)&gmp_rep, l) != 0;
#else
        return mpz_cmp_ui((mpz_srcptr)&gmp_rep, (uint64_t) l) != 0;
#endif
    }

    int32_t Integer::operator != (const int64_t l) const
    {
#if __GIVARO_SIZEOF_LONG < 8
        return this->operator != (Integer(l));
#else
        return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) != 0;
#endif
    }

    int32_t Integer::operator != (const uint64_t l) const
    {
#if __GIVARO_SIZEOF_LONG < 8
        return this->operator != (Integer(l));
#else
        return mpz_cmp_ui((mpz_srcptr)&gmp_rep, l) != 0;
#endif
    }

    //-------------------------------------------------inline comparaison operators
    int32_t operator != (double l, const Integer& n)
    {
        return n.operator != (l);
    }

    int32_t operator != (float l, const Integer& n)
    {
        return n.operator != (l);
    }

    int32_t operator != (int32_t l, const Integer& n)
    {
        return n.operator != (l);
    }

    int32_t operator != (int64_t l, const Integer& n)
    {
        return n.operator != (l);
    }

    int32_t operator != (uint64_t l, const Integer& n)
    {
        return n.operator != (l);
    }

    int32_t operator != (uint32_t l, const Integer& n)
    {
        return n.operator != (l);
    }

    // operator ==
    int32_t Integer::operator == (const Integer & l) const
    {
        return !this->operator!=(l);
    }

    int32_t Integer::operator == (const double l) const
    {
        return !this->operator!=(l);
    }

    int32_t Integer::operator == (const float l) const
    {
        return !this->operator!=(l);
    }

    int32_t Integer::operator == (const int32_t l) const
    {
        return !this->operator!=(l);
    }

    int32_t Integer::operator == (const uint32_t l) const
    {
        return !this->operator!=(l);
    }

    int32_t Integer::operator == (const int64_t l) const
    {
        return !this->operator!=(l);
    }

    int32_t Integer::operator == (const uint64_t l) const
    {
        return !this->operator!=(l);
    }

    //-------------------------------------------------inline comparaison operators
    int32_t operator == (double l, const Integer& n)
    {
        return n.operator == (l);
    }

    int32_t operator == (float l, const Integer& n)
    {
        return n.operator == (l);
    }

    int32_t operator == (int32_t l, const Integer& n)
    {
        return n.operator == (l);
    }

    int32_t operator == (int64_t l, const Integer& n)
    {
        return n.operator == (l);
    }

    int32_t operator == (uint64_t l, const Integer& n)
    {
        return n.operator == (l);
    }

    int32_t operator == (uint32_t l, const Integer& n)
    {
        return n.operator == (l);
    }

    // Operator >
    int32_t Integer::operator > (const Integer & l) const
    {
        return mpz_cmp((mpz_srcptr)&gmp_rep,  (mpz_srcptr)l.get_mpz_const()) > 0;
    }

    int32_t Integer::operator > (const double l) const
    {
        return mpz_cmp_d((mpz_srcptr)&gmp_rep,  l) > 0;
    }

    int32_t Integer::operator > (const float l) const
    {
        return mpz_cmp_d((mpz_srcptr)&gmp_rep, (float)  l) > 0;
    }

    int32_t Integer::operator > (const int32_t l) const
    {
        return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) > 0;
    }

    int32_t Integer::operator > (const uint32_t l) const
    {
#if __GIVARO_SIZEOF_LONG == 4
        return mpz_cmp_ui((mpz_srcptr)&gmp_rep, l) > 0;
#else
        return mpz_cmp_ui((mpz_srcptr)&gmp_rep, (uint64_t) l) > 0;
#endif
    }

    int32_t Integer::operator > (const int64_t l) const
    {
#if __GIVARO_SIZEOF_LONG == 8
        return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) > 0;
#else
        return this->operator > (Integer(l));
#endif
    }

    int32_t Integer::operator > (const uint64_t l) const
    {
#if __GIVARO_SIZEOF_LONG == 8
        return mpz_cmp_ui((mpz_srcptr)&gmp_rep, l) > 0;
#else
        return this->operator > (Integer(l));
#endif
    }

    //-------------------------------------------------inline comparaison operators
    int32_t operator > (double l, const Integer& n)
    {
        return n.operator < (l);
    }

    int32_t operator > (float l, const Integer& n)
    {
        return n.operator < (l);
    }

    int32_t operator > (int32_t l, const Integer& n)
    {
        return n.operator < (l);
    }

    int32_t operator > (int64_t l, const Integer& n)
    {
        return n.operator < (l);
    }

    int32_t operator > (uint64_t l, const Integer& n)
    {
        return n.operator < (l);
    }

    int32_t operator > (uint32_t l, const Integer& n)
    {
        return n.operator < (l);
    }

    // Operator <
    int32_t Integer::operator < (const Integer & l) const
    {
        return mpz_cmp((mpz_srcptr)&gmp_rep,  (mpz_srcptr)l.get_mpz_const()) < 0;
    }

    int32_t Integer::operator < (const double l) const
    {
        return mpz_cmp_d((mpz_srcptr)&gmp_rep,  l) < 0;
    }

    int32_t Integer::operator < (const float l) const
    {
        return mpz_cmp_d((mpz_srcptr)&gmp_rep, (float) l) < 0;
    }

    int32_t Integer::operator < (const uint32_t l) const
    {
#if __GIVARO_SIZEOF_LONG == 4
        return mpz_cmp_ui((mpz_srcptr)&gmp_rep, l) < 0;
#else
        return mpz_cmp_ui((mpz_srcptr)&gmp_rep, (uint64_t) l) < 0;
#endif
    }

    int32_t Integer::operator < (const uint64_t l) const
    {
#if __GIVARO_SIZEOF_LONG == 8
        return mpz_cmp_ui((mpz_srcptr)&gmp_rep, l) < 0;
#else
        return this->operator < (Integer(l));
#endif
    }

    int32_t Integer::operator < (const int32_t l) const
    {
        return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) < 0;
    }

    int32_t Integer::operator < (const int64_t l) const
    {
#if __GIVARO_SIZEOF_LONG == 8
        return mpz_cmp_si((mpz_srcptr)&gmp_rep, l) < 0;
#else
        return this->operator < (Integer(l));
#endif
    }

    //-------------------------------------------------inline comparaison operators
    int32_t operator < (double l, const Integer& n)
    {
        return n.operator > (l);
    }

    int32_t operator < (float l, const Integer& n)
    {
        return n.operator > (l);
    }

    int32_t operator < (int32_t l, const Integer& n)
    {
        return n.operator > (l);
    }

    int32_t operator < (int64_t l, const Integer& n)
    {
        return n.operator > (l);
    }

    int32_t operator < (uint64_t l, const Integer& n)
    {
        return n.operator > (l);
    }

    int32_t operator < (uint32_t l, const Integer& n)
    {
        return n.operator > (l);
    }

    // Operator >=
    int32_t Integer::operator >= (const Integer & l) const
    {
        return !this->operator<(l);
    }

    int32_t Integer::operator >= (const double l) const
    {
        return !this->operator<(l);
    }

    int32_t Integer::operator >= (const float l) const
    {
        return !this->operator<(l);
    }

    int32_t Integer::operator >= (const int32_t l) const
    {
        return !this->operator<(l);
    }

    int32_t Integer::operator >= (const uint32_t l) const
    {
        return !this->operator<(l);
    }

    int32_t Integer::operator >= (const int64_t l) const
    {
        return !this->operator<(l);
    }

    int32_t Integer::operator >= (const uint64_t l) const
    {
        return !this->operator<(l);
    }

    //-------------------------------------------------inline comparaison operators
    int32_t operator >= (double l, const Integer& n)
    {
        return n.operator <= (l);
    }

    int32_t operator >= (float l, const Integer& n)
    {
        return n.operator <= (l);
    }

    int32_t operator >= (int32_t l, const Integer& n)
    {
        return n.operator <= (l);
    }

    int32_t operator >= (int64_t l, const Integer& n)
    {
        return n.operator <= (l);
    }

    int32_t operator >= (uint64_t l, const Integer& n)
    {
        return n.operator <= (l);
    }

    int32_t operator >= (uint32_t l, const Integer& n)
    {
        return n.operator <= (l);
    }

    // Operator <=
    int32_t Integer::operator <= (const Integer & l) const
    {
        return !this->operator>(l);
    }

    int32_t Integer::operator <= (const double l) const
    {
        return !this->operator>(l);
    }

    int32_t Integer::operator <= (const float l) const
    {
        return !this->operator>(l);
    }

    int32_t Integer::operator <= (const uint32_t l) const
    {
        return !this->operator>(l);
    }

    int32_t Integer::operator <= (const uint64_t l) const
    {
        return !this->operator>(l);
    }

    int32_t Integer::operator <= (const int32_t l) const
    {
        return !this->operator>(l);
    }

    int32_t Integer::operator <= (const int64_t l) const
    {
        return !this->operator>(l);
    }

    //-------------------------------------------------inline comparaison operators
    int32_t operator <= (double l, const Integer& n)
    {
        return n.operator >= (l);
    }

    int32_t operator <= (float l, const Integer& n)
    {
        return n.operator >= (l);
    }

    int32_t operator <= (int32_t l, const Integer& n)
    {
        return n.operator >= (l);
    }

    int32_t operator <= (int64_t l, const Integer& n)
    {
        return n.operator >= (l);
    }

    int32_t operator <= (uint64_t l, const Integer& n)
    {
        return n.operator >= (l);
    }

    int32_t operator <= (uint32_t l, const Integer& n)
    {
        return n.operator >= (l);
    }


    // compare to 1 and 0
    int32_t isOne(const Integer& a)
    {
        return ! mpz_cmp_ui((mpz_srcptr)&(a.gmp_rep), 1U);
    }
    int32_t isMOne(const Integer& a)
    {
        return ! mpz_cmp_si((mpz_srcptr)&(a.gmp_rep), -1);
    }

    int32_t nonZero(const Integer& a)
    {
        return mpz_cmp_ui((mpz_srcptr)&(a.gmp_rep), 0U);
        // return (mpz_sgn((mpz_srcptr)&(a.gmp_rep)) != 0) ; // BB which one is faster ?
    }

    int32_t isZero(const Integer& a)
    {
        return ! mpz_cmp_ui((mpz_srcptr)&(a.gmp_rep), 0U);
        // return (mpz_sgn((mpz_srcptr)&(a.gmp_rep)) == 0) ; // BB which one is faster ?
    }
    int32_t isZero(const int16_t a)
    {
        return a ==0;
    }
    int32_t isZero(const int32_t a)
    {
        return a ==0;
    }
    int32_t isZero(const int64_t a)
    {
        return a ==0;
    }
    int32_t isZero(const uint16_t a)
    {
        return a ==0U;
    }
    int32_t isZero(const uint32_t a)
    {
        return a ==0U;
    }
    int32_t isZero(const uint64_t a)
    {
        return a ==0U;
    }

}

#endif // __GIVARO_gmpxx_gmpxx_int_compare_C

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
