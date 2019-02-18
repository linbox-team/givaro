// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Original authors (LinBox): B. Hovinen, JG Dumas, C. Pernet
// Imported and adapted by: A. Breust
// ==========================================================================

#ifndef __Givaro_field_gf2_INL
#define __Givaro_field_gf2_INL

namespace Givaro
{
    inline std::ostream& GF2::write (std::ostream& os) const
    {
	    return os << "GF2";
    }


    inline std::istream& GF2::read (std::istream& is)
    {
	    return is;
    }


    inline std::ostream& GF2::write (std::ostream &os, const Element& x) const
    {
	    return os << x;
    }


   inline  std::istream& GF2::read (std::istream& is, Element& x) const
    {
        is >> x;
        return is;
    }

    inline std::istream& GF2::read (std::istream& is, BitReference x) const
    {
        bool a;
	    is >> a;
	    x = a;
	    return is;
    }

    inline GF2::Element& GF2::add (Element& x, const Element& y, const Element& z) const
    {
	    return x = y ^ z;
    }

    inline GF2::BitReference GF2::add (BitReference x, const Element& y, const Element& z) const
    {
	    return x = y ^ z;
    }

    inline GF2::Element& GF2::sub (Element& x, const Element& y, const Element& z) const
    {
	    return x = y ^ z;
    }

    inline GF2::BitReference GF2::sub (BitReference x, const Element& y, const Element& z) const
    {
	    return x = y ^ z;
    }

    inline GF2::Element& GF2::mul (Element& x, const Element& y, const Element& z) const
    {
	    return x = y & z;
    }

    inline GF2::BitReference GF2::mul (BitReference x, const Element& y, const Element& z) const
    {
	    return x = y & z;
    }

    inline GF2::Element& GF2::div (Element& x, const Element& y, const Element& z ) const
    {
        GIVARO_ASSERT_MATHDIV0(z, "Error: division by zero in GF2::div");
        return x = y;
    }

    inline GF2::BitReference GF2::div (BitReference x, const Element& y, const Element& z ) const
    {
        GIVARO_ASSERT_MATHDIV0(z, "Error: division by zero in GF2::div");
        return x = y;
    }

    inline GF2::Element& GF2::neg (Element& x, const Element& y) const
    {
	    return x = y;
    }

    inline GF2::BitReference GF2::neg (BitReference x, const Element& y) const
    {
	    return x = y;
    }

    inline GF2::Element& GF2::inv (Element& x, const Element& y) const
    {
        GIVARO_ASSERT_MATHDIV0(y, "Error: division by zero in GF2::inv");
        return x = y;
    }
    inline GF2::BitReference GF2::inv (BitReference x, const Element& y) const
    {
        GIVARO_ASSERT_MATHDIV0(y, "Error: division by zero in GF2::inv");
        return x = y;
    }

    inline GF2::BitReference GF2::axpy (BitReference r, const Element& a, const Element& x, const Element& y) const
    {
	    return r = (a & x) ^ y;
    }

    inline GF2::Element& GF2::axpy (Element& r, const Element& a, const Element& x, const Element& y) const
    {
	    return r = (a & x) ^ y;
    }

    inline GF2::BitReference GF2::axmy (BitReference r, const Element& a, const Element& x, const Element& y) const
    {
	    return r = (a & x) ^ y;
    }

    inline GF2::Element& GF2::axmy (Element& r, const Element& a, const Element& x, const Element& y) const
    {
	    return r = (a & x) ^ y;
    }

    inline GF2::BitReference GF2::maxpy (BitReference r, const Element& a, const Element& x, const Element& y) const
    {
	    return r = (a & x) ^ y;
    }

    inline GF2::Element& GF2::maxpy (Element& r, const Element& a, const Element& x, const Element& y) const
    {
	    return r = (a & x) ^ y;
    }

    inline GF2::Element& GF2::addin (Element& x, const Element& y) const
    {
	    return x ^= y;
    }

    inline GF2::BitReference GF2::addin (BitReference x, const Element& y) const
    {
	    return x = x ^ y;
    }

    inline GF2::Element& GF2::subin (Element& x, const Element& y) const
    {
	    return x ^= y;
    }

    inline GF2::BitReference GF2::subin (BitReference x, const Element& y) const
    {
	    return x = x ^ y;
    }

    inline GF2::Element& GF2::mulin (Element& x, const Element& y) const
    {
	    return x &= y;
    }

    inline GF2::BitReference GF2::mulin (BitReference x, const Element& y) const
    {
	    return x = (bool)x & y;
    }

    inline GF2::Element& GF2::divin (Element& x, const Element& y ) const
    {
        GIVARO_ASSERT_MATHDIV0(y, "Error: division by zero in GF2::divin");
        return x;
    }

    inline GF2::BitReference GF2::divin (BitReference x, const Element& y ) const
    {
        GIVARO_ASSERT_MATHDIV0(y, "Error: division by zero in GF2::divin");
        return x;
    }

    inline GF2::Element& GF2::negin (Element& x) const
    {
	    return x;
    }

    inline GF2::BitReference GF2::negin (BitReference x) const
    {
	    return x;
    }

    inline GF2::Element& GF2::invin (Element& x) const
    {
        GIVARO_ASSERT_MATHDIV0(x, "Error: division by zero in GF2::invin");
        return x;
    }

    inline GF2::BitReference GF2::invin (BitReference x) const
    {
        GIVARO_ASSERT_MATHDIV0(x, "Error: division by zero in GF2::invin");
        return x;
    }

    inline GF2::Element& GF2::axpyin (Element& r, const Element& a, const Element& x) const
    {
	    return r ^= a & x;
    }

    inline GF2::BitReference GF2::axpyin (BitReference r, const Element& a, const Element& x) const
    {
	    return r = r ^ (a & x);
    }

    inline GF2::Element& GF2::axmyin (Element& r, const Element& a, const Element& x) const
    {
	    return r ^= a & x;
    }

    inline GF2::BitReference GF2::axmyin (BitReference r, const Element& a, const Element& x) const
    {
	    return r = r ^ (a & x);
    }

    inline GF2::Element& GF2::maxpyin (Element& r, const Element& a, const Element& x) const
    {
	    return r ^= a & x;
    }

    inline GF2::BitReference GF2::maxpyin (BitReference r, const Element& a, const Element& x) const
    {
	    return r = r ^ (a & x);
    }
}

#endif


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
