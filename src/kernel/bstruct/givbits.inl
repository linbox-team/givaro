// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/bstruct/givbits.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Author: T. Gautier
// $Id: givbits.inl,v 1.3 2011-02-02 16:23:55 briceboyer Exp $
// ==========================================================================

#ifndef __GIVARO_bits_INL
#define __GIVARO_bits_INL
namespace Givaro {
    // -- Copy operators
    inline
    Bits& Bits::copy( const Bits& src )
    { rep.copy( src.rep ); return *this; }

    inline
    Bits& Bits::logcopy( const Bits& src )
    { rep.copy( src.rep ); return *this; }

    inline Bits& Bits::operator= (const Bits& B) { return copy(B); }

    //-------------------------------------------------inline << operators
    inline std::ostream& operator<< (std::ostream& o, const Bits& a)
    { return a.print(o); }
} // namespace Givaro
#endif // __GIVARO_bits_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
