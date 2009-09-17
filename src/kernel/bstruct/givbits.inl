// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/bstruct/givbits.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.
// Author: T. Gautier
// $Id: givbits.inl,v 1.2 2009-09-17 14:28:22 jgdumas Exp $
// ==========================================================================

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

