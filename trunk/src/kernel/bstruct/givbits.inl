// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/bstruct/givbits.inl,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Author: T. Gautier
// $Id: givbits.inl,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
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

