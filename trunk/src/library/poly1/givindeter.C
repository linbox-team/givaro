// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givindeter.C,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givindeter.C,v 1.4 2008-07-01 15:40:33 jgdumas Exp $
// ==========================================================================
// Description:

#include <iostream>
#include <string.h>
#include "givaro/givindeter.h"


Indeter& Indeter::operator=( const Indeter& s )
{
  name = s.name;
  return *this;
}

int Indeter::compare(const Indeter& b)  const
{ 
  return name.compare(b.name); 
}

std::ostream& operator<< (std::ostream& o, const Indeter& X) 
{ 
//   return o << '[' << X.name.baseptr() << ']'; 
  return o << X.name ; 
}

 std::istream& operator>> (std::istream& s_in, Indeter& X) 
 { 
   return s_in>>X.name;
 }

