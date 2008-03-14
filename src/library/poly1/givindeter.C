// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givindeter.C,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givindeter.C,v 1.3 2008-03-14 21:32:15 pernet Exp $
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
   s_in>>X.name;
 }

