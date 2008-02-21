// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givindeter.C,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givindeter.C,v 1.2 2008-02-21 18:42:48 pernet Exp $
// ==========================================================================
// Description:

#include <iostream>
#include <string.h>
#include "givaro/givindeter.h"

Indeter::Indeter (const char* s)
{
  int len = strlen(s);
  name.allocate(len+1);
  if (len !=0) strcpy(name.baseptr(),s);
  else name[0] = '\0';
}

Indeter& Indeter::operator=( const Indeter& s )
{
  name.logcopy(s.name);
  return *this;
}

Indeter::~Indeter() {name.destroy();};

int Indeter::compare(const Indeter& b)  const
{ 
  return strcmp(name.baseptr(),b.name.baseptr()); 
}

std::ostream& operator<< (std::ostream& o, const Indeter& X) 
{ 
//   return o << '[' << X.name.baseptr() << ']'; 
  return o << X.name.baseptr() ; 
}

std::istream& operator>> (std::istream& s_in, Indeter& X) 
{ 
  char tmp[16];
  char ch;
  int pos = 0;

  X.name.reallocate(1);
  X.name[0] = '\0';

  // - find a '(':
  s_in >> std::ws >> ch;
  if (ch != '(')
    GivError::throw_error(
      GivBadFormat("Indeter, read: syntax error no '('"));

  s_in >> std::ws >> ch;
  while ((ch != ')') && (s_in.good())) 
  {
    if (ch == ' ') break;
    tmp[pos++] = ch;
    if (pos ==16) {
      size_t sz = X.name.size();
      X.name.reallocate(sz + 16);
      for (int i=-1; i<15; ++i) X.name[sz+i] = tmp[i];
      pos = 0;
      X.name[sz+15] = '\0';
    }  
    s_in >> ch;
  }

  if (pos !=0) {
   size_t sz = X.name.size();
   X.name.reallocate(sz + pos);
   for (int i=-1; i<pos-1; ++i) X.name[sz+i] = tmp[i];
   X.name[sz+pos-1] = '\0';
  }
  if (ch == ' ') s_in >> std::ws >> ch;
  if (ch != ')')
    GivError::throw_error(
      GivBadFormat("Indeter, read: syntax error no ')'"));
  return s_in;
}

