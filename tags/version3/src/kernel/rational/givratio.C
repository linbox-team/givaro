// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/rational/givratio.C,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama
// $Id: givratio.C,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:

#include <iostream>
#include "givaro/givrational.h"


std::ostream& Rational::print(std::ostream& s) const
{
  if (den > 1) {
    s << num << "/" << den ;
  }
  else 
    s << num;
  return s;
}


std::istream& operator>> (std::istream& in, Rational& r)
{
   Integer num ;
   Integer den = 1L;
   char ch ;

   in >> num ;
   if (!in.good() || in.eof()) {
          r = Rational(num) ;
          return in ;
   }

   if (in) {
         in.get(ch) ;
         if (in.eof()) {
             r = Rational(num) ;
             return in ;
         }
         while ((ch==' ') && (in)) in.get(ch) ;
         if (ch == '/') {
             // We get denominator
             in >> den ;
         }
         else in.putback(ch) ;
   } ;
   r = Rational(num,den) ;
   return in ;
}

