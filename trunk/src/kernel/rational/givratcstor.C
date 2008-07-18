// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/rational/givratcstor.C,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama
// $Id: givratcstor.C,v 1.3 2008-07-18 12:42:37 jgdumas Exp $
// ==========================================================================
// Description:

#include "givaro/givrational.h"
#include "givaro/giverror.h"
#include "givaro/givpower.h"
#include <math.h>
#include <float.h>
#include <iostream>

#if !defined(__MWERKS__) && (defined(__GNUC__) && (__GNUC__ == 2))
#define __GIVARO_OLD_SSTREAM__
#include <strstream>
#else
// -- new interface for string stream
#include <sstream>
#endif

  // -- Predefined cstes
const Rational Rational::zero = givNoInit();
const Rational Rational::one  = givNoInit();

Rational::ReduceFlag Rational::flags = Rational::Reduce ;
void Rational::SetReduce() { Rational::flags = Rational::Reduce ; }
void Rational::SetNoReduce() { Rational::flags = Rational::NoReduce ; }

// Explicit instanciation
template double power(double x, unsigned int p) ; 



        struct ieee {
#if     __BYTE_ORDER == __BIG_ENDIAN
            uint64 negative:1;
            uint64 exponent:11;
            uint64 mantissa:52;
#endif                          /* Big endian.  */
#if     __BYTE_ORDER == __LITTLE_ENDIAN
            uint64 mantissa:52;
            uint64 exponent:11;
            uint64 negative:1;
#endif                          /* Little endian.  */
        };

Rational::Rational(double x) {
    union { 
        uint64 l; 
        ieee u; 
        double d; 
    } t; // temp

    t.d = x;
    if (t.u.exponent == 0) {
            // Denormal numbers
        num = (x<0.?-t.u.mantissa:t.u.mantissa);
        den = 1;
        *this/=Rational(Integer(1)<<1074);
    } else {
        const long shift = 1075-t.u.exponent;
        t.u.exponent = 1076;
        if (shift > 0) {
            Integer tt( t.u.mantissa+4503599627370496ULL );
            num = (x<0.?-tt:tt);
            den = Integer(1)<<shift;
        } else {
            Integer tt( t.u.mantissa+4503599627370496ULL);
            tt <<=(-shift);
            num = (x<0.?-tt:tt);
            den = 1;
        }
    }
    if (Rational::flags == Rational::Reduce) reduce();
}


//   ------------------------------ Rational(Neutral n )
Rational::Rational(Neutral n ) : den(Integer::one)
{
   if (n == Neutral::zero) {
      num = Integer::zero;
   }
   else { // n = one 
      num = Integer::one;
   }
}

//   ------------------------------ Rational(int n)
Rational::Rational(int n ) : num(n), den(Integer::one) 
{ }


//   ------------------------------ Rational(long n)
Rational::Rational(long n ) : num(n), den(Integer::one) 
{ }

//   ------------------------------ Rational(long n, long d )
Rational::Rational(long n, long d )
{
  if (d == 0)
    {
      throw GivMathDivZero("[Rational::Rational]: null denominator of the rational.") ;
    }
  
  if (n == 0)
    {
      num = Integer::zero;
      den = Integer::one;
    }
  if (d > 0)
    {
      num = Integer(n);
      den = Integer(d);
    }
  else
    {
      num = Integer(-n);
      den = Integer(-d);
    }
    reduce();
}


//   ------------------------------ Rational(const char* s )
Rational::Rational(const char* s ) 
{
#ifdef __GIVARO_OLD_SSTREAM__
  std::istrstream input (s) ;
#else
  std::istringstream input (s) ;
#endif
  Rational r ;
  input >> r ;
  operator= (r) ;
}


//   ------------------------------ Rational(const Integer &n)
Rational::Rational(const Integer &n) : den(Integer::one)
{
  if (isZero(n))
    {
      num = Integer::zero;
    }
  else
    {
      num = n ;
    }
}

// ------------------------------ Rational(const Integer &n, const Integer &d)
// If red == 1 then the rational is reduced (gcd computation!)
Rational::Rational(const Integer &n, const Integer &d, int red)
{
  if (isZero(d))
    {
      throw GivMathDivZero( "[Rational::Rational]: null denominator of the rational.") ;
    }
  
  if (isZero(n))
    {
      num = Integer::zero;
      den = Integer::one;
    }
  if (sign(d) > 0)
    {
      num = n ;
      den = d ;
    }
  else
    {
      num = -n;
      den = -d;
    }
  if (red == 1) reduce();
}

//   ------------------------------ Rational(const Rational &n)
Rational::Rational(const Rational &r) : num(r.num), den(r.den) 
{ } 

Rational::Rational( givNoInit gi) : num(Integer::zero), den(Integer::one) 
{ } 

// ------ Initialization module:
//GivModule Rational::Module (Rational::Init,
//                            Rational::End, 
//                            InitAfter(Integer::Module),
//                            "Rational package") ; 
void Rational::Init(int* argc, char***argv)
{
  new ((Rational*)&Rational::zero.num) Integer(Integer::zero);
  new ((Rational*)&Rational::zero.den) Integer(Integer::one);
  new ((Rational*)&Rational::one.num) Integer(Integer::one);
  new ((Rational*)&Rational::one.den) Integer(Integer::one);
}

void Rational::End()
{}

