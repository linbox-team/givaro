// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/rational/givratcstor.C,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama
// $Id: givratcstor.C,v 1.2 2005-06-14 14:53:14 pernet Exp $
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

//!!!Integer or unsigned for n and d ?
//Might be better to use Rational(char *) for full prec. input, have have this 
//sufficently fast for numerous initialisations.
Rational::Rational(double x)
//snarfed from netlib (file frac.c)
//Result is already reduced
{
  throw GivError("Rational::Rational(double x) : undefined") ;
/*
  // The following values are used to test for too small
  // and too large values of v for an integer fraction to
  // represent.
  
  const double MAX = DBL_MAX;
  const double MIN = DBL_MIN;
  
  // we assume a base 2 representation for doubles.
  static const double error = 1.0 / power(2.0, DSIGNIF - 1);
  long n = 0, d = 1;
  int sgn = (x > 0) ? 1 : -1;
  if (x < 0)
    x = -x;
  
  int D, N, t;
  double epsilon, r , m;

  if (x < MIN || x > MAX || error < 0.0)
    {
      num = Integer(n);
      den = Integer(d);
      return;
    }
  
  d = D = 1;
  n = int(x);
  N = n + 1;
  r = 0.0 ;
  goto three;
  
 one:   
  if (r > 1.0)
    goto two;
  r = 1.0 / r;
 
 two:  
  N += n *int(r);
  D += d * int(r);
  n += N;
  d += D;
 
 three:  
  r = 0.0;
  if ( x * d == double(n) )
    goto four;
  r = (N - x * D) / (x * d - n);
  if (r > 1.0)
    goto four; 
  t = N;
  N = n;
  n = t;
  t = D;
  D = d;
  d = t;

 four:
  epsilon = fabs(1.0 - n / (x * d));
  if (epsilon <= error)
    goto six;
  m = 1.0;
  do {
    m *= 10.0;
  } while (m * epsilon < 1.0);
  epsilon = 1.0/m * (int(0.5 + m * epsilon));
 
 six:
  if (epsilon <= error)
    {
      num = Integer(n * sgn);
      den = Integer(d);
      return;
    }
  if (r != 0.0)
    goto one;
*/
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
// If red == 1 then the rational is reduce (gcd computation!)
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

