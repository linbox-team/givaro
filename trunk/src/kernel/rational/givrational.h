#ifndef _RATIONAL_H_
#define _RATIONAL_H_
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/rational/givrational.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: M. Samama, T. Gautier
// $Id: givrational.h,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:

#include "givaro/givinteger.h"

// ----------------------------------- Class Rational

class Rational {

public :
  // Cstor et dstor
  Rational(Neutral n = Neutral::zero) ;
  Rational(int n) ;
  Rational(long n) ;
  Rational(long n, long d ) ;
  Rational(double x) ;
  Rational(const char* s) ;
  Rational(const Integer& n) ;
  Rational(const Integer& n, const Integer& d, int red = 1 ) ;
    // Rational number reconstruction
  Rational(const Integer& f, const Integer& m, const Integer& k, bool recurs = false ) ;
  Rational(const Rational&) ;
  //~Rational();

  // Predefined cstes
  static const Rational zero ;
  static const Rational one ;

  // Logical and physical copies
  Rational& operator = (const Rational& );
  Rational& logcpy (const Rational& ) ; 
  Rational& copy (const Rational& ) ; 

//------------------Equalities and inequalities between rationals 
  friend int compare(const Rational& a, const Rational& b) ;
  friend int absCompare(const Rational& a, const Rational& b) ;


//----------------Elementary arithmetic between Rational 
  Rational operator + (const Rational& r) const ;
  Rational operator - (const Rational& r) const ;
  Rational operator -() const ;
  Rational operator +() const ;
  Rational operator * (const Rational& r) const ;
  Rational operator / (const Rational &r) const ;
  Rational& operator += (const Rational& r) ;
  Rational& operator -= (const Rational& r) ;
  Rational& operator *= (const Rational& r) ;
  Rational& operator /= (const Rational &r) ;

//-----------------------------------------Arithmetic functions
  friend const Rational pow(const Rational &r, const long l);
  
//-----------------------------------------Miscellaneous
  friend const Integer floor (const Rational &r) ;
  friend const Integer ceil  (const Rational &r) ;
  friend const Integer round (const Rational &r) ;
  friend const Integer trunc (const Rational &r) ;
inline friend const Rational abs  (const Rational &r) ;


  const Integer nume() const ;
  const Integer deno() const ;
inline friend unsigned long length (const Rational& r) ;
inline friend int sign   (const Rational& r) ;
inline friend int iszero (const Rational& r) ;
inline friend int isone  (const Rational& r) ;
inline friend int isinteger(const Rational& r);  

  std::ostream& print ( std::ostream& o ) const ;

inline Rational reduce(const Rational& R) const ;

static void SetReduce() ; 
static void SetNoReduce() ; 
  
protected: // Internal Representation : num/den 
  Integer num, den;

public:
enum ReduceFlag { Reduce = 0x1, NoReduce = 0x0 } ;
protected:
static ReduceFlag flags ;    // flags that indicates is reduction is done or not
                             // by default = Reduce 
  void reduce(void) ;
  // -- module initialization
static void Init(int* argc, char***argv);
static void End();
friend class GivModule;
    // Rational number reconstruction
    bool ratrecon(const Integer& f, const Integer& m, const Integer& k, bool recurs = false ) ;

public:
  // - exportation of the module
static GivModule Module; 
  // -- Cstor for Zero and One to delay initialization after the main
  Rational( givNoInit );
}; // ----------------------------------- End of Class Rationalional 

extern std::istream& operator>> (std::istream& in, Rational& r) ;

#include "givaro/givrational.inl"

//------------------------------------------------------ Class RationalDom
class RationalDom : public Rational {
public:
  typedef Rational Rep;

  // -- Cstor
  RationalDom() : one(1), zero(0) {}

  int operator==( const RationalDom& BC) const { return 1;}
  int operator!=( const RationalDom& BC) const { return 0;}

  // -- Constants
  const Rational one;
  const Rational zero;

  // -- assignement
  Rep& init( Rep& a ) const{ return a; }
  Rep& init  ( Rep& a, const Rep& b) const { return a = b ; }
  Rep& assign( Rep& a, const Rep& b) const { return a = b ; }

  // -- arithmetic operators
  Rep& mul( Rep& r, const Rep& a, const Rep& b ) const { return r = a * b; };
  Rep& div( Rep& r, const Rep& a, const Rep& b ) const { return r = a / b; };
  Rep& add( Rep& r, const Rep& a, const Rep& b ) const { return r = a + b; };
  Rep& sub( Rep& r, const Rep& a, const Rep& b ) const { return r = a - b; };

  Rep& mulin( Rep& r, const Rep& a) const { return r *= a; };
  Rep& divin( Rep& r, const Rep& a) const { return r /= a; };
  Rep& addin( Rep& r, const Rep& a) const { return r += a; };
  Rep& subin( Rep& r, const Rep& a) const { return r -= a; };

  Rep& axpy( Rep& r, const Rep& a, const Rep& b, const Rep& c ) const
  { return r = a * b + c; };
  Rep& axpyin( Rep& r, const Rep& a, const Rep& b ) const
  { return r += a * b; };
  Rep& axmy( Rep& r, const Rep& a, const Rep& b, const Rep& c ) const
  { return r = a * b - c; };
  Rep& axmyin( Rep& r, const Rep& a, const Rep& b ) const
  { return r -= a * b; };

  // -- unary methods
  Rep& neg( Rep& r, const Rep& a ) const { return r = -a; };
  Rep& inv( Rep& r, const Rep& a ) const { return r = Rational::one/a; };

  // - return n^l 
  Rep& pow(Rep& r, const Rep& n, const long l) const { return r = ::pow(n, l); }
  Rep& pow(Rep& r, const Rep& n, const int l) const { return r = ::pow(n, l); }


        // - Rational number reconstruction
    Rep& ratrecon(Rep& r, const Integer& f, const Integer& m, const Integer& k, bool recurs = false) {
        return r = Rational(f,m,k,recurs);
    }       
    Rep& ratrecon(Rep& r, const Integer& f, const Integer& m, bool recurs=true) {
        return r = Rational(f,m,::sqrt(m),recurs);
    }       


  // - Misc
  size_t length (const Rep& a) const { return ::length(a); }
  int sign    (const Rep& a) const { return ::sign(a); }
  int isone   (const Rep& a) const { return compare(a, one) ==0; }
  int iszero  (const Rep& a) const { return compare(a, zero) ==0; }
  int areEqual (const Rep& a, const Rep& b) const { return compare(a, b) ==0; }
  int areNEqual(const Rep& a, const Rep& b) const { return compare(a, b) !=0; }

  // -- IO
  // -- IO
  std::istream& read ( std::istream& i ) 
  { char ch;
    i >> std::ws >> ch; 
    if (ch != 'R') 
      GivError::throw_error(GivBadFormat("RationalDom::read: bad signature domain"));
    return i;
  }
  std::ostream& write( std::ostream& o ) const { return o << 'R'; }
  std::istream& read ( std::istream& i, Rep& n) const { return i >> n; }
  std::ostream& write( std::ostream& o, const Rep& n) const { return n.print(o); }
};


#endif
