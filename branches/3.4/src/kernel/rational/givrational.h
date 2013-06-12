// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/rational/givrational.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama, T. Gautier
// $Id: givrational.h,v 1.13 2011-02-02 16:23:56 bboyer Exp $
// ==========================================================================
#ifndef __GIVARO_rational_H
#define __GIVARO_rational_H
// #define __GIVARO_GMPplusplus_rational_H

#include "givaro/givinteger.h"
#include "givaro/givmodule.h"

namespace Givaro {

// ----------------------------------- Functions Rational

class Rational ;
int compare(const Rational& a, const Rational& b) ;
int absCompare(const Rational& a, const Rational& b) ;
const Rational pow(const Rational &r, const long l);
const Integer floor (const Rational &r) ;
const Integer ceil  (const Rational &r) ;
const Integer round (const Rational &r) ;
const Integer trunc (const Rational &r) ;
const Rational abs  (const Rational &r) ;
const Rational pow (const Rational& n, unsigned int l);
const Rational pow (const Rational& n, unsigned long l);
unsigned long length (const Rational& r) ;
int sign   (const Rational& r) ;
int isZero (const Rational& r) ;
int isOne  (const Rational& r) ;
int isInteger(const Rational& r);

class RationalDom;

// ----------------------------------- Class Rational

class Rational {

public :
        // Cstor et dstor
    Rational(Neutral n = Neutral::zero) ;
    Rational(int n) ;
    Rational(long n) ;
    Rational(unsigned long n) ;
    Rational(long n, long d ) ;
    Rational(unsigned long n, unsigned long d ) ;
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

    friend const Rational pow (const Rational& n, unsigned int l) {
        Rational r;
        r.num = ::Givaro::pow(n.num, l); r.den =  ::Givaro::pow(n.den, l);
        return r;
    }

    friend const Rational pow (const Rational& n, unsigned long l) {
        Rational r;
        r.num =  ::Givaro::pow(n.num, l); r.den =  ::Givaro::pow(n.den, l);
        return r;
    }

    const Integer nume() const ;
    const Integer deno() const ;
    inline friend unsigned long length (const Rational& r) ;
    inline friend int sign   (const Rational& r) ;
    inline friend int isZero (const Rational& r) ;
    inline friend int isOne  (const Rational& r) ;
    inline friend int isInteger(const Rational& r);

    std::ostream& print ( std::ostream& o ) const ;

    inline Rational reduce(const Rational& R) const ;

    static void SetReduce() ;
    static void SetNoReduce() ;

        // -- Cast operators
    operator short() const { return (short)(int) *this; }
    operator unsigned short() const { return (unsigned short) (unsigned int) *this; }
    operator unsigned char() const { return (unsigned char)(unsigned int) *this; }
    operator unsigned int() const { return (unsigned int) (this->num/this->den); }
    operator int() const  { return (int) (this->num/this->den); }
    operator signed char() const { return (signed char) (int) *this; }
    operator unsigned long() const { return (unsigned long) (this->num/this->den); }
    operator long() const { return (long) (this->num/this->den); }
#ifndef __GIVARO__DONOTUSE_longlong__
    operator unsigned long long() const { return (unsigned long long) (this->num/this->den); }
    operator long long() const { return (long long) (this->num/this->den); }
#endif
    operator std::string() const { return std::string(this->num)+'/'+std::string(this->den); }
    operator float() const { return ((float)this->num)/((float)this->den); }
    operator double() const { return ((double)this->num)/((double)this->den); }

protected: // Internal Representation : num/den
    Integer num, den;

public:
    enum ReduceFlag { Reduce = 0x1, NoReduce = 0x0 } ;
protected:
    static ReduceFlag flags ;    // flags that indicates is reduction is done or not
        // by default = Reduce
    Rational& reduce() ;
        // -- module initialization
    static void Init(int* argc, char***argv);
    static void End();
    friend class GivModule;
    friend class RationalDom;
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
class RationalDom  {
public:
    typedef Rational Element;
    typedef Rational Rep;

        // -- Cstor
    RationalDom() : one(1), zero(0) {}
    template<class X> RationalDom(const X& x) : one(1), zero(0) {}

    int operator==( const RationalDom& ) const { return 1;}
    int operator!=( const RationalDom& ) const { return 0;}

        // -- Constants
    const Rational one;
    const Rational zero;

    unsigned long characteristic() const { return 0UL; }
    Integer& characteristic(Integer& p) const { return p=characteristic();}

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
    Rep& maxpy( Rep& r, const Rep& a, const Rep& b, const Rep& c ) const
        { return r = c - a * b; };
    Rep& axmy( Rep& r, const Rep& a, const Rep& b, const Rep& c ) const
        { return r = a * b - c; };
    Rep& axmyin( Rep& r, const Rep& a, const Rep& b ) const
        { return r = a * b - r ; };
    Rep& maxpyin( Rep& r, const Rep& a, const Rep& b ) const
        { return r -= a * b; };

        // -- unary methods
    Rep& neg( Rep& r, const Rep& a ) const { return r = -a; };
    Rep& inv( Rep& r, const Rep& a ) const { r.num=a.den; r.den=a.num; return r; }
    Rep& negin( Rep& r ) const { r.num=-r.num; return r; }
    Rep& invin( Rep& r ) const { std::swap(r.num,r.den); return r; }

        // - return n^l
    Rep& pow(Rep& r, const Rep& n, const unsigned long l) const { return r =  ::Givaro::pow(n, l); }
    Rep& pow(Rep& r, const Rep& n, const unsigned int l) const { return r =  ::Givaro::pow(n, l); }


        // - Rational number reconstruction
    Rep& ratrecon(Rep& r, const Integer& f, const Integer& m, const Integer& k, bool recurs = false) const {
        return r = Rational(f,m,k,recurs);
    }
    Rep& ratrecon(Rep& r, const Integer& f, const Integer& m, bool recurs=true) const {
        return r = Rational(f,m, ::Givaro::sqrt(m),recurs);
    }


        // - Misc
    size_t length (const Rep& a) const { return  ::Givaro::length(a); }
    int sign    (const Rep& a) const { return  ::Givaro::sign(a); }
    int isOne   (const Rep& a) const { return compare(a, one) ==0; }
    int isZero  (const Rep& a) const { return compare(a, zero) ==0; }
    int areEqual (const Rep& a, const Rep& b) const { return compare(a, b) ==0; }
    int areNEqual(const Rep& a, const Rep& b) const { return compare(a, b) !=0; }
    template< class RandIter > Rep& random(RandIter& g, Rep& r, long s = 1) const { return r=Rational(Integer::random(s), Integer::nonzerorandom(s)); }
    template< class RandIter > Rep& random(RandIter& g, Rep& r, const Rep& b) const { Integer rnum,rden; Integer::random(rnum,b.nume()); Integer::nonzerorandom(rden,b.deno()); return r=Rational(rnum,rden); }
    template< class RandIter > Rep& nonzerorandom(RandIter& g, Rep& r, long s = 1) const { return r=Rational(Integer::nonzerorandom(s), Integer::nonzerorandom(s)); }
    template< class RandIter > Rep& nonzerorandom (RandIter& g,Rep& r, const Rep& b) const { Integer rnum,rden; Integer::nonzerorandom(rnum,b.nume()); Integer::nonzerorandom(rden,b.deno()); return r=Rational(rnum,rden); }


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

} //namespace Givaro

#endif // __GIVARO_rational_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
