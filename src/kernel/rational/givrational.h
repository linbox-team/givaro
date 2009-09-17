// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/rational/givrational.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.
// Authors: M. Samama, T. Gautier
// $Id: givrational.h,v 1.7 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
#ifndef _RATIONAL_H_
#define _RATIONAL_H_

#include "givaro/givinteger.h"
#include "givaro/givmodule.h"


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

    friend const Rational pow (const Rational& n, unsigned int l) { 
        Rational r;
        r.num = ::pow(n.num, l); r.den = ::pow(n.den, l); 
        return r;
    }
    
    friend const Rational pow (const Rational& n, unsigned long l) { 
        Rational r;
        r.num = ::pow(n.num, l); r.den = ::pow(n.den, l); 
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
    operator short() const { return (int) *this; }
    operator unsigned short() const { return (unsigned int) *this; }
    operator unsigned char() const { return (unsigned int) *this; }
    operator unsigned int() const { return (unsigned int) (this->num/this->den); }
    operator int() const  { return (int) (this->num/this->den); }
    operator signed char() const { return (int) *this; }
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
    typedef Rational Element;
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
    Rep& amxy( Rep& r, const Rep& a, const Rep& b, const Rep& c ) const
        { return r = a - b * c; };
    Rep& axmy( Rep& r, const Rep& a, const Rep& b, const Rep& c ) const
        { return r = a * b - c; };
    Rep& axmyin( Rep& r, const Rep& a, const Rep& b ) const
        { return r -= a * b; };

        // -- unary methods
    Rep& neg( Rep& r, const Rep& a ) const { return r = -a; };
    Rep& inv( Rep& r, const Rep& a ) const { return r = Rational::one/a; };

        // - return n^l 
    Rep& pow(Rep& r, const Rep& n, const unsigned long l) const { return r = ::pow(n, l); }
    Rep& pow(Rep& r, const Rep& n, const unsigned int l) const { return r = ::pow(n, l); }


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
    int isOne   (const Rep& a) const { return compare(a, one) ==0; }
    int isZero  (const Rep& a) const { return compare(a, zero) ==0; }
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
