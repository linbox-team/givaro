#ifndef _GIVARO_MONTG32_H_
#define _GIVARO_MONTG32_H_
// ==========================================================================
// author: JG Dumas (from P. Zimmermann's Montgomery implementation)
// $Id: givmontg32.h,v 1.5 2005-02-02 19:08:29 pernet Exp $
// ==========================================================================
//
#include "givbasictype.h"
#include "giverror.h"
#include "giv_randiter.h"
#include <math.h>


// ==========================================================================
// -- This class implement the standard arithmetic with Modulo elements:
//    Reduction is made through Montgomery's reduction
//    Representation of a is by storing (aB).
//
//    We must have (p-1)^2 + p*(B-1) < B^2,
//    i.e. p<=40504 for B=2^16
// ==========================================================================

#define B32 65536UL
#define MASK32 65535UL
#define HALF_BITS32 16

template<class TYPE> class Montgomery;


template<>
class Montgomery<Std32> {
public:
        // ----- Exported Types and constantes
    typedef uint32 Residu_t;                    // - type to store residue
    enum { size_rep = sizeof(Residu_t) };      // - size of the storage type
        // ----- Representation of element of the domain Montgomery
    typedef uint32 Rep;
    typedef uint32 element;
    typedef element Element;

        // ----- Constructor 
    Montgomery() : _p(0UL), _dp(0.0), zero(0UL), one(1UL) {}

    Montgomery( Residu_t p, int expo = 1)
        : _p(p), _Bp(B32%p), _B2p( (_Bp<<HALF_BITS32) % p), _B3p( (_B2p<<HALF_BITS32) % p), _nim( -Montgomery<Std32>::invext(_p,B32) ), _dp((double)p), _invdp(1.0/(double)p), zero(0UL), one( redcsal(_B2p) ) {}

    Montgomery( const Montgomery<Std32>& F)
        : _p(F._p), _Bp(F._Bp), _B2p( F._B2p), _B3p( F._B3p), _nim(F._nim),_dp(F._dp), _invdp(F._invdp), zero(0UL), one(F.one) { }


    int operator==( const Montgomery<Std32>& BC) const { return _p == BC._p;}
    int operator!=( const Montgomery<Std32>& BC) const { return _p != BC._p;}
 
    Montgomery<Std32>& operator=( const Montgomery<Std32>& F) { 
        this->_p = F._p; 
        this->_Bp = F._Bp;
        this->_B2p = F._B2p;
        this->_B3p = F._B3p;
        this->_nim = F._nim;
        this->_dp = F._dp;
        this->_invdp = F._invdp;
        return *this;
    }

        // ----- Access to the modulus 
    Residu_t residu() const;
    Residu_t size() const {return _p;}
    Rep access( const Rep a ) const { return a; }
    Residu_t characteristic() const { return _p; }
    Residu_t characteristic(Residu_t p) const { return p=_p; }
    Residu_t cardinality() const { return _p; }


        // ----- Access to the modulus 
    Rep& init( Rep& a ) const;
    Rep& init( Rep& r, const long a) const ;
    Rep& init( Rep& r, const unsigned long a) const ;
    Rep& init( Rep& a, const int i) const ;
    Rep& init( Rep& a, const unsigned int i) const ;
    Rep& init ( Rep& r, const Integer& residu ) const ;
    
        // Initialisation from double ( added for FFLAS usage) (C Pernet)
    Rep& init( Rep& a, const double i) const;
    Rep& init( Rep& a, const float i) const;

    unsigned long int& convert(unsigned long int& r, const Rep a) const { 
	    uint32 ur;
	    return r = (unsigned long)redc(ur,a);}

    uint32& convert(uint32& r, const Rep a) const { 
	    unsigned long ur;
	    return r = (uint32)convert(ur, a);
    }

     int32& convert(int32& r, const Rep a) const { 
	    unsigned long ur;
	    return r = (int32)convert(ur, a);
    }

   long int& convert(long int& r, const Rep a) const { 
	    unsigned long ur;
	    return r = (long int)convert(ur, a);
    }

    Integer& convert(Integer& i, const Rep a) const {
        unsigned long ur;
        return i = (Integer)convert(ur, a);
    }        

        // Conversion to double ( added for FFLAS usage) (C Pernet)
    float& convert(float& r, const Rep a ) const { 
	    unsigned long ur;
	    return r = (float)convert(ur, a); }
    double& convert(double& r, const Rep a ) const {
	    unsigned long ur;
	    return r = (double)convert(ur, a); }
    
        // ----- Misc methods
    int iszero( const Rep a ) const;
    int isone ( const Rep a ) const;
    int isZero( const Rep a ) const;
    int isOne ( const Rep a ) const;
    size_t length ( const Rep a ) const;

        // ----- Equality between two elements
    int areEqual(const  Rep& a, const Rep& b) const { 
	    return a==b;
    }

        // ----- Operations with reduction: r <- a op b mod p, r <- op a mod p
    Rep& mul (Rep& r, const Rep a, const Rep b) const;
    Rep& div (Rep& r, const Rep a, const Rep b) const;
    Rep& add (Rep& r, const Rep a, const Rep b) const;
    Rep& sub (Rep& r, const Rep a, const Rep b) const;
    Rep& neg (Rep& r, const Rep a) const;
    Rep& inv (Rep& r, const Rep a) const;

    Rep& mulin (Rep& r, const Rep a) const;
    Rep& divin (Rep& r, const Rep a) const;
    Rep& addin (Rep& r, const Rep a) const;
    Rep& subin (Rep& r, const Rep a) const;
    Rep& negin (Rep& r) const;
    Rep& invin (Rep& r) const;

        // -- axpy: r <- a * x + y mod p
    Rep& axpy  (Rep& r, const Rep a, const Rep b, const Rep c) const;
        // -- axpyin: r <- r + a * x mod p
    Rep& axpyin(Rep& r, const Rep a, const Rep b) const;
        // -- axmy: r <- a * x - y mod p
    Rep& axmy  (Rep& r, const Rep a, const Rep b, const Rep c) const;
        // -- axmyin: r <- r - a * x mod p
    Rep& axmyin(Rep& r, const Rep a, const Rep b) const;
        // -- Misc: r <- a mod p
    Rep& assign ( Rep& r, const Rep a) const;

        // ----- random generators
    template< class RandIter > Rep& random(RandIter&, Rep& r) const ;
    template< class RandIter > Rep& random(RandIter&, Rep& r, long s) const ;
    template< class RandIter > Rep& random(RandIter&, Rep& r, const Rep& b) const ;
    template< class RandIter > Rep& nonzerorandom(RandIter&, Rep& r) const ;
    template< class RandIter > Rep& nonzerorandom(RandIter&, Rep& r, long s) const ;
    template< class RandIter > Rep& nonzerorandom(RandIter&, Rep& r, const Rep& b) const ;
    
    typedef GIV_randIter< Montgomery<Std32> , Rep > randIter;

        // --- IO methods
    std::istream& read ( std::istream& s );
    std::ostream& write( std::ostream& s ) const;
    std::istream& read ( std::istream& s, Rep& a ) const;
    std::ostream& write( std::ostream& s, const Rep a ) const;

protected:
        // -- based on modular inverse, d = a*u + b*v
//   static const int32 gcdext ( int32& u, int32& v, const int32 a, const int32 b );
    int32& gcdext (int32& d, int32& u, int32& v, const int32 a, const int32 b ) const;
    int32& invext (int32& u, const int32 a, const int32 b ) const;
    int32 invext(const int32 a, const int32 b ) const;



    Element& redc(Element&, const Element) const ;
    Element redcal(const Element) const;
    Element redcsal(const Element) const;
    Element& redcin(Element&) const;
    Element& redcs(Element&, const Element) const;
    Element& redcsin(Element&) const;
    

protected:
        // -- data representation of the domain:
    Residu_t _p;
    Residu_t _Bp;
    Residu_t _B2p;
    Residu_t _B3p;
    Residu_t _nim;
    double _dp;
    double _invdp;
    

    static void Init();
    static void End();

public:
        // ----- Constantes 
    const Rep zero;
    const Rep one;
};


#include "givmontg32.inl"

#endif
