#ifndef _GIV_POLY1_DENSE_H_
#define _GIV_POLY1_DENSE_H_
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givpoly1dense.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givpoly1dense.h,v 1.16 2008-07-18 12:42:55 jgdumas Exp $
// ==========================================================================
// Description: univariate polynom over T
// - we assume that T is a ring (0,1,+,*) with:

#include <iostream>
#include <vector>
#include "givaro/givdegree.h"
#include "givaro/givindeter.h"
#include "givaro/givinteger.h"
#include "givaro/givrandom.h"


template < typename T, typename A=std::allocator<T> > 
class givvector : public std::vector<T,A> {
    typedef givvector<T,A>     Self_t;
public:
    givvector() : std::vector<T,A>() {}
    givvector(size_t s) : std::vector<T,A>(s) { }
    givvector(const Self_t& p, givNoCopy xxx) : std::vector<T,A>(p) {}
    givvector(const Self_t& p, givWithCopy xxx) : std::vector<T,A>(p) {}
    Self_t& reallocate (size_t s) { this->resize(s); return *this; }
    Self_t& logcopy(const Self_t& src) { return *this = src; }
    Self_t& copy(const Self_t& src) { return *this = src; }
    int areEqual(const Self_t& p) const { 
        return *this == p; 
    }
    int areNEqual(const Self_t& p) const { 
        return *this != p;
    }

};

//  -------------------------------------------- Class Poly1Dom<Domain>
template <class Domain>
class Poly1Dom<Domain,Dense> {
protected:  //  -- Representation 
    Domain 		_domain;  // -- subdomain
    Indeter		_x;	  // -- for I/O, if any
public :

        // -- Exported types
    typedef          Domain		       Domain_t;
    typedef typename Domain::Element           Type_t;

        // -- Self_t
    typedef          Poly1Dom<Domain,Dense>    Self_t;

        // -- The representation of a dense polynomial.
        // assuming that we have correct operator, especially size, allocate
        // , reallocate
        // - zero is Rep.size() ==0 or Rep.size() =1 && Rep[0] ==0
        // - Rep.size() is the degree + 1 if !=0
//     typedef          Array0<Type_t>            Storage_t;
    typedef          givvector<Type_t>            Storage_t;
    typedef          Storage_t                 Rep;
    typedef          Storage_t                 Element;

    Poly1Dom (const Domain& d, const Indeter& X = Indeter() );
    Poly1Dom (const Self_t&);

    int operator==( const Poly1Dom<Domain,Dense>& BC) const 
        { return _domain == BC._domain;}
    int operator!=( const Poly1Dom<Domain,Dense>& BC) const 
        { return _domain != BC._domain;}

        // -- Return the domain of the entries
    const Domain& subdomain() const { return _domain; }
        // -- Return the domain of the entries
    const Domain& getdomain() const { return _domain; }
    Domain& setdomain(const Domain& D) { return _domain = D; }
            

        // -- Constantes
     Rep zero;
     Rep one;
	
        // -- Init polynomial 
    Rep& init(Rep& a) const;
        // -- Init polynomial with value : F.init(p[0],cste)
    template<class XXX>
    Rep& init(Rep& p, const XXX &cste ) const;

        // -- Allocate a polynomial with deg+1 coefficients, each of them are 
        // set to zero, except the leading coef which is set to one.
    Rep& init (Rep& r, const Degree deg) const;

        // -- For polynomial = lcoeff X^deg 
    template<class XXX>
    Rep& init (Rep& p, const Degree deg , const XXX& lcoeff) const;

        //    F.assign(P[deg], lcoeff); 
    Rep& assign (Rep& p, const Degree deg , const Type_t& lcoeff) const;
        // -- Assign polynomial with field value : F.assign(p[0],cste)
    Rep& assign(Rep& p, const Type_t &cste ) const {
	    return assign(p, Degree(0), cste);
    }
        // -- Assignment p = q
    Rep& assign( Rep& p, const Rep& q) const;

        // -- Convert polynomials : F.convert(cste, p[0])
    template<class XXX>
    XXX& convert(XXX& p, const Rep &) const;

    template<class UU, template<class XX> class Vect>
    Vect<UU>& convert( Vect<UU>&, const Rep& P ) const ;

        // -- Dstor
    ~Poly1Dom ();

        // -- Comparaison operator
    int isZero  ( const Rep& P ) const;
    //    int isZero  ( const Rep& P ) const{return iszero(P);}
    int isOne   ( const Rep& P ) const;
    int areEqual ( const Rep& P, const Rep& Q ) const;
    int areNEqual( const Rep& P, const Rep& Q ) const;

        // -- Returns the leading coefficients
    Type_t& leadcoef(Type_t& c, const Rep& P) const;

        // -- Returns the i-th coefficients
    Type_t& getEntry(Type_t& c, const Degree& i, const Rep& P) const;

        // -- Returns the degree of polynomial
    Degree& degree(Degree& d, const Rep& P) const;

        // -- Returns the valuation of polynomial
    Degree& val(Degree& d, const Rep& P) const;
  
        // -- Compute the degree of P
    Rep& setdegree( Rep& P ) const;

        // -- Evaluation on one point.
    Type_t& eval(Type_t& pval, const Rep& P, const Type_t& val) const;

        // -- Returns the differentiate polynomial
    Rep& diff( Rep& P, const Rep& Q) const;

        // -- 
    std::istream& read ( std::istream& i );
    std::ostream& write( std::ostream& o ) const;
    std::istream& read ( std::istream& i, Rep& n) const;
    std::ostream& write( std::ostream& o, const Rep& n) const;

        // -- Arithmetics operators
    Rep& addin ( Rep& res, const Rep& u ) const;
    Rep& add ( Rep& res, const Rep& u, const Rep& v ) const;
    Rep& add ( Rep& res, const Rep& u, const Type_t& val ) const;
    Rep& add ( Rep& res, const Type_t& val, const Rep& v ) const;

    Rep& subin ( Rep& res, const Rep& u ) const;
    Rep& sub ( Rep& res, const Rep& u, const Rep& v ) const;
    Rep& sub ( Rep& res, const Rep& u, const Type_t& val ) const;
    Rep& sub ( Rep& res, const Type_t& val, const Rep& v ) const;

    Rep& negin ( Rep& res ) const;
    Rep& neg ( Rep& res, const Rep& u ) const;

    Rep& mulin ( Rep& q, const Rep& a ) const;
    Rep& mulin ( Rep& q, const Type_t& a ) const;
    Rep& mul   ( Rep& q, const Rep& a, const Rep& b ) const;
    Rep& mul   ( Rep& q, const Type_t& a, const Rep& b ) const;
    Rep& mul   ( Rep& q, const Rep& a, const Type_t& b ) const;

    Rep& shiftin ( Rep&, int ) const;
    Rep& shift   ( Rep&, const Rep&, int ) const;

    Rep& divin ( Rep& q, const Rep& a ) const;
    Rep& divin ( Rep& q, const Type_t& a ) const;
    Rep& div   ( Rep& q, const Rep& a, const Rep& b ) const;
    Rep& div   ( Rep& q, const Type_t& a, const Rep& b ) const;
    Rep& div   ( Rep& q, const Rep& a, const Type_t& b ) const;

    Rep& modin ( Rep& q, const Rep& a ) const;
    Rep& modin ( Rep& q, const Type_t& a ) const;
    Rep& mod   ( Rep& q, const Rep& a, const Rep& b ) const;
    Rep& mod   ( Rep& q, const Type_t& a, const Rep& b ) const;
    Rep& mod   ( Rep& q, const Rep& a, const Type_t& b ) const;


    Rep& axpy  (Rep& r, const Rep& a, const Rep& x, const Rep& y) const;
    Rep& axpy  (Rep& r, const Type_t& a, const Rep& x, const Rep& y) const;
    Rep& axpyin(Rep& r, const Rep& a, const Rep& x) const;
    Rep& axpyin(Rep& r, const Type_t& a, const Rep& x) const;
        // -- amxy: r <- c - a * b
    Rep& amxy  (Rep& r, const Rep& a, const Rep& b, const Rep& c) const;
    Rep& amxy  (Rep& r, const Type_t& a, const Rep& b, const Rep& c) const;
        // -- amxyin: r -= a*b
    Rep& amxyin(Rep& r, const Rep& a, const Rep& b) const;
    Rep& amxyin(Rep& r, const Type_t& a, const Rep& b) const;
        // -- axmy: r <- a * x - y
    Rep& axmy  (Rep& r, const Rep& a, const Rep& x, const Rep& y) const;
    Rep& axmy  (Rep& r, const Type_t& a, const Rep& x, const Rep& y) const;

        // A = q*B + r
    Rep& divmod( Rep& q, Rep& r, const Rep& a, const Rep& b ) const;

        // m*A = q*B + r
    Rep& pdivmod( Rep& q, Rep& r, Type_t& m, const Rep& a, const Rep& b ) const;
    Rep& pmod( Rep& r, Type_t& m, const Rep& a, const Rep& b ) const;
    Rep& pmod( Rep& r, const Rep& a, const Rep& b ) const;
    Rep& pdiv( Rep& q, Type_t& m, const Rep& a, const Rep& b ) const;
    Rep& pdiv( Rep& q, const Rep& a, const Rep& b ) const;


        // -- gcd D = gcd(P,Q) = P*U+Q*V;
    Rep& gcd ( Rep& D, const Rep& P, const Rep& Q) const;  
    Rep& gcd ( Rep& D, Rep& U, Rep& V, const Rep& P, const Rep& Q) const;  
    Rep& lcm ( Rep& D, const Rep& P, const Rep& Q) const;  
        // -- modular inverse of P : U P = 1 + V Q
    Rep& invmod ( Rep& U, const Rep& P, const Rep& Q) const;  

        // -- misc
        // -- W <-- P^n
    Rep& pow( Rep& W, const Rep& P, long n) const;
        // -- W <-- P^n [ U ]
    Rep& powmod( Rep& W, const Rep& P, IntegerDom::Element pwr, const Rep& U) const;
    template < class MyInt >
    Rep& powmod( Rep& W, const Rep& P, MyInt pwr, const Rep& U) const {
        return powmod(W, P, (IntegerDom::Element)pwr, U);
    }

        // -- W <-- P(X^b)
    Rep& power_compose( Rep& W, const Rep& P, long b) const;

        // -- n th cyclotomic polynomial
    Rep& cyclotomic( Rep& P, long n) const;


        // -- Random generators
        // -- Random dense polynomial of degree 0
    template< class RandIter > Rep& random(RandIter& g, Rep& r) const;
        // -- Random dense polynomial of size s
    template< class RandIter > Rep& random(RandIter& g, Rep& r, long s) const ;
        // -- Random dense polynomial of degree d
    template< class RandIter > Rep& random(RandIter& g, Rep& r, Degree s) const ;

    Rep& random(GivRandom& g, Rep& r, Degree s) const ;
        // -- Random dense polynomial with same size as b. 
    template< class RandIter > Rep& random(RandIter& g, Rep& r, const Rep& b) const;

    template< class RandIter > Rep& nonzerorandom(RandIter& g, Rep& r) const;
    template< class RandIter > Rep& nonzerorandom(RandIter& g, Rep& r, long s) const;
    template< class RandIter > Rep& nonzerorandom(RandIter& g, Rep& r, Degree s) const ;
    template< class RandIter > Rep& nonzerorandom(RandIter& g, Rep& r, const Rep& b) const;
    
        // -- Square free decomposition
     size_t& sqrfree(size_t& Nfact, Rep* Fact, const Rep& P) const;


}; //  ------------------------------- End Of The Class Poly1Dom<Type_t>

#endif
