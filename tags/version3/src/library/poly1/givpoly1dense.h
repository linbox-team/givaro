#ifndef _GIV_POLY1_DENSE_H_
#define _GIV_POLY1_DENSE_H_
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givpoly1dense.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givpoly1dense.h,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description: univariate polynom over T
// - we assume that T is a ring (0,1,+,*) with:

// #include "givarray0.h"
#include <iostream>
#include <vector>
#include "givaro/givdegree.h"
#include "givaro/givindeter.h"
#include "givaro/givinteger.h"



template < typename T > class givvector : public std::vector<T> {
    typedef givvector<T>     Self_t;
public:
    givvector() : std::vector<T>() {}
    givvector(size_t s) : std::vector<T>(s) { }
    givvector(const Self_t& p, givNoCopy xxx) : std::vector<T>(p) {}
    givvector(const Self_t& p, givWithCopy xxx) : std::vector<T>(p) {}
    Self_t& reallocate (size_t s) { resize(s); return *this; }
    Self_t& logcopy(const Self_t& src) { return *this = src; }
    Self_t& copy(const Self_t& src) { return *this = src; }
    int areEqual(const Self_t& p) const { 
        return *this == p; 
    }
    int areNEqual(const Self_t& p) const { 
        return *this != p;
    }
};

// template < typename T >
// std::ostream& operator<< (std::ostream& o, const givvector<T>& v) {
//     o << "Poly (s [v_1, ..,  v_s]) : ";
//     o << v.size();
//     for(size_t i=0; i<v.size(); ++i) {
//         o << " " << v[i];
//     }
// }    
    


//  -------------------------------------------- Class Poly1Dom<Domain>
template <class Domain>
class Poly1Dom<Domain,Dense> {
protected:  //  -- Representation 
    Domain& 		_domain;  // -- subdomain
    Indeter		_x;	  // -- for I/O, if any
public :

        // -- Exported types
    typedef          Domain		       Domain_t;
    typedef typename Domain::element               Type_t;

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
    typedef          Storage_t                 element;

    Poly1Dom (Domain& d, const Indeter& X = Indeter() );
    Poly1Dom (const Self_t&);

    int operator==( const Poly1Dom<Domain,Dense>& BC) const 
        { return _domain == BC._domain;}
    int operator!=( const Poly1Dom<Domain,Dense>& BC) const 
        { return _domain != BC._domain;}

        // -- Return the domain of the entries
    Domain& subdomain() const { return _domain; }
        // -- Return the domain of the entries
    Domain& getdomain() const { return _domain; }

        // -- Constantes
    const Rep zero;
    const Rep one;
	
        // -- Init polynomial 
    Rep& init(Rep& a) const;
    Rep& init(Rep& p, const Integer &cste ) const;
        // -- Init polynomial with field value : F.assign(p[0],cste)
    Rep& assign(Rep& p, const Type_t &cste ) const;

        // -- Allocate a polynomial with deg+1 coefficients, each of them are 
        // set to zero, except the leading coef which is set to one.
    Rep& init (Rep& r, const Degree deg) const;

        // -- For polynomial = lcoeff X^deg 
    Rep& init (Rep& p, const Degree deg , const Integer& lcoeff) const;
        //    F.assign(P[deg], lcoeff); 
    Rep& assign (Rep& p, const Degree deg , const Type_t& lcoeff) const;

        // -- Assignment p = q
    Rep& assign( Rep& p, const Rep& q) const;

        // -- Dstor
    ~Poly1Dom ();

        // -- Comparaison operator
    int iszero  ( const Rep& P ) const;
    int isone   ( const Rep& P ) const;
    int areEqual ( const Rep& P, const Rep& Q ) const;
    int areNEqual( const Rep& P, const Rep& Q ) const;

        // -- Returns the leading coefficients
    Type_t& leadcoef(Type_t& c, const Rep& P) const;

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

        // -- misc
        // -- W <-- P^n
    Rep& pow( Rep& W, const Rep& P, long n) const;
        // -- W <-- P^n [ U ]
    Rep& powmod( Rep& W, const Rep& P, IntegerDom::element pwr, const Rep& U) const;
    template < class MyInt >
    Rep& powmod( Rep& W, const Rep& P, MyInt pwr, const Rep& U) const {
        return powmod(W, P, (IntegerDom::element)pwr, U);
    }

        // -- W <-- P(X^b)
    Rep&  power_compose( Rep& W, const Rep& P, long b) const;

        // -- n th cyclotomic polynomial
    Rep& cyclotomic( Rep& P, long n) const;


        // -- Random generators
        // -- Random dense polynomial of degree 0
    template< class RandIter > Rep& random(RandIter& g, Rep& r) const;
        // -- Random dense polynomial of size s
    template< class RandIter > Rep& random(RandIter& g, Rep& r, long s) const ;
        // -- Random dense polynomial of degree d
    template< class RandIter > Rep& random(RandIter& g, Rep& r, Degree s) const ;
        // -- Random dense polynomial with same size as b. 
    template< class RandIter > Rep& random(RandIter& g, Rep& r, const Rep& b) const;

    template< class RandIter > Rep& nonzerorandom(RandIter& g, Rep& r) const;
    template< class RandIter > Rep& nonzerorandom(RandIter& g, Rep& r, long s) const;
    template< class RandIter > Rep& nonzerorandom(RandIter& g, Rep& r, Degree s) const ;
    template< class RandIter > Rep& nonzerorandom(RandIter& g, Rep& r, const Rep& b) const;
    
        // -- Square free decomposition
    void sqrfree(size_t& Nfact, Rep* Fact, const Rep& P) const;


}; //  ------------------------------- End Of The Class Poly1Dom<Type_t>

#endif