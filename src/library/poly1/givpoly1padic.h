// ================================================================= //
// (C) The Linbox Group 1999
// Time-stamp: <25 Sep 01 19:04:26 Jean-Guillaume.Dumas@imag.fr> 
// ================================================================= //

#ifndef _GIV_POLY1_P_ADIC_H_
#define _GIV_POLY1_P_ADIC_H_
#include <givaro/givinteger.h>
#include <givaro/givpoly1.h>
#include <math.h>
#include <iostream>

template<class Domain, class Tag> class Poly1PadicDom;



template<class Domain>
class Poly1PadicDom<Domain,Dense> : public Poly1Dom<Domain,Dense>, public IntegerDom {
    typedef Poly1Dom<Domain,Dense> Poly_t;
    
    typedef typename Poly_t::Rep      Rep;
public:
    typedef typename Poly_t::element  element;
    typedef typename Poly_t::element  pol_element;
    typedef typename IntegerDom::element              int_element;

    Poly1PadicDom (Domain& d, const Indeter& X) : Poly_t (d,X), IntegerDom() {}
    Poly1PadicDom (const Poly_t& P) : Poly_t (P), IntegerDom() {}
    Poly1PadicDom (const Poly_t& P, const IntegerDom& D) : Poly_t (P), IntegerDom(D) {}



    std::ostream& write( std::ostream& o, const pol_element& p) { 
        return Poly_t::write(o, p); 
    }
    
            


        // Horner like evaluation of the polynomial for p = _domain.size()
    template<class vect>
    IntegerDom::element& eval( IntegerDom::element& E, const vect& P) {
        typename vect::const_reverse_iterator pi = P.rbegin();
        IntegerDom::init(E, _domain.convert(*pi) );
        for (++pi;pi != P.rend();++pi) {
            IntegerDom::mulin(E, _domain.size() );
            IntegerDom::addin(E, _domain.convert(*pi) );
        }
        return E;
    }

        // Gain is 40% to 75% compared to Integers !!!
    template<class vect>
    unsigned long& eval( unsigned long& E, const vect& P) {
        typename vect::const_reverse_iterator pi = P.rbegin();
        _domain.convert(E,*pi);
        for (++pi;pi != P.rend();++pi) {
            E *= _domain.size();
            E += _domain.convert(*pi);
        }
        return E;
    }

    template<class unsignedinttype, class vect>
    unsignedinttype& eval( unsignedinttype& E, const vect& P) {
        typename vect::const_reverse_iterator pi = P.rbegin();
        _domain.convert(E,*pi);
        for (++pi;pi != P.rend();++pi) {
            E *= _domain.size();
            E += _domain.convert(*pi);
        }
        return E;
    }

    template<class elem, class vect>
    elem& evaldirect( elem& E, const vect& P) {
        typename vect::const_reverse_iterator pi = P.rbegin();
        E = elem(*pi);
        for (++pi;pi != P.rend();++pi) {
            E *= _domain.size();
            E += elem(*pi);
        }
        return E;
    }


        // Reconstruction of the radix decomposition mod _domain.size()
        // E into a polynomial P of degree n-1 or of lowest degree if n==0.
    template<class vect>
    vect& radix(vect& P, const IntegerDom::element& E, long n = 0) {
        if (n < 1) n = logp(E,_domain.size()) + 1;
        if (n == 1) {
            typename Domain::element e;
            return Poly_t::init(P, Degree(0), _domain.init(e, E) );
        }
        IntegerDom::element iq, ir;
        vect Q; 
            long t = (n+1)/2;
            IntegerDom::element q;
            IntegerDom::pow(q, _domain.size(), t);
            IntegerDom::divmod(iq, ir, E, q);
            radix(Q, iq, n-t);
            radix(P, ir, t);
            Degree dp; degree(dp,P); ++dp;
            for(long i=t; dp<i; --i)
                P.push_back(_domain.zero);
        P.insert(P.end(),Q.begin(),Q.end());
        return setdegree(P);
    } 


    template<class vect>
    vect& fastradixdirect(vect& P, const double& E, unsigned long n) {
        if (n <= 1) {
            P.resize(0);
            typedef typename vect::value_type elem;
            P.push_back( elem( E ) );
            return P;
        }
        double iq, ir;
        vect Q; 
            long t = (n+1)/2;
            double q = std::pow(double(_domain.size()), double(t));
            iq = floor( E / q );
            ir = E - iq*q;
            radixdirect(Q, iq, n-t);
            radixdirect(P, ir, t);
            Degree dp; degree(dp,P); ++dp;
            for(long i=t; dp<i; --i)
                P.push_back(_domain.zero);
        P.insert(P.end(),Q.begin(),Q.end());
        return setdegree(P);
    } 


    template<class vect>
    vect& radixdirect(vect& P, const double& E, unsigned long n) {
        P.resize(n);
        double r = E, s;
        for(typename vect::iterator pit=P.begin(); pit != P.end(); ++pit) {
            s = floor( r / _domain.size() );
            *pit = typename vect::value_type(r-s*_domain.size());
            r = s;
        }
        return P;
    } 
};

#endif 
