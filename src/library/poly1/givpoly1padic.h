// ================================================================= //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <28 Oct 15 04:14:34 Jean-Guillaume.Dumas@imag.fr>
// ================================================================= //

/** @file givpoly1padic.h
 * @ingroup poly1
 * @brief NO DOC.
 */

#ifndef __GIVARO_poly1_p_adic_H
#define __GIVARO_poly1_p_adic_H
#include <givaro/givinteger.h>
#include <givaro/givpoly1.h>
#include <cmath>
#include <iostream>

namespace Givaro {

    template<class Domain, class Tag=Dense> class Poly1PadicDom;

    //! Poly1 p-adic.
    template<class Domain>
    class Poly1PadicDom<Domain,Dense> : public Poly1Dom<Domain,Dense>, public IntegerDom {
        using Poly1Dom<Domain,Dense>::_domain;
        typedef Poly1Dom<Domain,Dense> Poly_t;

        typedef typename Poly_t::Rep      Rep;
    public:
        typedef typename Poly_t::Element  Element;
        typedef typename Poly_t::Element  pol_Element;
        typedef typename IntegerDom::Element              int_Element;

        Poly1PadicDom (Domain& d, const Indeter& X) : Poly_t (d,X), IntegerDom() {}
        Poly1PadicDom (const Poly_t& P) : Poly_t (P), IntegerDom() {}
        Poly1PadicDom (const Poly_t& P, const IntegerDom& D) : Poly_t (P), IntegerDom(D) {}

        using Poly_t::write;

        // Horner like evaluation of the polynomial for p = _domain.size()
        template<class vect>
        IntegerDom::Element& eval( IntegerDom::Element& E, const vect& P) {
            typename vect::const_reverse_iterator pi = P.rbegin();
            _domain.convert(E, *pi);
            for (++pi;pi != P.rend();++pi) {
                IntegerDom::mulin(E, _domain.size() );
                IntegerDom::Element epi;
                IntegerDom::addin(E, _domain.convert(epi, *pi) );
            }
            return E;
        }

        // Gain is 40% to 75% compared to Integers !!!
        template<class vect>
        uint64_t& eval( uint64_t& E, const vect& P) {
            typename vect::const_reverse_iterator pi = P.rbegin();
            _domain.convert(E,*pi);
            for (++pi;pi != P.rend();++pi) {
                E *= _domain.size();
                uint64_t epi;
                E += _domain.convert(epi, *pi);
            }
            return E;
        }

        template<class unsignedinttype, class vect>
        unsignedinttype& eval( unsignedinttype& E, const vect& P) {
            typename vect::const_reverse_iterator pi = P.rbegin();
            _domain.convert(E,*pi);
            for (++pi;pi != P.rend();++pi) {
                E *= _domain.size();
                unsignedinttype epi;
                E += _domain.convert(epi,*pi);
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
        // See e.g. [von zur Gathen, Gerhard 1999], Modern Computer Algebra
        // Algorithm 9.14
        template<class vect>
        vect& radix(vect& P, const IntegerDom::Element& E, int64_t n = 0) {
            if (n < 1) n = logp(E,_domain.size()) + 1;
            if (n == 1) {
                // Could also be
                // typename Domain::Element e;
                // Poly_t::assign(P, Degree(0), _domain.init(e, E) );
                // But Poly_t::init uses less temporaries
                return Poly_t::init(P, Degree(0), E );
            }
            IntegerDom::Element iq, ir;
            vect Q;
            int64_t t = (n+1)/2;
            IntegerDom::Element q;
            IntegerDom::pow(q, _domain.size(), t);
            IntegerDom::divmod(iq, ir, E, q);
            radix(Q, iq, n-t);
            radix(P, ir, t);
            Degree dp; this->degree(dp,P); ++dp;
            for(int64_t i=t; dp<i; --i)
                P.push_back(_domain.zero);
            P.insert(P.end(),Q.begin(),Q.end());
            return this->setdegree(P);
        }


        // vect is supposed to be a vector of doubles
        // Therefore there is no automatic conversion
        template<class vect>
        vect& fastradixdirect(vect& P, const double& E, uint64_t n) {
            if (n <= 1) {
                P.resize(0);
                typedef typename vect::value_type elem;
                P.push_back( static_cast<elem>(E) );
                return P;
            }
            double iq, ir;
            vect Q;
            int64_t t = (int64_t)(n+1)/2;
            double q = std::pow(double(_domain.size()), double(t));
            iq = floor( E / q );
            ir = E - iq*q;
            radixdirect(Q, iq, n-t);
            radixdirect(P, ir, t);
            Degree dp; degree(dp,P); ++dp;
            for(int64_t i=t; dp<i; --i)
                P.push_back(_domain.zero);
            P.insert(P.end(),Q.begin(),Q.end());
            return setdegree(P);
        }


        template<class vect, class TT>
        vect& radixdirect(vect& P, const TT& E, uint64_t n) {
            P.resize(n);
            TT r = E, s=r;
            for(typename vect::iterator pit=P.begin(); pit != P.end(); ++pit) {
                s /= _domain.size();
                *pit = typename vect::value_type(r-s*_domain.size());
                r = s;
            }
            return P;
        }

        template<class vect>
        vect& radixdirect(vect& P, const double& E, uint64_t n) {
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

} // Givaro
#endif // __GIVARO_poly1_p_adic_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
