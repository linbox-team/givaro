// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givinterp.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: JG Dumas
// $Id: givinterp.h,v 1.3 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================
/** @file givinterpgeom-multip.h
 * @ingroup poly1
 * @brief Interpolation at geometric points
 * @bib
 * - A Bostan and E Schost, <i>Polynomial evaluation and interpolation on special sets of points</i>,
 * Journal of Complexity 21(4): 420-446, 2005.
 */

#ifndef __GIVARO_multiple_interpolation_at_geometric_points_H
#define __GIVARO_multiple_interpolation_at_geometric_points_H

#include "givaro/givconfig.h"
#include "givaro/giverror.h"
#include "givaro/givpoly1.h"
#include <givaro/givtruncdomain.h>
#include <vector>

namespace Givaro {

    //! Newton (multip)
    template<class Domain, bool REDUCE = true>
    struct NewtonInterpGeomMultip : TruncDom<Domain>  {
        typedef std::vector< typename Domain::Element > Vect_t;
        typedef typename TruncDom<Domain>::Polynomial_t Polynomial;
        typedef Polynomial Element;
        typedef typename TruncDom<Domain>::Element Truncated;
        typedef typename TruncDom<Domain>::Type_t Type_t;
        typedef typename std::vector< Polynomial > VectPoly_t;

    private:
        Type_t _q;
        Type_t powerq;
        Type_t _ui;
        Polynomial _qi;
        Polynomial _mui;
        VectPoly_t _wi;
        Polynomial _g2;
        bool flip;
        Degree _deg;

    public:
        // Usage :
        // Cstor + initialize with first evaluation (at 1)
        // Then calls to operator(), will evaluate the blackbox at q^i
        // Finally calls to Newton would yield the coefficients in the Newton basis
        // While interpolator calls Newton, and then transforms to monomial basis

        NewtonInterpGeomMultip (const Domain& d, const Indeter& X = Indeter() ) : TruncDom<Domain>(d,X), powerq(d.one), _ui(d.one), flip(false), _deg(0) {
            d.generator(_q); // get a primitive root
        }

        template<typename BlackBox>
        void initialize (const BlackBox& bb) {
            _qi.resize(0);
            _mui.resize(0);
            _wi.resize(0);
            _g2.resize(0);
            this->_domain.assign(powerq, this->_domain.one);
            this->_domain.assign(_ui, this->_domain.one);
            _deg = 0;
            _qi.push_back(this->_domain.one);
            _mui.push_back(this->_domain.one);
            Vect_t v0;
            bb(v0, this->_domain.one);

            _wi.resize(v0.size());

            typename Vect_t::const_iterator iter_v0 = v0.begin();
            for(typename VectPoly_t::iterator iter_wi = _wi.begin();
                iter_wi != _wi.end(); ++iter_wi, ++iter_v0)
                iter_wi->push_back(*iter_v0);

            _g2.push_back(this->_domain.one);
            flip = true;
        }

        template<typename BlackBox>
        void operator() (const BlackBox& bb) {
            Type_t qi;
            this->_domain.mul(qi, _qi.back(), powerq);
            _qi.push_back(qi);

            this->_domain.mulin(powerq, _q);

            Type_t ui, mui;
            this->_domain.sub(ui, powerq, this->_domain.one);

            this->_domain.mul(mui, powerq, _mui.back() );
            this->_domain.divin(mui, ui );
            this->_domain.negin(mui);
            _mui.push_back(mui);

            this->_domain.mulin(_ui, ui );

            Vect_t vi;
            bb(vi, powerq);

            typename Vect_t::iterator iter_vi = vi.begin();
            typename VectPoly_t::iterator iter_wi = _wi.begin();
            for( ; iter_vi != vi.end(); ++iter_vi, ++iter_wi) {
                Type_t wi;
                this->_domain.div(wi, *iter_vi, _ui);
                iter_wi->push_back(wi);
            }

            Type_t gi;
            this->_domain.div(gi, qi, _ui);
            if (flip) this->_domain.negin(gi);
            _g2.push_back(gi);

            flip = !flip;
            ++_deg;
        }



        VectPoly_t& Newton(VectPoly_t& inter) {
            this->getpoldomain().setdegree(_g2);

            Truncated QU;
            this->assign(QU,_g2); // truncated

            inter.resize(_wi.size());
            typename VectPoly_t::iterator iter_inter = inter.begin();
            typename VectPoly_t::iterator iter_wi = _wi.begin();
            for( ; iter_wi != _wi.end(); ++iter_wi, ++iter_inter) {
                Truncated G,W;
                this->getpoldomain().setdegree(*iter_wi);
                this->assign(W, *iter_wi); // truncated

                // truncated multiplication
                this->mul(G, W, QU, 0, _deg);

                this->convert(*iter_inter, G); // trunc to polynomial

                for(size_t i=0; i<iter_inter->size(); ++i)
                    this->_domain.divin((*iter_inter)[i], _qi[i]);
            }

            return inter;
        }


        VectPoly_t& interpolator(VectPoly_t& inter) {
            this->Newton(inter);

            Polynomial _rev_mui;
            this->getpoldomain().setdegree(_mui);
            this->getpoldomain().reverse(_rev_mui, _mui);
            Truncated U;
            this->assign(U, _rev_mui); // truncated

            typename VectPoly_t::iterator iter_inter = inter.begin();
            for( ; iter_inter != inter.end(); ++iter_inter) {

                Type_t mvi;
                Polynomial mwi(_qi.size()), mzi(_qi.size());
                for(size_t i=0; i<iter_inter->size(); ++i) {
                    this->_domain.mul(mvi, (*iter_inter)[i],_qi[i]);
                    if (i & 1) this->_domain.negin(mvi);
                    this->_domain.div(mwi[i], mvi, _mui[i]);
                    this->_domain.div(mzi[i], _mui[i], _qi[i]);
                    if (i & 1) this->_domain.negin(mzi[i]);
                }

                this->getpoldomain().setdegree( mwi);

                Truncated G,W;
                this->assign(W, mwi); // truncated

                // Transposed multiplication (U has been reversed)
                this->mul(G, U, W, _deg, _deg * 2);
                this->divin(G,_deg);

                this->convert(*iter_inter, G); // trunc to polynomial

                for(size_t i=0; i<iter_inter->size(); ++i)
                    this->_domain.mulin( (*iter_inter)[i], mzi[i]);
            }

            return inter;
        }

    };

} // Givaro


#endif // __GIVARO_multiple_interpolation_at_geometric_points_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
