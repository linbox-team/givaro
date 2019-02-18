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

/** @file givinterp.h
 * @ingroup poly1
 * @brief NO DOC
 */

#ifndef __GIVARO_interpolation_H
#define __GIVARO_interpolation_H

#include "givaro/givconfig.h"
#include "givaro/giverror.h"
#include "givaro/givpoly1.h"

namespace Givaro {

    //! Interpolation
    template<class Domain, bool REDUCE = true>
    struct Interpolation : Poly1Dom<Domain,Dense>  {
        typedef std::vector< typename Domain::Element > Vect_t;
        typedef typename Poly1Dom<Domain, Dense>::Element Element;


        Interpolation (const Domain& d, const Indeter& X = Indeter() ) : Poly1Dom<Domain,Dense>(d,X), Pi(this->one) {}

        void operator() ( const typename Domain::Element& x, const typename Domain::Element& f) {
            DD.push_back(f);
            Element M;
            this->init(M);

            if (DD.size() > 1) {
                // Pi *= (X-x_{n});
                this->mul(M, Pi, Points.back());
                this->shiftin(Pi, 1); // Pi *= X;
                this->subin(Pi, M);

                // Update divided differences
                typename Domain::Element tmp;
                this->_domain.init(tmp);
                // Adds the last evaluation
                for( typename Vect_t::reverse_iterator prev = DD.rbegin(), next = DD.rbegin(), point = Points.rbegin(); ++next != DD.rend(); ++prev, ++point)
                    // [x_i, ..., x_{n+1}] = ([x_i, ..., x_n] - [x_{i+1}, ..., x_{n+1}])
                    //			 / ( x_i - x_{n+1} )
                    this->_domain.divin(
                                        this->_domain.subin(*next, *prev),
                                        this->_domain.sub(tmp, *point, x)
                                       );
            }
            Points.push_back(x);


            // inter += [x_0, ..., x_{n+1}] (X-x_0)...(X-x_n)
            this->mul(M, Pi, DD.front());
            this->addin(inter, M);

        }

        Element& interpolator() {
            return inter;
        }

    private:

        Element inter;
        Element Pi;
        Vect_t Points, DD;
    };

} // Givaro

#endif // __GIVARO_interpolation_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
