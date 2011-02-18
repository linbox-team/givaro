// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givinterp.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: JG Dumas
// $Id: givinterp.h,v 1.3 2011-02-02 16:23:56 bboyer Exp $
// ==========================================================================
// Description:
// Interpolation at geometric points
// see: Polynomial evaluation and interpolation on special sets of points,
// A. Bostan and E. Schost, Journal of Complexity 21(4): 420-446, 2005.

#ifndef __GIVARO_interpolation_at_geometric_points_H
#define __GIVARO_interpolation_at_geometric_points_H

#include "givaro/givconfig.h"
#include "givaro/giverror.h"
#include "givaro/givpoly1.h"
#include <givaro/givtruncdomain.h>

template<class Domain, bool REDUCE = true>
struct NewtonInterpGeom : TruncDom<Domain>  {
    typedef std::vector< typename Domain::Element > Vect_t;
    typedef typename TruncDom<Domain>::Polynomial_t Polynomial;
    typedef Polynomial Element;
    typedef typename TruncDom<Domain>::Element Truncated;
    typedef typename TruncDom<Domain>::Type_t Type_t;
//     using typename TruncDom<Domain>::Type_t;

    Type_t _q;
    Type_t powerq;
//     Vect_t _qi;
//     Vect_t _ui;
    Polynomial _qi;
    Polynomial _ui;
    Polynomial _mui;
    Polynomial _wi;
    Polynomial _g2;
    bool flip;
    Degree _deg;

    NewtonInterpGeom (const Domain& d, const Indeter& X = Indeter() ) : TruncDom<Domain>(d,X), powerq(d.one), flip(false), _deg(0) {
        d.generator(_q); // get a primitive root
    }

    template<typename BlackBox>
    void initialize (const BlackBox& bb) {
        _qi.resize(0);
        _ui.resize(0);
        _mui.resize(0);
        _wi.resize(0);
        _g2.resize(0);
        this->_domain.assign(powerq, this->_domain.one);
        _deg = 0;
        _qi.push_back(this->_domain.one);
        _ui.push_back(this->_domain.one);
        _mui.push_back(this->_domain.one);
        Type_t v0;
        bb(v0, this->_domain.one);
        _wi.push_back(v0);
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

        this->_domain.mulin(ui, _ui.back() );
        _ui.push_back(ui);
        
        Type_t vi;
        bb(vi, powerq);

        Type_t wi;
//         this->_domain.write(std::cout << "vi :", vi) << std::endl;
//         this->_domain.write(std::cout << "ui :", ui) << std::endl;
        this->_domain.div(wi, vi, ui);
//         this->_domain.write(std::cout << "wi :", wi) << std::endl;
        _wi.push_back(wi);
        
        Type_t gi;
//         this->_domain.write(std::cout << "qi :", qi) << std::endl;
//         this->_domain.write(std::cout << "ui :", ui) << std::endl;
        this->_domain.div(gi, qi, ui);
//         this->_domain.write(std::cout << "gi :", gi) << std::endl;
        if (flip) this->_domain.negin(gi);
//         this->_domain.write(std::cout << "gi :", gi) << std::endl;
        _g2.push_back(gi);

        flip = !flip;
        ++_deg;
    }
    


    Polynomial& Newton(Polynomial& inter) {
        this->getpoldomain().setdegree(_wi);
        this->getpoldomain().setdegree(_g2);

        Truncated G,W,QU;
        this->assign(W, _wi); // truncated
        this->assign(QU,_g2); // truncated

//         this->write(std::cout << "W : ", W) << std::endl;
//         this->write(std::cout << "G2: ", QU) << std::endl;
        
        

        this->mul(G, W, QU, 0, _deg);
//         this->write(std::cout << "G : ", G) << std::endl;
        

        this->convert(inter, G); // trunc to polynomial

        for(size_t i=0; i<inter.size(); ++i) {
//             this->_domain.write(std::cout << "BEF i[i]: ", inter[i]) << std::endl;
//             this->_domain.write(std::cout << "    q[i]: ", _qi[i]) << std::endl;
            this->_domain.divin(inter[i], _qi[i]);
//             this->_domain.write(std::cout << "AFT i[i]: ", inter[i]) << std::endl;
            
        }
            

//         this->getpoldomain().write(std::cout << "qi: ", _qi) << std::endl;
//         this->getpoldomain().write(std::cout << "ui: ", _ui) << std::endl;
//         this->getpoldomain().write(std::cout << "wi: ", _wi) << std::endl;
        

        return inter;
    }


    Polynomial& interpolator(Polynomial& inter) {
        this->Newton(inter);

        Polynomial mvi(_qi.size()), mwi(_qi.size()), mzi(_qi.size());
        for(size_t i=0; i<inter.size(); ++i) {
            this->_domain.mul(mvi[i],inter[i],_qi[i]);
            if (i & 1) this->_domain.negin(mvi[i]);
            this->_domain.div(mwi[i], mvi[i], _mui[i]);
            this->_domain.div(mzi[i], _mui[i], _qi[i]);
            if (i & 1) this->_domain.negin(mzi[i]);
        }

        this->getpoldomain().setdegree(_mui);
        this->getpoldomain().setdegree( mwi);

        this->getpoldomain().reversein(_mui);

        Truncated G,U,W;
        this->assign(U, _mui); // truncated
        this->assign(W, mwi); // truncated

        this->mul(G, U, W, _deg, _deg * 2);

        this->divin(G,_deg);

        this->convert(inter, G); // trunc to polynomial
        
        for(size_t i=0; i<inter.size(); ++i) {
            this->_domain.mulin(inter[i], mzi[i]);            
        }
        
        return inter;
    }
    
};



#endif // __GIVARO_interpolation_at_geometric_points_H
