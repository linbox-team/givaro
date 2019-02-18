// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givpoly1axpy.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: J-G. Dumas
// $Id: givpoly1axpy.inl,v 1.6 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================

#ifndef __GIVARO_poly1axpy_INL
#define __GIVARO_poly1axpy_INL

namespace Givaro {

    // axpy, axmy, maxpy
    // J.G.D. 16.11.2006
    // A lot can be done to optimize those
    // Except for axpy, axpyin, maxpy with a a Type_t,
    // all of them use a temporary vector where
    // a temporary value only would be sufficient.

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::axpy  (Rep& r, const Rep& a, const Rep& x, const Rep& y) const
    {
        return this->addin( this->mul(r,a,x), y );
    }


    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::axpy  (Rep& r, const Type_t& a, const Rep& x, const Rep& y) const
    {
        typename Rep::const_iterator ix = x.begin(), iy = y.begin();
        if (y.size() > x.size()) {
            r.resize(y.size());
            typename Rep::iterator ir = r.begin();
            for( ; ix != x.end(); ++ir, ++ix, ++iy)
                this->_domain.axpy(*ir, a, *ix, *iy);
            for( ; ir != r.end(); ++ir, ++iy)
                this->_domain.assign(*ir, *iy);
        } else {
            r.resize(x.size());
            typename Rep::iterator ir = r.begin();
            for( ; iy != y.end(); ++ir, ++ix, ++iy)
                this->_domain.axpy(*ir, a, *ix, *iy);
            for( ; ir != r.end(); ++ir, ++ix)
                this->_domain.mul(*ir, a, *ix);
        }
        return r;
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::axpyin(Rep& r, const Rep& a, const Rep& x) const
    {
        Rep tmp; this->init(tmp);
        this->assign(tmp,r);
        return this->axpy(r,a,x,tmp);
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::axpyin(Rep& r, const Type_t& a, const Rep& x) const{
        typename Rep::const_iterator ix = x.begin();
        if (x.size() > r.size()) {
            for(typename Rep::iterator ir = r.begin() ; ir != r.end(); ++ir, ++ix)
                this->_domain.axpyin(*ir, a, *ix);
            Type_t tmp;
            for( ; ix != x.end(); ++ix)
                r.push_back( this->_domain.mul(tmp, a, *ix) );
        } else {
            for(typename Rep::iterator ir = r.begin() ; ix != x.end(); ++ir, ++ix)
                this->_domain.axpyin(*ir, a, *ix);
        }
        return r;
    }
    // -- maxpy: r <- c - a * b
    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::maxpy  (Rep& r, const Rep& a, const Rep& b, const Rep& c) const{
        Rep tmp; this->init(tmp);
        return this->sub(r,c,this->mul(tmp,a,b));
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::maxpy  (Rep& r, const Type_t& a, const Rep& b, const Rep& c) const{
        size_t sC = c.size();
        size_t sB = b.size();
        size_t sR = r.size();
        if (sB == 0) { r.copy(c); return r; }
        if (sC == 0) { return this->negin( this->mul(r,a,b) ); }
        size_t i, max = sC < sB ? sB : sC;
        if (sR != max) r.resize(max);
        if (sC < sB)
        {
            for (i=0; i<sC; ++i) _domain.maxpy(r[i], a, b[i], c[i]);
            for (; i<sB; ++i) _domain.negin( _domain.mul(r[i], a, b[i]) );
        }
        else {
            for (i=0; i<sB; ++i) _domain.maxpy(r[i], a, b[i], c[i]);
            for (; i<sC; ++i) _domain.assign(r[i], c[i]);
        }
        return r;
    }

    // -- maxpyin: r -= a*b
    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::maxpyin(Rep& r, const Rep& a, const Rep& b) const{
        Rep tmp; this->init(tmp);
        return this->subin(r, this->mul(tmp,a,b));
    }
    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::maxpyin(Rep& r, const Type_t& a, const Rep& b) const{
        Rep tmp; this->init(tmp);
        return this->subin(r, this->mul(tmp,a,b));
    }
    // -- axmy: r <- a * x - y
    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::axmy  (Rep& r, const Rep& a, const Rep& x, const Rep& y) const{
        return this->subin(this->mul(r, a, x),y);
    }
    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::axmy  (Rep& r, const Type_t& a, const Rep& x, const Rep& y) const{
        return this->subin(this->mul(r, a, x),y);
    }
    // -- axmyin: r <- a * x - r
    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::axmyin (Rep& r, const Rep& a, const Rep& x) const
    {
        this->maxpyin(r, a, x);
        return this->negin(r);
    }
    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::axmyin (Rep& r, const Type_t& a, const Rep& x) const
    {
        this->maxpyin(r, a, x);
        return this->negin(r);
    }

} // Givaro

#endif // __GIVARO_poly1axpy_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
