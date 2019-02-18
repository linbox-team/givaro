// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givpoly1cstor.inl,v $
// Copyright(c)'1994-2017 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givpoly1cstor.inl,v 1.15 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================
#ifndef __GIVARO_poly1_cstor_INL
#define __GIVARO_poly1_cstor_INL

namespace Givaro {
    template<class Domain>
    inline Poly1Dom<Domain,Dense>::Poly1Dom(const Domain& d, const Indeter& X ) :
        _domain(d), _x(X)
        , zero(1,d.zero), one(1,d.one)
        , mOne(1,d.mOne)
    {}

    template<class Domain>
    inline Poly1Dom<Domain,Dense>::Poly1Dom(const Self_t& P) :
        _domain(P._domain), _x(P._x)
        ,zero(P.zero), one(P.one),mOne(P.mOne)
    {}

    template<class Domain>
    inline Poly1Dom<Domain,Dense>::~Poly1Dom()
    {
    }


    template<class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::init( Rep& P) const
    { P.resize(0); return P; }

    /*
       template<class Domain>
       inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::init( Rep& P, const Rep& Q ) const
       {
       Degree degQ;
       degree(degQ,Q);
       if (degQ <0) {
       P.resize(0);
       return P;
       }
       P.resize(++degQ);
       for (int i=0; degQ>i; ++i)
       _domain.init(P[i], Q[i]);
       return P;
       }
       */

    template<class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::assign( Rep& P, const Rep& Q ) const
    {
        Degree degQ;
        degree(degQ,Q);
        if (degQ <0) {
            P.resize(0);
            return P;
        }
        P.resize((size_t)++degQ);
        // degQ >=0
        for (size_t i=0; (size_t)degQ.value()>i; ++i)
            _domain.assign(P[i], Q[i]);
        return P;
    }


    template<class Domain> template<class XXX>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::init( Rep& P, const XXX& Val ) const
    {
        P.resize(1);
        _domain.init(P[0], Val);
        return P;
    }

    template<class Domain> template<class XXX>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::init( Rep& P, const std::initializer_list<XXX>& L ) const
    {

        P.resize(L.size());
        auto il(L.begin());
        for(auto &ip : P) {
            _domain.init(ip,*il);
            ++il;
        }
        return P;
    }


    template<class Domain>
    inline typename Poly1Dom<Domain,Dense>::Type_t& Poly1Dom<Domain,Dense>::assign(typename Poly1Dom<Domain,Dense>::Type_t & Val, const typename Poly1Dom<Domain,Dense>::Rep& P ) const
    {
        if (P.size())
            return _domain.assign(Val, P[0]);
        else
            return _domain.assign(Val, _domain.zero);
    }

    template<class Domain> template<class XXX>
    inline XXX& Poly1Dom<Domain,Dense>::convert( XXX& Val, const typename Poly1Dom<Domain,Dense>::Rep& P ) const
    {
        if (P.size())
            return _domain.convert(Val, P[0]);
        else
            return _domain.convert(Val, _domain.zero);
    }

    template<class Domain> template<class UU, template<class XX> class Vect>
    inline Vect<UU>& Poly1Dom<Domain,Dense>::convert( Vect<UU>& Val, const typename Poly1Dom<Domain,Dense>::Rep& P ) const
    {
        Val.resize( P.size() );
        typename Vect<UU>::iterator vit = Val.begin();
        typename Rep::const_iterator        pit = P.begin();
        for ( ; pit != P.end(); ++pit, ++vit)
            _domain.convert(*vit, *pit);
        return Val;
    }



    template<class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::init( Rep& P, const Degree deg ) const
    {
        P.resize((size_t)value(deg+1));

        size_t sz = P.size();
        for (size_t i=0; i<sz-1; ++i)
            _domain.assign(P[i], _domain.zero);
        _domain.assign(P[sz-1], _domain.one);
        return P;
    }

    template<class Domain> template<class XXX>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::init
    ( Rep& P, const Degree d, const XXX& Val ) const
    {
        long deg = value(d);
        P.resize((size_t)deg+1);
        assert (deg>=0);
        for (size_t i=0; i<(size_t)deg; ++i)
            _domain.assign(P[i], _domain.zero);
        _domain.init(P[(size_t)deg], Val);

        if (_domain.isZero(P[(size_t)deg])) {
            P.resize(0);
        }
        return P;
    }


    template<class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::assign
    ( Rep& P, const Degree d, const Type_t& lcoeff ) const
    {
        long deg = value(d);
        if (_domain.isZero(lcoeff)) {
            P.resize(0);
        }
        else {
            P.resize((size_t)deg+1);
            assert (deg>=0);
            for (size_t i=0; i<(size_t)deg; ++i)
                _domain.assign(P[i], _domain.zero);
            _domain.assign(P[(size_t)deg], lcoeff);
        }
        return P;
    }

} // Givaro
#endif // __GIVARO_poly1_cstor_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
