// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givpoly1addsub.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givpoly1addsub.inl,v 1.4 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================

#ifndef __GIVARO_poly_addsub_INL
#define __GIVARO_poly_addsub_INL

namespace Givaro {
    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::addin (Rep& R, const Rep& P) const
    {
        //     this->write(this->write(std::cout, R) << " += ", P) << std::endl;

        size_t i;
        size_t sP = P.size();
        size_t sR = R.size();
        if (sP == 0) return R;
        if (sR == 0) { return assign(R,P); }
        //   if (sR == 0) { R.copy(P); return R; }
        //   if (sR == sP){ _supportdomain.addin(R,P); return; }
        if (sR < sP) {
            Rep tmp; tmp = P;
            for (i=0; i<sR; ++i) _domain.addin(tmp[i], R[i]);
            R = tmp;
        }
        else {
            for (i=0; i<sP; ++i) _domain.addin(R[i], P[i]);
        }
        return R;
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::add(Rep& R, const Rep& P, const Rep& Q) const
    {
        size_t sP = P.size();
        size_t sQ = Q.size();
        size_t sR = R.size();
        if (sP == 0) { R = Q; return R; }
        if (sQ == 0) { R = P; return R; }
        // JGD 04.11.1999
        //   if (sP == sQ) {
        //     R.resize(sP);
        //     _supportdomain.add(R,P,Q);
        //     return;
        //   }
        size_t i, max = sP < sQ ? sQ : sP;
        if (sR != max) R.resize(max);
        if (sP < sQ)
        {
            for (i=0; i<sP; ++i) _domain.add(R[i], P[i], Q[i]);
            //     for (; i<sQ; ++i) _domain.assign(R[i], Q[i]);
            // JGD 05.11.1999
            for (; i<sQ; ++i) R[i] = Q[i];
        }
        else {
            for (i=0; i<sQ; ++i) _domain.add(R[i], P[i], Q[i]);

            //     for (; i<sP; ++i) _domain.assign(R[i], P[i]);
            // JGD 05.11.1999
            for (; i<sP; ++i) R[i] = P[i];
        }
        return R;
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::add
    (Rep& R, const Rep& P, const Type_t& Val) const
    {
        size_t sP = P.size();
        if (sP == 0)  {
            R.resize(1);
            _domain.assign(R[0],Val);
        }
        else {
            assign(R, P);
            _domain.add(R[0],P[0],Val);
        }
        return R;
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::add
    (Rep& R, const Type_t& Val, const Rep& P) const
    {
        size_t sP = P.size();
        if (sP == 0)  {
            R.resize(1);
            _domain.assign(R[0],Val);
        }
        else {
            assign(R, P);
            _domain.add(R[0],Val, P[0]);
        }
        return R;
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::addin
    (Rep& R, const Type_t& Val) const
    {
        size_t sR = R.size();
        if (sR == 0)  {
            R.resize(1);
            _domain.assign(R[0],Val);
        } else
            _domain.addin(R[0],Val);
        return R;
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::subin (
                                                                                Rep& R, const typename Rep::iterator Rbeg,
                                                                                const Rep& P, const typename Rep::const_iterator Pbeg, const typename Rep::const_iterator Pend) const
    {
        // PRECONDITION: NO reallocation, R MUST be of larger degree than P
        typename Rep::iterator ri=Rbeg;
        typename Rep::const_iterator pi=Pbeg;
        for( ; pi != Pend; ++pi, ++ri) _domain.subin(*ri,*pi);
        return R;
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::subin (
                                                                                Rep& R,
                                                                                const Rep& P, const typename Rep::const_iterator Pbeg, const typename Rep::const_iterator Pend) const
    {
        // PRECONDITION: P of larger degree than R
        size_t sP = (size_t)(Pend-Pbeg);
        size_t sR = R.size();
        R.resize(sP);
        typename Rep::iterator ri=R.begin(), Rend=ri+sR;
        typename Rep::const_iterator pi=Pbeg;
        for (; ri!=Rend; ++ri, ++pi) _domain.subin(*ri, *pi);
        for (; pi != Pend; ++ri, ++pi) _domain.neg(*ri, *pi);
        return setdegree(R);
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::subin (
                                                                                Rep& R, const typename Rep::iterator Rbeg, const typename Rep::iterator Rend,
                                                                                const Rep& P, const typename Rep::const_iterator Pbeg, const typename Rep::const_iterator Pend) const{
        size_t sP = (size_t) (Pend-Pbeg);
        size_t sR = (size_t)(Rend-Rbeg);
        if (sP == 0) return R;
        if (sR < sP) {
            return subin(R, P, Pbeg, Pend);
        }
        else {
            return subin(R, Rbeg, P, Pbeg, Pend);
        }
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::subin (Rep& R, const Rep& P) const
    {
        size_t sP = P.size();
        size_t sR = R.size();
        if (sP == 0) return R;
        if (sR == 0) { return neg(R,P); }
        if (sR < sP)
            return setdegree( subin(R, P, P.begin(), P.end()) );
        else
            return setdegree( subin(R, R.begin(), P, P.begin(), P.end()) );
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::sub(Rep& R, const Rep& P, const Rep& Q) const
    {
        size_t sP = P.size();
        size_t sQ = Q.size();
        size_t sR = R.size();
        if (sQ == 0) { R = P; return R; }
        if (sP == 0) { return neg(R,Q); }
        //   if (sP == sQ) {
        //     R.resize(sP);
        //     _supportdomain.sub(R,P,Q);
        //     return;
        //   }
        size_t i, max = sP < sQ ? sQ : sP;
        if (sR != max) R.resize(max);
        if (sP < sQ)
        {
            for (i=0; i<sP; ++i) _domain.sub(R[i], P[i], Q[i]);
            for (; i<sQ; ++i) _domain.neg(R[i], Q[i]);
        }
        else {
            for (i=0; i<sQ; ++i) _domain.sub(R[i], P[i], Q[i]);
            for (; i<sP; ++i) _domain.assign(R[i], P[i]);
        }
        return R;
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::sub
    (Rep& R, const Rep& P, const Type_t& Val) const
    {
        size_t sP = P.size();
        if (sP == 0)  {
            R.resize(1);
            _domain.neg(R[0],Val);
        }
        else {
            assign(R, P);
            _domain.sub(R[0],P[0],Val);
        }
        return R;
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::subin
    (Rep& R, const Type_t& Val) const
    {
        size_t sR = R.size();
        if (sR == 0)  {
            R.resize(1);
            _domain.neg(R[0],Val);
        } else
            _domain.subin(R[0],Val);
        return R;
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::sub
    (Rep& R, const Type_t& Val, const Rep& P) const
    {
        size_t sP = P.size();
        if (sP == 0)  {
            R.resize(1);
            _domain.neg(R[0],Val);
        }
        else {
            neg(R, P);
            _domain.add(R[0],Val, P[0]);
        }
        return R;
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::negin (Rep& R ) const
    {
        //     _supportdomain.negin(R);
        size_t sR = R.size();
        if (sR == 0) { return R; }
        for (size_t i=0; i<sR; ++i) _domain.negin(R[i]);
        return R;
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::neg (Rep& R, const Rep& P ) const
    {
        //   _supportdomain.neg(R,P);
        size_t sP = P.size();
        R.resize(sP);
        if (sP == 0) { return R; }
        for (size_t i=0; i<sP; ++i) _domain.neg(R[i],P[i]);
        return R;
    }

} // Givaro

#endif // __GIVARO_poly_addsub_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
