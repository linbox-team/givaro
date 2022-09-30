// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givpoly1misc.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givpoly1misc.inl,v 1.22 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================
// Description:

#include <algorithm>

#ifndef __GIVARO_poly_misc_INL
#define __GIVARO_poly_misc_INL

namespace Givaro {
    template<class Domain>
    inline int Poly1Dom<Domain,Dense>::isZero (const Rep& P) const
    {
        setDegree(const_cast<Rep&>(P));
        if (P.size() ==0) return 1;
        if (P.size() ==1) return _domain.isZero(P[0]);
        else return 0;
    }

    template<class Domain>
    inline int Poly1Dom<Domain,Dense>::isOne ( const Rep& P ) const
    {
        setDegree(const_cast<Rep&>(P));
        if (P.size() ==1)
            return _domain.isOne(P[0]);
        else
            return 0;
    }

    template<class Domain>
    inline int Poly1Dom<Domain,Dense>::isMOne ( const Rep& P ) const
    {
        setDegree(const_cast<Rep&>(P));
        if (P.size() ==1)
            return _domain.isMOne(P[0]);
        else
            return 0;
    }

    template<class Domain>
    inline bool Poly1Dom<Domain,Dense>::isUnit ( const Rep& P ) const
    {
        setDegree(const_cast<Rep&>(P));
        if (P.size() ==1)
            return _domain.isUnit(P[0]);
        else
            return 0;
    }

    template<class Domain>
    inline int Poly1Dom<Domain,Dense>::areEqual (const Rep& P, const Rep& Q) const
    {
        setDegree(const_cast<Rep&>(P));
        setDegree(const_cast<Rep&>(Q));
        if (P.size() != Q.size()) return 0;
        for( typename Element::const_iterator pit = P.begin(), qit = Q.begin();
             pit != P.end();
             ++pit, ++qit)
            if ( !_domain.areEqual(*pit, *qit) ) return 0;
        return 1;
    }


    template<class Domain>
    inline int Poly1Dom<Domain,Dense>::areNEqual (const Rep& P, const Rep& Q) const
    {
        setDegree(const_cast<Rep&>(P));
        setDegree(const_cast<Rep&>(Q));
        if (P.size() != Q.size()) return 1;
        for( typename Element::const_iterator pit = P.begin(), qit = Q.begin();
             pit != P.end();
             ++pit, ++qit)
            if ( !_domain.areEqual(*pit, *qit) ) return 1;
        return 0;
    }


    // -- Compute the degree of P
    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::setdegree ( Rep& P ) const
    {
        int sz = (int) (P.size() - 1);
        if (P.size() <= 0) {P.resize(0); return P;}
        if (_domain.isZero(P[(size_t)sz]) ==0) {
            return P;
        }
        for (int j=sz; j--; )
            if (_domain.isZero(P[(size_t)j]) ==0) {
                P.resize((size_t)j+1);
                return P;
            }
        P.resize(0);
        return P;
    }

    template <class Domain>
    inline Degree& Poly1Dom<Domain,Dense>::degree(Degree& deg, const Rep& P) const
    {
        int sz = (int) P.size();
        if (sz ==0) {
            return deg = Degree::deginfty;
        }
        if (_domain.isZero(P[(size_t)sz-1])) {
            setDegree(const_cast<Rep&>(P));
            sz = (int) P.size();
        }
        return deg = (Degree) (sz-1);
    }

    template <class Domain>
    inline Degree Poly1Dom<Domain,Dense>::degree(const Rep& P) const
    {
        Degree d; return degree(d,P);
    }


    template<class Domain>
    inline typename Poly1Dom<Domain,Dense>::Type_t& Poly1Dom<Domain,Dense>::leadcoef (Type_t& c, const Rep& P) const
    {
        Degree dP;
        degree(dP, P);
        if (dP == Degree::deginfty)
            return _domain.assign(c, _domain.zero);
        else
            return _domain.assign(c, P[(size_t)dP.value()]);
    }


    // -- Returns the i-th coefficients
    template<class Domain>
    inline typename Poly1Dom<Domain,Dense>::Type_t& Poly1Dom<Domain,Dense>::getEntry (Type_t& c, const Degree& i, const Rep& P) const
    {
        Degree dP;
        degree(dP, P);
        if (dP < i) return _domain.assign(c, _domain.zero);
        else return _domain.assign(c, P[i.value()]);
    }

    template<class Domain>
    inline typename Poly1Dom<Domain,Dense>::Type_t
    Poly1Dom<Domain,Dense>::setEntry(Rep &P, const Type_t&c, const Degree&i) const
    {
        Degree dP;
        degree(dP, P);

        if (_domain.isZero(c)) {
            if (dP < i) { /* nothing happens */
                return c ;
            }
            else if (dP == i) { /* degree is killed */
                _domain.assign(P[i.value()], c);
                setDegree(P);
                return c ;
            }
            else { /* element is killed */
                return _domain.assign(P[i.value()], c);
            }
        }
        /*  c != 0 */
        if (dP < i) {
            P.resize(i.value()+1);
        }
        return _domain.assign(P[i.value()], c);
    }




    template <class Domain>
    inline Degree& Poly1Dom<Domain,Dense>::val(Degree& d, const Rep& P) const
    {
        size_t sz = P.size();
        if (sz ==0) {
            return d = Degree::deginfty;
        }
        if (!_domain.isZero(P[0])) {
            return d = 0;
        }
        for (size_t i=1; i<sz; ++i)
        {
            if (!_domain.isZero(P[i])) {
                // return d = (Degree)i;
                return d = (uint64_t)i;
            }
        }
        return d=0;
    }


    // Horner's scheme for evaluation
    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Type_t& Poly1Dom<Domain,Dense>::eval (Type_t& res, const Rep& P, const Type_t& Val) const
    {
        typename Domain::Element tmp;
        _domain.init(tmp,0U);
        Degree dP ; degree(dP, P);
        if (dP == Degree::deginfty) _domain.assign(res, _domain.zero);
        else {
            _domain.assign(res, P[(size_t)dP.value()]);
            for (int i = (int)dP.value(); i--; )
                _domain.assign(res,_domain.axpy(tmp, res, Val, P[(size_t)i]));
        }
        return res;
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::diff(Rep& P, const Rep& Q) const
    {
        Degree dQ;
        degree(dQ, Q);
        if ((dQ == Degree::deginfty) || (dQ == 0)) {
            P.resize(0);
            return P;
        }
        P.resize((size_t)dQ.value());
        Type_t cste; _domain.assign(cste, _domain.zero);
        for (int i=0; dQ>i; ++i) {
            _domain.add(cste, cste, _domain.one);
            _domain.mul(P[(size_t)i], Q[(size_t)i+1], cste);
        }
        return P;
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::pow( Rep& W, const Rep& P, uint64_t p) const
    {
        // TODO: manage a negative exponent ...
        Rep puiss2;
        assign(puiss2, P); // -- P**(2^k)
        Rep tmp;
        assign(W,one);
        while (p != 0) {
            if (p & 0x1) {
                mul( tmp, W, puiss2);
                assign(W,tmp);
            }
            if ((p >>= 1) != 0) {
                mul( tmp, puiss2, puiss2);
                assign(puiss2,tmp);
            }
        }
        return W;
    }


    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::powmod( Rep& W, const Rep& P, IntegerDom::Element pwr, const Rep& U) const
    {
        IntegerDom ID;
        // ID.write(cerr << "\n----------- POWMOD -----------\n pwr: ", pwr) << endl;
        // write(cerr << "P: ",P) << endl;
        // write(cerr << "U: ",U) << endl;
        Rep puiss, tmp;
        mod(puiss, P, U);
        assign(W,one);

        Integer n(pwr);
        if (n<0) {
            std::cerr << "Powering with negative exponent not implemented" << std::endl;
            n = -n;
        }
        while(n>0) {
            if (n & 1U) {
                mulin(W,puiss);
                modin(W,U);
            }
            sqr(tmp,puiss);
            mod(puiss,tmp,U);
            n >>= 1;
        }

        // write(cerr << "W: ", W) << "\n----------- END POWMOD -----------" <<  endl;
        return setDegree(W);
    }


    // -- Random dense polynomial of degree 0
    template <class Domain> template<class RandomIterator>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::random(RandomIterator& g, Rep& r) const
    {
        return random(g, r,Degree(0));
    }

    // -- Random dense polynomial of size s
    template <class Domain> template<class RandomIterator>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::random(RandomIterator& g, Rep& r, uint64_t s) const
    {
        return random(g, r,Degree(s-1));
    }


    // -- Random dense polynomial of degree d
    template <class Domain> template<class RandomIterator>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::random(RandomIterator& g, typename Poly1Dom<Domain,Dense>::Rep& r, Degree d) const
    {
        r.resize((size_t)d.value()+1);
        _domain.nonzerorandom(g, r[(size_t)d.value()]);
        for (int i=(int)d.value(); i--;)
            _domain.random(g,r[(size_t)i]);
        return r;
    }


#if 0
    template <class Domain> template<class RandomIterator>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::random(RandomIterator& g, Rep& r, Degree d) const
    {
        r.resize((uint64_t)d.value()+1);
        while (_domain.isZero(_domain.init(r[d.value()], g()))) {};
        for (int i=d.value(); i--;)
            _domain.init(r[(uint64_t)i],g());
        _domain.nonzerorandom(g, r[d.value()]);
        for (int i=d.value(); i--;)
            _domain.random(g,r[(uint64_t)i]);
        return r;
    }
#endif

    // -- Random dense polynomial with same size as b.
    template <class Domain> template<class RandomIterator>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::random(RandomIterator& g, Rep& r, const Rep& b) const
    {
        return random(g, r, b.size());
    }


    template <class Domain> template<class RandomIterator>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::nonzerorandom(RandomIterator& g, Rep& r) const
    {
        return random(g, r);
    }


    template <class Domain> template<class RandomIterator>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::nonzerorandom(RandomIterator& g, Rep& r, uint64_t s) const
    {
        return random(g, r, s);
    }

    template <class Domain> template<class RandomIterator>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::nonzerorandom(RandomIterator& g, Rep& r, Degree d) const
    {
        return random(g, r, d);
    }

    template <class Domain> template<class RandomIterator>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::nonzerorandom(RandomIterator& g, Rep& r, const Rep& b) const
    {
        return random(g, r, b);
    }




    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::reverse( Rep& P, const Rep& Q) const {

        P.resize(Q.size());
        std::reverse_copy(Q.begin(), Q.end(), P.begin());
        this->setDegree(P);
        return P;
    }



    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::reversein( Rep& P) const
    {
        std::reverse(P.begin(), P.end());
        this->setDegree(P);
        return P;
    }

} // Givaro

#endif // __GIVARO_poly_misc_INL

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
