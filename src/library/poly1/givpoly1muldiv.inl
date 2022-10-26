// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givpoly1muldiv.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givpoly1muldiv.inl,v 1.14 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================
#ifndef __GIVARO_poly1_muldiv_INL
#define __GIVARO_poly1_muldiv_INL
#include "givaro/givpower.h"
#include "givaro/giverror.h"

namespace Givaro {

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::shiftin ( Rep& R, int s) const
    {
        R.insert(R.begin(), s, this->_domain.zero );
        return R;
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::shift ( Rep& R, const Rep& a, int s) const
    {
        R = a;
        return R.shiftin(R, s);
    }


    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::mulin( Rep& R, const Type_t& u ) const
    {
        for(typename Rep::iterator ri = R.begin();ri!=R.end();++ri)
            _domain.mulin(*ri, u);
        return R;

        //  return _supportdomain.mulin(R,u);
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::mulin( Rep& R, const Rep& P ) const
    {
        size_t sR = R.size();
        size_t sP = P.size();
        Rep tmp(sR+sP);
        mul(tmp, R, P);
        //   R.logcopy(tmp);
        //   return R;
        return assign(R,tmp);
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::sqr( Rep& R, const Rep& P) const
    {
        const size_t sP = P.size();
        if (sP ==0) { R.resize(0); return R; }
        size_t sR = sP<<1;
        if (R.size() != --sR) R.resize(sR);

        // Generic square handler
        return sqr(R, R.begin(), R.end(), P, P.begin(), P.end());
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::mul( Rep& R, const Rep& P, const Rep& Q ) const
    {
        size_t sR = R.size();
        size_t sP = P.size();
        size_t sQ = Q.size();
        if ((sQ ==0) || (sP ==0)) { R.resize(0); return R; }
        if (sR != sP+sQ-1) R.resize(sR = sP+sQ-1);

        // Generic multiplication handler
        // Can use e.g. Karatsuba multiplication
        mul(R, R.begin(), R.end(),
            P, P.begin(), P.end(),
            Q, Q.begin(), Q.end());

        return setdegree(R);
    }

    // Compute truncated mul: only the coefficients inside the degree interval, included
    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::mul( Rep& R, const Rep& P, const Rep& Q, const Degree& Val, const Degree& deg) const
    {
        size_t sR = R.size();
        size_t sP = P.size();
        size_t sQ = Q.size();
        if ((sQ ==0) || (sP ==0)) { R.resize(0); return R; }
        size_t newS = (size_t)value(deg-Val)+1;
        if (sR != newS) R.resize(sR = newS);
        for(typename Rep::iterator ri=R.begin(); ri!= R.end(); ++ri)
            *ri = _domain.zero;

        for(size_t i=0; i<sR; ++i) {
            long k=(long)i+Val.value();
            size_t j=0;
            if (static_cast<size_t>(k)>=sQ) {
                j=static_cast<size_t>(k);
                k=long(sQ-1);
                j-=(size_t)k;
            }
            for( ; (j<sP) && (k>=0); ++j,--k) {
                _domain.axpyin(R[(size_t)i],P[(size_t)j],Q[(size_t)k]);
            }
        }
        return setdegree(R);
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::mul
    ( Rep& R, const Rep& P, const Type_t& u ) const
    {
        typename Rep::const_iterator ip = P.begin();
        R.resize(P.size());
        for(typename Rep::iterator ir = R.begin(); ir != R.end(); ++ir, ++ip)
            this->_domain.mul(*ir, *ip, u);
        return R;
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::mul
    ( Rep& R, const Type_t& u, const Rep& P ) const
    {
        return this->mul(R,P,u);
    }
} // Givaro

//#include <typeinfo>

namespace Givaro {


    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::modpowxin ( Rep& A, const Degree& l) const
    {
        A.resize(l.value());
        return setdegree(A); // A mod X^l
    }


    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::modpowx ( Rep& Am, const Rep& A, const Degree& l) const
    {
        assign(Am, A);
        return modpowxin(Am, l); // A mod X^l
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::newtoninviter ( Rep& G, Rep& S, Rep& Am, const Rep& A, const Degree& i) const
    {
            // Precondition i>=0
        sqr(S, G);						// G^2
        addin(G, G);					// 2G
        Am.resize(i.value());			// Compute only up to deg i
        mul(Am, Am.begin(), Am.end(),	// A * G^2
            A, A.begin(), A.begin() + std::min(static_cast<std::ptrdiff_t>(i.value()), A.end()-A.begin()),
            S, S.begin(), S.end());
        return subin(G, Am);			// 2G-AG^2 = G(2-AG), N-R iteration
    }


    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::invmodpowx ( Rep& G, const Rep& A, const Degree& l) const
    {
			// Precondition A is invertible
            // Precondition l>=0
        Rep S, Am; init(S); init(Am);
        S.reserve(l.value()); Am.reserve(l.value());

        assign(G, one);
        getdomain().inv(G[0],A[0]);				// Precondition A is invertible

        for(Degree i(2); i<l; i<<=1) {
            newtoninviter(G, S, Am, A, i);		// 2G-AG^2 mod X^i
        }

        return newtoninviter(G, S, Am, A, l);	// 2G-AG^2 mod X^l
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::divin(Rep& R, const Type_t& u) const
    {
#ifdef __GIVARO_DEBUG
        if (_domain.isZero(u)) GivError::throw_error(GivMathDivZero("[Poly1Dom<D>::divin]"));
#endif
        size_t sz =R.size();
        for (unsigned int i=0; i<sz; ++i)
            _domain.divin(R[i],u);
        return setdegree(R);
    }


    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::div(Rep& R, const Rep& P, const Type_t& u) const
    {
#ifdef __GIVARO_DEBUG
        if (_domain.isZero(u)) GivError::throw_error(GivMathDivZero("[Poly1Dom<D>::div]"));
#endif
        size_t sP =P.size();
        R.resize((size_t)sP);
        for (unsigned int i=0; i<sP; ++i)
            _domain.div(R[i],P[i],u);
        return setdegree(R);
    }


    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::div(Rep& R, const Type_t& u, const Rep& P) const
    {
#ifdef __GIVARO_DEBUG
        if (isZero(P)) GivError::throw_error(GivMathDivZero("[Poly1Dom<D>::div]"));
#endif
        if (_domain.isZero(u)) { return assign(R,zero);}
        size_t sP =P.size();
        if (sP >1) { R.resize(0); return R; }
        size_t sR =R.size();
        if (sR !=1) R.resize(1);
        _domain.div(R[0], u, P[0]);
        return setdegree(R);
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::div(Rep& Q, const Rep& A, const Rep& B) const
    {
        Degree degB; degree(degB, B);
#ifdef __GIVARO_DEBUG
        if (degB == Degree::deginfty)
            GivError::throw_error(GivMathDivZero("[Poly1Dom<D>::div]"));
#endif
        Degree degA; degree(degA, A);
        if ( degA < degB ) {
            return assign(Q, zero);
        }
        if (degB == 0) // cste
        {
            return div(Q, A, B[0]);
        }

            // Fast division: via multiplications
            // A = B Q + R, thus degQ = degA-degB
            // Thus X^a A = (X^b B)(X^q Q) + X^(q+1) (X^{q-1} R)
            // Thus rev(A) = rev(B) rev(Q) + X^(q+1) rev(R)
            // Thus rev(Q) = rev(A) rev(B)^{-1} mod X^(q+1)
            // Thus let degX = q+1 = a-b+1 = degA-degB+1
        Degree degX(degA); degX -= degB; ++degX;

        Rep T, S; init(T); init(S); S.reserve(degX.value()+1);
        reverse(T, B);
        invmodpowx(S, T, degX);	// rev(B)^{-1} mod X^l

        reverse(T, A);

        Q.resize(degX.value());		// deg(Q) = l - 1

        mul(Q, Q.begin(), Q.end(),	// rev(Q) = rev(A) rev(B)^{-1} mod X^l
            S, S.begin(), S.end(),
            T, T.begin(), T.end());

        return reversein(Q);
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::divin(Rep& Q, const Rep& A) const
    {
        Rep B;
        div(B,Q,A);
        return assign(Q,B);
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::invin (Rep& R ) const
    {
        Rep P; init(P);
        inv(P,R);
        return assign(R,P);
    }

    // Requires isUnit(P).
    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::inv (Rep& R, const Rep& P ) const
    {
        return div(R,one,P);
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::modin(Rep& R, const Type_t& u) const
    {
#ifdef __GIVARO_DEBUG
        if (_domain.isZero(u)) GivError::throw_error(GivMathDivZero("[Poly1Dom<D>::modin]"));
#endif
        R.resize(0);
        return R;
    }


    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::mod(Rep& R, const Rep& P, const Type_t& u) const
    {
#ifdef __GIVARO_DEBUG
        if (_domain.isZero(u)) GivError::throw_error(GivMathDivZero("[Poly1Dom<D>::mod]"));
#endif
        R.resize(0);
        return R;
    }


    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::mod(Rep& R, const Type_t& u, const Rep& P) const
    {
#ifdef __GIVARO_DEBUG
        if (isZero(P)) GivError::throw_error(GivMathDivZero("[Poly1Dom<D>::mod]"));
#endif
        size_t sP =P.size();
        if (sP >1) {
            R.resize(1);
            _domain.assign(R[0], u);
            return R;
        }
        R.resize(0); // else deg(R)<deg(P)=0 implies R=0
        return R;
    }


    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::modin(Rep& A, const Rep& B) const
    {
        // In place remainder
        // A is written with next remainder in
        // the division algorithm written at the end.
        // Last step is erasing of the first values.
        //     write(std::cerr << "Rem(", A) << " ,";
        //     write(std::cerr, B) << ", X) mod " << _domain.size();
        long i = (long)(A.size()-B.size());
        if (i >= 0) {
            typedef typename Rep::value_type TT;
            TT l;
            typename Rep::reverse_iterator ai,aai;
            typename Rep::const_reverse_iterator bi;
            for (; i>=0; --i) {
                ai = A.rbegin();
                bi = B.rbegin();
                _domain.div(l,*ai,*bi);
                aai = A.rbegin();
                for(++bi,++ai;bi!=B.rend();++bi,++ai,--i) {
                    _domain.maxpy(*aai,l,*bi,*ai);
                    if (! _domain.isZero(*aai)) break;
                }
                if (bi!=B.rend())
                    for(++bi,++ai,++aai;bi!=B.rend();++bi,++ai,++aai)
                        _domain.maxpy(*aai,l,*bi,*ai);
                for(;ai!=A.rend();++ai,++aai)
                    *aai = *ai;
                *aai = _domain.zero;
            }
            //         write(std::cerr << " = ", A) << ";" << std::endl;
            A.erase(A.begin(), A.begin()+(ssize_t)(A.size()-B.size()-(size_t)i));
        }
        //     write(std::cerr << " = ", setdegree(A)) << ";" << std::endl;
        return setdegree(A);
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::mod(Rep& R, const Rep& A, const Rep& B) const
    {
        Rep Q;
// write(std::cerr, A) << " - (";
// write(std::cerr, B) << ") * (";
        divmod(Q,R,A,B);
// write(std::cerr, Q) << ") + (";
// write(std::cerr, R) << ") mod " << characteristic() << ';' << std::endl;
        return R;
    }

    // #include <typeinfo>

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::divmod( Rep& Q, Rep& R, const Rep& A,  const Rep& B) const
    // returns Q such that A = B Q + R
    {
        div(Q,A,B);
        return maxpy(R,Q,B,A); // R <-- A - Q * B
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::divmodin( Rep& Q, Rep& R, const Rep& B) const
    // returns Q such that R = B Q + newR
    {
        div(Q, R, B);
        return maxpyin(R, Q, B);
    }


    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::pdivmod
    ( Rep& Q, Rep& R, Type_t& m, const Rep& A, const Rep& B) const
    // returns Q ...
    {
        Degree degB; degree(degB, B);
#ifdef __GIVARO_DEBUG
        if (degB == Degree::deginfty)
            GivError::throw_error(GivMathDivZero("[Poly1Dom<D>::div]"));
#endif
        Degree degA; degree(degA, A);
        if (degA == Degree::deginfty) {
            assign(R, zero);
            _domain.assign(m, _domain.one);
            return assign(Q, zero);
        }
        if (degB == 0) // cste
        {
            assign(R, zero);
            _domain.assign(m, B[0]);
            return assign(Q, A);
        }
        if (degA ==0)
        {
            assign(R, zero);
            _domain.assign(m, _domain.one);
            return assign(Q, zero);
        }
        if (degB > degA) {
            assign(R, A);
            _domain.assign(m, _domain.one);
            return assign(Q, zero);
        }

        long degQuo = value(degA-degB);
        long degRem = value(degA);
        Q.resize((size_t)degQuo+1);
        assign(R,A);

        Type_t tmp, lB;
        _domain.assign(lB, B[degB.value()]);
        _domain.assign(m, _domain.one);
        long i,j;
        for (i=degQuo; i>=0; --i)
        {
            // == ld X^ (degRem-degQ)
            _domain.assign(Q[degQuo], R[degRem]);

            // rem <- lB*rem - lQ*x^(degRem-degB)*B
            for (j=0; j<degQuo; j++)
                _domain.mulin (R[j], lB);
            for (j=0; degB>j; j++)
            {
                _domain.mulin(R[j+degQuo], lB);
                _domain.maxpyin(R[j+degQuo], Q[degQuo], B[j]);
            }
            _domain.assign(R[degRem],_domain.zero); degQuo--; degRem--;
            _domain.mulin(m, lB);
        }
        R.resize((size_t)degRem+1);
        setdegree(R);
        return setdegree(Q);
        //  Poly1Dom<Domain,Dense>::Rep U,V;
        //  assign(U,A);
        //  mulin(U,m);
        //  write(std::cout << "m*A:", U) << std::endl;
        //  mul(U,Q,B);
        //  write(std::cout << "Q*B:", U) << std::endl;
    }



    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::pmod
    ( Rep& R, Type_t& m, const Rep& A, const Rep& B) const
    {
        Degree degB; degree(degB, B);
#ifdef __GIVARO_DEBUG
        if (degB == Degree::deginfty)
            GivError::throw_error(GivMathDivZero("[Poly1Dom<D>::div]"));
#endif
        Degree degA; degree(degA, A);
        if (degA == Degree::deginfty) {
            _domain.assign(m, _domain.one);
            return assign(R, zero);
        }
        if (degB == 0) // cste
        {
            _domain.assign(m, B[0]);
            return assign(R, zero);
        }
        if (degA ==0)
        {
            _domain.assign(m, _domain.one);
            return assign(R, zero);
        }
        if (degB > degA) {
            _domain.assign(m, _domain.one);
            return assign(R, A);
        }

        Degree degR = degA;
        assign(R,A);

        Type_t tmp, lB;
        _domain.assign(lB, B[degB.value()]);
        //write(std::cout << "B:", B) << std::endl;
        //_domain.write(std::cout << "lB:", lB) << "^" << degA-degB+1 << std::endl;
        //   _domain.pow(m, lB, degA.value()-degB.value()+1);
        dom_power(m, lB, degA.value()-degB.value()+1,_domain);
        //_domain.write(std::cout << "m:", m) << std::endl;
        for (; degB<= degR; )
        {
            long d = degR.value()-degB.value();
            // R <- lB*R - lR*x^(degR-degB)*B
            for (long j=0; degB>j; j++)
            {
                _domain.mulin (R[j+d], lB);
                _domain.maxpyin(R[j+d], R[degR.value()], B[j]);
            }
            for (long j=0; j<d; ++j)
                _domain.mulin (R[j], lB);
            _domain.assign(R[degR.value()],_domain.zero);
            degree(degR, R);
        }
        R.resize((size_t)degR.value()+1);
        return setdegree(R);
    }

} // Givaro
#endif // __GIVARO_poly1_muldiv_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
