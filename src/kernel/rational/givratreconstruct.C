// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/rational/givratreconstruct.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: Jean-Guillaume Dumas
// $Id: givratreconstruct.C,v 1.4 2010-04-12 15:54:39 jgdumas Exp $
// ==========================================================================
// Description:
#include "givaro/givrational.h"
#include <iostream>

namespace Givaro {

    // Let f,m \in Z
    // Find num, den such that
    // either num = den f mod m     (1)
    // or     num den^-1 = f mod m  (2)
    // verifying |num|<k and 0 <= den <= m/k
    //
    // List of remarks :
    // * it seems that forcereduce tells whether we solve (1) or (2)
    //   -> change variable name ?
    // * Rename, refactor e.g.
    //     ratrecon(7 args) ~= RationalReconstruction(7 args)
    //     rename ratrecon(5 args) or remove from interface ?
    // * What does RationalReconstruction(7 args) do
    //   that ratrecon(7 args) do not do
    // * Should the recursive aspect bool recurs be done by ratrecon ?

    bool Rational::ratrecon(
        Integer& num, Integer& den, const Integer& f, const Integer& m,
        const Integer& k, bool forcereduce, bool recurs)
    {

#ifdef GIVARO_RATRECON_DEBUG
        std::clog << "RatRecon : " << f << " mod " << m
                  << ", num bounded by: " << k
                  << ", reducedfraction: " << forcereduce
                  << ", expandboundwhenfailure: " << recurs
                  << std::endl;
#endif


        Integer r0, t0, r1, t1, q, u;
        r0=m;
        t0=0;
        r1=f;
        if (f<0) r1+= m;
        t1=1;
        while(r1>=k)
        {
#ifdef GIVARO_RATRECON_DEBUG
			std::clog << "r0:=" << r0 << ';' << std::endl
                      << "r1:=" << r1 << ';' << std::endl
                      << "q:=" << q << ';' << std::endl
                      << "t0:=" << t0 << ';' << std::endl
                      << "t1:=" << t1 << ';' << std::endl;
#endif
            q = r0;
            q /= r1;    // r0/r1

            u = r1;
            r1 = r0;	// r1 <-- r0
            r0 = u;	    // r0 <-- r1
            Integer::maxpyin(r1,u,q);

            u = t1;
            t1 = t0;	// r1 <-- r0
            t0 = u;	    // r0 <-- r1
            Integer::maxpyin(t1,u,q);
        }

        // [GG, MCA, 1999] Theorem 5.26
        // (i)
        if (t1 < 0) {
            num = -r1;
            den = -t1;
        } else {
            num = r1;
            den = t1;
        }

        if (forcereduce) {

                // (ii)
            if (gcd(num,den) != 1) {
#ifdef GIVARO_RATRECON_DEBUG
            std::clog << "num:=" << num << ';' << std::endl
                      << "den:=" << den << ';' << std::endl
                      << "g:=" << gcd(num,den) << ';'
                      << std::endl;
            std::clog << "r0:=" << r0 << ';' << std::endl
                      << "t0:=" << t0 << ';' << std::endl
                      << "k:=" << k << ';'
                      << std::endl;
#endif

                if (num == 0) {
                    if ((f % m) == 0) {
                        return true;
                    } else {
                        if (!recurs)
                            std::cerr
                                << "*** Error *** There exists no rational reconstruction of "
                                << f
                                << " modulo "
                                << m
                                << " with |numerator| < "
                                << k
                                << std::endl
                                << "*** Error *** But "
                                << num
                                << " = "
                                << den
                                << " * "
                                << f
                                << " modulo "
                                << m
                                << std::endl;
                        return false;
                    }
                }

                q = r0;
                q += r1;
                q -= k;
                q /= r1;

#ifdef GIVARO_RATRECON_DEBUG
            std::clog << "q:=" << q  << ';' << std::endl;
#endif

                r0 -= q * r1;
                t0 -= q * t1;

                if (t0 < 0) {
                    num = -r0;
                    den = -t0;
                } else {
                    num = r0;
                    den = t0;
                }

                if (t0 > m/k) {
                    if (!recurs)
                        std::cerr
                            << "*** Error *** No rational reconstruction of "
                            << f
                            << " modulo "
                            << m
                            << " with denominator <= "
                            << (m/k)
                            << std::endl;
                }
                if (gcd(num,den) != 1) {
                    if (!recurs)
                        std::cerr
                            << "*** Error *** There exists no rational reconstruction of "
                            << f
                            << " modulo "
                            << m
                            << " with |numerator| < "
                            << k
                            << std::endl
                            << "*** Error *** But "
                            << num
                            << " = "
                            << den
                            << " * "
                            << f
                            << " modulo "
                            << m
                            << std::endl;
                    return false;
                }
            }
        }
#ifdef GIVARO_RATRECON_DEBUG
        std::clog << "RatRecon End " << std::endl;
#endif
        return true;
    }

    bool Rational::ratrecon
    (const Integer& f, const Integer& m, const Integer& k,
     bool forcereduce, bool recurs) {
        return Rational::ratrecon(this->num, this->den,
                                  f,m,k,forcereduce,recurs);
    }

    Rational::Rational
    (const Integer& f, const Integer& m, const Integer& k,
     bool recurs ) {
        bool res = this->ratrecon(f,m,k,Rational::flags,recurs);
        if (recurs)
            for( Integer newk = k + 1; (!res) && (newk<f) ; newk<<=1)
                res = this->ratrecon(f,m,newk,Rational::flags,true);
    }

    bool Rational::RationalReconstruction
    (Integer& a, Integer& b, const Integer& f, const Integer& m,
     const Integer& k, bool forcereduce, bool recursive)
    {
        bool res(true);
        Integer x(f);
        if (x<0) {
            if ((-x)>m)
                x %= m;
            if (x<0)
                x += m;
        }
        else {
            if (x>m)
                x %= m;
        }

        if (x == 0) {
            a = 0;
            b = 1;
        }
        else {
            res = ratrecon(a,b,x,m,k, forcereduce, recursive);
            if (recursive)
                for( Integer newk = k + 1; (!res) && (newk<f) ; newk<<=1)
                    res = ratrecon(a,b,x,m,newk,forcereduce, true);
        }
        return res;
    }

    bool Rational::RationalReconstruction
    (Integer& a, Integer& b, const Integer& x, const Integer& m) {
        return ratrecon(a, b, x, m, Givaro::sqrt(m), true, true);
    }
    bool Rational::RationalReconstruction
    (Integer& a, Integer& b, const Integer& x, const Integer& m,
     const Integer& a_bound, const Integer& b_bound) {
        Integer bound = x/b_bound;
        ratrecon(a,b,x,m,
                 (bound>a_bound?bound:a_bound),
                 true, false);
        return b <= b_bound;
    }


} // namespace Givaro

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
