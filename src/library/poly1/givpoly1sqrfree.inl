// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givpoly1sqrfree.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givpoly1sqrfree.inl,v 1.9 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================
// Description:
#ifndef __GIVARO_poly1_sqrfree_INL
#define __GIVARO_poly1_sqrfree_INL

namespace Givaro {

    template <class Domain>
    size_t& Poly1Dom<Domain,Dense>::sqrfree(size_t& Nfact, Rep* Fact, const Rep& P) const
    {
        if (Nfact ==0) return Nfact;
        GIVARO_ASSERT( (Fact !=0), "nul pointer");

        unsigned long count = 0;
        Rep A,B,C,D,W,Z,Y;
        Type_t lc;
        leadcoef(lc, P);
        init(Fact[count], 0, lc);
        // write(cout << "P:", P) << endl;
        // _domain.write(cout << "lc(P):", lc) << endl;
        assign(A, P);
        diff(B, A);
        // write(cout << "P':", B) << endl;
        gcd(D, A, B);
        // write(cout << "Gcd(P,P'):", D) << endl;
        div(A, P, lc);
        // write(cout << "A/lc:", A) << endl;
        leadcoef(lc, D);
        // _domain.write(cout << "lc(Gcd):", lc) << endl;
        div(C, D, lc);
        // write(cout << "Gcd/lc:", C) << endl;

        if (areEqual(C, one)) {
            assign(W, A);
        }
        else {
            div(W, A, C);
            // write(cout << "W:", W) << endl;
            div(Y, B, C);
            // write(cout << "Y:", Y) << endl;
            diff(Z, W);
            // write(cout << "W':", Z) << endl;
            sub(Z, Y, Z);
            while (!isZero(Z))
            {
                gcd(Fact[count], W, Z);
                // write(cout << "L" << count << ":", Fact[count]) << endl;

                div(C, W, Fact[count]); assign(W, C);
                div(Y, Z, Fact[count]);
                diff(Z, W);
                sub(Z, Y, Z);
                if (++count > Nfact) return Nfact;
            }
        }
        assign(Fact[count], W);
        //write(cout << "L" << count << ":", Fact[count]) << endl;
        return Nfact = ++count;
    }

} // Givaro

#endif // __GIVARO_poly1_sqrfree_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
