// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givpoly1gcd.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Author: J-L. Roch, T. Gautier, J-G. Dumas
// $Id: givpoly1gcd.inl,v 1.11 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================
#ifndef __GIVARO_poly1_gcd_INL
#define __GIVARO_poly1_gcd_INL

namespace Givaro {

    // computes d = unitary gcd(P,Q) and u,v such that :
    //   u*P+v*Q = d
    // u and v are the unique polynomisals such that :
    //   deg(u) <= deg(Q)-deg(d) and deg(v) <= deg(P) - deg(d)
    // friend Poly1<T> gcd (const Poly1<T> &P, const Poly1<T> &Q );
    // computes d = unitary gcd(P,Q)
    // ===========================================================================
    // A coder par PRS : cf Geddes page 283
    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::gcd ( Rep& G, const Rep& P, const Rep& Q) const
    {

        Rep U;
        Degree degU, degG;
        degree(degU,P);
        degree(degG,Q);
        if ((degU < 0) || (degG == 0)) return assign(G, Q);
        if ((degG < 0) || (degU == 0)) return assign(G, P);

        if (degU >= degG) {
            assign(U, P);
            assign(G, Q);
        }
        else {
            assign(U, Q);
            assign(G, P);
        }
        // -- PRS (U,G) using pmod:
        Type_t g;
        _domain.assign(g, _domain.one);
        Degree degR;
        Rep R;
        do {
            // write(cout << "gcd: U:", U) << endl;
            // write(cout << "gcd: G:", G) << endl;
            mod( R, U, G);
            setdegree(R);
            degree(degR, R);
            // write(cout << "mod: R:", R) << endl;
            if (degR < 0) break;
            assign(U,G);
            assign(G,R);
        } while (1);

        degree(degG, G);
        if (degG <= 0) assign(G,_domain.one);
        // write(cout << "gcd: G:", G) << endl;
        return G;
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::gcd ( Rep& F, Rep& S0, Rep& T0, const Rep& A, const Rep& B) const
    {
        Type_t r0, r1, tt; _domain.init(r0);_domain.init(r1);_domain.init(tt);
        Rep G; init(G);
        Degree degF, degG;
        degree(degF,A); degree(degG,B);
        if ((degF < 0) || (degG == 0)) {
            _domain.inv( tt, leadcoef(r0,B)); assign(T0, tt);
            assign(S0, zero);
            assign(F, B);
            return mulin(F,tt);
        }
        if ((degG < 0) || (degF == 0)) {
            _domain.inv( tt, leadcoef(r0,A)); assign(S0, tt);
            assign(T0, zero);
            assign(F, A);
            return mulin(F,tt);
        }


        //   if (degF >= degG) {
        assign(F, A);
        assign(G, B);
        //   }
        //   else {
        //     assign(F, B);
        //     assign(G, A);
        //   }

        leadcoef(r0,  F);
        leadcoef(r1, G);

        divin(F,r0);
        divin(G,r1);

        Rep S1,R1,T1,Q,TMP, TMP2;
        init(S1);init(R1);init(T1);init(Q);init(TMP);init(TMP2);


        _domain.inv(tt,r0); assign(S0, tt);
        assign(S1,zero);
        assign(T0,zero);
        _domain.inv(tt,r1); assign(T1, Degree(0), tt);

        while ( ! isZero(G) ) {
            divmod(Q,R1,F,G);
            leadcoef(r1, R1); if (_domain.isZero(r1)) _domain.assign(r1,_domain.one);
            assign(F,G);
            div(G,R1,r1);
            mul(TMP,Q,S1); sub(TMP2,S0,TMP); assign(S0,S1); div(S1,TMP2,r1);
            mul(TMP,Q,T1); sub(TMP2,T0,TMP); assign(T0,T1); div(T1,TMP2,r1);
        }

        return F;
    }


    // #include <typeinfo>

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::invmod ( Rep& S0, const Rep& A, const Rep& B) const
    {
        //     std::cerr << "BEG invmod of " << typeid(*this).name() << std::endl;
        //     std::cerr << "with domain " << typeid(_domain).name() << std::endl;
        //        write(std::cout << "A:=", A) << ';' << std::endl;
        //        write(std::cout << "B:=", B) << ';' << std::endl;
        Type_t r0, r1, tt;
        Rep F, G;
        Degree degF, degG;
        degree(degF,A); degree(degG,B);
        if ((degF <= 0) || (degG <= 0) ) {
            return assign(S0, Degree(0), _domain.inv( tt, leadcoef(r0,A)));
        }

        assign(F, A);
        assign(G, B);

        leadcoef(r0, F);
        leadcoef(r1, G);
        //       getdomain().write(std::cout << "r0: ", r0) << std::endl;
        //       getdomain().write(std::cout << "r1: ", r1) << std::endl;

        divin(F,r0);
        divin(G,r1);
        //        write(std::cout << "F1: ", F) << std::endl;
        //        write(std::cout << "G1: ", G) << std::endl;

        Rep S1,R1,Q,TMP, TMP2;

        assign(S0, Degree(0), _domain.inv(tt,r0) );
        assign(S1, zero);

        while ( ! isZero(G) ) {
            //        write(std::cout << "F: ", F) << std::endl;
            //        write(std::cout << "G: ", G) << std::endl;

            divmod(Q,R1,F,G);
            //       write(std::cout << "Q: ", Q) << std::endl;
            //       write(std::cout << "R: ", R1) << std::endl;

            leadcoef(r1, R1); if (_domain.isZero(r1)) _domain.assign(r1,_domain.one);
            //       getdomain().write(std::cout << "l: ", r1) << std::endl;

            assign(F,G);
            div(G,R1,r1);
            //       write(std::cout << "Fn: ", F) << std::endl;
            //       write(std::cout << "Gn: ", G) << std::endl;
            mul(TMP,Q,S1); sub(TMP2,S0,TMP); assign(S0,S1); div(S1,TMP2,r1);
            //       write(std::cout << "S: ", S1) << std::endl;
        }

        //       write(std::cout << "S: ", S0) << std::endl;
        //     std::cerr << "END invmod of " << typeid(*this).name() << std::endl;
        return S0;
    }

    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::invmodunit ( Rep& S0, const Rep& A, const Rep& B) const
    {
        //     std::cerr << "BEG invmodunit of " << typeid(*this).name() << std::endl;
        Rep F, G;
        Degree degF, degG;
        degree(degF,A); degree(degG,B);
        if ((degF <= 0) || (degG <= 0) ) {
            return assign(S0, Degree(0), _domain.one);
        }

        //   if (degF >= degG) {
        assign(F, A);
        assign(G, B);
        //   }
        //   else {
        //     assign(F, B);
        //     assign(G, A);
        //   }

        Rep S1,R1,Q,TMP, TMP2;

        assign(S0,one);
        assign(S1,zero);

        while ( ! isZero(G) ) {
            //       write(std::cout << "F: ", F) << std::endl;
            //       write(std::cout << "G: ", G) << std::endl;

            divmod(Q,R1,F,G);
            //       write(std::cout << "Q: ", Q) << std::endl;
            //       write(std::cout << "R: ", R1) << std::endl;

            assign(F,G);
            assign(G,R1);
            //       write(std::cout << "Fn: ", F) << std::endl;
            //       write(std::cout << "Gn: ", G) << std::endl;
            mul(TMP,Q,S1); sub(TMP2,S0,TMP); assign(S0,S1); assign(S1,TMP2);
            //       write(std::cout << "S: ", S1) << std::endl;
        }
        //     std::cerr << "END invmodunit of " << typeid(*this).name() << std::endl;

        return S0;
    }


    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::lcm ( Rep& F, const Rep& A, const Rep& B) const
    {
        //     write(write(std::cerr << "A: ", A) << ", B: ", B) << std::endl;

        Rep G, S0, T0;
        Degree degA, degB, degG;
        degree(degA,A); degree(degB,B);
        if (degA < 0) return assign(F, Degree(0), _domain.zero);
        if (degB < 0) return assign(F, Degree(0), _domain.zero);
        if (degB == 0) return assign(F, A);
        if (degA == 0) return assign(F, B);

        if (degA >= degB) {
            assign(F, A);
            assign(G, B);
        }
        else {
            assign(F, B);
            assign(G, A);
        }

        Type_t r0, r1, tt;
        leadcoef(r0, F);
        leadcoef(r1, G);

        divin(F,r0);
        divin(G,r1);

        Rep S1,R1,T1,Q,TMP, TMP2;

        assign(S0, 0, _domain.inv(tt,r0) );
        assign(S1,zero);
        assign(T0,zero);
        assign(T1, 0, _domain.inv(tt,r1) );

        while ( ! isZero(G) ) {
            divmod(Q,R1,F,G);
            leadcoef(r1, R1); if (_domain.isZero(r1)) _domain.assign(r1,_domain.one);
            assign(F,G);
            div(G,R1,r1);
            mul(TMP,Q,S1); sub(TMP2,S0,TMP); assign(S0,S1); div(S1,TMP2,r1);
            mul(TMP,Q,T1); sub(TMP2,T0,TMP); assign(T0,T1); div(T1,TMP2,r1);


            //   std::cerr << "Degree Q : " << degree(degF,Q) << std::endl;
            //   std::cerr << "Degree R1 : " << degree(degF,R1) << std::endl;
            //   std::cerr << "Degree F : " << degree(degF,F) << std::endl;
            //   std::cerr << "Degree G : " << degree(degF,G) << std::endl;
            //   std::cerr << "Degree S0 : " << degree(degF,S0) << std::endl;
            //   std::cerr << "Degree T0 : " << degree(degF,T0) << std::endl;
            //   std::cerr << "Degree S1 : " << degree(degF,S1) << std::endl;
            //   std::cerr << "Degree T1 : " << degree(degF,T1) << std::endl;

            //   write( write( write( write( write( std::cerr, G ) << " = ((", T1) << ") * (", A) << ")) + ((", S1) << ") * (", B) << "))" << std::endl;

        }

        degree(degG, G);

        //     write(write(std::cerr << "S1: ", S1) << ", T1: ", T1) << std::endl;

        if ( degG <= 0) {
            // If normalisation is needed
            //       if (degA >= degB) {
            //           return mul(G, S1, A);
            //       } else {
            //           return mul(G, T1, B);
            //       }
            //       leadcoef(tt, G);
            //       return div(F,G,tt);
            if (degA >= degB) {
                return mul(F, S1, A);
            } else {
                return mul(F, T1, B);
            }
        } else {
            return mul(F, A, B);
        }
    }


} // GIVARO

#endif // __GIVARO_poly1_gcd_INL

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
