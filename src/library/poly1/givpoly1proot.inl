// =================================================================== //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Givaro / Athapascan-1
// Irreducible polynomial finder
// Primitive root finder
// Time-stamp: <10 Jul 07 15:31:18 Jean-Guillaume.Dumas@imag.fr>
// =================================================================== //
#ifndef __GIVARO_poly_primitive_root_INL
#define __GIVARO_poly_primitive_root_INL

// Reasonible cyclotomic polynomial size
#define CYCLO_DEGREE_BOUND 1000
#define CYCLO_TIMES_FACTOR 8

#include <stdlib.h>
#include <list>
#include <vector>
#include <givaro/givinteger.h>
#include <givaro/givintnumtheo.h>
#include <givaro/givdegree.h>
#include <givaro/givpoly1factor.h>
#include <givaro/givpoly1cyclo.inl>

namespace Givaro {

    ////////////////////////////////////////////
    // BINOM
    ////////////////////////////////////////////

    template<class Domain, class Tag, class RandomIterator >
    template<class Residue>
    inline bool Poly1FactorDom<Domain,Tag, RandomIterator>::find_irred_binomial (Element& R, Degree n, Residue MOD) const
    {
        for(Residue a=0; a<MOD; ++a) {
            _domain.assign(R[0],(Type_t)a);
            if (is_irreducible(R))
                return true;
        }
        return false;
    }

    template<class Domain, class Tag, class RandomIterator >
    // template<>
    inline bool Poly1FactorDom<Domain,Tag, RandomIterator>::find_irred_binomial (Element& R, Degree n, bool MOD) const
    {
        return false;
    }

    template<class Domain, class Tag, class RandomIterator >
    template<class Residue>
    inline bool Poly1FactorDom<Domain,Tag, RandomIterator>::find_irred_binomial (Element& R, Degree n, Residue MOD, Element IXE) const
    {
        for(Residu_t a=0; a<MOD; ++a) {
            _domain.assign(R[0],static_cast<Element_t>(a));
            if (is_irreducible(R) && (is_prim_root(IXE,R) ))
                return true;
        }

        return false;
    }

    template<class Domain, class Tag, class RandomIterator >
    // template<>
    inline bool Poly1FactorDom<Domain,Tag, RandomIterator>::find_irred_binomial (Element& R, Degree n, bool MOD, Element IXE) const
    {
        return false;
    }


    template<class Domain, class Tag, class RandomIterator >
    template<class Residue>
    inline bool Poly1FactorDom<Domain,Tag, RandomIterator>::find_irred_binomial2 (Element& R, Degree n, Residue MOD, Element IXE) const
    {
        for(Residu_t a=0; a<MOD; ++a) {
            _domain.assign(R[0],(Element_t)a);
            if (is_irreducible2(R) && (is_prim_root(IXE,R) ))
                return true;
        }

        return false;
    }

    template<class Domain, class Tag, class RandomIterator >
    // template<>
    inline bool Poly1FactorDom<Domain,Tag, RandomIterator>::find_irred_binomial2 (Element& R, Degree n, bool MOD, Element IXE) const
    {
        return false;
    }

    ////////////////////////////////////////////
    // TRINOM
    ////////////////////////////////////////////

    template<class Domain, class Tag, class RandomIterator >
    template<class Residue>
    inline bool Poly1FactorDom<Domain,Tag, RandomIterator>::find_irred_trinomial (Element& R, Degree n, Residue MOD) const
    {

        for(long d=1;d<=(n.value()/2);++d) {
            for(Residu_t b=0; b<MOD; ++b) {
                _domain.assign(R[(size_t)d],(Type_t)b);
                for(Residu_t a=1; a<MOD; ++a) {
                    _domain.assign(R[0],(Type_t)a);
                    if (is_irreducible(R))
                        return true;
                }
            }
            // _domain.assign(R[0],_domain.zero);
            // JGD 21.10.02
            _domain.assign(R[(size_t)d],_domain.zero);
        }
        return false ;
    }

    template<class Domain, class Tag, class RandomIterator >
    // template<>
    inline bool Poly1FactorDom<Domain,Tag, RandomIterator>::find_irred_trinomial (Element& R, Degree n, bool MOD) const
    {

        _domain.assign(R[0],_domain.one);
        for(long d=1;d<=(n.value()/2);++d) {
            _domain.assign(R[(size_t)d],_domain.one);
            if (is_irreducible(R))
                return true;
            _domain.assign(R[(size_t)d],_domain.zero);
        }
        return false ;
    }

    template<class Domain, class Tag, class RandomIterator >
    template<class Residue>
    inline bool Poly1FactorDom<Domain,Tag, RandomIterator>::find_irred_trinomial (Element& R, Degree n, Residue MOD, Element IXE) const
    {
        for(long d=2;d<=(n.value()/2);++d) {
            for(Residu_t b=0; b<MOD; ++b) {
                _domain.assign(R[(size_t)d],(Element_t)b);
                for(Residu_t a=1; a<MOD; ++a) {
                    _domain.assign(R[0],(Element_t)a);
                    if (is_irreducible(R) && (is_prim_root(IXE,R) ))
                        return true;
                }
            }
            // _domain.assign(R[0],_domain.zero);
            // JGD 21.10.02
            _domain.assign(R[(size_t)d],_domain.zero);
        }
        return false;
    }

    template<class Domain, class Tag, class RandomIterator >
    // template<>
    inline bool Poly1FactorDom<Domain,Tag, RandomIterator>::find_irred_trinomial (Element& R, Degree n, bool MOD, Element IXE) const
    {

        _domain.assign(R[0],_domain.one);
        for(long d=2;d<=(n.value()/2);++d) {
            _domain.assign(R[(size_t)d],_domain.one);
            if (is_irreducible(R) && (is_prim_root(IXE,R) ))
                return true;
            _domain.assign(R[(size_t)d],_domain.zero);
        }
        return false ;
    }

    template<class Domain, class Tag, class RandomIterator >
    template<class Residue>
    inline bool Poly1FactorDom<Domain,Tag, RandomIterator>::find_irred_trinomial2 (Element& R, Degree n, Residue MOD, Element IXE) const
    {
        for(long d=2;d<n.value();++d) {
            for(Residu_t b=0; b<MOD; ++b) {
                _domain.assign(R[(size_t)d],(Element_t)b);
                for(Residu_t a=1; a<MOD; ++a) {
                    _domain.assign(R[0],(Element_t)a);
                    if (is_irreducible2(R) && (is_prim_root(IXE,R) ))
                        return true;
                }
            }
            _domain.assign(R[0],_domain.zero);
        }
        return false;

    }

    template<class Domain, class Tag, class RandomIterator >
    // template<>
    inline bool Poly1FactorDom<Domain,Tag, RandomIterator>::find_irred_trinomial2 (Element& R, Degree n, bool MOD, Element IXE) const
    {

        _domain.assign(R[0],_domain.one);
        for(long d=2;d<=n.value();++d) {
            _domain.assign(R[(size_t)d],_domain.one);
            if (is_irreducible2(R) && (is_prim_root(IXE,R) ))
                return true;
            _domain.assign(R[(size_t)d],_domain.zero);
        }
        return false ;
    }

    ////////////////////////////////////////////
    // RANDOM
    ////////////////////////////////////////////
    template<class Domain, class Tag, class RandomIterator >
    template<class Residue>
    inline bool Poly1FactorDom<Domain,Tag, RandomIterator>::find_irred_randomial (Element& R, Degree n, Residue NUM) const
    {
#ifdef INFLOOPDEBUG
        int no_inf_loop =(int) n.value()/2+5 ;
#endif

        do {
            this->random( (RandomIterator&)_g, R, n); // must cast away const
            _domain.assign(R[(size_t)n.value()],_domain.one);
            for(Residu_t a=0; a<NUM; ++a) {
                _domain.random(_g, R[0]);
                if (is_irreducible(R))
                    return true;
            }
        }
#ifdef INFLOOPDEBUG
        while(--no_inf_loop);
#else
        while(1);
#endif

        return false;
    }

    template<class Domain, class Tag, class RandomIterator >
    // template<>
    inline bool Poly1FactorDom<Domain,Tag, RandomIterator>::find_irred_randomial (Element& R, Degree n, bool tag) const
    {
#ifdef INFLOOPDEBUG
        int no_inf_loop = n.value()/2+5 ;
#endif

        do {
            this->random( (RandomIterator&)_g, R, n); // must cast away const
            _domain.assign(R[(size_t)n.value()],_domain.one);
            _domain.assign(R[0],_domain.one);
            if (is_irreducible(R))
                return true;
        }
#ifdef INFLOOPDEBUG
        while(--no_inf_loop);
#else
        while(1);
#endif

        return false;
    }

    template<class Domain, class Tag, class RandomIterator >
    template<class Residue>
    inline bool Poly1FactorDom<Domain,Tag, RandomIterator>::find_irred_randomial (Element& R, Degree n, Residue NUM, Element IXE) const
    {
#ifdef INFLOOPDEBUG
        int no_inf_loop = (int)n.value()/2+5 ;
#endif
        do {
            this->random( (RandomIterator&)_g, R, n); // must cast away const
            _domain.assign(R[(size_t)n.value()],_domain.one);
            for(Residu_t a=0; a<NUM; ++a) {
                _domain.random(_g, R[0]);
                if (is_irreducible(R) && (is_prim_root(IXE,R) ))
                    return true;
            }
        }
#ifdef INFLOOPDEBUG
        while(--no_inf_loop);
#else
        while(1);
#endif

        return false;
    }

    template<class Domain, class Tag, class RandomIterator >
    // template<>
    inline bool Poly1FactorDom<Domain,Tag, RandomIterator>::find_irred_randomial (Element& R, Degree n, bool tag, Element IXE) const
    {
#ifdef INFLOOPDEBUG
        int no_inf_loop = n.value()/2+5 ;
#endif

        do {
            this->random( (RandomIterator&)_g, R, n); // must cast away const
            _domain.assign(R[(size_t)n.value()],_domain.one);
            _domain.assign(R[0],_domain.one);
            if (is_irreducible(R) && (is_prim_root(IXE,R) ))
                return true;
        }
#ifdef INFLOOPDEBUG
        while(--no_inf_loop);
#else
        while(1);
#endif

        return false;
    }

    template<class Domain, class Tag, class RandomIterator >
    template<class Residue>
    inline bool Poly1FactorDom<Domain,Tag, RandomIterator>::find_irred_randomial2 (Element& R, Degree n, Residue NUM, Element IXE) const
    {
#ifdef INFLOOPDEBUG
        int no_inf_loop = n.value()/2+5 ;
#endif
        do {
            this->random( (RandomIterator&)_g, R, n); // must cast away const
            _domain.assign(R[(size_t)n.value()],_domain.one);
            for(Residu_t a=0; a<NUM; ++a) {
                _domain.random(_g, R[0]);
                if (is_irreducible2(R) && (is_prim_root(IXE,R) ))
                    return true;
            }
        }
#ifdef INFLOOPDEBUG
        while(--no_inf_loop);
#else
        while(1);
#endif

        return false;
    }

    template<class Domain, class Tag, class RandomIterator >
    // template<>
    inline bool Poly1FactorDom<Domain,Tag, RandomIterator>::find_irred_randomial2 (Element& R, Degree n, bool tag, Element IXE) const
    {
#ifdef INFLOOPDEBUG
        int no_inf_loop = n.value()/2+5 ;
#endif
        do {
            this->random( (RandomIterator&)_g, R, n); // must cast away const
            _domain.assign(R[(size_t)n.value()],_domain.one);
            _domain.assign(R[0],_domain.one);
            if (is_irreducible2(R) && (is_prim_root(IXE,R) ))
                return true;
        }
#ifdef INFLOOPDEBUG
        while(--no_inf_loop);
#else
        while(1);
#endif
        return false;

    }

    ////////////////////////////////////////////


    // ---------------------------------------------------------------
    // Monic irreducible polynomial of degree n over Z/pZ
    // having 2, 3 nonzero terms or dividing a cyclotomic polynomial
    // of degree < CYCLO_DEGREE_BOUND or a random one.
    // ---------------------------------------------------------------

    template<class Domain, class Tag, class RandomIterator >
    inline typename Poly1FactorDom<Domain,Tag, RandomIterator>::Element& Poly1FactorDom<Domain,Tag, RandomIterator>::creux_random_irreducible (Element& R, Degree n) const
    {
        this->init(R, n);

        Residu_t MOD = _domain.residu();


        // Search for an irreducible BINOMIAL : X^n + a
        // WARNING : Here we may have X^n + x,
        // where a = representation of x, and sometimes a != x.
        if (find_irred_binomial(R,n,MOD))
            return R ;

        // Search for an irreducible TRINOMIAL : X^n + b*X^i + a
        // Precondition : n >= 2
        assert(n.value()>=2);
        // WARNING : same warning as for the binomial.
        // JGD 21.10.02
        // for(Residu_t d=2;d<n.value();++d)

        if (find_irred_trinomial(R,n,MOD))
            return R ;

        // Search for a monic irreducible Polynomial
        // with random Elements

        if (find_irred_randomial(R,n,MOD))
            return R ;
        else
            throw "could not find a random polynomial" ;

    }

    template<class Domain, class Tag, class RandomIterator >
    inline typename Poly1FactorDom<Domain,Tag, RandomIterator>::Element& Poly1FactorDom<Domain,Tag, RandomIterator>::random_irreducible (Element& R, Degree n) const
    {
        // Search for a monic irreducible Polynomial
        // with random Elements
        this->init(R, n);
        Residu_t MOD = _domain.residu();

        if (find_irred_randomial(R,n,MOD))
            return R ;
        else
            throw "could not find a random polynomial" ;

    }


    // ---------------------------------------------------------------
    // Monic irreducible polynomial of degree n over Z/pZ
    // having 2, 3 nonzero terms or or a random one,
    // with X as a primitive root.
    // ---------------------------------------------------------------
    template<class Domain, class Tag, class RandomIterator >
    inline typename Poly1FactorDom<Domain,Tag, RandomIterator>::Element& Poly1FactorDom<Domain,Tag, RandomIterator>::ixe_irreducible (Element& R, Degree n) const
    {
        this->init(R, n);
        Element IXE;
        this->init(IXE,Degree(1));
        Residu_t MOD = _domain.residu();

        // Search for an irreducible BINOMIAL : X^n + a
        // WARNING : Here we may have X^n + x,
        // where a = representation of x, and sometimes a != x.

        if (find_irred_binomial(R,n,MOD,IXE))
            return R ;

        // Search for an irreducible TRINOMIAL : X^n + b*X^i + a
        // Precondition : n >= 2
        assert(n.value()>=2);
        // WARNING : same warning as for the binomial.
        // // JGD 21.10.02
        // for(unsigned long d=2;d<n.value();++d)

        if (find_irred_trinomial(R,n,MOD,IXE))
            return R ;


        // Search for a monic irreducible Polynomial
        // with random Elements
        if (find_irred_randomial(R,n,MOD,IXE))
            return R ;
        else
            throw "could not find a random polynomial" ;


    }

    template<class Domain, class Tag, class RandomIterator >
    inline typename Poly1FactorDom<Domain,Tag, RandomIterator>::Element& Poly1FactorDom<Domain,Tag, RandomIterator>::ixe_irreducible2 (Element& R, Degree n) const
    {
        this->init(R, n);
        Element IXE;
        this->init(IXE,Degree(1));
        Residu_t MOD = _domain.residu();

        // Search for an irreducible BINOMIAL : X^n + a
        // WARNING : Here we may have X^n + x,
        // where a = representation of x, and sometimes a != x.

        if (find_irred_binomial2(R,n,MOD,IXE))
            return R ;

        // Search for an irreducible TRINOMIAL : X^n + b*X^i + a
        // Precondition : n >= 2
        assert(n.value()>=2);
        // WARNING : same warning as for the binomial.

        if (find_irred_trinomial2(R,n,MOD,IXE))
            return R ;



        // Search for a monic irreducible Polynomial
        // with random Elements

        if (find_irred_randomial2(R,n,MOD,IXE))
            return R ;
        else
            throw "could not find a random polynomial" ;


    }
    // ---------------------------------------------------------------
    // Irreducibility tests
    // ---------------------------------------------------------------

    template<class Domain, class Tag, class RandomIterator>
    inline bool Poly1FactorDom<Domain,Tag, RandomIterator>::is_irreducible2( const Rep& P
                                                                             , Residu_t MOD ) const
    {
        // Square free ?
        Rep W,D; this->gcd(W,diff(D,P),P);
        Degree d, dP;
        if (this->degree(d,W) > 0) return 0;
        IntFactorDom<> FD;

        int64_t n = this->degree(dP,P).value();
        IntFactorDom<>::Rep qn;

        FD.pow( qn, IntFactorDom<>::Rep(MOD), n);
        Rep Unit, G1; this->init(Unit, Degree(1));
        this->powmod(G1, Unit, qn, P);
        if (this->degree(d, sub(D,G1,Unit)) >= 0) return 0;

        std::vector<IntFactorDom<>::Rep> Lp; std::vector<uint64_t> Le;
        FD.set(Lp, Le, n );
        for( std::vector<IntFactorDom<>::Rep>::const_iterator p = Lp.begin(); p != Lp.end(); ++p) {
            int64_t ttmp;
            FD.pow( qn, IntFactorDom<>::Rep(MOD), n/FD.convert(ttmp,*p) );
            this->powmod(G1, Unit, qn, P);
            if (this->degree(d, sub(D,G1,Unit)) < 0) return 0;
        }

        return 1;
    }




    // ---------------------------------------------------------------
    // Primitive Root over Z/pZ / F
    // returns 1 if P is a generator.
    // ---------------------------------------------------------------

    template<class Domain, class Tag, class RandomIterator>
    bool Poly1FactorDom<Domain,Tag, RandomIterator>::is_prim_root( const Rep& P, const Rep& F)  const
    {
        bool isproot = 0;
        Rep A, G;
        this->mod(A,P,F);
        Degree d;
        if ( this->degree(d, this->gcd(G,A,F)) == 0) {
            Residu_t MOD = _domain.residu();
            IntFactorDom<> FD;
            IntFactorDom<>::Element IMOD( MOD ), q, qp;
            this->degree(d,F);
            //         FD.pow(q ,IMOD, d.value());
            //         FD.sub(qp, q, FD.one);
            FD.subin( FD.pow(qp ,IMOD, d.value()) , FD.one);
            std::list< IntFactorDom<>::Element > L;
            FD.set(L, qp);
            L.sort();
            std::list< IntFactorDom<>::Element >::iterator li = L.begin();
            isproot = 1;
            for(;(li != L.end()) && isproot; ++li) {
                isproot = ( ! this->isOne(this->powmod(G, A, FD.div(q, qp , *li), F) ) );
            }
        }
        return isproot;
    }

    template<class Domain, class Tag, class RandomIterator>
    inline typename IntegerDom::Element Poly1FactorDom<Domain,Tag, RandomIterator>::order( const Rep& P, const Rep& F)  const
    {
        bool isproot = 0;
        Rep A, G; mod(A,P,F);
        Degree d;
        if ( this->degree(d, this->gcd(G,A,F)) == 0) {
            Residu_t MOD = _domain.residu();
            IntFactorDom<> FD;
            IntFactorDom<>::Element IMOD( MOD ), g, gg, tt, qp;
            this->degree(d,F);
            //         FD.pow(q ,IMOD, d.value());
            //         FD.sub(qp, q, FD.one);
            FD.subin( FD.pow(qp ,IMOD, d.value()) , FD.one);
            std::list< IntFactorDom<>::Element > L;
            FD.set(L, qp);
            L.sort();
            std::list< IntFactorDom<>::Element >::iterator li = L.begin();
            isproot = 1;
            for(;(li != L.end()) && isproot; ++li)
                isproot = ( ! this->isOne(this->powmod(G, A, FD.div(g, qp , *li), F) ) );

            if (isproot)
                return qp;
            else {
                for(--li;li!=L.end();++li)
                    while ( FD.isZero(FD.mod(tt,g,*li)) && (this->isOne(this->powmod(G, A, FD.div(gg,g,*li), F))))
                        g.copy(gg);
                return g;
            }
        }
        IntegerDom ID;
        return ID.zero;
    }

    template<class Domain, class Tag, class RandomIterator >
    inline typename Poly1FactorDom<Domain,Tag, RandomIterator>::Rep& Poly1FactorDom<Domain,Tag, RandomIterator>::give_prim_root(Rep& R, const Rep& F)  const
    {
        Degree n; this->degree(n,F);
        Residu_t MOD = _domain.residu();
        //    this->write(std::cout << "Give Pr: ", F) << std::endl;


        // Search for a primitive BINOMIAL : X^i + a
        for(Degree di=1;di<n;++di) {
            this->init(R, di);
            //         for(Residu_t a=MOD; a--; )
            for(Residu_t a=0; a<MOD;++a ) {
                _domain.assign(R[0],(Element_t)a);
                if (is_prim_root(R,F))
                    return R;
            }
        }
        // Search for a primitive TRINOMIAL : X^i + b*X^j + a
        for(Degree di=2;di<n;++di) {
            this->init(R, di);
            for(Degree dj=1;dj<di;++dj)
                //             for(Residu_t b=MOD; b--;)
                for(Residu_t b=0; b<MOD;++b) {
                    _domain.assign(R[(size_t)dj.value()],(Element_t)b);
                    //                 for(Residu_t a=MOD; a--;)
                    for(Residu_t a=0; a<MOD;++a ) {
                        _domain.assign(R[0],(Element_t)a);
                        if (is_prim_root(R,F))
                            return R;
                    }
                }
        }

        // Search for a primitive Polynomial
        // with random Elements
        do {
            this->random( (RandomIterator&)_g, R, n); // must cast away const
            _domain.assign(R[(size_t)n.value()],_domain.one);
            for(Residu_t a=0; a<MOD; ++a) {
                _domain.assign(R[0],(Element_t)a);
                if (is_prim_root(R,F))
                    return R;
            }
        } while(1);
    }


    template<class Domain, class Tag, class RandomIterator >
    inline typename Poly1FactorDom<Domain,Tag, RandomIterator>::Rep& Poly1FactorDom<Domain,Tag, RandomIterator>::give_random_prim_root(Rep& R, const Rep& F)  const
    {
        Degree n; this->degree(n,F);
        Residu_t MOD = _domain.residu();

        // Search for a primitive Polynomial
        // with random Elements
        do {
            this->random( (RandomIterator&)_g, R, n); // must cast away const
            _domain.assign(R[(size_t)n.value()],_domain.one);
            for(Residu_t a=0; a<MOD; ++a) {
                _domain.assign(R[0],(Element_t)a);
                if (is_prim_root(R,F))
                    return R;
            }
        } while(1);
    }


    template<class Domain, class Tag, class RandomIterator >
    inline typename Poly1FactorDom<Domain,Tag, RandomIterator>::Rep& Poly1FactorDom<Domain,Tag, RandomIterator>::random_prim_root(Rep& P, Rep& R, Degree n)  const
    {
        // P is irreducible
        // R is a primitive root. i.e R generates (Z_p)/P.
        // returns R
        return give_prim_root(R, random_irreducible(P,n));
    }

} // Givaro

#endif // __GIVARO_poly_primitive_root_INL

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
