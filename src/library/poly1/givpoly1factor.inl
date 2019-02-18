// ================================================================= //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Givaro / Athapascan-1
// Irreducibily test
// Factorisations de  Polynomes dans Fp[X] :
//      Distinct Degree
//      Cantor-Zassenhaus
//      Berlekamp : in LinBox
// Time-stamp: <26 Feb 08 13:37:26 Jean-Guillaume.Dumas@imag.fr>
// ================================================================= //
#ifndef __GIVARO_poly1_facto_INL
#define __GIVARO_poly1_facto_INL
#include <givaro/givpower.h>
#include <givaro/givcaster.h>

//!@todo use NTL if available ?


namespace Givaro {

    // ---------------------------------------------------------------
    // Splits a polynomial into prime factors of same degree
    // ---------------------------------------------------------------

    template<class Domain, class Tag, class RandomIterator>
    template< template<class, class> class Container, template <class> class Alloc >
    inline void Poly1FactorDom<Domain,Tag, RandomIterator>::SplitFactor(
                                                                        Container< Rep, Alloc<Rep> > & L
                                                                        , const Rep& G
                                                                        , Degree d
                                                                        , Residu_t MOD) const
    {
        Degree dG; this->degree(dG,G);
        if (dG == d)
            L.push_back(G);
        else {
            int splitted = 0;
            while (! splitted) {
                Rep G1, G2;
                this->gcd(G1, G, this->random(_g, G2, dG-1) );
                Degree dG1; this->degree(dG1,G1);
                // write(std::cerr << "SF rd: ", G2) << std::endl;
                // write(std::cerr << "SF G1: ", G1) << std::endl;
                if ( dG1 != dG) {
                    if (dG1 > 0 ) {
                        splitted = 1;
                        SplitFactor ( L, G1, d, MOD) ;
                        this->div(G2, G, G1);
                        SplitFactor ( L, G2, d, MOD) ;
                        return ;
                    }

                    Integer iMOD; Caster(iMOD, MOD);
                    Integer pp = (power(iMOD, d.value()) - 1)/2;
                    // std::cerr << "pp: " << pp << std::endl;
                    Rep tp;
                    this->gcd(G1, G,
                              this->subin( this->powmod(tp, G2, pp, G),
                                           _domain.one)
                             );

                    this->degree(dG1,G1);
                    // write(std::cerr << "SF t2: ", tp2) << std::endl;
                    // write(std::cerr << "SF G1: ", G1) << std::endl;
                    if ( ( dG1 != dG) && (dG1 > 0 ) ) {
                        splitted = 1 ;
                        SplitFactor ( L, G1, d, MOD) ;
                        this->div( G2, G, G1);
                        SplitFactor ( L, G2, d, MOD) ;
                    }
                }
            }
        }
    }


    template<class Domain, class Tag, class RandomIterator>
    inline typename Poly1FactorDom<Domain,Tag, RandomIterator>::Rep& Poly1FactorDom<Domain,Tag, RandomIterator>::SplitFactor(
                                                                                                                             Rep& G1
                                                                                                                             , const Rep& G
                                                                                                                             , Degree d
                                                                                                                             , Residu_t MOD)  const  {
        Degree dG;this->degree(dG,G);
        if (dG == d)
            return G1 = G ;
        else {
            while (1) {
                Rep tmp;
                this->gcd(G1, G, this->random(_g, tmp, d));
                // write(std::cerr << "SF rd: ", tmp) << std::endl;
                // write(std::cerr << "SF G1: ", G1) << std::endl;
                Degree dG1; this->degree(dG1,G1);
                if ( dG1 != dG) {
                    if (dG1 > 0 ) {
                        return G1;
                    }
                    Integer iMOD; Caster(iMOD, MOD);
                    Integer pp = (power(iMOD, d.value()) - 1)/2;
                    Rep tp, tp2, G2;
                    this->gcd(G2,G, this->sub(tp2, this->powmod(tp, tmp, pp, G) , _domain.one) );
                    Degree dG2; this->degree(dG2,G2);
                    // write(std::cerr << "SF t2: ", tp2) << std::endl;
                    // write(std::cerr << "SF G2: ", G2) << std::endl;
                    if ( dG2 != dG) {
                        if ( dG2 > 0 ) {
                            return G1 = G2;
                        }
                        // UNNECESSARY : ANYTHING FOUND BY G3 WOULD HAVE THE COFACTOR IN G2
                        Rep G3; this->gcd(G3, G, this->add(tp2,tp,_domain.one) );
                        Degree dG3; this->degree(dG3,G3);
                        // write(std::cerr << "SF t3: ", tp2) << std::endl;
                        // write(std::cerr << "SF G3: ", G3) << std::endl;
                        if (( dG3 != dG) && (dG3 > 0 )) {
                            return G1 = G3;
                        }
                    }
                }
            }
        }
    }


    // ---------------------------------------------------------------
    // Splits a polynomial into divisors of homogenous prime factors
    // ---------------------------------------------------------------

    template<class Domain, class Tag, class RandomIterator>
    template< template<class, class> class Container, template <class> class Alloc >
    inline void Poly1FactorDom<Domain,Tag, RandomIterator>::DistinctDegreeFactor(
                                                                                 Container< Rep, Alloc<Rep> > & L
                                                                                 , const Rep& f
                                                                                 , Residu_t MOD)  const  {

        // write(std::cerr << "DD in: ", f) << std::endl;
        Rep W, D, P = f;
        Degree dP;
        Rep Unit, G1;
        this->init(Unit, Degree(1));
        W = Unit;
        this->degree(dP,P); Degree dPo = (dP/2);
        for(Degree dp = 1; dp <= dPo; ++dp) {
            // std::cerr << "DD degree: " << dp << std::endl;
            D = W;
            this->powmod(W, D, MOD, P);
            this->gcd (G1,this->sub(D,W,Unit), P) ;
            Degree dG1; this->degree(dG1,G1);
            // write(std::cerr << "DD found: ", G1) << ", of degree " << dG1 << std::endl;
            if ( dG1 > 0 ) {
                SplitFactor (L, G1, dp, MOD);
                this->divin(P,G1);
            }
        }
        this->degree(dP,P);
        if (dP > 0)
            L.push_back(P);
        // write(std::cerr << "DD: ", P) << std::endl;
    }

    // ---------------------------------------------------------------
    // Cantor-Zassenhaus Polynomial factorization over Z/pZ
    // ---------------------------------------------------------------

    template<class Domain, class Tag, class RandomIterator>
    template< template<class, class> class Container, template <class> class Alloc>
    inline void
    Poly1FactorDom<Domain,Tag, RandomIterator>::CZfactor( Container< Rep, Alloc<Rep> > & Lf,
                                                          Container< uint64_t, Alloc<uint64_t> > & Le,
                                                          const Rep& P,
                                                          Residu_t MOD)  const
    {
        // write(std::cerr << "CZ in: ", P) << std::endl;
        Degree dp;
        this->degree(dp,P);
        size_t nb=(size_t)dp.value()+1;
        Rep * g = new Rep[nb];
        this->sqrfree(nb,g,P);
        // std::cerr << "CZ sqrfree: " << nb << std::endl;
        for(size_t i = 0; i<nb;++i) {
            size_t this_multiplicity = Lf.size();
            DistinctDegreeFactor(Lf, g[i], MOD) ;
            Le.resize(Lf.size());
            for( ; this_multiplicity < Lf.size(); ++this_multiplicity)
                Le[this_multiplicity] = i+1;
            // std::cerr << "multiplicities";
            // for (typename Container< uint64_t, Alloc<uint64_t> >::const_iterator e=Le.begin(); e!=Le.end(); ++e)
            // std::cerr << " " << *e;
            // std::cerr << std::endl;

        }
        ::delete [] g;
    }



    // ---------------------------------------------------------------
    // Irreducibility tests
    // ---------------------------------------------------------------

    template<class Domain, class Tag, class RandomIterator>
    inline bool Poly1FactorDom<Domain,Tag, RandomIterator>::is_irreducible( const Rep& P
                                                                            , Residu_t MOD ) const
    {
        Rep W,D;
        this->gcd(W,this->diff(D,P),P);
        Degree d, dP;
        if (this->degree(d,W) > 0) return 0;
        // Distinct degree free ?
        Rep Unit, G1;
        this->init(Unit, Degree(1));
        W = Unit;
        this->degree(dP,P); Degree dPo = (dP/2);
        for(Degree dp = 1; dp <= dPo; ++dp) {
            D = W;
            this->powmod(W, D, MOD, P);
            this->gcd (G1, this->sub(D,W,Unit), P) ;
            if ( this->degree(d,G1) > 0 ) return 0;
        }
        return 1;
    }


    // ---------------------------------------------------------------
    // Gives one non-trivial factor of P if P is reducible
    // returns P otherwise
    // ---------------------------------------------------------------

    template<class Domain, class Tag, class RandomIterator>
    inline typename Poly1FactorDom<Domain,Tag, RandomIterator>::Rep& Poly1FactorDom<Domain,Tag, RandomIterator>::factor(
                                                                                                                        Rep& W
                                                                                                                        , const Rep& P
                                                                                                                        , Residu_t MOD)  const
    {
        // write(cerr << "In factor P:", P) << endl;
        // Square free ?
        Rep D; this->gcd(W,diff(D,P),P);
        Degree d, dP;
        // write(cerr << "In factor P':", D) << "(deg: " << degree(d,D) << ")" << endl;
        // write(cerr << "In factor P^P':", W) << "(deg: " << degree(d,W) << ")" << endl;

        if (this->degree(d,W) > 0) return W;
        // Distinct degree free ?
        Rep Unit, G1; init(Unit, Degree(1));
        // write(cerr << "In factor U:", Unit) << endl;
        W.copy(Unit);
        this->degree(dP,P); Degree dPo = (dP/2);
        for(Degree dp = 1; dp <= dPo; ++dp) {
            // write(cerr << "In factor W:(deg: " << degree(d,W) << "):", W) << endl;
            this->powmod(W, D.copy(W), MOD, P);
            this->gcd (G1, sub(D,W,Unit), P) ;
            Degree dG1; this->degree(dG1,G1);
            if ( dG1 > 0 ) {
                if (dG1 < dP)
                    return W.copy(G1);
                else
                    return SplitFactor(W,G1,dp,MOD);
            }
        }
        return W.copy(P);
    }

} // Givaro


#endif // __GIVARO_poly1_facto_INL

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
