// ================================================================= //
// Givaro / Athapascan-1
// Irreducibily test
// Factorisations de  Polynomes dans Fp[X] :
//      Distinct Degree
//      Cantor-Zassenhaus
//      Berlekamp : in LinBox
// Time-stamp: <06 Jan 05 18:01:43 Jean-Guillaume.Dumas@imag.fr> 
// ================================================================= //
#ifndef _GIV_POLY1_FACTO_INL_
#define _GIV_POLY1_FACTO_INL_
#include <givaro/givpower.h>
#include <givaro/givtimer.h>

// ---------------------------------------------------------------
// Splits a polynomial into prime factors of same degree
// ---------------------------------------------------------------

template<class Domain, class Tag, class RandIter>
template< template<class> class Container > 
inline void Poly1FactorDom<Domain,Tag, RandIter>::SplitFactor( 
    Container< Rep > & L
    , const Rep& G
    , Degree d
    , typename Domain::Residu_t MOD) const {
    Degree dG;degree(dG,G);
    if (dG == d)
        L.push_back(G);
    else {
        int splitted = 0;
        while (! splitted) {
            Rep tmp, G1;
            gcd(G1, G, random(_g, tmp, d));
            Degree dG1; degree(dG1,G1);
            if ( dG1 != dG) {
                if (dG1 > 0 ) {
                    splitted = 1;
                    SplitFactor ( L, G1, d, MOD) ;
                }

                typename Domain::Residu_t pp = (power(MOD, d.value()) - 1)/2;
                Rep tp, tp2, G2;
                gcd(G2,G, sub(tp2, powmod(tp, tmp, pp, G) , one) );
                Degree dG2; degree(dG2,G2);
                if ( dG2 != dG) {
                    if ( dG2 > 0 ) {
                        splitted = 1 ;
                        SplitFactor ( L, G2, d, MOD) ;
                    }
                    Rep G3; gcd(G3, G, add(tp2,tp,one) );
                    Degree dG3; degree(dG3,G3);
                    if (( dG3 != dG) && (dG3 > 0 )) {
                        splitted = 1 ;
                        SplitFactor ( L, G3, d, MOD) ;
                    }
                }
            }
        }
    }
}


template<class Domain, class Tag, class RandIter>
inline typename Poly1FactorDom<Domain,Tag, RandIter>::Rep& Poly1FactorDom<Domain,Tag, RandIter>::SplitFactor(
    Rep& G1
    , const Rep& G
    , Degree d
    , typename Domain::Residu_t MOD)  const  {
    Degree dG;degree(dG,G);
    if (dG == d)
        return G1.copy(G) ;
    else {
        while (1) {
            Rep tmp;
            gcd(G1, G, random(_g, tmp, d));
            Degree dG1; degree(dG1,G1);
            if ( dG1 != dG) {
                if (dG1 > 0 ) {
                    return G1;
                }
                typename Domain::Residu_t pp = (power(MOD, d.value()) - 1)/2;
                Rep tp, tp2, G2;
                gcd(G2,G, sub(tp2, powmod(tp, tmp, pp, G) , one) );
                Degree dG2; degree(dG2,G2);
                if ( dG2 != dG) {
                   if ( dG2 > 0 ) {
                        return G1.copy(G2);
                    }
                    Rep G3; gcd(G3, G, add(tp2,tp,one) );
                    Degree dG3; degree(dG3,G3);
                    if (( dG3 != dG) && (dG3 > 0 )) {
                        return G1.copy(G3);
                    }
                }
            }
        }
    }
}


// ---------------------------------------------------------------
// Splits a polynomial into divisors of homogenous prime factors
// ---------------------------------------------------------------

template<class Domain, class Tag, class RandIter>
template< template<class> class Container > 
inline void Poly1FactorDom<Domain,Tag, RandIter>::DistinctDegreeFactor(
    Container< Rep > & L
    , const Rep& f
    , typename Domain::Residu_t MOD)  const  {
    // srand48(BaseTimer::seed());
// write(std::cerr << "DD in: ", f) << std::endl;
    Rep W, D, P = f;
    Degree dP;
    Rep Unit, G1; init(Unit, Degree(1), _domain.one);
    W.copy(Unit);
    degree(dP,P); Degree dPo = (dP/2);
    for(Degree dp = 1; dp <= dPo; ++dp) {
        powmod(W, D.copy(W), MOD, P);
        gcd (G1,sub(D,W,Unit), P) ;
        Degree dG1; degree(dG1,G1);
        if ( dG1 > 0 ) {
            SplitFactor (L, G1, dp, MOD);
            divin(P,G1);
        }
    }
    degree(dP,P);    
    if (dP > 0)
        L.push_back(P);
// write(std::cerr << "DD: ", P) << std::endl;
}

// ---------------------------------------------------------------
// Cantor-Zassenhaus Polynomial factorization over Z/pZ
// ---------------------------------------------------------------

template<class Domain, class Tag, class RandIter>
template< template<class> class Container > 
inline void Poly1FactorDom<Domain,Tag, RandIter>::CZfactor( 
    Container< Rep > & Lf
    , Container< unsigned long > & Le
    , const Rep & P
    , typename Domain::Residu_t MOD ) const  {
// write(std::cerr << "CZ in: ", P) << std::endl;
    Degree dp; degree(dp,P);
    size_t nb=dp.value()+1; 
    Rep * g = new Rep[nb];
    sqrfree(nb,g,P);
// std::cerr << "CZ sqrfree: " << nb << std::endl;
    for(size_t i = 0; i<nb;++i) {
        size_t this_multiplicity = Lf.size();
        DistinctDegreeFactor(Lf, g[i], MOD) ;
        Le.resize(Lf.size());
        for( ; this_multiplicity < Lf.size(); ++this_multiplicity)
            Le[this_multiplicity] = i+1;
// std::cerr << "multiplicities";
// for (typename Container< unsigned long >::const_iterator e=Le.begin(); e!=Le.end(); ++e)
// std::cerr << " " << *e;
// std::cerr << std::endl;
        
    }
    ::delete [] g;
}



// ---------------------------------------------------------------
// Irreducibility tests
// ---------------------------------------------------------------

template<class Domain, class Tag, class RandIter>
inline bool Poly1FactorDom<Domain,Tag, RandIter>::is_irreducible(
    const Rep& P
    , typename Domain::Residu_t MOD ) const  {
        // Square free ?
    Rep W,D; gcd(W,diff(D,P),P);
    Degree d, dP;
    if (degree(d,W) > 0) return 0;
        // Distinct degree free ?
    Rep Unit, G1; init(Unit, Degree(1), _domain.one);
    W.copy(Unit);
    degree(dP,P); Degree dPo = (dP/2);
    for(Degree dp = 1; dp <= dPo; ++dp) {
        powmod(W, D.copy(W), MOD, P);
        gcd (G1, sub(D,W,Unit), P) ;
        if ( degree(d,G1) > 0 ) return 0;
    }
    return 1;
}


// ---------------------------------------------------------------
// Gives one non-trivial factor of P if P is reducible
// returns P otherwise
// ---------------------------------------------------------------

template<class Domain, class Tag, class RandIter>
inline typename Poly1FactorDom<Domain,Tag, RandIter>::Rep& Poly1FactorDom<Domain,Tag, RandIter>::factor(
    Rep& W
    , const Rep& P
    , typename Domain::Residu_t MOD)  const {
// write(cerr << "In factor P:", P) << endl;
        // Square free ?
    Rep D; gcd(W,diff(D,P),P);
    Degree d, dP;
// write(cerr << "In factor P':", D) << "(deg: " << degree(d,D) << ")" << endl;
// write(cerr << "In factor P^P':", W) << "(deg: " << degree(d,W) << ")" << endl;
 
    if (degree(d,W) > 0) return W;
        // Distinct degree free ?
    Rep Unit, G1; init(Unit, Degree(1), _domain.one);
// write(cerr << "In factor U:", Unit) << endl;
    W.copy(Unit);
    degree(dP,P); Degree dPo = (dP/2);
    for(Degree dp = 1; dp <= dPo; ++dp) {
// write(cerr << "In factor W:(deg: " << degree(d,W) << "):", W) << endl;
        powmod(W, D.copy(W), MOD, P);
        gcd (G1, sub(D,W,Unit), P) ;
        Degree dG1; degree(dG1,G1);
        if ( dG1 > 0 )
            if (dG1 < dP)
                return W.copy(G1);
            else
                return SplitFactor(W,G1,dp,MOD);
    }
    return W.copy(P);
}


#endif 
