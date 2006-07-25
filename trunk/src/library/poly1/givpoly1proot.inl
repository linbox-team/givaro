// =================================================================== //
// Givaro / Athapascan-1
// Irreducible polynomial finder
// Primitive root finder
// Time-stamp: <21 Jul 06 10:15:04 Jean-Guillaume.Dumas@imag.fr> 
// =================================================================== //
#ifndef _GIVARO_POLY_PRIMITIVE_ROOT_
#define _GIVARO_POLY_PRIMITIVE_ROOT_

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

// ---------------------------------------------------------------
// Monic irreducible polynomial of degree n over Z/pZ
// having 2, 3 nonzero terms or dividing a cyclotomic polynomial
// of degree < CYCLO_DEGREE_BOUND or a random one.
// ---------------------------------------------------------------

template<class Domain, class Tag, class RandIter >
inline typename Poly1FactorDom<Domain,Tag, RandIter>::Element& Poly1FactorDom<Domain,Tag, RandIter>::creux_random_irreducible (Element& R, Degree n) const {
    init(R, n, _domain.one);
    Residu_t MOD = _domain.residu();

        // Search for an irreducible BINOMIAL : X^n + a
        // WARNING : Here we may have X^n + x, 
        // where a = representation of x, and sometimes a != x.
    for(Residu_t a=0; a<MOD; ++a) {
        _domain.assign(R[0],a);
        if (is_irreducible(R))
            return R;
    }
        // Search for an irreducible TRINOMIAL : X^n + b*X^i + a
        // Precondition : n >= 2
        // WARNING : same warning as for the binomial.
	// JGD 21.10.02
    // for(Residu_t d=2;d<n.value();++d) {
    for(long d=2;d<=(n.value()/2);++d) {
        for(Residu_t b=0; b<MOD; ++b) {
            _domain.assign(R[d],b);
            for(Residu_t a=1; a<MOD; ++a) {
                _domain.assign(R[0],a);
                if (is_irreducible(R))
                    return R;
            }
        }
	// _domain.assign(R[0],_domain.zero);
	// JGD 21.10.02
        _domain.assign(R[d],_domain.zero);
    }

        // Search for a monic irreducible Polynomial
        // with random Elements
    do {
        this->random( (RandIter&)_g, R, n); // must cast away const 
        _domain.assign(R[n.value()],_domain.one);
        for(Residu_t a=0; a<MOD; ++a) {
            _domain.assign(R[0],a);
            if (is_irreducible(R))
                return R;
        }
    } while(1);
}

template<class Domain, class Tag, class RandIter >
inline typename Poly1FactorDom<Domain,Tag, RandIter>::Element& Poly1FactorDom<Domain,Tag, RandIter>::random_irreducible (Element& R, Degree n) const {
        // Search for a monic irreducible Polynomial
        // with random Elements
    init(R, n, _domain.one);
    Residu_t MOD = _domain.residu();

    do {
        this->random( (RandIter&)_g, R, n); // must cast away const 
        _domain.assign(R[n.value()],_domain.one);
        for(Residu_t a=0; a<MOD; ++a) {
            _domain.assign(R[0],a);
            if (is_irreducible(R))
                return R;
        }
    } while(1);
}


// ---------------------------------------------------------------
// Monic irreducible polynomial of degree n over Z/pZ
// having 2, 3 nonzero terms or or a random one,
// with X as a primitive root.
// ---------------------------------------------------------------
template<class Domain, class Tag, class RandIter >
inline typename Poly1FactorDom<Domain,Tag, RandIter>::Element& Poly1FactorDom<Domain,Tag, RandIter>::ixe_irreducible (Element& R, Degree n) const {
    init(R, n, _domain.one);
    Element IXE;
    init(IXE,Degree(1),_domain.one);
    Residu_t MOD = _domain.residu();

        // Search for an irreducible BINOMIAL : X^n + a
        // WARNING : Here we may have X^n + x, 
        // where a = representation of x, and sometimes a != x.
    for(Residu_t a=0; a<MOD; ++a) {
        _domain.assign(R[0],a);
        if (is_irreducible(R) && (is_prim_root(IXE,R) ))
            return R;
    }
        // Search for an irreducible TRINOMIAL : X^n + b*X^i + a
        // Precondition : n >= 2
        // WARNING : same warning as for the binomial.
    // // JGD 21.10.02
    // for(unsigned long d=2;d<n.value();++d) {
    for(long d=2;d<=(n.value()/2);++d) {
        for(Residu_t b=0; b<MOD; ++b) {
            _domain.assign(R[d],b);
            for(Residu_t a=1; a<MOD; ++a) {
                _domain.assign(R[0],a);
                if (is_irreducible(R) && (is_prim_root(IXE,R) ))
                    return R;
            }
        }
        // _domain.assign(R[0],_domain.zero);
        // JGD 21.10.02
        _domain.assign(R[d],_domain.zero);
    }

        // Search for a monic irreducible Polynomial
        // with random Elements
    do {
        this->random( (RandIter&)_g, R, n); // must cast away const 
        _domain.assign(R[n.value()],_domain.one);
        for(Residu_t a=0; a<MOD; ++a) {
            _domain.assign(R[0],a);
            if (is_irreducible(R) && (is_prim_root(IXE,R) ))
                return R;
        }
    } while(1);
}

template<class Domain, class Tag, class RandIter >
inline typename Poly1FactorDom<Domain,Tag, RandIter>::Element& Poly1FactorDom<Domain,Tag, RandIter>::ixe_irreducible2 (Element& R, Degree n) const {
    init(R, n, _domain.one);
    Element IXE;
    init(IXE,Degree(1),_domain.one);
    Residu_t MOD = _domain.residu();

        // Search for an irreducible BINOMIAL : X^n + a
        // WARNING : Here we may have X^n + x, 
        // where a = representation of x, and sometimes a != x.
    for(Residu_t a=0; a<MOD; ++a) {
        _domain.assign(R[0],a);
        if (is_irreducible2(R) && (is_prim_root(IXE,R) ))
            return R;
    }
        // Search for an irreducible TRINOMIAL : X^n + b*X^i + a
        // Precondition : n >= 2
        // WARNING : same warning as for the binomial.
    for(Residu_t d=2;d<n.value();++d) {
        for(Residu_t b=0; b<MOD; ++b) {
            _domain.assign(R[d],b);
            for(Residu_t a=1; a<MOD; ++a) {
                _domain.assign(R[0],a);
                if (is_irreducible2(R) && (is_prim_root(IXE,R) ))
                    return R;
            }
        }
        _domain.assign(R[0],_domain.zero);
    }

        // Search for a monic irreducible Polynomial
        // with random Elements
    do {
        this->random( (RandIter&)_g, R, n); // must cast away const 
        _domain.assign(R[n.value()],_domain.one);
        for(Residu_t a=0; a<MOD; ++a) {
            _domain.assign(R[0],a);
            if (is_irreducible2(R) && (is_prim_root(IXE,R) ))
                return R;
        }
    } while(1);
}
// ---------------------------------------------------------------
// Irreducibility tests
// ---------------------------------------------------------------

template<class Domain, class Tag, class RandIter>
inline bool Poly1FactorDom<Domain,Tag, RandIter>::is_irreducible2(
    const Rep& P
    , Residu_t MOD ) const  {
        // Square free ?
    Rep W,D; this->gcd(W,diff(D,P),P);
    Degree d, dP;
    if (degree(d,W) > 0) return 0;
    IntFactorDom<> FD; 
    
    long n = degree(dP,P).value();
    IntFactorDom<>::Rep qn;

    FD.pow( qn, IntFactorDom<>::Rep(MOD), n);
    Rep Unit, G1; init(Unit, Degree(1), _domain.one);
    this->powmod(G1, Unit, qn, P);
    if (degree(d, sub(D,G1,Unit)) >= 0) return 0;    

    std::vector<IntFactorDom<>::Rep> Lp; std::vector<unsigned long> Le;
    FD.set(Lp, Le, n );
    for( std::vector<IntFactorDom<>::Rep>::const_iterator p = Lp.begin(); p != Lp.end(); ++p) {
        long ttmp; 
        FD.pow( qn, IntFactorDom<>::Rep(MOD), n/FD.convert(ttmp,*p) );
        this->powmod(G1, Unit, qn, P);
        if (degree(d, sub(D,G1,Unit)) < 0) return 0;    
    }
    
    return 1;
}




// ---------------------------------------------------------------
// Primitive Root over Z/pZ / F
// returns 1 if P is a generator. 
// ---------------------------------------------------------------

template<class Domain, class Tag, class RandIter>
bool Poly1FactorDom<Domain,Tag, RandIter>::is_prim_root( const Rep& P, const Rep& F)  const {
    bool isproot = 0;
    Rep A, G; mod(A,P,F);
    Degree d;
    if ( degree(d, this->gcd(G,A,F)) == 0) {
        Residu_t MOD = _domain.residu();
        IntFactorDom<> FD;
        IntFactorDom<>::Element IMOD( MOD ), q, qp;
        degree(d,F);
//         FD.pow(q ,IMOD, d.value());
//         FD.sub(qp, q, FD.one);
        FD.subin( FD.pow(qp ,IMOD, d.value()) , FD.one);
        std::list< IntFactorDom<>::Element > L;
        FD.set(L, qp);
        L.sort();
        std::list< IntFactorDom<>::Element >::iterator li = L.begin();
        isproot = 1;
        for(;(li != L.end()) && isproot; ++li)
            isproot = ( ! this->isOne(this->powmod(G, A, FD.div(q, qp , *li), F) ) );
    }
    return isproot;
}

template<class Domain, class Tag, class RandIter>
inline typename IntegerDom::Element Poly1FactorDom<Domain,Tag, RandIter>::order( const Rep& P, const Rep& F)  const {
    bool isproot = 0;
    Rep A, G; mod(A,P,F);
    Degree d;
    if ( degree(d, this->gcd(G,A,F)) == 0) {
        Residu_t MOD = _domain.residu();
        IntFactorDom<> FD;
        IntFactorDom<>::Element IMOD( MOD ), g, gg, tt, qp;
        degree(d,F);
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

template<class Domain, class Tag, class RandIter >
inline typename Poly1FactorDom<Domain,Tag, RandIter>::Rep& Poly1FactorDom<Domain,Tag, RandIter>::give_prim_root(Rep& R, const Rep& F)  const {
    Degree n; degree(n,F);
    Residu_t MOD = _domain.residu();
//    this->write(std::cout << "Give Pr: ", F) << std::endl;
    
    
        // Search for a primitive BINOMIAL : X^i + a
    for(Degree di=1;di<n;++di) {
        init(R, di, _domain.one);
//         for(Residu_t a=MOD; a--; ) {
        for(Residu_t a=0; a<MOD;++a ) {
            _domain.assign(R[0],a);
            if (is_prim_root(R,F))
                return R;
        }
    }
        // Search for a primitive TRINOMIAL : X^i + b*X^j + a
    for(Degree di=2;di<n;++di) {
        init(R, di, _domain.one);
        for(Degree dj=1;dj<di;++dj)
//             for(Residu_t b=MOD; b--;) {
            for(Residu_t b=0; b<MOD;++b) {
                _domain.assign(R[dj.value()],b);
//                 for(Residu_t a=MOD; a--;) {
                for(Residu_t a=0; a<MOD;++a ) {
                    _domain.assign(R[0],a);
                    if (is_prim_root(R,F))
                        return R;
                }
            }
    }

        // Search for a primitive Polynomial
        // with random Elements
    do {
        this->random( (RandIter&)_g, R, n); // must cast away const 
        _domain.assign(R[n.value()],_domain.one);
        for(Residu_t a=0; a<MOD; ++a) {
            _domain.assign(R[0],a);
            if (is_prim_root(R,F))
                return R;
        }
    } while(1);
}


template<class Domain, class Tag, class RandIter >
inline typename Poly1FactorDom<Domain,Tag, RandIter>::Rep& Poly1FactorDom<Domain,Tag, RandIter>::give_random_prim_root(Rep& R, const Rep& F)  const {
    Degree n; degree(n,F);
    Residu_t MOD = _domain.residu();
    
        // Search for a primitive Polynomial
        // with random Elements
    do {
        this->random( (RandIter&)_g, R, n); // must cast away const 
        _domain.assign(R[n.value()],_domain.one);
        for(Residu_t a=0; a<MOD; ++a) {
            _domain.assign(R[0],a);
            if (is_prim_root(R,F))
                return R;
        }
    } while(1);
}


template<class Domain, class Tag, class RandIter > 
inline typename Poly1FactorDom<Domain,Tag, RandIter>::Rep& Poly1FactorDom<Domain,Tag, RandIter>::random_prim_root(Rep& P, Rep& R, Degree n)  const {
        // P is irreducible
        // R is a primitive root. i.e R generates (Z_p)/P. 
        // returns R
    return give_prim_root(R, random_irreducible(P,n));
}


#endif 
