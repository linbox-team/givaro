// =================================================================== //
// Givaro : Euler's phi function
//          Primitive roots.
// Needs list structures : stl ones for instance
// Time-stamp: <17 Apr 03 14:33:45 Jean-Guillaume.Dumas@imag.fr> 
// =================================================================== //
#include "givintnumtheo.h"
#include <list>
#include <vector>

#include "givintrns.h"
#include "givpower.h"


// =================================================================== //
// Euler's phi function
// =================================================================== //
template<class RandIter>
typename IntNumTheoDom<RandIter>::Rep& IntNumTheoDom<RandIter>::phi(Rep& res, const Rep& n) const {
    if (isleq(n,1)) return res=n;
    if (isleq(n,3)) return sub(res,n,one);
    std::list<Rep> Lf;
    set(Lf,n);
    return phi(res,Lf,n);
}


template<class RandIter>
template< template<class> class Container> typename IntNumTheoDom<RandIter>::Rep& IntNumTheoDom<RandIter>::phi(Rep& res, const Container<Rep>& Lf, const Rep& n) const {
    if (isleq(n,1)) return res=n;
    if (isleq(n,3)) return sub(res,n,one);
    res = n; Rep t,m;
    for(typename Container<Rep>::const_iterator f=Lf.begin(); f!=Lf.end(); ++f) 
        mul(res, divexact(t,res,*f), sub(m, *f, one));
    return res;
}

// =================================================================== //
// Möbius function
// =================================================================== //
template<class RandIter>
template< template<class> class Container> short IntNumTheoDom<RandIter>::mobius(const Container<unsigned long>& lpow) const {
    if (lpow.size()) {
        short mob = 1;
        for(typename Container<unsigned long>::const_iterator i=lpow.begin();i != lpow.end(); ++i) {
            if (*i > 1) {
                 return 0;
            } else
                mob = -mob;
        }
        return mob;
    } else
        return 1;
}
    
template<class RandIter>
short IntNumTheoDom<RandIter>::mobius(const Rep& a) const {
    std::list< Rep> lr;
    std::list<unsigned long> lp;
    set(lr, lp, a);
    return mobius(lp);
}


// =================================================================== //
// Primitive Root
// =================================================================== //

template<class RandIter>
typename IntNumTheoDom<RandIter>::Rep& IntNumTheoDom<RandIter>::prim_root(Rep& A, unsigned long& runs, const Rep& n) const {
        // n must be in {2,4,p^m,2p^m} where p is an odd prime
        // else infinite loop

    if (isleq(n,4)) return sub(A,n,one);
    if (iszero(mod(A,n,4))) return A=Rep(zero);
    Rep p,ismod2, q, no2, root; 
    if (iszero(mod(ismod2,n,2))) divexact(no2,n,2); else no2=n;
    p=no2;
    int k = 1; 
    while (! isprime(p) ) {
        sqrt(root, p);
        while (mul(q,root,root) == p) {
            p = root;
            sqrt(root,p);
        }
        if (! isprime(p) ) {
            q=p;
            while( p == q ) factor(p, q);
            divin(q,p);
            if (q < p) p = q;
        }
    }
    if (iszero(ismod2)) mul(q,p,2); else q=p;
    for(;q != n;++k,q*=p);
    Rep phin, tmp; 
    phi(phin,p);
    std::list<Rep> Lf;
    set(Lf,phin);
    typename std::list<Rep>::iterator f;
    for(f=Lf.begin();f!=Lf.end();++f)
            div(*f,phin,*f);
    int found; runs = 0;
    A=2;
    found = ++runs;
    for(f=Lf.begin();(f!=Lf.end() && found);f++)
        found = (! isone( powmod(tmp,A,*f,p)) );
    if (! found) {
        A=3;
        found = ++runs;
        for(f=Lf.begin();(f!=Lf.end() && found);f++)
            found = (! isone( powmod(tmp,A,*f,p)) );
    }
    if (! found) {
        A=5;
        found = ++runs;
        for(f=Lf.begin();(f!=Lf.end() && found);f++)
            found = (! isone( powmod(tmp,A,*f,p)) );
    }
    if (! found) {
        A=6;
        found = ++runs;
        for(f=Lf.begin();(f!=Lf.end() && found);f++)
            found = (! isone( powmod(tmp,A,*f,p)) );
    }
    while (! found) {
       do {
            random(_g, A, p);
            addin( modin(A,sub(tmp,p,7)) , 7);
        } while ( ! isone(gcd(tmp,A,p)) );
        found = ++runs;
        for(f=Lf.begin();(f!=Lf.end() && found);f++)
            found = (! isone( powmod(tmp,A,*f,p)) );
    }
    if (k == 1) {
        if (iszero(ismod2) && iszero(mod(ismod2, A, 2)))
            return A+=p;
        else
            return A;
    } else {
        if (! is_prim_root(A,no2))
            A+=p;
        if (iszero(ismod2) && iszero(mod(ismod2, A, 2)))
            return A+=no2;
        else
            return A;
    }
}

template<class RandIter>
typename IntNumTheoDom<RandIter>::Rep& IntNumTheoDom<RandIter>::prim_root(Rep& A, const Rep& n) const { unsigned long runs; return prim_root(A, runs, n); }


template<class RandIter>
typename IntNumTheoDom<RandIter>::Rep& IntNumTheoDom<RandIter>::prim_root_of_prime(Rep& A, const Rep& n) const { 
    
    std::vector<Rep> Lf;
    Rep phin; sub(phin,n,one);
    set(Lf,phin);
    return prim_root_of_prime(A, Lf, phin, n);
}


inline Integer& ppin(Integer& res, const Integer& prime) {
    IntegerDom I;
    Integer tmp;
    while( I.iszero(I.mod(tmp, res, prime)) ) {
        I.divexact(res, tmp = res, prime);
    }
    return res;
}

    
/// Add Jacobi for quadratic nonresidue    

template<class RandIter>
template<class Array>
typename IntNumTheoDom<RandIter>::Rep& IntNumTheoDom<RandIter>::prim_root_of_prime(Rep& A, const Array& aLf, const Rep& phin, const Rep& n) const { 
    
    Rep tmp, expo, temp;
    A = one;

    Rep prime(2), Aorder=one;
    std::vector<Rep> Lf = aLf, newLf, oldLf; 
    newLf.reserve(Lf.size()); 
    oldLf.reserve(Lf.size()); 

    Rep primeorder;
    
    for(bool exemp = true; exemp; nextprimein(prime) ) {
        A = prime;
        primeorder = phin;
        for(typename Array::const_iterator f = Lf.begin(); f != Lf.end(); ++f) {
            powmod(tmp, prime, div(expo, primeorder, *f), n);
            if (isone(tmp)) {
                newLf.push_back(*f);
                while (iszero(mod(tmp,expo,*f)) && isone( powmod(tmp, prime, div(temp, expo, *f), n) ) ) { expo = temp; }
                primeorder = expo;
//                 std::cerr << "2 Order (Div): " << primeorder << std::endl;
            } else {
                oldLf.push_back(*f);
                exemp = false;
//                 std::cerr << "2 Order : " << primeorder << std::endl;
            }
        }
    }
    
    Aorder = primeorder;
    Lf = newLf;

        // Now we have A with order > 1, we need to add other primes

//     std::cerr << "Prime : 2" << std::endl;
//     std::cerr << "Root : " << A << std::endl;
//     std::cerr << "Order : " << Aorder << std::endl;
    
    for ( ; islt(Aorder,phin); nextprimein(prime) ) {
        newLf.resize(0); oldLf.resize(0);

        for(typename Array::const_iterator f = Lf.begin(); f != Lf.end(); ++f) {
            powmod(tmp, prime, div(expo, phin, *f), n);
            if (isone(tmp)) {
                newLf.push_back(*f);
            } else {
                oldLf.push_back(*f);
            }
        }

        if (oldLf.size() > 0) {
            Rep g = phin;

            for(typename Array::const_iterator f = oldLf.begin(); f != oldLf.end(); ++f) {
                ppin(g, *f);
                ppin(Aorder, *f);
            }

            powmod(tmp, prime, g, n);

            modin( mulin(A, tmp), n );

            mulin(Aorder, div(tmp, phin, g));

            Lf = newLf;
        }

//         std::cerr << "Prime : " << prime << std::endl;
//         std::cerr << "Root : " << A << std::endl;
//         std::cerr << "Order : " << Aorder << std::endl;

    }

    
    return A;
}

    
        
template<class RandIter>
typename IntNumTheoDom<RandIter>::Rep& IntNumTheoDom<RandIter>::lowest_prim_root(Rep& A, const Rep& n) const {
        // n must be in {2,4,p^m,2p^m} where p is an odd prime
        // else returns zero
    if (isleq(n,4)) return sub(A,n,one);
    if (iszero(mod(A,n,4))) return A=zero;
    Rep phin, tmp; 
    phi(phin,n);
    std::list<Rep> Lf;
    set(Lf,phin);
    typename std::list<Rep>::iterator f;
    for(f=Lf.begin();f!=Lf.end();++f)
            div(*f,phin,*f);
    int found=0;
    for(A = 2;(isleq(A,n) && (! found));addin(A,1)) {
        if (isone(gcd(tmp,A,n))) {
            found = 1;
            for(f=Lf.begin();(f!=Lf.end() && found);f++)
                found = (! isone( powmod(tmp,A,*f,n)) );
        }
    }
    if (isleq(A,n))
        return subin(A,1);
    else
        return A=zero; 
}

template<class RandIter>
bool IntNumTheoDom<RandIter>::is_prim_root(const Rep& p, const Rep& n) const {
        // returns 0 if failed
    bool found=false;
    Rep phin, tmp; 
    phi(phin,n);
    std::list<Rep> Lf;
    set(Lf,phin);
    typename std::list<Rep>::iterator f=Lf.begin();
    Rep A; mod(A,p,n);
    if (isone(gcd(tmp,A,n))) {
        found = true;
        for(;(f!=Lf.end() && found);f++) {
//             found = ( powmod(A,phin / (*f),n) != 1);
            found = (! isone( powmod(tmp,A, div(tmp,phin,*f),n)) );
        }
    }
    return found;
}

template<class RandIter>
bool IntNumTheoDom<RandIter>::isorder(const Rep& g, const Rep& p, const Rep& n) const {
        // returns 1 if p is of order g in Z/nZ
    Rep tmp;
    return (isone( pow(tmp, p, Integer2long(g)) ) && areEqual( g, order(tmp,p,n) ) );
}

template<class RandIter>
typename IntNumTheoDom<RandIter>::Rep& IntNumTheoDom<RandIter>::order(Rep& g, const Rep& p, const Rep& n) const {
        // returns 0 if failed
    Rep A; mod(A,p,n);
    if (iszero(A))
	return g = zero;
    if (isone(A))
	return g = one;
    int noprimroot=0;
    Rep phin,gg,tmp;
    phi(phin,n);
    std::list<Rep> Lf;
    set(Lf,phin);
    Lf.sort();
    typename std::list<Rep>::iterator f=Lf.begin();
    if (isone(gcd(tmp,A,n))) {
        noprimroot = 0;
        for(;f!=Lf.end();++f)
            if ( noprimroot = isone(powmod(tmp,A, div(g,phin,*f),n)) )
                break;
        if (noprimroot) {
            for(;f!=Lf.end();++f)
                while (iszero(mod(tmp,g,*f)) && isone( powmod(tmp,A, div(gg,g,*f),n) ) )
                    g = gg;
            return g;
        } else
            return g=phin;
    }
    return g=zero;
}

template<class RandIter>
typename IntNumTheoDom<RandIter>::Rep& IntNumTheoDom<RandIter>::prim_inv(Rep& A, const Rep& n) const {
    if (isleq(n,4)) return sub(A,n,this->one);
    if (areEqual(n,8)) return init(A,3);
    return prim_base(A, n);
}

template<class RandIter>
typename IntNumTheoDom<RandIter>::Rep& IntNumTheoDom<RandIter>::prim_elem(Rep& A, const Rep& n) const {
    if (isleq(n,4)) { 
        Rep tmp; 
        return this->sub(A,n,this->one); 
    }
    
    if (areEqual(n,8)) return init(A,2);
    return prim_base(A, n);
}

template<class RandIter>
typename IntNumTheoDom<RandIter>::Rep& IntNumTheoDom<RandIter>::prim_base(Rep& A, const Rep& m) const {
        // Prerequisite : m > 4, and m != 8.
    std::vector<Rep> Lp; std::vector<unsigned long> Le;
    unsigned long nbf = set(Lp, Le, m);
    std::vector<Rep> Pe(nbf); 
    std::vector<Rep> Ra(nbf);
    
    typename std::vector<Rep>::const_iterator p=Lp.begin();
    typename std::vector<unsigned long>::const_iterator e=Le.begin() ;
    typename std::vector<Rep>::iterator pe = Pe.begin();
    typename std::vector<Rep>::iterator a = Ra.begin() ;
    for( ;p!=Lp.end();++p, ++e, ++pe, ++a) {
        dom_power( *pe, *p, *e, *this);
        if (areEqual(*p,2))
            init(*a, 3);
        else
            prim_root(*a, *pe);
    }
    
    IntRNSsystem<std::vector> RNs( Pe );
    RNs.RnsToRing( A, Ra );
    return A;
}


template<class RandIter>
typename IntNumTheoDom<RandIter>::Rep& IntNumTheoDom<RandIter>::zeta_primpow(Rep & z, const Rep& p, const unsigned long e) const {
        // Prerequisite : p prime.
    if (areEqual(p, 2)) {
        if (e<=3) return init(z,e);
        return dom_power(z, p, e-2, *this);
    } else {
        Rep tmp;
        return mulin( dom_power(z, p, e-1, *this), sub(tmp, p, one) );
    }
}

template<class RandIter>
typename IntNumTheoDom<RandIter>::Rep& IntNumTheoDom<RandIter>::zeta_inv_primpow(Rep & z, const Rep& p, const unsigned long e) const {
        // Prerequisite : p prime.
    if (areEqual(p, 2)) {
        if (e<=2) return init(z,e);
        if (e==3) return init(z,2);
        return dom_power(z, p, e-2, *this);
    } else {
        Rep tmp;
        return mulin( dom_power(z, p, e-1, *this), sub(tmp, p, one) );
    }
}



    
template<class RandIter>
typename IntNumTheoDom<RandIter>::Rep& IntNumTheoDom<RandIter>::zeta_inv(Rep & z, const Rep& m) const {
        if (areEqual(m,2)) return init(z,1);
        if (areEqual(m,3) || areEqual(m,4) || areEqual(m,8) ) return init(z,2);
        return zeta_base(z, m);
}

template<class RandIter>
typename IntNumTheoDom<RandIter>::Rep& IntNumTheoDom<RandIter>::zeta(Rep & z, const Rep& m) const {
        if (areEqual(m,2)) return init(z,1);
        if (areEqual(m,3) || areEqual(m,4)) return init(z,2);
        if (areEqual(m,8) ) return init(z,3);
        return zeta_base(z, m);
}


template<class RandIter>
typename IntNumTheoDom<RandIter>::Rep& IntNumTheoDom<RandIter>::zeta_base(Rep & z, const Rep& m) const {
        // Prerequisite: m > 4, and m != 8.
        std::vector<Rep> Lp; std::vector<unsigned long> Le;
        unsigned long nbf = set(Lp, Le, m);

        zeta_inv_primpow(z, Lp.front(), Le.front() );

        if (nbf == 1) return z;

        typename std::vector<Rep>::const_iterator p=Lp.begin();
        typename std::vector<unsigned long>::const_iterator e=Le.begin() ;
        for( ++p, ++e; p != Lp.end(); ++p, ++e) {
            Rep tmp;
            zeta_inv_primpow(tmp, *p, *e);
            Rep g;
            gcd(g, z, tmp);
            mulin(z, divin(tmp, g));
        }
        
        return z;
}
