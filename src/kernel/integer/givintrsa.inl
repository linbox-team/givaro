// =================================================================== //
// Givaro : Prime numbers
//              RSA public-key cipher codes
// Time-stamp: <16 Sep 03 19:10:58 Jean-Guillaume.Dumas@imag.fr> 
// =================================================================== //

#ifndef _GIVARO_RSA_Public_KEY_
#define _GIVARO_RSA_Public_KEY_

#include <iostream>
#include "givinteger.h"
#include "givintrsa.h"

// =================================================================== //
// log[10]
// =================================================================== //

template<class RandIter>
long IntRSADom<RandIter>::log(const element& n, const long b = 10) const {
    long res = 0;
    for(element p = n; p>=b; ++res, divin(p,b) ) {}
    return res;
}   
    

// =================================================================== //
// Text conversions
// =================================================================== //

                
template<class RandIter>
std::ostream& IntRSADom<RandIter>::ecriture_str(std::ostream& o, const element& n) const {

    element p = n, a, b;
    for(long i = _lm; --i>=0; ) {
//	    std::cerr << "pi: " << p << std::endl;
	a = p >> (i<<3);
//	    std::cerr << "a: " << a << std::endl;
    	o << char( int(a) );
	p -= a << (i<<3);
//	    std::cerr << "pf: " << p << std::endl;
    }
    return o;
}


template<class RandIter>
std::ostream& IntRSADom<RandIter>::ecriture_Int(std::ostream& o, const element& p) const {
    return o << p << std::endl;
}



    

template<class RandIter>
std::ostream& IntRSADom<RandIter>::encipher(std::ostream& o, std::istream& in) const {
    unsigned char x;
    element res,r;
    
    do { 
        res = 0;
        for(int i=0; i<_lm; ++i) {
            x = in.get();
            if (in.eof())
                break;
            res <<= 8;
            res += x;
        }
        ecriture_Int(o, powmod(r, res,_k,_m) );
    } while (! in.eof());
    
    return o;
}           

template<class RandIter>
std::ostream& IntRSADom<RandIter>::decipher(std::ostream& o, std::istream& in) {
    double length = _lm * 2.4082399653118495617; // _lm * 8*log[10](2)
    char * tmp = new char[(long)length+2];
    element r;

    if (_fast_impl) {
        element p, q, pr, k, phi, t, delt ;
        pr = _k*_u-1;
        k = pr/_m;
        ++k;
        phi = pr/ k;
        t = _m-phi+1;
        t >>= 1;
        sqrt(delt, t*t-_m);
        p = t-delt;
        q = t+delt;
//         std::cerr << "p: " << p << ", q: " << q << std::endl;
        element gd, a, b, c, r1, r2;
        gcd(gd, a, b, q, p);
        a *= q;
        a %= _m;
        b *= p;
        b %= _m;

        do {
            in >> tmp;
            if (in.eof()) break;
            powmod(r, Integer(tmp)%p, _u, p);
            powmod(r2, Integer(tmp)%q, _u, q);
            r *= a;
            r %= _m ;
            r2 *= b;
            r2 %= _m;
            r += r2;
            if (r > _m) r -= _m;
            else if (r < zero) r+=_m;
//             powmod(r, Integer(tmp),_u,_m);
            ecriture_str(o, r);
        } while (! in.eof());
        
    } else {
        do {
            in >> tmp;
            if (in.eof()) break;
            ecriture_str(o,powmod(r, Integer(tmp),_u,_m));
        } while (! in.eof());
    }
    
    return o;
}        




template<class RandIter>
typename IntRSADom<RandIter>::element& IntRSADom<RandIter>::strong_prime(random_generator& g, long psize, element& p) const {
    element q,t,r,s;

    if (psize > 9) {
        random(g,t,(psize>>1)-2);
        t += Integer(1)<<( (psize>>1)-2 );
        random(g,s,(psize>>1)-2);
        s += Integer(1)<<( (psize>>1)-2 );
    } else {
        random(g,t,3); 
        random(g,s,3); 
    }
    nextprimein( t );
    nextprimein( s );
    
    

    r = t<<1;
    ++r;
    while( ! isprime(r,4) ) {
        r += t<<1;
    }




        // q = 2(s^(r-2) mod r)s - 1
    powmod(q, s, r-2, r);
    q <<= 1;
    q *= s;
    --q;

    
        // 2rs
    r *= s;
    r <<= 1;
    
    p = q+r;
    while( ! isprime(p,4) ) {
        p += r;
    }

    return p;
}

    


// =================================================================== //
// Keys generation
// public keys are m and k, the secret key is u.
// ciphering is computing	: x^k mod m
// deciphering is computing	: b^u mod m
// since for any x, x^(k.u) = x mod m
// =================================================================== //
template<class RandIter>
void IntRSADom<RandIter>::keys_gen(random_generator& g, long psize, long qsize, element& m, element& k, element& u) const {
    element p,q, d, l;

    strong_prime(g, psize, p);
    do  strong_prime(g, qsize, q); while (q == p);
    

    element phim; mul(phim, sub(d,p,one), sub(l,q,one));
    mul(m, p, q);

    element v, gd;

    if (_fast_impl) {
        mod(k,SIMPLE_EXPONENT, phim);
        gcd(gd,u,v,k,phim);
    } else {
        do {
            random(g,k,phim);
        } while (gcd(gd,u,v,k,phim) != 1);
    }
    modin(u,phim);
    if ( islt(u,zero) ) addin(u,phim);
}

// =================================================================== //
// Breaking codes
// =================================================================== //
template<class RandIter>
typename IntRSADom<RandIter>::element& IntRSADom<RandIter>::point_break(element& u) {
    if ( iszero(_u) ) {
        element p,v,d, pm;
        factor(p, _m);
        mul(pm, sub(v,p,one), subin( div(d,_m,p), one ) );
        gcd(d,_u,v,_k,pm);
        if (islt(_u,zero)) addin(_u, pm);
    }
    return u = _u;
}

#endif
