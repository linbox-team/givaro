// =================================================================== //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.
// Givaro : Prime numbers
//              RSA public-key cipher codes
//              ECB mode (UNSECURE !!!)
// Time-stamp: <07 May 09 13:44:00 Jean-Guillaume.Dumas@imag.fr> 
// =================================================================== //

#ifndef _GIVARO_RSA_Public_KEY_
#define _GIVARO_RSA_Public_KEY_

#include <iostream>
#include <givaro/givinteger.h>
#include <givaro/givrandom.h>
#include <givaro/givintrsa.h>

// =================================================================== //
// log[10]
// =================================================================== //

template<class RandIter>
long IntRSADom<RandIter>::log(const Element& n, const long b = 10) const {
    long res = 0;
    for(Element p = n; p>=b; ++res, this->divin(p,b) ) {}
    return res;
}   
    

// =================================================================== //
// Text conversions
// =================================================================== //

                
template<class RandIter>
std::ostream& IntRSADom<RandIter>::ecriture_str(std::ostream& o, const Element& n) const {

    Element p = n, a, b;
    long i = _lm-1;
        // First char is ignored as it is zero or this->random enabling CBC
    a = p >> (i<<3);
    p -= a << (i<<3);
    for(--i; i>=0; --i) {
	a = p >> (i<<3);
    	o << char( int(a) );
	p -= a << (i<<3);
    }
    return o;
}

template<class RandIter>
std::ostream& IntRSADom<RandIter>::ecriture_str_last(std::ostream& o, const Element& n) const {
    Element p = n, a, b;
    long i = _lm-1, nbzeroes(0);
        // First char is ignored as it is zero or this->random enabling CBC
    a = p >> (i<<3);
    p -= a << (i<<3);
    for(--i; i>=1; --i) {
	a = p >> (i<<3);
            // Treatment of trailing zeroes
        if ( char(int(a)) ) {
            for(int j =0; j <nbzeroes; ++j)
                o << char(0);
            o << char( int(a) );
            nbzeroes = 0;
            p -= a << (i<<3);
        } else {
            ++nbzeroes;
        }
    }
    return o << char(int(p));
}


template<class RandIter>
std::ostream& IntRSADom<RandIter>::ecriture_Int(std::ostream& o, const Element& p) const {
    return o << p << std::endl;
}



    
// CBC mode enciphering
template<class RandIter>
std::ostream& IntRSADom<RandIter>::encipher(std::ostream& o, std::istream& in) const {
    RandIter generator;
    unsigned char x;
    Element res,r;
    o << generator.seed() << std::endl;
    Element ancien(generator());
    int imax = (_lm-1)<<3;
    do { 
        res = 0;
        for(int i=0; i<_lm-1; ++i) {
            x = in.get();
            if (in.eof()) {
                    // Adding zeroes at end of file
                res <<= ( 8*(_lm-1-i) );
                break;
            }
            res <<= 8;
            res += x;
        }

            // Padding randomly to enable CBC
            // indeed, decryption will return (res^ancien) % _n
        while ( (res^ancien) > _n) {
            res += ((Integer)((unsigned char)generator()) << imax);
        }

        powmod(r, res^ancien,_e,_n);
        ecriture_Int(o, r );
        ancien = r;
    } while (! in.eof());
    return o;
}           

// CBC mode deciphering
template<class RandIter>
std::ostream& IntRSADom<RandIter>::decipher(std::ostream& o, std::istream& in) {
    double length = _lm * 2.4082399653118495617; // _lm * 8*log[10](2)
    char * tmp = new char[(long)length+2];
    Element r;
    unsigned long seed; in >> seed;
    GivRandom generator(seed);

    if (_fast_impl) {
        Element p, q, pr, k, phi, t, delt ;
        pr = _e*_d;
	--pr;
        k = pr/_n;
	do {
            ++k;
            phi = pr/ k;
        } while (k*phi != pr);
        t = _n-phi; ++t;
        t >>= 1;
        sqrt(delt, t*t-_n);
        p = t-delt;
        q = t+delt;
        Element gd, a, b, c, r1, r2;
        this->gcd(gd, a, b, q, p);

        b *= p;
        b %= _n;

        Element ancien( generator() );
        
        in >> tmp; 
        do {
            powmod(r, Integer(tmp)%p, _d, p);
            powmod(r2, Integer(tmp)%q, _d, q);
	    // Chinese reconstruction
	    r2 -= r;
            r2 *= b;
            r2 %= _n;
            r2 += r;
	    // Must be positive
	    if (r2 > _n) r2 -= _n;
	    else if (r2 < IntFactorDom<RandIter>::zero) r2+=_n;
            r2 ^= ancien;
            ancien = Integer(tmp);
            in >> tmp;
            if (in.eof()) {
                    // Treatment of trailing zeroes
                ecriture_str_last(o, r2);
                break;
            } else
                ecriture_str(o, r2);
        } while (! in.eof());
    } else {      
        Element ancien( generator() );
        in >> tmp;
        do {
            powmod(r, Integer(tmp),_d,_n);
            r ^= ancien;
            ancien = Integer(tmp);
            in >> tmp;
            if (in.eof()) {
                ecriture_str_last(o, r);
                break;
            } else
                ecriture_str(o, r^ancien);
        } while (! in.eof());
    }
    
    return o;
}        




template<class RandIter>
typename IntRSADom<RandIter>::Element& IntRSADom<RandIter>::strong_prime(random_generator& g, long psize, Element& p) const {
    Element q,t,r,s;

    if (psize > 9) {
        this->random(g,t,(psize>>1)-2);
        t += Integer(1)<<( (psize>>1)-2 );
        this->random(g,s,(psize>>1)-2);
        s += Integer(1)<<( (psize>>1)-2 );
    } else {
        this->random(g,t,3); 
        this->random(g,s,3); 
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
void IntRSADom<RandIter>::keys_gen(random_generator& g, long psize, long qsize, Element& m, Element& k, Element& u) const {
    Element p, q;
    keys_gen(g,psize,qsize,m,k,u,p,q);
}

template<class RandIter>
void IntRSADom<RandIter>::keys_gen(random_generator& g, long psize, long qsize, Element& m, Element& k, Element& u, Element& p, Element& q) const {
    Element d, l;

    strong_prime(g, psize, p);
    do  strong_prime(g, qsize, q); while (q == p);
    

    Element phim; mul(phim, sub(d,p,IntFactorDom<RandIter>::one), sub(l,q,IntFactorDom<RandIter>::one));
    mul(m, p, q);

    Element v, gd;

    if (_fast_impl) {
        mod(k,SIMPLE_EXPONENT, phim);
        this->gcd(gd,u,v,k,phim);
    } else {
        do {
            this->random(g,k,phim);
        } while (this->gcd(gd,u,v,k,phim) != 1);
    }
    modin(u,phim);
    if ( islt(u,IntFactorDom<RandIter>::zero) ) addin(u,phim);
}

// =================================================================== //
// Breaking codes
// =================================================================== //
template<class RandIter>
typename IntRSADom<RandIter>::Element& IntRSADom<RandIter>::point_break(Element& u) {
    if ( isZero(_d) ) {
        Element p,v,d, pm;
        factor(p, _n);
        mul(pm, sub(v,p,IntFactorDom<RandIter>::one), subin( this->div(d,_n,p), IntFactorDom<RandIter>::one ) );
        this->gcd(d,_d,v,_e,pm);
        if (islt(_d,IntFactorDom<RandIter>::zero)) addin(_d, pm);
    }
    return u = _d;
}

#endif
