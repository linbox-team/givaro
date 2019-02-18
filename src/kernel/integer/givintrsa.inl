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

#ifndef __GIVARO_rsa_public_key_INL
#define __GIVARO_rsa_public_key_INL

#include <iostream>
#include <givaro/givinteger.h>
#include <givaro/givrandom.h>
#include <givaro/givintrsa.h>

namespace Givaro {

    // =================================================================== //
    // log[10]
    // =================================================================== //

    template<class MyRandIter>
    int64_t IntRSADom<MyRandIter>::log(const Element& n, const int64_t b ) const
    {
        int64_t res = 0;
        for(Element p = n; p>=b; ++res, this->divin(p,b) ) {}
        return res;
    }


    // =================================================================== //
    // Text conversions
    // =================================================================== //


    template<class MyRandIter>
    std::ostream& IntRSADom<MyRandIter>::ecriture_str(std::ostream& o, const Element& n) const
    {

        Element p = n, a, b;
        int64_t i = _lm-1;
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

    template<class MyRandIter>
    std::ostream& IntRSADom<MyRandIter>::ecriture_str_last(std::ostream& o, const Element& n) const
    {
        Element p = n, a, b;
        int64_t i = _lm-1, nbzeroes(0);
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


    template<class MyRandIter>
    std::ostream& IntRSADom<MyRandIter>::ecriture_Int(std::ostream& o, const Element& p) const
    {
        return o << p << std::endl;
    }




    // CBC mode enciphering
    template<class MyRandIter>
    std::ostream& IntRSADom<MyRandIter>::encipher(std::ostream& o, std::istream& in) const
    {
        MyRandIter generator;
        unsigned char x;
        Element res,r;
        o << generator.seed() << std::endl;
        Element ancien(generator());
        int64_t imax = (_lm-1)<<3;
        do {
            res = 0;
            for(int i=0; i<_lm-1; ++i) {
                x = (unsigned char) in.get();
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
    template<class MyRandIter>
    std::ostream& IntRSADom<MyRandIter>::decipher(std::ostream& o, std::istream& in)
    {
        Element r;
        uint64_t seed; in >> seed;
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

            Element ancien( generator() ), itmp;

            in >> itmp;
            do {
                powmod(r, itmp%p, _d, p);
                powmod(r2, itmp%q, _d, q);
                // Chinese reconstruction
                r2 -= r;
                r2 *= b;
                r2 %= _n;
                r2 += r;
                // Must be positive
                if (r2 > _n) r2 -= _n;
                else if (r2 < IntFactorDom<MyRandIter>::zero) r2+=_n;
                r2 ^= ancien;
                ancien = itmp;
                in >> itmp;
                if (in.eof()) {
                    // Treatment of trailing zeroes
                    ecriture_str_last(o, r2);
                    break;
                } else
                    ecriture_str(o, r2);
            } while (! in.eof());
        } else {
            Element ancien( generator() ), itmp;
            in >> itmp;
            do {
                powmod(r, itmp,_d,_n);
                r ^= ancien;
                ancien = itmp;
                in >> itmp;
                if (in.eof()) {
                    ecriture_str_last(o, r);
                    break;
                } else
                    ecriture_str(o, r^ancien);
            } while (! in.eof());
        }

        return o;
    }




    template<class MyRandIter>
    typename IntRSADom<MyRandIter>::Element& IntRSADom<MyRandIter>::strong_prime(random_generator& g, int64_t psize, Element& p) const
    {
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
        this->nextprimein( t );
        this->nextprimein( s );



        r = t<<1;
        ++r;
        while( ! this->isprime(r,4) ) {
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
        while( ! this->isprime(p,4) ) {
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
    template<class MyRandIter>
    void IntRSADom<MyRandIter>::keys_gen(random_generator& g, int64_t psize, int64_t qsize, Element& m, Element& k, Element& u) const
    {
        Element p, q;
        keys_gen(g,psize,qsize,m,k,u,p,q);
    }

    template<class MyRandIter>
    void IntRSADom<MyRandIter>::keys_gen(random_generator& g, int64_t psize, int64_t qsize, Element& m, Element& k, Element& u, Element& p, Element& q) const
    {
        Element d, l;

        strong_prime(g, psize, p);
        do  strong_prime(g, qsize, q); while (q == p);


        Element phim;
        this->mul(phim,
                  this-> sub(d,p,IntFactorDom<MyRandIter>::one),
                  this-> sub(l,q,IntFactorDom<MyRandIter>::one));
        this->mul(m, p, q);

        Element v, gd;

        if (_fast_impl) {
            this->mod(k,SIMPLE_EXPONENT, phim);
            this->gcd(gd,u,v,k,phim);
        } else {
            do {
                this->random(g,k,phim);
            } while (this->gcd(gd,u,v,k,phim) != 1);
        }
        this->modin(u,phim);
        if ( this->islt(u,IntFactorDom<MyRandIter>::zero) )
            this->addin(u,phim);
    }

    // =================================================================== //
    // Breaking codes
    // =================================================================== //
    template<class MyRandIter>
    typename IntRSADom<MyRandIter>::Element& IntRSADom<MyRandIter>::point_break(Element& u)
    {
        if ( isZero(_d) ) {
            Element p,v,d, pm;
            this->factor(p, _n);
            this->mul(pm, this->sub(v,p,IntFactorDom<MyRandIter>::one),
                      this->subin( this->div(d,_n,p), IntFactorDom<MyRandIter>::one ) );
            this->gcd(d,_d,v,_e,pm);
            if (this->islt(_d,IntFactorDom<MyRandIter>::zero))
                this->	addin(_d, pm);
        }
        return u = _d;
    }

} // namespace Givaro

#endif // __GIVARO_rsa_public_key_INL

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
