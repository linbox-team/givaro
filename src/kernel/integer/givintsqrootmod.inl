// ============================================================= //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <12 Mar 20 15:44:49 Jean-Guillaume.Dumas@imag.fr>
// Givaro : Modular square roots
// Author : Yanis Linge, Jean-Guillaume Dumas
// ============================================================= //

#ifndef __GIVARO_sqrootmod_INL
#define __GIVARO_sqrootmod_INL

#include <givaro/givtimer.h>

namespace Givaro {


        // ======================================================== //
        // Modular Square root functions
        // ======================================================== //
    template <class MyRandIter> inline typename IntSqrtModDom<MyRandIter>::Rep &
    IntSqrtModDom<MyRandIter>::sqrootmod (Rep & x, const Rep & a, const Rep & n) const {
        std::vector < Rep > Lf;
        std::vector < uint64_t > Le;
        Father_t::set (Lf, Le, n);

        typename std::vector < Rep >::const_iterator Lf_iter = Lf.begin ();
        typename std::vector < uint64_t >::const_iterator Le_iter = Le.begin ();

        std::vector < Rep > roots;
        Rep tmp;

            // Build prime powers
        std::vector < Rep > Pe (Lf.size ());
        typename std::vector < Rep >::iterator Pe_iter = Pe.begin ();
        for (; Pe_iter != Pe.end (); ++Pe_iter, ++Lf_iter, ++Le_iter)
            dom_power (*Pe_iter, *Lf_iter, (long)*Le_iter, *this);

        Lf_iter = Lf.begin ();
        Le_iter = Le.begin ();
        Pe_iter = Pe.begin ();

            // roots mod powers of primes
        for (; Lf_iter != Lf.end (); ++Lf_iter, ++Le_iter, ++Pe_iter){
                // root mod a power of 2
            if (*Lf_iter == 2U){
                roots.push_back (
                    this->sqrootmodpoweroftwo (tmp, a, *Le_iter, *Pe_iter));
            } else {
                roots.push_back (
                    this->sqrootmodprimepower (tmp, a, *Lf_iter, *Le_iter, *Pe_iter));
            }
        }

            // Chinese Remaindering
        IntRNSsystem < std::vector, std::allocator > RNs (Pe);

        RNs.RnsToRing (x, roots);
        x = (x<0?-x:x);
        return x;
    }




    template <class MyRandIter> inline typename IntSqrtModDom<MyRandIter>::Rep &
    IntSqrtModDom<MyRandIter>::sqrootmodprime (Rep & x,
                                               const Rep & a,
                                               const Rep & p) const {
//         std::cerr << "p:= " << p << ';' << std::endl;
//         std::cerr << "a:= " << a << ';' << std::endl;
        Rep amp (a); Integer::modin(amp,p);
        if (amp == 0U || amp == 1) return x = amp;

        if (legendre (amp, p) == -1){
            std::cerr << amp << " is not a quadratic residue mod " << p << " (SQRMP)" << std::endl;
            return x = -1;
        }

        if ((p & 3U) == 3U) {			// If p = 3 mod 4
            Rep ppu (p); ++ppu; ppu >>= 2U;	// ppu = (p+1)/4
            return powmod (x, amp, ppu, p);		// powmod (x,a,(p+1)/4,p);
        }

            // O. Atkin
        if ((p & 7U) == 5U) {			// If p = 5 mod 8
            Rep tmp;
            Rep puis (p); puis -= 1; puis >>= 2U;// puis = (p-1)/4
            powmod (tmp, amp, puis, p);

            if (tmp == 1) {
                puis = p; puis += 3U; puis >>= 3U;// puis = (p+3)/8
                return powmod (x, amp, puis, p);
            }
            puis = p; puis -= 5U; puis >>= 3U;	// puis = (p-5)/8

            Rep a4 (amp); a4 <<= 2;
            powmod (x, a4, puis, p);
            x *= amp; x <<= 1;			// 2a(4a)^{(p-5)/8}
            return x %= p;
        }

        size_t l = (size_t) ceil (logtwo (p) - 1);

            // S. Mueller
        if ((p & 15U) == 9U) {			// If p = 9 mod 16
            Rep i (amp); i <<= 1;
            Rep puis (p); puis -= 1; puis >>= 2U;// puis = (p-1)/4
            powmod (x, i, puis, p);			// (2a)^{(p-1}/4} is +1 or -1
            if (x != 1) x = -1;

            Rep d; while (legendre (Rep::nonzerorandom (d, l), p) == x) ;
            puis = p; puis -= 9U; puis >>= 4U;	// puis = (p-9)/16
            i *= d; i *= d;
            powmod(x, i, puis, p);			// (2d^2a)^{(p-9)/16}
            i *= x; i%=p; i *= x; i%=p;		// i=2d^2x^2a ; i^2 = -1
            --i;
            x *= d; x%=p; x *= i; x%=p; x *= amp;	// xda(i-1)
            return x %= p;				// +/- x is a root
        }

            // Tonelli and Shanks
            // [H. Cohen, Algorithm 1.5.1, p33,
            // A course in computational algebraic number theory]
        Rep p1 (p); --p1;
        Rep q (p1);
        int64_t e (0);
        for( ; (q & 1U) == 0; ++e) q >>= 1;

            // now we have e and q such that : p-1=q*2^e with q odd
            // we need a non quadratic element : g
        Rep g; while (legendre (Rep::nonzerorandom (g, l), p) != -1) ;

        Rep z;
        powmod (z, g, q, p);	// z = g^q mod p
            //Initialize
        Rep y (z);
        Rep tmp (q); tmp -= 1; tmp >>= 1;
        powmod (x, amp, tmp, p);	// a^{(q-1)/2} mod p
        Rep b (x);
        b *= x; b *= amp; b %= p;	// ax^2
        x *= amp; x %= p;		// ax

            // Find exponent
        int64_t m(1), r(e);
        Rep b2k, t, puis(r);
        while (b != 1){
            b2k = b;
            for(m = 0; b2k != 1; ++m) {
                b2k *= b2k; b2k %= p;
            } // m smallest such that b^{2^m} is 1 mod p
            if (m == r){
                std::cerr << amp << " is not a quadratic residue mod " << p << " (NoExp)" << std::endl;
                return x = -1;
            }
            int64_t lpuis = r; lpuis -= m; --lpuis;
            puis = 1; puis <<= lpuis;	// 2^{m-r-1}
            powmod (t, y, puis, p);		// t = y^{ 2^{m-r-1} } mod p
            y = t; y *= t; y %= p;		// y = t^2 mod p
            r = m;				// r = m
            x *= t; x %= p;			// x = xt mod p
            b *= y; b %= p;			// b = by mod p
        }
        return x;
    }

    template <class MyRandIter> inline typename IntSqrtModDom<MyRandIter>::Rep &
    IntSqrtModDom<MyRandIter>::sqrootmodprimepower (Rep & x,
                                                    const Rep & a,
                                                    const Rep & p,
                                                    const uint64_t k,
                                                    const Rep & pk) const{

        Rep tmpa(a); Integer::modin(tmpa,pk);
        if(tmpa==0) return x=0;
        if(tmpa==1) return x=1;
        if (k == 1) return sqrootmodprime (x, tmpa, p);

        if ((tmpa%p)==0){
            Rep b(tmpa);
            uint64_t t=0;
            for( ; (b%p) == 0; ++t) b/=p; // a = b p^t and p does not divide b

            if((t&1)==0){
                Rep sqrtb;
                sqrootmodprimepower(sqrtb,b,p,k,pk);
                powmod(x,p,(t>>1),pk);
                x*=sqrtb;
                return x%=pk;
            }
            else{
                std::cerr <<tmpa << "is not a quadratic residue mod " << pk << " (SQRMPP)" << std::endl;
                return x=-1;
            }
        }
            //linear version
        if (k < 3 ) return sqrootlinear (x, a, p, k);
        else{
                //quadratic version
            uint64_t kdivtwo(k>>1);
            if ((k & 1) == 1){ // kdivtwo = (k-1)/2
                Rep sqpkdivp; pow(sqpkdivp,p,kdivtwo);

                    //x1^2 = a mod (p^((k-1)/2))
                sqrootmodprimepower (x, a, p, kdivtwo, sqpkdivp);
                if (x == -1) return x;

                    //x0^2 = a mod (p^(k-1))
                sqroothensellift (x, a, p, kdivtwo, sqpkdivp);
                if (x == -1) return x;

                Rep pkdivp (pk); pkdivp /= p;
                    //x2^2 = a mod (p^k)
                return sqrootonemorelift (x, a, p, k-1, ((pkdivp)));
            } else { // kdivtwo = k/2
                Rep sqpk; pow(sqpk,p,kdivtwo);

                    //x1^2 = a mod (p^(k/2))
                sqrootmodprimepower(x, a, p, kdivtwo, sqpk);
                if (x == -1) return x = -1;

                    //x0^2 = a mod (p^k)
                return sqroothensellift (x, a, p, kdivtwo, sqpk);
            }
        }
        return x;
    }

    template <class MyRandIter> inline typename IntSqrtModDom<MyRandIter>::Rep &
    IntSqrtModDom<MyRandIter>::sqrootmodpoweroftwo (Rep & x,
                                                    const Rep & a,
                                                    const uint64_t k,
                                                    const Rep & pk) const {
        Rep tmpa (a); Integer::modin(tmpa,pk);
        x = 0;
            //first cases k = 1,2,3
        if (k == 1) return x = tmpa;
        if (k == 2) {
            if (tmpa == 0) return x = 0;
            if (tmpa == 1) return x = 1;
            else {
                std::cerr << tmpa << "is not a quadratic residue mod " << pk << " (SQRM2k)" << std::endl;
                return x = -1;
            }
        }
        if (k == 3) {
            if (tmpa == 0) return x = 0;
            if (tmpa == 1) return x = 1;
            if (tmpa == 4) return x = 2;
            else{
                std::cerr << tmpa << " is not a quadratic residue mod " << pk << " (SQRM2k, k=3)" << std::endl;
                return x = -1;
            }
        }
            // General case k >= 4
        if(tmpa==0) return x=0;
        if(tmpa==1) return x=1;
        if ((tmpa & 1U)==0){

            Rep b(tmpa);
            uint64_t t=0;
            for( ; (b & 1U) == 0; ++t) b>>=1; // a = b p^t and p does not divide b

            if ((t & 1U)==0) {
                Rep sqrtpt(1); sqrtpt<<=(t>>1);
                sqrootmodpoweroftwo(x,b,k,pk);
                x <<= (t>>1); // x <-- x * 2^{t/2}
                return x%=pk;
            } else {
                std::cerr << tmpa  << "is not a quadratic residue mod " << pk << " (SQRM2k, k>=4)" << std::endl;
                return x=-1;
            }
        }


            //linear version
        if (k < 29) return sqroottwolinear (x, a, k);
        else {

            Rep un (1);
            uint64_t kdivtwoplusone(k);
            kdivtwoplusone >>= 1; ++kdivtwoplusone;
                // is k/2+1 if k is even, (k-1)/2+1 otherwise

            Rep pkmulttwo (pk); pkmulttwo <<= 1;
            Rep pkdivtwo (pk); pkdivtwo >>= 1;

            if ((k & 1) == 0){
                    //if k is even
                Rep sqrt_pk_mult_two (2); sqrt_pk_mult_two <<= (k>>1);
                    //x0^2=a mod (2^{k/2+1})
                sqrootmodpoweroftwo (x, tmpa, kdivtwoplusone, (sqrt_pk_mult_two));
                if (x == -1) return x;

                    //x1^2=a mod (2^k)
                return sqrootmodtwolift (x, tmpa, kdivtwoplusone, (sqrt_pk_mult_two));
            } else {
                    //if k is odd
                Rep sqrt_pkdivtwo_mult_two (2); sqrt_pkdivtwo_mult_two <<= (k>>1);
                    //x0^2=a mod (2^{k/2+1})
                sqrootmodpoweroftwo (x, tmpa,kdivtwoplusone, (sqrt_pkdivtwo_mult_two));
                if (x == -1) return x;

                    //x1^2=a mod (p^{k-1})
                sqrootmodtwolift (x, tmpa, kdivtwoplusone, (sqrt_pkdivtwo_mult_two));
                if (x == -1) return x;

                Rep u(tmpa);
                Integer::maxpyin(u,x,x); u %= pk;
                    //if x is a square root of a mod p^k
                if (u == 0) return x;

                    //if x is not square root of a mod p^k
                    //x + (p^{k-2}) is a square root of a mod p^k
                return x += pk>>2;
            }
        }
        return x;
    }

    template <class MyRandIter> inline typename IntSqrtModDom<MyRandIter>::Rep &
    IntSqrtModDom<MyRandIter>::sqrootlinear (Rep & x,
                                             const Rep & a,
                                             const Rep & p,
                                             const uint64_t k) const {
        sqrootmodprime(x,a,p);
        Rep pk(p);
        for(uint64_t i=1;i<k;i++){
            sqrootonemorelift(x,a,p,i,pk);
            pk *= p;
        }
        return x;
    }

    template <class MyRandIter> inline typename IntSqrtModDom<MyRandIter>::Rep &
    IntSqrtModDom<MyRandIter>::sqroottwolinear (Rep & x,
                                                const Rep & a,
                                                const uint64_t k) const {
            //first cases k = 1,2,3
        sqrootmodpoweroftwo(x, a, 3, 8);
        if (x == -1 || k<4) return x;

        Rep pk(16);
        Rep pk2(4);
        for(uint64_t i=4;i<=k;i++){
            if(((x*x)%pk)!=(a%pk)){
                x+=pk2;
            }
            pk2=pk;
            pk2>>=1;
            pk<<=1;
        }
        return x;
    }


    template <class MyRandIter> inline typename IntSqrtModDom<MyRandIter>::Rep &
    IntSqrtModDom<MyRandIter>::sqroothensellift (Rep & x,
                                                 const Rep & a,
                                                 const Rep & p,
                                                 const uint64_t k,
                                                 const Rep & pk) const {
            //we have a square root of a mod p^k : x0
            //x = x0 + h*p^k mod p^{2k}
            //with h = ((((a-x0^2) mod p^{2k})/p^k)*(2x0)^{-1} mod p^k) mod p^(2k)
            //is a square root of a mod p^{2k}
        Rep u(a);
        Integer::maxpyin(u,x,x);
        if(u == 0) return x;

        u /= pk;
            //    u %= pk;
            //u=(a-x0^2)/p^k

        Rep h(x<<1);
        this->invin (h, pk);
        h *= u; h %= pk;
            // h = (a-x0^2)/(2*x0*p^k) modulo pk

        return Integer::axpyin(x,h,pk);
    }

    template <class MyRandIter> inline typename IntSqrtModDom<MyRandIter>::Rep &
    IntSqrtModDom<MyRandIter>::sqrootonemorelift (Rep & x0,
                                                  const Rep & a,
                                                  const Rep & p,
                                                  const uint64_t k,
                                                  const Rep & pk) const {
        Rep u(a);
        Integer::maxpyin(u,x0,x0);
        u /= pk; u %= p;
        if (u == 0) return x0;

            //u=(a-x0^2)/p^k

        Rep h(x0<<1);
        this->invin (h, p);
        h *= u; h %= p;
            // h = (a-x0^2)/(2*x0*p^k) modulo p

        return Integer::axpyin(x0,h,pk);
    }

    template <class MyRandIter> inline typename IntSqrtModDom<MyRandIter>::Rep &
    IntSqrtModDom<MyRandIter>::sqrootmodtwolift (Rep & x,
                                                 const Rep & a,
                                                 const uint64_t k,
                                                 const Rep & pk) const {
            //we have a square root of a mod 2^k : x0 and we have
            //x = x0 + h*2^{k-1}
            //with h = ((((a-x0^2)mod 2^{2k-2})/2^k)*x0^{-1}mod 2^{k-1}) mod 2^{k-1}
            //is a square root of a mod 2^{2k-2}
        Rep u(a);
        Integer::maxpyin(u,x,x);
        u /= pk;

        Rep pk1(pk); pk1 >>= 1;
        u %= pk1;
        if (u == 0) return x;

        Rep h(x);
        invin(h,pk1);
        h *= u; h %= pk1;

        return Integer::axpyin(x,h,pk1);
    }



// =================================================================== //
// Brillhart decomposition as a sum of squares
// =================================================================== //
    template <class MyRandIter> inline void
    IntSqrtModDom<MyRandIter>::Brillhart(Rep& a, Rep& b, const Rep& p) const {
        GIVARO_REQUIRE(this->isprime(p),"isprime");
        GIVARO_STATE(Rep t);
        GIVARO_REQUIRE(this->isOne(this->mod(t,p,4U)),"prime = 1 mod 4");

        Integer x,s,q,r(p+this->mOne);
        sqrootmodprime(x,r,p);

        Givaro::sqrt(s,p);
        r=p;
        b=x>(p>>1)?p-x:x;

        Integer::mod(a,r,b);
        if (! this->isOne(a)) {
            while(a>s) {
                r=b;
                b=a;
                Integer::mod(a,r,b);
            }

            r=b;
            b=a;
            Integer::mod(a,r,b);
        } // otherwise x^2+1 = p is already correct in a and b

        GIVARO_ENSURE(this->isZero(a*a+b*b-p),"prime as sum of squares");
    }



// =================================================================== //
// Modular decomposition as a sum of squares
// =================================================================== //


        // Fast under ERH
    template <class MyRandIter> inline void
    IntSqrtModDom<MyRandIter>::sumofsquaresmodprime
    (Rep& a, Rep& b, const Rep& k, const Rep& p) const {
        sumofsquaresmodprimeDeterministic(a,b,k,p);
    }

	// Warning:
    //  the least quadratic non-residue is deterministic,
    //  but sub calls may not be deterministic.
    //  for instance sqrootmodprime is Las Vegas
    template <class MyRandIter> inline void
    IntSqrtModDom<MyRandIter>::sumofsquaresmodprimeDeterministic
    (Rep& a, Rep& b, const Rep& k, const Rep& p) const {
        GIVARO_REQUIRE(this->isprime(p),"isprime");

        Integer r(k);
        Integer::modin(r,p);
        if (this->isZero(r)) {
            a=this->zero;
            b=this->zero;
        } else {
            if (this->isOne(legendre(r,p))) {
                b=this->zero;
                this->sqrootmodprime(a,r,p);
            } else {
                Integer s(r); --s; // r-1 is a square ?
                if (this->isOne(legendre(s,p))) {
                        // then 1+(sqrt(s))^2=r
                    a=this->one;
                    this->sqrootmodprime(b,s,p);
                } else {
                    Integer lsnqr(2);
                        // Under ERH, least quad. non-residue
                        // should be lower than 3/2log^2(p)
                        // [Th. 6.35, Primality Tests on Commutator Curves,
                        //  U. Tubingen PhD 2001, Sebastian Wedeniwski]
                    for( ; legendre(lsnqr,p) != -1; ++lsnqr);

                    sumofsquaresmodprimewithnonresidue(a,b,r,lsnqr,p);
                }
            }
        }

        GIVARO_ENSURE(this->isZero( (a*a+b*b-k)%p ),"modular sum of squares");
    }

    template <class MyRandIter> inline void
    IntSqrtModDom<MyRandIter>::sumofsquaresmodprimeMonteCarlo
    (Rep& a, Rep& b, const Rep& k, const Rep& p) const {
        GIVARO_REQUIRE(this->isprime(p),"isprime");

        Integer r(k);
        Integer::modin(r,p);
        if (this->isZero(r)) {
            a=this->zero;
            b=this->zero;
        } else {
            if (this->isOne(legendre(r,p))) {
                b=this->zero;
                this->sqrootmodprime(a,r,p);
            } else {
                Integer s(r); --s; // r-1 is a square ?
                if (this->isOne(legendre(s,p))) {
                        // then 1+(sqrt(s))^2=r
                    a=this->one;
                    this->sqrootmodprime(b,s,p);
                } else {
                    Integer s,t;
                    while (
                        legendre (Integer::nonzerorandom (s, p.bitsize()), p)
                        != -1) {};
                        // Now s is not a residue
                    t=s;
                    for(--t ; legendre(t,p) == -1; --t);
                        // Now t is a quadratic residue and t+1 is not
                    sumofsquaresmodprimewithnonresidue(a,b,r,++t,p);
                }
            }
        }

        GIVARO_ENSURE(this->isZero( (a*a+b*b-k)%p ),"modular sum of squares");
    }

    template <class MyRandIter> inline void
    IntSqrtModDom<MyRandIter>::sumofsquaresmodprimewithnonresidue
    (Rep& a, Rep& b, const Rep& k, const Rep& s, const Rep& p) const {
        GIVARO_REQUIRE(this->isprime(p),"isprime");
        GIVARO_REQUIRE(legendre(s,p) == -1, "non-residue");
        GIVARO_REQUIRE(legendre(s-1,p) == 1, "quadratic residue");


        this->sqrootmodprime(b,s-1,p); // s = 1+b^2

        Integer r(k), il; Givaro::inv(il, s, p);
        r *= il;
        Integer::modin(r,p); // r/s mod p

        this->sqrootmodprime(a,r,p); // k/s = a^2

            // Now k = a^2(1+b^2)
        b *= a;
        b %= p;

        GIVARO_ENSURE(this->isZero( (a*a+b*b-k)%p ),"modular sum of squares");
    }

        // Unconditonal
    template <class MyRandIter> inline void
    IntSqrtModDom<MyRandIter>::sumofsquaresmodprimeNoERH
    (Rep& a, Rep& b, const Rep& k, const Rep& p) const {
        GIVARO_REQUIRE(this->isprime(p),"isprime");

        Integer r(k);
        Integer::modin(r,p);
        if (this->isZero(r)) {
            a=this->zero;
            b=this->zero;
        } else {
            if (this->isOne(legendre(r,p))) {
                b=this->zero;
                sqrootmodprime(a,r,p);
            } else {

            Integer t(1),h(p);
            h <<= 2U;	// h is 4p

            this->mod(r,p,4U);
            if (r == 1U)
                r -= (k%4); // Warning k can be negative
            else
                r += (k%4); // Warning k can be negative
            r *= p;
            r += k;
                // now r is k mod p and 1 mod 4

            while(! this->isprime(r)) {
                r += h;
            }   // r is prime and still k mod p and 1 mod 4

            Brillhart(a,b,r);
        }

        this->modin(a,p);
        this->modin(b,p);
        Integer half(p>>1);
        if (a>half) this->sub(a,p,a);
        if (b>half) this->sub(b,p,b);
        }

        GIVARO_ENSURE(this->isZero( (a*a+b*b-k)%p ),"modular sum of squares");
    }




} // namespace Givaro

#endif // __GIVARO_sqrootmod_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
