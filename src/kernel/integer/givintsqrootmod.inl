// ============================================================= //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <20 Jan 11 11:02:30 Jean-Guillaume.Dumas@imag.fr>
// Givaro : Modular square roots
// Author : Yanis Linge
// ============================================================= //

#ifndef __GIVARO_sqrootmod_INL
#define __GIVARO_sqrootmod_INL

#include <givaro/givtimer.h>

template <class RandIter> inline typename IntSqrtModDom<RandIter>::Rep &
IntSqrtModDom<RandIter>::sqrootmodprime (Rep & x,
                                         const Rep & a,
                                         const Rep & p) const {
//         std::cerr << "p:= " << p << ';' << std::endl;
//         std::cerr << "a:= " << a << ';' << std::endl;
    Rep amp (a); amp %=p;
    if (amp == 0UL || amp == 1UL) return x = amp;

    if (legendre (amp, p) == -1){
        std::cerr << amp << " is not a quadratic residue mod " << p << std::endl;
        return x = -1;
    }

    if ((p & 3UL) == 3UL) {			// If p = 3 mod 4
        Rep ppu (p); ++ppu; ppu >>= 2UL;	// ppu = (p+1)/4
        return powmod (x, amp, ppu, p);		// powmod (x,a,(p+1)/4,p);
    }

        // O. Atkin
    if ((p & 7UL) == 5UL) {			// If p = 5 mod 8
        Rep tmp;
        Rep puis (p); puis -= 1UL; puis >>= 2UL;// puis = (p-1)/4
        powmod (tmp, amp, puis, p);

        if (tmp == 1UL) {
            puis = p; puis += 3UL; puis >>= 3UL;// puis = (p+3)/8
            return powmod (x, amp, puis, p);
        }
        puis = p; puis -= 5UL; puis >>= 3UL; 	// puis = (p-5)/8

        Rep a4 (amp); a4 <<= 2;
        powmod (x, a4, puis, p);
        x *= amp; x <<= 1; 			// 2a(4a)^{(p-5)/8}
        return x %= p;
    }

    size_t l = (size_t) ceil (logtwo (p) - 1);

        // S. Mueller
    if ((p & 15UL) == 9UL) {			// If p = 9 mod 16
        Rep i (amp); i <<= 1UL;
        Rep puis (p); puis -= 1UL; puis >>= 2UL;// puis = (p-1)/4
        powmod (x, i, puis, p);			// (2a)^{(p-1}/4} is +1 or -1
 	if (x != 1UL) x = -1L;

        Rep d; while (legendre (Rep::nonzerorandom (d, l), p) == x) ;
        puis = p; puis -= 9UL; puis >>= 4UL; 	// puis = (p-9)/16
        i *= d; i *= d;
        powmod(x, i, puis, p);			// (2d^2a)^{(p-9)/16}
        i *= x; i *= x; --i;			// i=(2d^2x^2a-1) ; i^2 = -1

        x *= d; x *= amp; x *= i;		// xda(i-1)
        return x %= p;				// +/- x is a root
    }

        // Tonelli and Shanks
	// [H. Cohen, Algorithm 1.5.1, p33,
        // A course in computational algebraic number theory]
    Rep p1 (p); --p1;
    Rep q (p1);
    long e (0);
    for( ; (q & 1UL) == 0; ++e) q >>= 1;

	// now we have e and q such that : p-1=q*2^e with q odd
	// we need a non quadratic element : g
    Rep g; while (legendre (Rep::nonzerorandom (g, l), p) != -1) ;

    Rep z;
    powmod (z, g, q, p);	// z = g^q mod p
	//Initialize
    Rep y (z);
    Rep tmp (q); tmp -= 1; tmp >>= 1;
    powmod (x, amp, tmp, p); 	// a^{(q-1)/2} mod p
    Rep b (x);
    b *= x; b *= amp; b %= p;	// ax^2
    x *= amp; x %= p;		// ax

        // Find exponent
    long m(1), r(e);
    Rep b2k, t, puis(r);
    while (b != 1){
        b2k = b;
        for(m = 0; b2k != 1; ++m) {
            b2k *= b2k; b2k %= p;
        } // m smallest such that b^{2^m} is 1 mod p
        if (m == r){
            std::cerr << amp << " is not a quadratic residu mod " << p << std::endl;
            return x = -1;
        }
        long lpuis = r; lpuis -= m; --lpuis;
        puis = 1; puis <<= lpuis;	// 2^{m-r-1}
        powmod (t, y, puis, p);		// t = y^{ 2^{m-r-1} } mod p
        y = t; y *= t; y %= p;		// y = t^2 mod p
        r = m;				// r = m
        x *= t; x %= p;			// x = xt mod p
        b *= y; b %= p;			// b = by mod p
    }
    return x;
}

template <class RandIter> inline typename IntSqrtModDom<RandIter>::Rep &
IntSqrtModDom<RandIter>::sqrootmodprimepower (Rep & x,
                                              const Rep & a,
                                              const Rep & p,
                                              const unsigned long k,
                                              const Rep & pk) const{

    Rep tmpa(a); tmpa%=pk;
    if(tmpa==0) return x=0;
    if(tmpa==1) return x=1;
    if (k == 1) return sqrootmodprime (x, tmpa, p);

    if ((tmpa%p)==0){
        Rep b(tmpa);
        unsigned long t=0;
        for( ; (b%p) == 0; ++t) b/=p; // a = b p^t and p does not divide b

        if((t&1UL)==0){
            Rep sqrtb;
            sqrootmodprimepower(sqrtb,b,p,k,pk);
            powmod(x,p,(t>>1),pk);
            x*=sqrtb;
            return x%=pk;
        }
        else{
            std::cerr <<tmpa << "is not a quadratic residu mod " << pk << std::endl;
            return x=-1;
        }
    }
	//linear version
    if (k < 3 ) return sqrootlinear (x, a, p, k);
    else{
	//quadratic version
        unsigned long kdivtwo(k>>1);
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

template <class RandIter> inline typename IntSqrtModDom<RandIter>::Rep &
IntSqrtModDom<RandIter>::sqrootmodpoweroftwo (Rep & x,
                                              const Rep & a,
                                              const unsigned long k,
                                              const Rep & pk) const {
    Rep tmpa (a); tmpa %= pk;
    x = 0;
        //first cases k = 1,2,3
    if (k == 1) return x = tmpa;
    if (k == 2) {
        if (tmpa == 0) return x = 0;
        if (tmpa == 1) return x = 1;
        else {
            std::cerr << tmpa << "is not a quadratic residu mod " << pk << std::endl;
            return x = -1;
        }
    }
    if (k == 3) {
        if (tmpa == 0) return x = 0;
        if (tmpa == 1) return x = 1;
        if (tmpa == 4) return x = 2;
        else{
            std::cerr << tmpa << " is not a quadratic residu mod " << pk << " (case k = 3)" << std::endl;
            return x = -1;
        }
    }
        // General case k >= 4
    if(tmpa==0) return x=0;
    if(tmpa==1) return x=1;
    if ((tmpa & 1UL)==0){

        Rep b(tmpa);
        unsigned long t=0;
        for( ; (b & 1UL) == 0; ++t) b>>=1; // a = b p^t and p does not divide b

        if ((t & 1UL)==0) {
            Rep sqrtpt(1); sqrtpt<<=(t>>1);
            sqrootmodpoweroftwo(x,b,k,pk);
            x <<= (t>>1); // x <-- x * 2^{t/2}
            return x%=pk;
        } else {
            std::cerr << tmpa  << "is not a quadratic residu mod " << pk << std::endl;
            return x=-1;
        }
    }


        //linear version
    if (k < 29) return sqroottwolinear (x, a, k);
    else {

        Rep un (1);
        unsigned long kdivtwoplusone(k);
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

template <class RandIter> inline typename IntSqrtModDom<RandIter>::Rep &
IntSqrtModDom<RandIter>::sqrootlinear (Rep & x,
                                       const Rep & a,
                                       const Rep & p,
                                       const unsigned long k) const {
    sqrootmodprime(x,a,p);
    Rep pk(p);
    for(unsigned long i=1;i<k;i++){
        sqrootonemorelift(x,a,p,i,pk);
        pk *= p;
    }
    return x;
}

template <class RandIter> inline typename IntSqrtModDom<RandIter>::Rep &
IntSqrtModDom<RandIter>::sqroottwolinear (Rep & x,
                                          const Rep & a,
                                          const unsigned long k) const {
        //first cases k = 1,2,3
    sqrootmodpoweroftwo(x, a, 3, 8);
    if (x == -1 || k<4) return x;

    Rep pk(16);
    Rep pk2(4);
    for(unsigned long i=4;i<=k;i++){
        if(((x*x)%pk)!=(a%pk)){
            x+=pk2;
        }
        pk2=pk;
        pk2>>=1;
        pk<<=1;
    }
    return x;
}


template <class RandIter> inline typename IntSqrtModDom<RandIter>::Rep &
IntSqrtModDom<RandIter>::sqroothensellift (Rep & x,
                                           const Rep & a,
                                           const Rep & p,
                                           const unsigned long k,
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

template <class RandIter> inline typename IntSqrtModDom<RandIter>::Rep &
IntSqrtModDom<RandIter>::sqrootonemorelift (Rep & x0,
                                            const Rep & a,
                                            const Rep & p,
                                            const unsigned long k,
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

template <class RandIter> inline typename IntSqrtModDom<RandIter>::Rep &
IntSqrtModDom<RandIter>::sqrootmodtwolift (Rep & x,
                                           const Rep & a,
                                           const unsigned long k,
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

#endif // __GIVARO_sqrootmod_INL
