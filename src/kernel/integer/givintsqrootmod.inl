// ============================================================= //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.
// Time-stamp: <17 Sep 09 17:09:44 Jean-Guillaume.Dumas@imag.fr> 
// Givaro : Modular square roots
// Author : Yanis Linge
// ============================================================= //

#include <givaro/givtimer.h>

template <class RandIter> inline typename IntSqrtModDom<RandIter>::Rep & 
IntSqrtModDom<RandIter>::sqrootmodprime (Rep & x, 
                                         const Rep & a, 
                                         const Rep & p) const {
    Integer res;
    Rep amp (a);
    amp %=p;
    if (amp == 0 || amp == 1) return x = amp;

    if (legendre (amp, p) == -1){
        std::cerr << amp << " is not a quadratic residue mod " << p << std::endl;
        return x = -1;
    }

	// p = 3 mod 4: x=a^{(p+1)/4 } mod p
    if ((p & 3UL) == 3UL) {
        Rep ppu (p);
        ++ppu;
        ppu >>= 2;			// ppu = (p+1)/4
        return powmod (x, amp, ppu, p);	// powmod (x,a,(p+1)/4,p);
    }

    size_t l = (size_t) ceil (logtwo (p) - 1);

    if ((p & 7UL) == 5UL) {
        Rep tmp;
        Rep puis (p);
        puis -= 1;
        puis >>= 2;

        powmod (tmp, amp, puis, p);
        if (tmp == 1UL) {
            puis = p;
            puis += 3;
            puis >>= 3;
            return powmod (x, amp, puis, p);
        } 
        puis = p;
        puis -= 5;
        puis >>= 3;
        Rep a4 (amp);
        a4 *= 4;
        powmod (x, a4, puis, p);
        x *= amp;
        x <<= 1;
        return x %= p;
    }

	// H. Cohen version
    Rep p1 (p);
    --p1;
    Rep q (p1);
    Rep e (0);
    while ((q & 1UL) == 0){
        q >>= 1;
        ++e;
    }
	// now we have e and q such like : p-1=q*2^e with q odd
	// we need a non quadratic element : g 
    Rep g; while (legendre (Rep::random (g, l), p) != -1);

    Rep z;
    powmod (z, g, q, p);
	//Initialize
    Rep y (z);
    Rep r (e);
    Rep tmp (q);
    tmp -= 1;
    tmp >>= 1;
    powmod (x, amp, tmp, p);
    Rep b (x);
    b *= x;
    b *= amp;
    b %= p;
    x *= amp;
    x %= p;

    Rep m;
    Rep twopuism;
    Rep bpuis2puism;
    Rep t;
    Rep puis (r);
    while ((b % p) != 1){
        m = 1;
        twopuism = 2;
        Rep bpuis2puism;
        powmod (bpuis2puism, b, twopuism, p);
        while (bpuis2puism != 1){
            ++m;
            pow (twopuism, 2, m);
            powmod (bpuis2puism, b, twopuism, p);
        }
        if (m == r){
            x = -1;
            std::cerr << amp << " is not a quadratic residu mod " << p << std::endl;
        }
        puis = r;
        puis -= m;
        --puis;
        pow (puis, 2, puis);
        powmod (t, y, puis, p);
        powmod (y, t, 2, p);
        r = m;
        x *= t;
        x %= p;
        b *= y;
        b %= p;
    }
    return x;
}

template <class RandIter> inline typename IntSqrtModDom<RandIter>::Rep &
IntSqrtModDom<RandIter>::sqrootmodprimepower (Rep & x, 
                                              const Rep & a, 
                                              const Rep & p, 
                                              const unsigned long k, 
                                              const Rep & pk) const{
    Rep tmpa(a);
    tmpa%=pk;
    if(tmpa==0) return x=0;
    if(tmpa==1) return x=1;
    if (k == 1) return sqrootmodprime (x, tmpa, p);

    if ((tmpa%p)==0){
        Rep b(tmpa);
        unsigned long t=0;
        while((b%p) == 0){
            b/=p;
            ++t;
        }
        if((t%2)==0){
            Rep sqrtb;
            sqrootmodprimepower(sqrtb,b,p,k,pk);
            Rep sqrtpt;
            powmod(sqrtpt,p,(t/2),pk);
            x=sqrtpt;
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
        Rep x0;
        Rep x1 (0);
        unsigned long kdivtwo = k;
        kdivtwo >>= 1;
        Rep pkdivp (pk);
        pkdivp /= p;
        unsigned long kmone = k;
        --kmone;
        if ((k & 1) == 1){
                //x1^2 = a mod (p^((k-1)/2))
            sqrootmodprimepower (x1, a, p, kdivtwo, sqrt (pkdivp));	
            if (x1 == -1) return x = -1;

                //x0^2 = a mod (p^(k-1))
            sqroothensellift (x0, x1, a, p, kdivtwo, sqrt (pkdivp));
            if (x0 == -1) return x = -1;

                //x2^2 = a mod (p^k)
            return sqrootonemorelift (x, x0, a, p, kmone, ((pkdivp)));	
        } else {
            	//x1^2 = a mod (p^(k/2))
            sqrootmodprimepower(x1, a, p, kdivtwo, sqrt (pk));	
            if (x1 == -1) return x = -1;
            
                //x0^2 = a mod (p^k)
            return sqroothensellift (x, x1, a, p, kdivtwo, sqrt (pk));
        }
    }
    return x;
}

template <class RandIter> inline typename IntSqrtModDom<RandIter>::Rep &
IntSqrtModDom<RandIter>::sqrootmodpoweroftwo (Rep & x, 
                                              const Rep & a, 
                                              const unsigned long k,
                                              const Rep & pk) const {
    Rep tmpa (a);
    tmpa %= pk;
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
            std::cerr << tmpa << " is not a quadratic residu mod " << pk << " cas k = 3" << std::endl;
            return x = -1;
        }
    }
    if(tmpa==0) return x=0;
    if(tmpa==1) return x=1;
    if ((tmpa & 1UL)==0){

        Rep b(tmpa);
        unsigned long t=0;
        while ((b &1UL)==0) {
            b>>=1;
            ++t;
        }
        if ((t & 1UL)==0) {
                //powmod(sqrtpt,p,(t/2),pk);
            Rep sqrtpt(1); sqrtpt<<=(t/2);
            sqrootmodpoweroftwo(x,b,k,pk);
            x*=sqrtpt;
            return x%=pk;
        } else {
            std::cerr << tmpa  << "is not a quadratic residu mod " << pk << std::endl;
            return x=-1;
        }
    }


        //linear version
    if (k < 29) return sqroottwolinear (x, a, k);
    else {
        if ((tmpa & 1UL) == 0){
                //if a is even we search b odd and t such that a=b*2^t
            Rep b (tmpa);
            unsigned long t = 0;
            while ((b & 1UL) == 0){
                b >>= 1;
                ++t;
            }
            if (((t & 1UL) == 1) || ((b & 7UL) != 1)){
                    //if t is odd
                std::cerr << tmpa << " is not a quadratic residu mod " << pk << " cas t impaire" << std::endl;
                return x = -1;
            }
            else{
                // u^2 = b mod p^k 
                sqrootmodpoweroftwo(x, b, k, pk);	
                unsigned long tdivtwo = t;
                tdivtwo >>= 1;
                    // x = u*2^(t/2)
                return x <<= tdivtwo;
            }
        }
        else{
            Rep x1;
            Rep x0;
            Rep un (1);
            unsigned long kdivtwoplusone = k;
            kdivtwoplusone >>= 1;
            ++kdivtwoplusone;
            Rep pkmulttwo (pk);
            pkmulttwo <<= 1;
            Rep pkdivtwo (pk);
            pkdivtwo >>= 1;
            if ((k & 1) == 0){
                    //if k is even 
                Rep sqrt_pk_mult_two (1);
                sqrt_pk_mult_two = sqrt (pk);
                sqrt_pk_mult_two <<= 1;
                    //x0^2=a mod (2^{k/2+1})
                sqrootmodpoweroftwo (x0, tmpa, kdivtwoplusone, (sqrt_pk_mult_two));
                if (x0 == -1) return x = -1;

                    //x1^2=a mod (2^k)
                return sqrootmodtwolift (x, x0, tmpa, kdivtwoplusone, (sqrt_pk_mult_two));
            } else {
                    //if k is odd
                Rep sqrt_pkdivtwo_mult_two (1);
                sqrt_pkdivtwo_mult_two = sqrt (pkdivtwo);
                sqrt_pkdivtwo_mult_two <<= 1;
                    //x0^2=a mod (2^{k/2+1})
                sqrootmodpoweroftwo (x0, tmpa,kdivtwoplusone, (sqrt_pkdivtwo_mult_two));
                if (x0 == -1) return x = -1;
                
                    //x1^2=a mod (p^{k-1})
                sqrootmodtwolift (x1, x0, tmpa, kdivtwoplusone, (sqrt_pkdivtwo_mult_two));
                if (x1 == -1) return x = -1;

                if ((x1 * x1) % pk == tmpa % pk){
                        //if x1 is a square root of a mod p^k
                    return x = x1;
                }
                    //if x is not square root of a mod p^k
                    //x1 + (p^{k-2}) is a square root of a mod p^k
                x = pk;
                x >>= 2;
                return x += x1;
            }
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
    Rep x1;
    for(unsigned long i=1;i<k;i++){	
        sqrootonemorelift(x1,x,a,p,i,pk);
        x=x1;
        pk*=p;
    }
    return x;
}    

template <class RandIter> inline typename IntSqrtModDom<RandIter>::Rep &
IntSqrtModDom<RandIter>::sqrootquad (Rep & x, 
                                     const Rep & a, 
                                     const Rep & p, 
                                     const unsigned long k, 
                                     const Rep & pk) const {
    sqrootmodprime (x, a, p);
    Rep pn (p);
    unsigned long ktmp = 1;
    Rep x1;
    for (unsigned long i = 1; pn <= pk; i++){
        sqroothensellift (x1, x, a, p, ktmp, pn);
        ktmp <<= 1;
        pn *= pn;
        x = x1;
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
IntSqrtModDom<RandIter>::sqroottwoquad (Rep & x, 
                                        const Rep & a,
                                        const unsigned long k, 
                                        const Rep & pk) const {
    sqrootmodpoweroftwo (x, a, 3, 8);

    Rep pn (8);
    Rep x1 (x);
    for (unsigned long i = 4; pn <= pk; i++){
        sqrootmodtwolift (x1, x, a, i, pn);
        x = x1;
        pn *= pn;
        pn >>= 2;
    }

    return x;
}

template <class RandIter> inline typename IntSqrtModDom<RandIter>::Rep &
IntSqrtModDom<RandIter>::sqroothensellift (Rep & x, 
                                           const Rep & x0, 
                                           const Rep & a,
                                           const Rep & p, 
                                           const unsigned long k, 
                                           const Rep & pk) const {
//we have a square root of a mod p^k : x0 
//x = x0 + h*p^k mod p^{2k} 
//with h = ((((a-x0^2) mod p^{2k})/p^k)*(2x0)^{-1} mod p^k) mod p^(2k)
//is a square root of a mod p^{2k}

    Rep rtmp (x0);
    x = 0;
    Rep p2k (pk);
    p2k *= pk;
    rtmp %= pk;
    Rep tmp (rtmp);
    tmp *= tmp;
    tmp -= a;
    tmp %= p2k;

	//if tmp = 0 x0 is already a square root of a mod p^{2k}
    if (tmp == 0) return x = rtmp;

    Rep h (0);
    inv (h, 2 * rtmp, pk);
    h *= -tmp;
    h /= pk;
    h %= p2k;
    x += h;
    x *= pk;
    x += rtmp;
    return x %= p2k;
}

template <class RandIter> inline typename IntSqrtModDom<RandIter>::Rep &
IntSqrtModDom<RandIter>::sqrootonemorelift (Rep & x, 
                                            const Rep & x0, 
                                            const Rep & a, 
                                            const Rep & p, 
                                            const unsigned long k, 
                                            const Rep & pk) const {
    Rep p2 (p);
    p2 *= p;
    Rep pkpone (pk);
    pkpone *= p;
    Rep x0multtwo (2);
    x0multtwo *= x0;

//p2=p^2 pkpone = p^{k+1} and x0multtwo = 2x0

    Rep u (x0);
    u *= x0;
    u -= a;
    if (u == 0) return x = x0;

    u /= pk;
    u %= p;
//u=(x0^2-a)/p^k

    Rep h;
    inv (h, x0multtwo, p);
    h *= -u;
// h = (x0^2-a)/(2*x0*p^k)

    x = h;
    x *= pk;
    x += x0;
    return x %= pkpone;
}

template <class RandIter> inline typename IntSqrtModDom<RandIter>::Rep &
IntSqrtModDom<RandIter>::sqrootmodtwolift (Rep & x, 
                                           const Rep & x0, 
                                           const Rep & a,
                                           const unsigned long k, 
                                           const Rep & pk) const {
//we have a square root of a mod 2^k : x0 and we have
//x = x0 + h*2^{k-1}
//with h = ((((a-x0^2)mod 2^{2k-2})/2^k)*x0^{-1}mod 2^{k-1}) mod 2^{k-1}
//is a square root of a mod 2^{2k-2}

    Rep rtmp (x0);
    x = 0;
    Rep p2k2 (pk);
    p2k2 *= pk;
    p2k2 <<= 2;
    rtmp %= pk;
    Rep tmp (rtmp);
    tmp *= tmp;
    tmp -= a;
    tmp %= p2k2;

    if (tmp == 0) return x = rtmp;

    Rep h (0);
    inv (h, rtmp, pk / 2);
    h *= -tmp;
    h /= pk;
    h %= (pk / 2);
    x += h;
    x *= (pk / 2);
    x += rtmp;
    return x %= p2k2;
}
