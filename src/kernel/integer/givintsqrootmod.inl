#include <givaro/givtimer.h>

template < class RandIter> 
inline typename  IntSqrtModDom<RandIter>::Rep & IntSqrtModDom<RandIter>::sqrootmodprime 
(Rep & x, const Rep & a, const Rep & p) const {

Integer res;
	Rep amp (a);
	amp %=p;
	if (amp == 0 || amp == 1){
		x = amp;
		return x;
	}

	if (legendre (amp, p) == -1){
		std::cerr << amp << " is not a quadratic residue mod " << p << std::endl;
		x = -1;
		return x;
		}

	if ((p & Rep (3)) == 3){

//if p = 3 mod 4 x=a^{(p+1)/4 } mod p
		Rep ppu (p);
		++ppu;
		ppu >>= 2;		// ppu = (p+1)/4

		powmod (x, amp, ppu, p);	// powmod (x,a,(p+1)/4,p);
		if (x < 0){
			x = -x;
		}
		return x;
	}

	size_t l = (size_t) ceil (logtwo (p) - 1);
    Rep::seeding (BaseTimer::seed ());
    Rep g;

	while (legendre (Rep::random (g, l), p) != -1);
//we need a non quadratic element : g 
    if ((p & Rep (7)) == 5){
		Rep tmp;
		Rep puis (p);
		puis -= 1;
		puis >>= 2;

		powmod (tmp, amp, puis, p);
		if (tmp == Rep (1)){
			puis = p;
			puis += 3;
			puis >>= 3;
			powmod (x, amp, puis, p);
			if (x < 0){
				x = -x;
			}
			return x;
		}
		else{
			puis = p;
			puis -= 5;
			puis >>= 3;
			Rep a4 (amp);
			a4 *= 4;
			powmod (x, a4, puis, p);
			x *= amp;
			x <<= 1;
			x %= p;
			if (x < 0){
				x = -x;
			}
			return x;
		}
	}

//Cohen version
	Rep p1 (p);
    --p1;
    Rep q (p1);
    Rep e (0);
    while ((q & Rep (1)) == 0){
		q >>= 1;
		++e;
	}
//we have e and q such like : p-1=q*2^e with q odd

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
    if (x < 0)
		x = -x;
		return x;
}

template < class RandIter> 
inline typename  IntSqrtModDom<RandIter>::Rep & IntSqrtModDom<RandIter>::sqrootmodprimepower (Rep & x, const Rep & a, const Rep & p, const unsigned long k, const Rep & pk) const{

	Rep tmpa(a);
	tmpa%=pk;
	if(tmpa==0){
		x=0;
		return x;
	}
    
	if (k == 1){
		sqrootmodprime (x, tmpa, p);
		return x;
	}


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
			x%=pk;
			if (x<0)
				x=-x;
			return x;
		}
		else{
			std::cerr <<tmpa << "is not a quadratic residu mod " << pk << std::endl;
			x=-1;
			return x;
		}
	}
	if (k == 1){
		sqrootmodprime (x, a, p);
		return x;
	}
    if (k < 3 ){
//linear version
		sqrootlinear (x, a, p, k);
		return x;
	}
    else{
//quad version
		Rep x0;
		Rep x1 (0);
		unsigned long kdivtwo = k;
		kdivtwo >>= 1;
		Rep pkdivp (pk);
		pkdivp /= p;
		unsigned long kmone = k;
		--kmone;
		if ((k & 1) == 1){
			sqrootmodprimepower (x1, a, p, kdivtwo, sqrt (pkdivp));	//x1^2 = a mod (p^((k-1)/2))
			if (x1 == -1){
				x = -1;
				return x;
			}
			sqroothensellift (x0, x1, a, p, kdivtwo, sqrt (pkdivp));	//x0^2 = a mod (p^(k-1))
			if (x0 == -1){
				x = -1;
				return x;
			}
			Rep x2;
			sqrootonemorelift (x2, x0, a, p, kmone, ((pkdivp)));	//x2^2 = a mod (p^k)
			if (x2 == -1){
				x = -1;
				return x;
			}
			if (x2 < 0)
				x2 = -x2;
			x = x2;
			return x;
		}
		else{
			sqrootmodprimepower(x1, a, p, kdivtwo, sqrt (pk));	//x1^2 = a mod (p^(k/2)1)
			if (x1 == -1){
				x = -1;
				return x;
			}
			sqroothensellift (x0, x1, a, p, kdivtwo, sqrt (pk));	//x0^2 = a mod (p^k)
			if (x0 == -1){
				x = -1;
				return x;
			}
			if (x0 < 0)
				x0 = -x0;
			x = x0;
			return x;
		}
	}
	return x;
}

template < class RandIter> 
inline typename  IntSqrtModDom<RandIter>::Rep & IntSqrtModDom<RandIter>::sqrootmodpoweroftwo (Rep & x, const Rep & a, const unsigned long k,const Rep & pk)const{

	Rep tmpa (a);
    tmpa %= pk;
    x = 0;
    //first cases k = 1,2,3
    if (k == 1){
		x = tmpa;
		return x;
	}
    if (k == 2){
		if (tmpa == 0){
			x = 0;
			return x;
	    }
		if (tmpa == 1){
			x = 1;
			return x;
	    }
		else{
			std::cerr << tmpa << "is not a quadratic residu mod " << pk << std::endl;
			x = -1;
			return x;
		}
	}
    if (k == 3){
		if (tmpa == 0){
			x = 0;
			return x;
		}
		if (tmpa == 1){
			x = 1;
			return x;
	    }
		if (tmpa == 0){
			x = 0;
			return x;
		}
		if (tmpa == 1){
			x = 1;
			return x;
	    }
		if (tmpa == 4){
			x = 2;
			return x;
	    }
		else{
			std::cerr << tmpa << " is not a quadratic residu mod " << pk << " cas k = 3" << std::endl;
			x = -1;
			return x;
		}
	}

	if ((tmpa&Rep(1))==0){
		Rep b(tmpa);
		unsigned long t=0;
		while((b&Rep(1)) == 0){
			b>>=1;
			++t;
		}
		if((t%2)==0){
			Rep sqrtb;
			sqrootmodpoweroftwo(sqrtb,b,k,pk);
			Rep sqrtpt(1);
			//powmod(sqrtpt,p,(t/2),pk);
			sqrtpt<<=(t/2);
			x=sqrtpt;
			x*=sqrtb;
			x%=pk;
			if (x<0)
				x=-x;
			return x;
		}
		else{
			std::cerr << tmpa << "is not a quadratic residu mod " << pk << std::endl;
			x=-1;
		    return x;
		}
	}


    if (k < 29){
	   //linear version
		sqroottwolinear (x, a, k);
		return x;
	}
    else{
		if ((tmpa & Rep (1)) == 0){
	       //if a is even we search b odd and t such that a=b*2^t
			Rep b (tmpa);
			unsigned long t = 0;
			while ((b & Rep (1)) == 0){
				b >>= 1;
				++t;
			}
			if (((t & 1) == 1) || ((b & Rep (7)) != 1)){
		   //if t is odd
				std::cerr << tmpa << " is not a quadratic residu mod " << pk << " cas t impaire" << std::endl;
				x = -1;
				return x;
			}
			else{
				Rep u (0);
				sqrootmodpoweroftwo(u, b, k, pk);	// u^2 = b mod p^k 
				x = u;
				unsigned long tdivtwo = t;
				tdivtwo >>= 1;
				x <<= tdivtwo;	// x = u*2^(t/2)
				if (x < 0)
					x = -x;
				return x;
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
				sqrootmodpoweroftwo (x0, tmpa, kdivtwoplusone, (sqrt_pk_mult_two));	//x0^2=a mod (2^{k/2+1})
				if (x0 == -1){
					x = -1;
					return x;
				}
				sqrootmodtwolift (x1, x0, tmpa, kdivtwoplusone, (sqrt_pk_mult_two));	//x1^2=a mod (2^k)
				if (x1 == -1){
					x = -1;
					return x;
				}
				if (x < 0)
					x = -x;
				x = x1;
				return x;
			}
			else{
		   //if k is odd
				Rep sqrt_pkdivtwo_mult_two (1);
				sqrt_pkdivtwo_mult_two = sqrt (pkdivtwo);
				sqrt_pkdivtwo_mult_two <<= 1;
				sqrootmodpoweroftwo (x0, tmpa,kdivtwoplusone, (sqrt_pkdivtwo_mult_two));	//x0^2=a mod (2^{k/2+1})
				if (x0 == -1){
					x = -1;
					return x;
				}
				sqrootmodtwolift (x1, x0, tmpa, kdivtwoplusone, (sqrt_pkdivtwo_mult_two));	//x1^2=a mod (p^{k-1})
				if (x1 == -1){
					x = -1;
			return x;
				}
				if ((x1 * x1) % pk == tmpa % pk){
		       //if x1 is a square root of a mod p^k
					x = x1;
					if (x < 0)
						x = -x;
					return x;
				}
		   //if x is not square root of a mod p^k
		   //x1 + (p^{k-2}) is a square root of a mod p^k
				x = pk;
				x >>= 2;
				x += x1;
				if (x < 0)
					x = -x;
				return x;
			}
		}
	}
    return x;
}

template < class RandIter> 
inline typename  IntSqrtModDom<RandIter>::Rep & IntSqrtModDom<RandIter>::sqrootlinear (Rep & x, const Rep & a, const Rep & p, const unsigned long k)const{

	sqrootmodprime(x,a,p);
	Rep pk(p);
	Rep x1;
	for(unsigned long i=1;i<k;i++){

					
			sqrootonemorelift(x1,x,a,p,i,pk);
			x=x1;
			if(x<0)
				x=-x;
		pk*=p;
	}
	return x;
}    

template < class RandIter> 
inline typename  IntSqrtModDom<RandIter>::Rep & IntSqrtModDom<RandIter>::sqrootquad (Rep & x, const Rep & a, const Rep & p, const unsigned long k, const Rep & pk) const{

	sqrootmodprime (x, a, p);
    Rep pn (p);
    unsigned long ktmp = 1;
    Rep x1;
    for (unsigned long i = 1; pn <= pk; i++){
		sqroothensellift (x1, x, a, p, ktmp, pn);
		ktmp <<= 1;
		pn *= pn;
		x = x1;
		if (x < 0)
			x = -x;
	}
    return x;
}

template < class RandIter> 
inline typename  IntSqrtModDom<RandIter>::Rep & IntSqrtModDom<RandIter>::sqroottwolinear (Rep & x, const Rep & a, const unsigned long k) const{

       //first cases k = 1,2,3

	sqrootmodprimepower(x, a, 3, 8, 1);
    if (x == -1)
		return x;

	if (k<4)
		return x;
	Rep pk(16);
	Rep pk2(4);
	//if(x==0){
	// }
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


template < class RandIter> 
inline typename  IntSqrtModDom<RandIter>::Rep & IntSqrtModDom<RandIter>::sqroottwoquad (Rep & x, const Rep & a, const unsigned long k, const Rep & pk) const{

	sqrootmodpoweroftwo (x, a, 3, 8);

     Rep pn (8);
     Rep x1 (x);
     for (unsigned long i = 4; pn <= pk; i++){
		sqrootmodtwolift (x1, x, a, i, pn);
		x = x1;
		if (x < 0)
			x = -x;
		pn *= pn;
		pn >>= 2;
	}

	return x;
}

template < class RandIter> 
inline typename  IntSqrtModDom<RandIter>::Rep & IntSqrtModDom<RandIter>::sqroothensellift (Rep & x1, const Rep & x0, const Rep & a,const Rep & p, const unsigned long k, const Rep & pk) const{

//we have a square root of a mod p^k : x0 
//x1 = x0 + h*p^k mod p^{2k} 
//with h = ((((a-x0^2) mod p^{2k})/p^k)*(2x0)^{-1} mod p^k) mod p^(2k)
//is a square root of a mod p^{2k}

	Rep rtmp (x0);
    x1 = 0;
    Rep p2k (pk);
    p2k *= pk;
    rtmp %= pk;
    Rep tmp (rtmp);
    tmp *= tmp;
    tmp -= a;
    tmp %= p2k;

    if (tmp == 0){
//if tmp = 0 x0 is already a square root of a mod p^{2k}
		x1 = rtmp;
		if (x1 < 0)
	    x1 = -x1;
		return x1;
	}

    Rep h (0);
    inv (h, 2 * rtmp, pk);
    h *= -tmp;
    h /= pk;
    h %= p2k;
    x1 += h;
    x1 *= pk;
    x1 += rtmp;
    x1 %= p2k;
    if (x1 < 0)
		x1 = -x1;
    return x1;

}

template < class RandIter> 
inline typename  IntSqrtModDom<RandIter>::Rep & IntSqrtModDom<RandIter>::sqrootonemorelift (Rep & x1, const Rep & x0, const Rep & a, const Rep & p, const unsigned long k, const Rep & pk) const{

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
    if (u == 0){
		x1 = x0;
		if (x1 < 0)
			x1 = -x1;
		return x1;
	}

	u /= pk;
    u %= p;
//u=(x0^2-a)/p^k

	Rep h;
    inv (h, x0multtwo, p);
    h *= -u;
// h = (x0^2-a)/(2*x0*p^k)

    x1 = h;
    x1 *= pk;
    x1 += x0;
    x1 %= pkpone;
    if (x1 < 0)
		x1 = -x1;

	return x1;
}

template < class RandIter> 
inline typename  IntSqrtModDom<RandIter>::Rep & IntSqrtModDom<RandIter>::sqrootmodtwolift (Rep & x1, const Rep & x0, const Rep & a,const unsigned long k, const Rep & pk) const{

//we have a square root of a mod 2^k : x0 and we have
//x1 = x0 + h*2^{k-1}
//with h = ((((a-x0^2)mod 2^{2k-2})/2^k)*x0^{-1}mod 2^{k-1}) mod 2^{k-1}
//is a square root of a mod 2^{2k-2}

	Rep rtmp (x0);
    x1 = 0;
    Rep p2k2 (pk);
    p2k2 *= pk;
    p2k2 <<= 2;
    rtmp %= pk;
    Rep tmp (rtmp);
    tmp *= tmp;
    tmp -= a;
    tmp %= p2k2;

    if (tmp == 0){
		x1 = rtmp;
		if (x1 < 0)
			x1 = -x1;
		return x1;
	}

    Rep h (0);
    inv (h, rtmp, pk / 2);
    h *= -tmp;
    h /= pk;
    h %= (pk / 2);
    x1 += h;
    x1 *= (pk / 2);
    x1 += rtmp;
    x1 %= p2k2;
    if (x1 < 0)
		x1 = -x1;
    return x1;


}
