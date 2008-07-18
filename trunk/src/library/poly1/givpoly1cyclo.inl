#ifndef _GIV_POLY1_CYCLO_INL_
#define _GIV_POLY1_CYCLO_INL_
// =============================================================== //
// Givaro / Athapascan-1
// Cyclotomic polynomials
// Time-stamp: <15 Jul 08 10:40:47 Jean-Guillaume.Dumas@imag.fr> 
// =============================================================== //
#include <givaro/givintfactor.h>

// ---------------------------------------------------------------
// Composition by a power of X
// ---------------------------------------------------------------
template<class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::power_compose(Rep& W, const Rep& P, long b) const
{
  Degree dp; degree(dp, P);
  Type_t lc;
  leadcoef(lc, P);
  init( W, b*dp.value(), lc); // all coeffs to zero ...
  for(long i=0;i<dp.value();++i) {
      _domain.assign(W[i*b], P[i]);
  }
  return setdegree(W);
}




// ---------------------------------------------------------------
// n th Cyclotomic polynomial 
// ---------------------------------------------------------------
template<class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::cyclotomic( Rep& P, long n) const
{
  // P must provide an indeterminate
//   Integer In(n);
    IntFactorDom<> IF;
    if (n <= 1) {
        init(P, Degree(1), _domain.one);
        _domain.assign(P[0], _domain.one);
        _domain.negin(P[0]);
        return setdegree(P);
    } else if (IF.isprime(n)) {
        init(P, Degree(n-1), _domain.one);
        for(size_t i=n-1;i--;)
            _domain.assign(P[i], _domain.one); 
        return setdegree(P);
    } 
    else {
        long q,f;
        q = n / 2;
        f = n % 2;
        if (f) {
            IntFactorDom<>::Rep If, Iq;
            IF.divexact(Iq, n, IF.factor(If,n) );
            IF.convert(f, If);
            IF.convert(q, Iq);
            Rep inter;
            cyclotomic(inter,q);
            if (q % f) {
                Rep res;
                power_compose(res, inter, f);
                return div(P, res, inter); 
            } else
                return power_compose(P,inter,f);
        } else { 
            if (q%2) {
                    // q odd
                cyclotomic(P,q);
                Degree di;
                degree(di, P);
                for(int i=1;i<=di.value();i+=2)
                    _domain.negin(P[i]);
                return setdegree(P);
            } else {
                    // q even
                Rep inter;
                cyclotomic(inter,q);
                Degree di;
                degree(di, inter);
                return power_compose(P, inter, 2);
            }
        }
    }
}

#endif
