// ==========================================================================
// file: givpowers.h 
// Time-stamp: <26 Feb 08 13:36:38 Jean-Guillaume.Dumas@imag.fr>
// (c) 1999 Givaro Team
// ==========================================================================

#ifndef _GIVARO_POWER_H_
#define _GIVARO_POWER_H_

// -------------------------------------------------------------
// Powering
// -------------------------------------------------------------

template<class TT, class UU> TT power(const TT n, const UU l) {
  if (l == 0) return 1 ;

  unsigned long p = l ;
  short is_assg = 0 ;

  TT res = TT(1) ;
  TT puiss  = n ;

  while (p != 0) {
      if (p & 0x1) {
          if (is_assg) 
              res *= puiss ;
          else { 
              is_assg = 1 ; 
              res = puiss ; 
          }    
          if ((p >>= 1) != 0) puiss = puiss * puiss ;
      }
  }
  return res ;
}
/*
#include <givaro/givinteger.h>
template<> Integer power(const Integer n, const long l) { return pow(n,l); }
template<> Integer power(const Integer n, const unsigned long l) { return pow(n,l); }
template<> Integer power(const Integer n, const int l) { return pow(n,l); }
template<> Integer power(const Integer n, const unsigned int l) { return pow(n,l); }
*/

template<class D, class TT> TT& dom_power(TT& res, const TT& n, long l, const D& F) {
  if (l == 0) return res = F.one ;

  unsigned long p = l ;
  short is_assg = 0 ;

  TT puiss = n, tmp ;
  res = F.one ;

  while (p != 0) {
      if (p & 0x1) {
          if (is_assg) {
              F.mulin(res,puiss) ;
          } else { 
              is_assg = 1 ; 
              res = puiss ; 
          }
          if ((p >>= 1) != 0) { F.mul(tmp,puiss,puiss) ; puiss = tmp; }
      }
  } 
  return res ;
}


/*
#include <math.h>

template<> double power<double>(const double a, const double e) {
   return pow(a,e);
}
*/

#endif
