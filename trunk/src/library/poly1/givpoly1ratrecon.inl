// ===============================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Author: J-G. Dumas
// Time-stamp: <05 Oct 09 14:07:29 Jean-Guillaume.Dumas@imag.fr>
// Description: generic rational fraction reconstruction
// ===============================================================
#ifndef __GIVARO_poly1_ratrecon_INL
#define __GIVARO_poly1_ratrecon_INL

template <class Domain>
bool Poly1Dom<Domain,Dense>::ratrecon(typename Poly1Dom<Domain,Dense>::Rep& N, typename Poly1Dom<Domain,Dense>::Rep& D, const typename Poly1Dom<Domain,Dense>::Rep& P, const typename Poly1Dom<Domain,Dense>::Rep& M, const Degree& dk) const {

  Degree degU, degV;
  this->degree(degU,P); this->degree(degV,M);
  if ((degU < dk) || (degV == 0)) { this->assign(N,P); this->assign(D,one); return true; }
  if ((degV < 0) || (degU == 0)) { this->assign(N,one); this->assign(D,one); return false; }

  typename Poly1Dom<Domain,Dense>::Rep U, V;
  this->assign(U, P);
  this->assign(V, M);

  Degree degN;
  typename Poly1Dom<Domain,Dense>::Rep Q, TMP, TMP2, D0;
  this->assign(D0,this->zero);
  this->assign(D,this->one);
  do {

      this->divmod(Q,N,V,U);

      this->assign(V,U);
      this->assign(U,N);
      this->maxpy(TMP2,Q,D,D0);
      this->assign(D0,D);
      this->assign(D,TMP2);
      this->degree(degN, N);
      if (degN <= dk) break;
  } while (degN>=0);

  typename Poly1Dom<Domain,Dense>::Rep G;
  Degree degG;
  if ( degree(degG, gcd(G, N, D)) > 0)
      return false;

  Type_t r; leadcoef(r, D);
  if (! _domain.isOne(r)) {
      this->divin(D, r);
      this->divin(N, r);
  }
  return true;
}
#endif // __GIVARO_poly1_ratrecon_INL
