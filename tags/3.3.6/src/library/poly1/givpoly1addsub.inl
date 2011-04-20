// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givpoly1addsub.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givpoly1addsub.inl,v 1.4 2011-02-02 16:23:56 bboyer Exp $
// ==========================================================================

#ifndef __GIVARO_poly_addsub_INL
#define __GIVARO_poly_addsub_INL

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::addin (Rep& R, const Rep& P) const
{
//     this->write(this->write(std::cout, R) << " += ", P) << std::endl;

  size_t i;
  size_t sP = P.size();
  size_t sR = R.size();
  if (sP == 0) return R;
  if (sR == 0) { return assign(R,P); }
//   if (sR == 0) { R.copy(P); return R; }
//   if (sR == sP){ _supportdomain.addin(R,P); return; }
  if (sR < sP) {
    Rep tmp; tmp.copy(P);
    for (i=0; i<sR; ++i) _domain.addin(tmp[i], R[i]);
    R.logcopy(tmp);
  }
  else {
    for (i=0; i<sP; ++i) _domain.addin(R[i], P[i]);
  }
  return R;
}

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::add(Rep& R, const Rep& P, const Rep& Q) const
{
  size_t sP = P.size();
  size_t sQ = Q.size();
  size_t sR = R.size();
  if (sP == 0) { R.copy(Q); return R; }
  if (sQ == 0) { R.copy(P); return R; }
// JGD 04.11.1999
//   if (sP == sQ) {
//     R.reallocate(sP);
//     _supportdomain.add(R,P,Q);
//     return;
//   }
  size_t i, max = sP < sQ ? sQ : sP;
  if (sR != max) R.reallocate(max);
  if (sP < sQ)
  {
    for (i=0; i<sP; ++i) _domain.add(R[i], P[i], Q[i]);
//     for (; i<sQ; ++i) _domain.assign(R[i], Q[i]);
// JGD 05.11.1999
    for (; i<sQ; ++i) R[i] = Q[i];
  }
  else {
    for (i=0; i<sQ; ++i) _domain.add(R[i], P[i], Q[i]);

//     for (; i<sP; ++i) _domain.assign(R[i], P[i]);
// JGD 05.11.1999
    for (; i<sP; ++i) R[i] = P[i];
  }
  return R;
}

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::add
 (Rep& R, const Rep& P, const Type_t& val) const
{
  size_t sP = P.size();
  if (sP == 0)  {
    R.reallocate(1);
    _domain.assign(R[0],val);
  }
  else {
    assign(R, P);
    _domain.add(R[0],P[0],val);
  }
return R;
}

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::add
 (Rep& R, const Type_t& val, const Rep& P) const
{
  size_t sP = P.size();
  if (sP == 0)  {
    R.reallocate(1);
    _domain.assign(R[0],val);
  }
  else {
    assign(R, P);
    _domain.add(R[0],val, P[0]);
  }
return R;
}

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::subin (Rep& R, const Rep& P) const
{
  size_t sP = P.size();
  size_t sR = R.size();
  if (sP == 0) return R;
  if (sR == 0) { return neg(R,P); }
//   if (sR == sP){ _supportdomain.subin(R,P); return; }
  if (sR < sP) {
    size_t i;
    Rep tmp; tmp.reallocate(sP);
    for (i=0; i<sR; ++i) _domain.sub(tmp[i], R[i], P[i]);
    for (; i<sP; ++i) _domain.neg(tmp[i], P[i]);
    R.logcopy(tmp);
  }
  else {
    for (size_t i=0; i<sP; ++i) _domain.subin(R[i], P[i]);
  }
return R;
}

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::sub(Rep& R, const Rep& P, const Rep& Q) const
{
  size_t sP = P.size();
  size_t sQ = Q.size();
  size_t sR = R.size();
  if (sQ == 0) { R.copy(P); return R; }
  if (sP == 0) { return neg(R,Q); }
//   if (sP == sQ) {
//     R.reallocate(sP);
//     _supportdomain.sub(R,P,Q);
//     return;
//   }
  size_t i, max = sP < sQ ? sQ : sP;
  if (sR != max) R.reallocate(max);
  if (sP < sQ)
  {
    for (i=0; i<sP; ++i) _domain.sub(R[i], P[i], Q[i]);
    for (; i<sQ; ++i) _domain.neg(R[i], Q[i]);
  }
  else {
    for (i=0; i<sQ; ++i) _domain.sub(R[i], P[i], Q[i]);
    for (; i<sP; ++i) _domain.assign(R[i], P[i]);
  }
  return R;
}

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::sub
 (Rep& R, const Rep& P, const Type_t& val) const
{
  size_t sP = P.size();
  if (sP == 0)  {
    R.reallocate(1);
    _domain.neg(R[0],val);
  }
  else {
    assign(R, P);
    _domain.sub(R[0],P[0],val);
  }
  return R;
}

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::sub
 (Rep& R, const Type_t& val, const Rep& P) const
{
  size_t sP = P.size();
  if (sP == 0)  {
    R.reallocate(1);
    _domain.neg(R[0],val);
  }
  else {
    neg(R, P);
    _domain.add(R[0],val, P[0]);
  }
  return R;
}

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::negin (Rep& R ) const
{
//     _supportdomain.negin(R);
    size_t sR = R.size();
    if (sR == 0) { return R; }
    for (size_t i=0; i<sR; ++i) _domain.negin(R[i]);
    return R;
}

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::neg (Rep& R, const Rep& P ) const
{
//   _supportdomain.neg(R,P);
  size_t sP = P.size();
  R.reallocate(sP);
  if (sP == 0) { return R; }
  for (size_t i=0; i<sP; ++i) _domain.neg(R[i],P[i]);
  return R;
}


#endif // __GIVARO_poly_addsub_INL
