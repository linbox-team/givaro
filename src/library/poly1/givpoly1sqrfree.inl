// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givpoly1sqrfree.inl,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givpoly1sqrfree.inl,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:


/** Sqrfree decomposition.
Decompose P such that: P = Fact[0]^0 * Fact[1]^1 * ... * Fact[P.degree()]^(P.degree()), 
with Fact[0] the leading coefficient.
The array Fact must be allocated before calling the function.
The size of Fact must be degP+1 is all factors should be computed. 
For more readeable version of the algorithm, see Geddes, p342.
@param Nfact [in] the size of Fact
@param Fact  [in] an array of dimension Nfact
@param Nfact [out] is the number of factor in the sqrfree decomposition
@param Fact  [out] contains at most Nfact factors of the decomposition.
*/
template <class Domain>
void Poly1Dom<Domain,Dense>::sqrfree(size_t& Nfact, Rep* Fact, const Rep& P) const
{
  if (Nfact ==0) return;
  GIVARO_ASSERT( (Fact !=0), "nul pointer");

  long count = 1;
  Rep A,B,C,W,Z,Y;
  Type_t lc;
  leadcoef(lc, P);
  init(Fact[0], 0, lc);
// write(cout << "P:", P) << endl;
// _domain.write(cout << "lc:", lc) << endl;
  assign(A, P);
  diff(B, A);
// write(cout << "P':", B) << endl;
  gcd(C, A, B);
// write(cout << "Gcd(P,P'):", C) << endl;
  div(A, P, lc);
// write(cout << "A/lc:", A) << endl;
  if (areEqual(C, one)) {
    assign(W, A);
  }
  else {
    div(W, A, C);
    div(Y, B, C);
    diff(Z, W);
    sub(Z, Y, Z);
    while (!iszero(Z)) 
    {
      gcd(Fact[count], W, Z);
      div(C, W, Fact[count]); assign(W, C);
      div(Y, Z, Fact[count]);
      diff(Z, W);
      sub(Z, Y, Z); 
      count++;
      if (count >= Nfact) return;
    } 
  }
  assign(Fact[count], W);
  Nfact = 1+ count;
  return;
} 

