// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givpoly1sqrfree.inl,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givpoly1sqrfree.inl,v 1.2 2005-01-06 17:10:50 jgdumas Exp $
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
size_t& Poly1Dom<Domain,Dense>::sqrfree(size_t& Nfact, Rep* Fact, const Rep& P) const
{
  if (Nfact ==0) return Nfact;
  GIVARO_ASSERT( (Fact !=0), "nul pointer");

  unsigned long count = 0;
  Rep A,B,C,W,Z,Y;
  Type_t lc;
  leadcoef(lc, P);
  init(Fact[count], 0, lc);
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
// write(cout << "W:", W) << endl;
    div(Y, B, C);
// write(cout << "Y:", Y) << endl;
    diff(Z, W);
// write(cout << "W':", Z) << endl;
    sub(Z, Y, Z);
    while (!iszero(Z)) 
    {
      gcd(Fact[count], W, Z);
// write(cout << "L" << count << ":", Fact[count]) << endl;

      div(C, W, Fact[count]); assign(W, C);
      div(Y, Z, Fact[count]);
      diff(Z, W);
      sub(Z, Y, Z); 
      if (++count > Nfact) return Nfact;
    } 
  }
  assign(Fact[count], W);
// write(cout << "L" << count << ":", Fact[count]) << endl;
  return Nfact = ++count;
} 

