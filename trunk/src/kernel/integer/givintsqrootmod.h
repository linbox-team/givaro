// =================================================================== //
// Givaro : Modular square roots
//          Behaviour when a is not a quadratic  residue: 
//
//
// 
// Author : Yanis Linge
// Time-stamp: <11 Aug 09 15:36:49 Jean-Guillaume.Dumas@imag.fr> 
// =================================================================== //

#ifndef _GIVARO_SQRTMOD_
#define _GIVARO_SQRTMOD_

#include <iostream>
#include "givaro/givinteger.h"
#include "givaro/givintprime.h"
#include "givaro/givintfactor.h"
#include "givaro/givintrns.h"
#include "givaro/givrandom.h"
#include "givaro/givpower.h"
#include <math.h>

template < class RandIter = GivRandom > class IntSqrtModDom:public IntFactorDom < RandIter >
{
public:
	typedef IntFactorDom < RandIter > Father_t;
  	typedef typename IntFactorDom < RandIter >::Rep Rep;

	IntSqrtModDom (RandIter g = RandIter ())
		:IntFactorDom < RandIter > (g){}

  // =================================================================== //
  // Modular Square root functions
  // =================================================================== //
  Rep & sqrootmod (Rep & x, const Rep & a, const Rep & n) const{
    std::vector < Rep > Lf;
    std::vector < unsigned long > Le;
    Father_t::set (Lf, Le, n);

    typename std::vector < Rep >::const_iterator Lf_iter = Lf.begin ();
    typename std::vector < unsigned long >::const_iterator Le_iter = Le.begin ();

    std::vector < Rep > roots;
    Rep tmp;

    // Build prime powers
    std::vector < Rep > Pe (Lf.size ());
    typename std::vector < Rep >::iterator Pe_iter = Pe.begin ();
    for (; Pe_iter != Pe.end (); ++Pe_iter, ++Lf_iter, ++Le_iter)
      dom_power (*Pe_iter, *Lf_iter, *Le_iter, *this);

    Lf_iter = Lf.begin ();
    Le_iter = Le.begin ();
    Pe_iter = Pe.begin ();
    // root mod a power of 2
    if (*Lf_iter == Integer (2)){
	roots.push_back (this->sqrootmodpoweroftwo (tmp, a, *Le_iter, *Pe_iter));
	++Le_iter;
	++Lf_iter;
        ++Pe_iter;
      }

    // roots mod powers of odd primes
    for (; Lf_iter != Lf.end (); ++Lf_iter, ++Le_iter, ++Pe_iter){
	roots.push_back (this->
			 sqrootmodprimepower (tmp, a, *Lf_iter, *Le_iter, *Pe_iter));
      }

    // Chinese Remaindering
    IntRNSsystem < std::vector, std::allocator > RNs (Pe);

    RNs.RnsToRing (x, roots);
    return x;
  }

  // =================================================================== //
  // Modular Square root sub-functions
  // =================================================================== //

  // p is supposed to be prime
  Rep & sqrootmodprime (Rep & x, const Rep & a, const Rep & p) const;

  // p is supposed to be prime, modulo is taken mod p^k
  Rep & sqrootmodprimepower (Rep & x, const Rep & a, const Rep & p, const unsigned long k, const Rep &pk) const;

  // modulo is taken mod 2^k
  Rep & sqrootmodpoweroftwo (Rep & x, const Rep & a,const unsigned long k, const Rep & pk) const;

  //p is supposed to be prime and odd
  //use only onemorelift
  Rep & sqrootlinear (Rep & x, const Rep & a,const Rep & p,const unsigned long k) const;

  //p is supposed to be prime and odd
  //use only hensellift
  Rep & sqrootquad (Rep & x, const Rep & a,const Rep & p,const unsigned long k, const Rep & pk) const;

  //using only onemorelift
  Rep & sqroottwolinear(Rep & x, const Rep & a,const unsigned long k) const;

  //using only hensellift
  Rep & sqroottwoquad(Rep & x, const Rep & a,const unsigned long k, const Rep & pk) const;

 
protected:

  // =================================================================== //
  // Liftings
  // =================================================================== //

  // PRECONDITION: x0^2 = a [p^k]
  // RETURNS: x1 s.t. x1^2 = a [p^{2k}]
  Rep & sqroothensellift (Rep & x1, const Rep & x0, const Rep & a, const Rep & p, const unsigned long k, const Rep & pk) const;

  // PRECONDITION: x0^2 = a [p^k]
  // RETURNS: x1 s.t. x1^2 = a [p^{k+1}]
  Rep & sqrootonemorelift (Rep & x1, const Rep & x0, const Rep & a, const Rep & p, const unsigned long k, const Rep & pk) const;

  // PRECONDITION: x0^2 = a [2^k], with k>=3, a and x0 are odd
  // RETURNS: x1 s.t. x1^2 = a [p^{2k-2}]
  Rep & sqrootmodtwolift (Rep & x1, const Rep & x0, const Rep & a, const unsigned long k, const Rep & pk) const;

};



#include "givaro/givintsqrootmod.inl"

#endif
