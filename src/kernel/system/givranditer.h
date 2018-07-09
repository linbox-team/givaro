// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: Giorgi Pascal <pascal.giorgi@ens-lyon.fr>
// ==========================================================================

/** @file givranditer.h
 * @ingroup zpz
 * @brief NO DOC
 * Givaro ring Elements generator
 */

#ifndef __GIVARO_randiter_H
#define __GIVARO_randiter_H

#include "givaro/givconfig.h"
#include "givaro/givrandom.h"
#include "givaro/givtimer.h"

// For ModularBalancedRandIter
#include <sys/time.h>
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE
#endif

#include <stdlib.h>
#include <limits>

namespace Givaro {

  /** Random ring Element generator.
   *   This class defines a ring Element generator for all givaro ring (Gfq and Zpz)
   *   throught a template argument as a ring.
   *   The random generator used is the givrandom.
   */
  template <class Ring , class Type>
    class GIV_randIter
  {
  public:

    /** @name Common Object Interface.
     * These methods are required of all LinBox random ring Element generators.
     */
    //@{

    /** Ring Element type.
     * The ring Element must contain a default constructor,
     * a copy constructor, a destructor, and an assignment operator.
     */
    typedef typename Ring::Element Element;

    /** Constructor from ring, sampling size, and seed.
     * The random ring Element iterator works in the ring F, is seeded
     * by seed, and it returns any one Element with probability no more
     * than 1/min(size, F.cardinality()).
     * A sampling size of zero means to sample from the entire ring.
     * A seed of zero means to use some arbitrary seed for the generator.
     * This implementation sets the sampling size to be no more than the
     * cardinality of the ring.
     * @param F LinBox ring archetype object in which to do arithmetic
     * @param size constant integer reference of sample size from which to
     *             sample (default = 0)
     * @param seed constant integer reference from which to seed random number
     *             generator (default = 0)
     */
  GIV_randIter(const  Ring& F,
	       const size_t size = 0,
	       const uint64_t seed = 0)
    : _givrand(seed), _ring(F)
    {}

    /** Copy constructor.
     * Constructs ALP_randIter object by copying the random ring
     * Element generator.
     * This is required to allow generator objects to be passed by value
     * into functions.
     * In this implementation, this means copying the random ring Element
     * generator to which R._randIter_ptr points.
     * @param  R ALP_randIter object.
     */
#ifdef _GIVARO_RAND_LEGACY
  GIV_randIter(GIV_randIter&) = default;
#endif
  GIV_randIter(GIV_randIter&&) noexcept = default;

    /** Assignment operator.
     * Assigns ALP_randIter object R to generator.
     * In this implementation, this means copying the generator to
     * which R._randIter_ptr points.
     * @param  R ALP_randIter object.
     */
    GIV_randIter& operator=(GIV_randIter&&) noexcept = default;
    GIV_randIter& operator=(const GIV_randIter&) = delete;

    /** Random ring Element creator with assignement.
     * This returns a random ring Element from the information supplied
     * at the creation of the generator.
     * @return random ring Element
     */
      Element& operator()(Element& elt)
      {
	    return ring().random (_givrand, elt);
      } 
      Element& random(Element& elt)
      {
          return this->operator()(elt);
      } 
      Element operator()()
      {
          Element tmp;
          return this->operator()(tmp);
      } 
      Element random()
      {
          return this->operator()();
      } 

    const Ring& ring() const { return _ring; }



    //@} Common Object Iterface

    /** @name Implementation-Specific Methods.
     */
    //@{

//     /// Default constructor
//   GIV_randIter(void) : _size(0), _givrand(), _ring() {}

    //@}

  private:

    /// Random generator
    GivRandom _givrand;

    /// Ring
    const Ring& _ring;

  }; //  class GIV_randIter

  /** Random ring Element generator.
   *   This class defines a ring Element generator for all givaro modular rings (Gfq and Modular)
   *   throught a template argument as a ring.
   *   The random generator used is the givrandom.
   */
  template <class Ring>
    class ModularRandIter
    {
    public:

      /** @name Common Object Interface.
       * These methods are required of all LinBox random ring Element generators.
       */
      //@{

      /** Ring Element type.
       * The ring Element must contain a default constructor,
       * a copy constructor, a destructor, and an assignment operator.
       */
      typedef typename Ring::Element Element;

      /** Constructor from ring, sampling size, and seed.
       * The random ring Element iterator works in the ring F, is seeded
       * by seed, and it returns any one Element with probability no more
       * than 1/min(F.cardinality()).
       * size has no meaning in ModularRandIter.
       * A seed of zero means to use some arbitrary seed for the generator.
       * @param F LinBox ring archetype object in which to do arithmetic
       * @param seed constant integer reference from which to seed random number
       *             generator (default = 0)
       */
    ModularRandIter(const  Ring& F, const size_t& size = 0, const uint64_t& seed = 0)
      : _givrand( seed ), _ring(F) {}

      /** Copy constructor.
       * Constructs ALP_randIter object by copying the random ring
       * Element generator.
       * This is required to allow generator objects to be passed by value
       * into functions.
       * In this implementation, this means copying the random ring Element
       * generator to which R._randIter_ptr points.
       * @param  R ALP_randIter object.
       */
    ModularRandIter(ModularRandIter&&) noexcept = default;
#ifdef _GIVARO_RAND_LEGACY
    ModularRandIter(ModularRandIter&) = default;
#endif

      /** Assignment operator.
       * Assigns ALP_randIter object R to generator.
       * In this implementation, this means copying the generator to
       * which R._randIter_ptr points.
       * @param  R ALP_randIter object.
       */
      ModularRandIter& operator= (ModularRandIter&&) noexcept = default;
      ModularRandIter& operator= (const ModularRandIter&) = delete;

      /** Random ring Element creator with assignement.
       * This returns a random ring Element from the information supplied
       * at the creation of the generator.
       * @return random ring Element
       */
        Element& operator()(Element& elt)
	{
                // Create new random Elements
            return ring().random(_givrand, elt);
	}

        Element& random(Element& elt)
	{
            return this->operator()(elt);
            
	}

        Element operator()()
        {
            Element tmp;
            return this->operator()(tmp);
        }
        
        Element random()
        {
            return this->operator()();
        }
        
      //@} Common Object Iterface

      /** @name Implementation-Specific Methods.
       */
      //@{

//       /// Default constructor
//     ModularRandIter(void) : _givrand(), _ring() {}

      //@}

        const Ring& ring() const { return _ring; }

    private:

      /// Random generator
      GivRandom _givrand;

      /// Ring
      const Ring& _ring;

    }; //  class ModularRandIter

  /** UnparametricRandIter
   *
   * Imported from FFLAS-FFPACK UnparametricRandIter - AB 2015-01-12
   **/

  template <class Ring>
    class GeneralRingRandIter {
  public:
    typedef typename Ring::Element Element;

      GeneralRingRandIter(const Ring &F, const size_t& size = 0, uint64_t seed = 0) : _F(F), _size(size), _givrand(seed)
    {}
      GeneralRingRandIter(GeneralRingRandIter&&) noexcept = default;
      GeneralRingRandIter& operator=(GeneralRingRandIter&&) noexcept = default;
#ifdef _GIVARO_RAND_LEGACY
      GeneralRingRandIter(GeneralRingRandIter&) = default;
      GeneralRingRandIter& operator=(GeneralRingRandIter&) = default;
#endif

      Element& operator() (Element& a)
      {
          return ring().init(a, uint64_t( (_size == 0?_givrand():_givrand()% (1_ui64<<_size))));
      }
      Element& random (Element& a)
      {
          return this->operator()(a);
      }
      Element operator() ()
      {
          Element a; return this->operator()(a);
      }
      Element random ()
      {
          return this->operator()();
      }

      const Ring& ring() const { return _F; }

  private:
    const Ring& _F;
    size_t _size; 
            /// Random generator
    GivRandom _givrand;

  };


  /** Random iterator for nonzero random numbers
   *
   * Wraps around an existing random iterator and ensures that the output
   * is entirely nonzero numbers.
   * Imported from FFLAS-FFPACK NonzeroRandIter - AB 2015-01-12
   **/
  template <class Ring, class RandIter = typename Ring::RandIter>
    class GeneralRingNonZeroRandIter
    {
    public:
    typedef typename Ring::Element Element;

    GeneralRingNonZeroRandIter(RandIter& r) : _r(r) {}

    Element& operator()(Element &a)
    {
      do _r.random(a); while ( ring().isZero(a));
      return a;
    }
    Element& random(Element &a)
    {
        return this->operator()(a);
    }

    Element operator()()
    {
        Element a; return this->operator()(a);
    }
    Element random()
    {
        return this->operator()();
    }

    const Ring& ring() const { return _r.ring(); }

    private:
    RandIter& _r;
    };

} // namespace Givaro

#endif // __GIVARO_randiter_H
