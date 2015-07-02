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
 * Givaro field Elements generator
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

  /** Random field Element generator.
   *   This class defines a field Element generator for all givaro field (Gfq and Zpz)
   *   throught a template argument as a field.
   *   The random generator used is the givrandom.
   */
  template <class Field , class Type>
    class GIV_randIter
  {
  public:

    /** @name Common Object Interface.
     * These methods are required of all LinBox random field Element generators.
     */
    //@{

    /** Field Element type.
     * The field Element must contain a default constructor,
     * a copy constructor, a destructor, and an assignment operator.
     */
    typedef typename Field::Element Element;

    /** Constructor from field, sampling size, and seed.
     * The random field Element iterator works in the field F, is seeded
     * by seed, and it returns any one Element with probability no more
     * than 1/min(size, F.cardinality()).
     * A sampling size of zero means to sample from the entire field.
     * A seed of zero means to use some arbitrary seed for the generator.
     * This implementation sets the sampling size to be no more than the
     * cardinality of the field.
     * @param F LinBox field archetype object in which to do arithmetic
     * @param size constant integer reference of sample size from which to
     *             sample (default = 0)
     * @param seed constant integer reference from which to seed random number
     *             generator (default = 0)
     */
  GIV_randIter(const  Field& F,
	       const size_t size = 0,
	       const uint64_t seed = 0)
    : _size(size), _givrand( GivRandom(seed) ), _field(F)
    {

      size_t cardinality    = size_t( F.size() );
      if ((_size > cardinality) || (_size == 0) )
	_size = cardinality;
    }

    /** Copy constructor.
     * Constructs ALP_randIter object by copying the random field
     * Element generator.
     * This is required to allow generator objects to be passed by value
     * into functions.
     * In this implementation, this means copying the random field Element
     * generator to which R._randIter_ptr points.
     * @param  R ALP_randIter object.
     */
  GIV_randIter(const GIV_randIter& R)
    : _size(R._size), _givrand(R._givrand) , _field(R._field) {}

    /** Destructor.
     * This destructs the random field Element generator object.
     * In this implementation, this destroys the generator by deleting
     * the random generator object to which _randIter_ptr points.
     */
    ~GIV_randIter(void) {}

    /** Assignment operator.
     * Assigns ALP_randIter object R to generator.
     * In this implementation, this means copying the generator to
     * which R._randIter_ptr points.
     * @param  R ALP_randIter object.
     */
    GIV_randIter<Field,Type>& operator=(const GIV_randIter<Field,Type>& R)
      {
	if (this != &R) // guard against self-assignment
	  {
	    _size = R._size;
	    _givrand = R._givrand;
	    _field = R._field;
	  }

	return *this;
      }

    /** Random field Element creator with assignement.
     * This returns a random field Element from the information supplied
     * at the creation of the generator.
     * @return random field Element
     */
      Element& operator()(Element& elt)
      {
	    return _field.random (_givrand, elt);
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

//     /// Default constructor
//   GIV_randIter(void) : _size(0), _givrand(), _field() {}

    //@}

  private:

    /// Sampling size
    size_t _size;

    /// Random generator
    GivRandom _givrand;

    /// Field
    const Field& _field;

  }; //  class GIV_randIter

  /** Random field Element generator.
   *   This class defines a field Element generator for all givaro modular rings (Gfq and Modular)
   *   throught a template argument as a field.
   *   The random generator used is the givrandom.
   */
  template <class Field>
    class ModularRandIter
    {
    public:

      /** @name Common Object Interface.
       * These methods are required of all LinBox random field Element generators.
       */
      //@{

      /** Field Element type.
       * The field Element must contain a default constructor,
       * a copy constructor, a destructor, and an assignment operator.
       */
      typedef typename Field::Element Element;

      /** Constructor from field, sampling size, and seed.
       * The random field Element iterator works in the field F, is seeded
       * by seed, and it returns any one Element with probability no more
       * than 1/min(size, F.cardinality()).
       * A sampling size of zero means to sample from the entire field.
       * A seed of zero means to use some arbitrary seed for the generator.
       * This implementation sets the sampling size to be no more than the
       * cardinality of the field.
       * @param F LinBox field archetype object in which to do arithmetic
       * @param seed constant integer reference from which to seed random number
       *             generator (default = 0)
       */
    ModularRandIter(const  Field& F, const size_t& size = 0, const size_t& seed = 0)
      : _givrand( GivRandom(seed) ), _field(F) {}

      /** Copy constructor.
       * Constructs ALP_randIter object by copying the random field
       * Element generator.
       * This is required to allow generator objects to be passed by value
       * into functions.
       * In this implementation, this means copying the random field Element
       * generator to which R._randIter_ptr points.
       * @param  R ALP_randIter object.
       */
    ModularRandIter(const ModularRandIter& R)
      : _givrand(R._givrand) , _field(R._field) {}

      /** Destructor.
       * This destructs the random field Element generator object.
       * In this implementation, this destroys the generator by deleting
       * the random generator object to which _randIter_ptr points.
       */
      ~ModularRandIter(void) {}

      /** Assignment operator.
       * Assigns ALP_randIter object R to generator.
       * In this implementation, this means copying the generator to
       * which R._randIter_ptr points.
       * @param  R ALP_randIter object.
       */
      ModularRandIter<Field>& operator=(const ModularRandIter<Field>& R)
	{
	  // guard against self-assignment
	  if (this != &R)
	    {
	      _givrand = R._givrand;
	      _field = R._field;
	    }
	  return *this;
	}

      /** Random field Element creator with assignement.
       * This returns a random field Element from the information supplied
       * at the creation of the generator.
       * @return random field Element
       */
        Element& operator()(Element& elt) 
	{
                // Create new random Elements
            return _field.random(_givrand, elt);
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
//     ModularRandIter(void) : _givrand(), _field() {}

      //@}

    private:

      /// Random generator
      GivRandom _givrand;

      /// Field
      const Field& _field;

    }; //  class ModularRandIter

  /** UnparametricRandIter
   *
   * Imported from FFLAS-FFPACK UnparametricRandIter - AB 2015-01-12
   **/

  template <class Ring>
    class GeneralRingRandIter {
  public:
    typedef typename Ring::Element Element;

      GeneralRingRandIter(const Ring &F, const size_t& size = 0, size_t seed = 0) : _F(F), _size(size), _givrand( seed==0? uint64_t(BaseTimer::seed()) : seed)
    {}
      GeneralRingRandIter(const GeneralRingRandIter<Ring> &R) : _F(R._F), _size(R._size) {}
      ~GeneralRingRandIter() {}

      Element& operator() (Element& a) const
      {
          return _F.init(a, uint64_t( (_size == 0?_givrand():_givrand()%_size)));
      }
      Element& random (Element& a) const
      {
          return this->operator()(a);
      }
      Element operator() () const
      {
          Element a; return this->operator()(a);
      }
      Element random () const
      {
          return this->operator()();
      }

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

    GeneralRingNonZeroRandIter(const Ring &F, RandIter &r) : _F(F), _r(r) {}
    GeneralRingNonZeroRandIter(const GeneralRingNonZeroRandIter& R) : _F(R._F), _r(R._r) {}
    ~GeneralRingNonZeroRandIter() {}

    Element& operator()(Element &a) 
    {
      do _r.random(a); while (_F.isZero(a));
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

    private:
    const Ring&     _F;
    RandIter& _r;
    };

} // namespace Givaro

#endif // __GIVARO_randiter_H
