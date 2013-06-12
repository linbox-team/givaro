//==================================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Author : Giorgi Pascal   pascal.giorgi@ens-lyon.fr
//==================================================================

/** @file giv_randiter.h
 * @ingroup zpz
 * @brief NO DOC
 * Givaro field Elements generator
 */

#ifndef __GIVARO_randiter_H
#define __GIVARO_randiter_H


#include "givaro/givinteger.h"
#include "givaro/givrandom.h"
#include "givaro/givzpztypes.h"

namespace Givaro {

template<class TAG> class ZpzDom ;

// -- Tag for arithmetic:
struct Std16 /*{public: typedef  int16_t type;}*/ ; // -- standard arithmetic over 16bits representations.
struct Std32 /*{public: typedef int32_t type;}*/ ; // -- standard arithmetic over 32bits representations.

struct Log16 ; // -- log arithmetic over 16bits representations.

template<> class ZpzDom<Std16>;
template<> class ZpzDom<Std32>;
template<> class ZpzDom<Log16>;
template<class TT> class GFqDom;

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
		 const Type& size = 0,
		 const Type& seed = 0)
            : _size(size), _givrand( GivRandom(seed) ), _field(F)
      {

	Type cardinality    = Type( F.size() );
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

    /** Random field Element creator.
     * This returns a random field Element from the information supplied
     * at the creation of the generator.
     * @return random field Element
     */
    Element& operator() (void)
    {
      // Create new random Elements
      long tmp = static_cast<long>((double (_givrand()) / double(_GIVRAN_MODULO_)) * double(_size));
      Element* x=new Element;
      _field.assign(*x , tmp);
      return *x;

    } // Element& operator() (void)


    /** Random field Element creator with assignement.
     * This returns a random field Element from the information supplied
     * at the creation of the generator.
     * @return random field Element
     */
    Element& random(Element& elt) const
      {
	// Create new random Elements
	//atroce
	//long tmp = static_cast<long>((double (_givrand()) / double(_GIVRAN_MODULO_)) * double(_size));
	//_field.assign(elt , tmp);
	_field.random (_givrand, elt);
	return elt;

      } // Element& random(Element& )



    //@} Common Object Iterface

    /** @name Implementation-Specific Methods.
     */
    //@{

    /// Default constructor
    GIV_randIter(void) : _size(0), _givrand(), _field() {}

    //@}

 private:

   /// Sampling size
      Type _size;

   /// Random generator
      GivRandom _givrand;

   /// Field
      Field _field;

  }; //  class GIV_randIter

} // namespace Givaro

#endif // __GIVARO_randiter_H
