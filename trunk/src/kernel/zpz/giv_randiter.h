/* Givaro field elements generator
 * Author : Giorgi Pascal   pascal.giorgi@ens-lyon.fr
 */

#ifndef _GIV_RANDITER_
#define _GIV_RANDITER_


#include "givaro/givinteger.h"
#include "givaro/givrandom.h"



template<class TAG>
class ZpzDom ;

// -- Tag for arithmetic:
class Std16 /*{public: typedef  int16 type;}*/ ; // -- standard arithmetic over 16bits representations.
class Std32 /*{public: typedef int32 type;}*/ ; // -- standard arithmetic over 32bits representations.

class Log16 ; // -- log arithmetic over 16bits representations.

template<> class ZpzDom<Std16>;
template<> class ZpzDom<Std32>;
template<> class ZpzDom<Log16>;
template<class TT> class GFqDom;

/** Random field element generator.
    This class defines a field element generator for all givaro field (Gfq and Zpz)
    throught a template argument as a field.
    The random generator used is the givrandom.
   */


template <class Field , class Type> class GIV_randIter
{ 
  public:
    
    /** @name Common Object Interface.
     * These methods are required of all LinBox random field element generators.
     */
    //@{
   
    /** Field element type.
     * The field element must contain a default constructor, 
     * a copy constructor, a destructor, and an assignment operator.
     */
    typedef typename Field::element element;    

    /** Constructor from field, sampling size, and seed.
     * The random field element iterator works in the field F, is seeded
     * by seed, and it returns any one element with probability no more
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
      : _size(size)  
      {	
	_field=F;
	GivRandom tmp(seed);
	_givrand=tmp;

	Type cardinality    = Type( F.size() );
	if ((_size > cardinality) || (_size == 0) )
	  _size = cardinality; 
      }

    /** Copy constructor.
     * Constructs ALP_randIter object by copying the random field
     * element generator.
     * This is required to allow generator objects to be passed by value
     * into functions.
     * In this implementation, this means copying the random field element
     * generator to which R._randIter_ptr points.
     * @param  R ALP_randIter object.
     */
    GIV_randIter(const GIV_randIter& R)
      : _size(R._size), _givrand(R._givrand) , _field(R._field) {}

    /** Destructor.
     * This destructs the random field element generator object.
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
 
    /** Random field element creator.
     * This returns a random field element from the information supplied
     * at the creation of the generator.
     * @return random field element
     */
    element& operator() (void)
    {
      // Create new random elements     
      long tmp = static_cast<long>((double (_givrand()) / double(_GIVRAN_MODULO_)) * double(_size));
      element* x=new element;
      _field.assign(*x , tmp);
      return *x;

    } // element& operator() (void)
	
    
    /** Random field element creator with assignement.
     * This returns a random field element from the information supplied
     * at the creation of the generator.
     * @return random field element
     */	
    element& random(element& elt)
      {
      // Create new random elements     
      long tmp = static_cast<long>((double (_givrand()) / double(_GIVRAN_MODULO_)) * double(_size));
      _field.assign(elt , tmp);
      return elt;
      
      } // element& random(element& )
      


    //@} Common Object Iterface
   
    /** @name Implementation-Specific Methods.
     * These methods are not required of all 
     * \Ref{LinBox Random field element generators}
     * and are included only for this implementation of the archetype.
     */
    //@{

    /// Default constructor
    GIV_randIter(void) : _size(0) {GivRandom tmp();_givrand=tmp;Field f(); _field=f;}
    
    //@}

 private:

   /// Sampling size
      Type _size;

   /// Random generator
      GivRandom _givrand;

   /// Field 	
      Field _field;  

  }; //  class GIV_randIter


#endif
