#include <gmp.h>
#include <givaro/givgfq.h>
#include <givaro/givpoly1.h>
#include <givaro/givpoly1factor.h>
#include "givaro/givtablelimits.h"


template<class Rt> Rt FF_EXPONENT_MAX(const Rt p, const Rt e = 1) {
    Rt f = 0;
    for(Rt i = p; (i < (Rt)FF_TABLE_MAX) && (f < e); ++f, i*=p) 
        ;
    return f;
}

#define NEED_POLYNOMIAL_REPRESENTATION(p,e) ((e) > FF_EXPONENT_MAX((p),(e)))

#define EXTENSION(q,expo) ( NEED_POLYNOMIAL_REPRESENTATION((q),(expo)) ? Extension<>((q), (expo)) : GFqDom<long>((q), (expo)) )

template<class BFT = GFqDom<long>  >
class Extension {
public:
    typedef Extension<BFT> Self_t;
    typedef BFT BaseField_t;
    typedef typename BFT::Residu_t Residu_t;
    typedef Poly1FactorDom< BFT, Dense > Pol_t;

    typedef typename BFT::element BFelement; 
    typedef typename Pol_t::element Polelement;

protected:

    BaseField_t _bF;
    Pol_t _pD;
    Polelement _irred;
    const Residu_t _characteristic; 
    const Residu_t _exponent;
    const Residu_t _extension_order;
    const Integer _cardinality;

public:

    bool extension_type () const { return true; }

    typedef Polelement Element;
    typedef Polelement element;


    Extension ( const Residu_t p, const Residu_t e = 1)
            : _bF(p, FF_EXPONENT_MAX(p,e) ), _pD( _bF, "ý"  ), _characteristic( p ), _exponent ( e ) , _extension_order( e-FF_EXPONENT_MAX(p,e)+1 ) , _cardinality( pow(Integer(p),(unsigned long)(e)) ) {
/*     cerr << "Pol Cstor" << endl; */
        unsigned long basedegree = FF_EXPONENT_MAX(p,e) ;
        if (basedegree >= e) 
            cerr << "ERROR : Try a direct extension field GFDom instead of a polynomial extension" << endl;
        else
            _pD.creux_random_irreducible( _irred, _extension_order );
    }

    Extension ( const BaseField_t& bF, const Residu_t ex = 1)
            : _bF( bF ), _pD( _bF, "ý"  ), _characteristic( bF.characteristic() ), _exponent( ex + bF.exponent() ), _extension_order( ex ), _cardinality( pow( Integer(bF.cardinality()), (unsigned long)(ex) ) ) {
        _pD.creux_random_irreducible( _irred, ex);
    }

    Extension ( const Self_t& eF)
            : _bF( eF._bF ), _pD( eF._pD ), _irred( eF._irred ), _characteristic( eF._characteristic ), _exponent( eF._exponent ) , _extension_order( eF._extension_order ), _cardinality( eF._cardinality ) { }

    Self_t & operator=(const Self_t& eF) {
        if (this != &eF) {
            _bF = eF._bF;
            _pD = eF._pD;
            _irred = eF._irred;
            _characteristic = eF._characteristic;
            _exponent = eF._exponent;
            _extension_order = eF._extension_order;
            _cardinality = eF._cardinality;
        }
        return *this;
    }
    
    template<class XXX>
    Polelement& init( Polelement& e, const XXX& i) const { 
        return _pD.modin( _pD.init(e, i), _irred) ; 
    }
    
    Polelement& assign( Polelement& e, const Polelement& a) const { 
        return _pD.assign(e, a) ; 
    }

    template<class XXX>
    XXX& convert( XXX& i, const Polelement& e) const { 
        return _pD.convert( i, e) ; 
    }

    Polelement& add (Polelement& r, const Polelement& a, const Polelement& b) const {
        return _pD.add( r, a, b);
    }
  
    Polelement& sub (Polelement& r, const Polelement& a, const Polelement& b) const {
        return _pD.sub( r, a, b);
    }
  
    Polelement& neg (Polelement& r, const Polelement& a) const {
        return _pD.neg( r, a );
    }
  
    Polelement& mul (Polelement& r, const Polelement& a, const Polelement& b) const {
        return _pD.modin( _pD.mul( r, a, b), _irred );
    }
  
    Polelement& inv (Polelement& r, const Polelement& a) const {
//          _pD.write(_pD.write(_pD.write( std::cerr << "(", _pD.invmod( r, a, _irred)) << ") * (", a) << ")   == 1 + V * (", _irred) << std::endl;
         return  _pD.invmod( r, a, _irred);
     }
  
    Polelement& div (Polelement& r, const Polelement& a, const Polelement& b) const {
        return _pD.modin( _pD.mulin( inv(r, b), a), _irred );
    }
    
    Polelement& axpy (Polelement& r, const Polelement& a, const Polelement& b, const Polelement& c) const {
        return _pD.modin( _pD.addin(_pD.mul( r, a, b), c), _irred );
    }
  
    Polelement& addin(Polelement& r, const Polelement& b) const {
        return _pD.addin( r, b);
    }
  
    Polelement& subin(Polelement& r, const Polelement& b) const {
        return _pD.subin( r, b);
    }
  
    Polelement& negin(Polelement& r) const {
        return _pD.negin( r );
    }
  
    Polelement& mulin(Polelement& r, const Polelement& b) const {
        return _pD.modin( _pD.mulin( r, b), _irred );
    }
  
    Polelement& invin(Polelement& r) const {
         Polelement a(r);
         return _pD.invmod( r, a, _irred);
     }
  
    Polelement& divin(Polelement& r, const Polelement& b) const {
        Polelement tmp;
        inv(tmp,b);
        return _pD.modin( _pD.mulin( r, tmp), _irred );
    }
    
    Polelement& axpyin(Polelement& r, const Polelement& b, const Polelement& c) const {
        Polelement tmp; _pD.mul(tmp,b,c);
        return _pD.modin( _pD.addin( r, tmp), _irred );
    }
  
    bool areEqual (const Polelement& b, const Polelement& c) const {
        return _pD.areEqual( b, c) ;
    }
            
    bool isZero (const Polelement& b) const {
        return _pD.isZero(b) ;
    }
            
    bool isOne (const Polelement& b) const {
        return _pD.isone(b) ;
    }
            
    
    Integer &cardinality (Integer &c) const 
        { return c=_cardinality; }

    Integer &characteristic (Integer &c) const
        { return c=_characteristic; }

    Residu_t characteristic() const {
        return _characteristic;
    }

    Residu_t exponent() const {
        return _exponent;
    }

    Residu_t order() const {
        return _extension_order;
    }


    const BaseField_t& base_field() const {
        return _bF;
    }
    
            
    const Pol_t&  polynomial_domain() const {
        return _pD;
    }
    
            

    std::ostream&  write( std::ostream& o ) const {
        return _pD.write( _pD.write(o) << "/(", _irred) << ")";
    }
    

    std::istream& read ( std::istream& s, Polelement& a ) const { 
        _pD.read( s, a); 
        _pD.modin( a, _irred); 
        return s; 
    }

    std::ostream& write( std::ostream& o, const Polelement& R) const {
        return _pD.write( o, R );
    }


    std::istream&  read( std::istream& o ) const {
	std::cerr << "READ Extension, NOT YET IMPLEMENTED" << std::endl;
        return o;
    }
    
};


template <class ExtensionField, class Type> 
class GIV_ExtensionrandIter
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
    typedef typename ExtensionField::Polelement element;    

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
    GIV_ExtensionrandIter(const  ExtensionField& F,
		 const Type& size = 0, 
		 const Type& seed = 0)
            : _size(size), _givrand( GivRandom(seed) ), _field(F)
      {	
          Type charact    = Type( F.characteristic() );
          if ((_size > charact) || (_size == 0) )
              _size = charact; 
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
    GIV_ExtensionrandIter(const GIV_ExtensionrandIter& R)
      : _size(R._size), _givrand(R._givrand) , _field(R._field) {}

    /** Destructor.
     * This destructs the random field element generator object.
     * In this implementation, this destroys the generator by deleting 
     * the random generator object to which _randIter_ptr points.
     */
    ~GIV_ExtensionrandIter(void) {}
    
    /** Assignment operator.
     * Assigns ALP_randIter object R to generator.
     * In this implementation, this means copying the generator to
     * which R._randIter_ptr points.
     * @param  R ALP_randIter object.
     */
    GIV_ExtensionrandIter<ExtensionField,Type>& operator=(const GIV_ExtensionrandIter<ExtensionField,Type>& R)
    {
      if (this != &R) // guard against self-assignment
      {
	_size = R._size;
	_givrand = R._givrand;
	_field = R._field;
      }

      return *this;
    }
 
    /** Random field element creator with assignement.
     * This returns a random field element from the information supplied
     * at the creation of the generator.
     * @return random field element
     */	
    element& random(element& elt) const
      {
      // Create new random elements     
          elt.resize(_field.order());
          for(typename element::iterator it = elt.begin(); it != elt.end() ; ++ it) {
              long tmp = static_cast<long>((double (_givrand()) / double(_GIVRAN_MODULO_)) * double(_size));
              (_field.base_field()).init(*it , tmp);
          }
          return elt;
      } // element& random(element& )
      
    /** Random field element creator with assignement.
     * This returns a random field element from the information supplied
     * at the creation of the generator.
     * @return random field element
     */	
    element& operator()(element& elt) const 
        {
            return this->random(elt);
        }

    /** Random field element creator.
     * This returns a random field element from the information supplied
     * at the creation of the generator.
     * @return random field element
     */
    element& operator() (void)
    {
      element* x=new element;
      return this->random(*x);

    } // element& operator() (void)
	
    

    //@} Common Object Iterface
   
    /** @name Implementation-Specific Methods.
     * These methods are not required of all 
     * \Ref{LinBox Random field element generators}
     * and are included only for this implementation of the archetype.
     */
    //@{

    /// Default constructor
    GIV_ExtensionrandIter(void) : _size(0) {GivRandom tmp();_givrand=tmp;ExtensionField f(); _field=f;}
    
    //@}

 private:

   /// Sampling size
      Type _size;

   /// Random generator
      GivRandom _givrand;

   /// ExtensionField 	
      ExtensionField _field;  

  }; //  class GIV_ExtensionrandIter
