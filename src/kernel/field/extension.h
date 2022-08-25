// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*! @file extension.h
 * @ingroup zpz
 * @brief NO DOX
 */

#ifndef __GIVARO_extension_H
#define __GIVARO_extension_H

#include <gmp.h>
#include <givaro/gfq.h>
#include <givaro/givconfig.h>
#include <givaro/givpoly1.h>
#include <givaro/givpoly1factor.h>
#include <givaro/givpoly1padic.h>
#include "givaro/givtablelimits.h"

namespace Givaro {

	// ----- Forward declaration
	template <class ExtensionField, class Type>
    class GIV_ExtensionrandIter;

	//! XXX
    template<class Rt>
	Rt FF_EXPONENT_MAX(const Rt p, const Rt maxe = _GIVARO_FF_MAXEXPONENT_)
	{
	Rt f = 0;
	for(Rt i = p; (i < (Rt)_GIVARO_FF_TABLE_MAX) && (f < maxe); ++f, i*=p) ;
	return f;
    }

	//! XXX
    template<class Rt>
	Rt FF_SUBEXPONENT_MAX(const Rt p, const Rt e)
	{
	Rt f = FF_EXPONENT_MAX(p,e);
	for( ; f > 1; --f)
            if ((e % f) == 0) break;
	return f;
    }

#define NEED_POLYNOMIAL_REPRESENTATION(p,e) ((e) > FF_SUBEXPONENT_MAX((p),(e)))

#define EXTENSION(q,expo) ( NEED_POLYNOMIAL_REPRESENTATION((q),(expo)) ? Extension<>((q), (expo)) : GFqDom<int64_t>((q), (expo)) )


	//! XXX
    template<typename Field>
	int64_t Exponent_Trait(const Field& F)
	{
        return 1;
    }


	//! XXX
    template<>
	inline int64_t Exponent_Trait(const GFqDom<int64_t>& F)
	{
        return F.exponent();
    }

    template<typename BaseField> class Extension;

	//! XXX
    template<typename BaseField>
    int64_t Exponent_Trait(const Extension<BaseField>& F)
	{
        return F.exponent();
    }

//! Extension
    template<class BFT = GFqDom<int64_t>  >
    class Extension {
    public:
	typedef          Extension<BFT>                            Self_t;
	typedef          BFT                                  BaseField_t;
	typedef typename BFT::Element                           BFElement;
	typedef typename BFT::Residu_t                           Residu_t;

	typedef          Poly1FactorDom< BFT, Dense >               Pol_t;
	typedef typename Pol_t::Element                        PolElement;

    protected:

	BaseField_t           _bF;
	Pol_t                 _pD;
	PolElement         _irred;
	Residu_t  _characteristic;
	Residu_t _extension_order;
	Residu_t        _exponent;
	Integer      _cardinality;

    public:

	bool extension_type () const
            {
	       	return true;
            }

	typedef PolElement Element;
	typedef Element* Element_ptr ;
	typedef const Element* ConstElement_ptr;



	Element zero;
	Element one;
	Element mOne;

	Extension() {}


	Extension ( const Residu_t p, const Residu_t e = 1, const Indeter Y="Y") :
		_bF(p, FF_SUBEXPONENT_MAX(p,e) ), _pD( _bF, Y  ), _characteristic( p )
            , _extension_order( e/FF_SUBEXPONENT_MAX(p,e) ), _exponent ( e )
            , _cardinality( pow(Integer(p),e) ), zero (_pD.zero)
            , one (_pD.one), mOne(_pD.mOne)
            {
                    /*     cerr << "Pol Cstor" << endl; */
		int64_t basedegree = FF_SUBEXPONENT_MAX(p,e) ;
		if (basedegree >= (int64_t)e) {
                    std::cerr << "WARNING : Try a direct extension field GFDom instead of a polynomial extension" << std::endl;
                    _bF = BaseField_t(p, 1);
                    _pD = Pol_t(_bF, Y);
                    _extension_order = _exponent;
		}
		_pD.creux_random_irreducible( _irred, (int64_t)_extension_order );
            }



	Extension ( const BaseField_t& bF, const Residu_t ex = 1, const Indeter Y="Y") :
	       	_bF( bF )
            , _pD( _bF, Y  )
            , _characteristic(  (Residu_t) bF.characteristic() )
            , _extension_order( (Residu_t)( ex ) )
            , _exponent(        (Residu_t)(ex * (Residu_t)Exponent_Trait(bF)) )
            , _cardinality(     (Integer) pow( Integer(bF.cardinality()) , (uint64_t)(ex) ) )
            , zero(             (Element)(_pD.zero))
            , one (             (Element)(_pD.one))
            , mOne (             (Element)(_pD.mOne))
            {
				Degree eo ((int64_t)_extension_order);
		if (_cardinality < (1<<20) )
                    _pD.creux_random_irreducible( _irred, eo);
		else
                    _pD.random_irreducible( _irred,  eo);
            }


        Extension ( const Pol_t& polydomain, const PolElement& Irred) :
                _bF( polydomain.getdomain() )
            , _pD( polydomain )
            , _irred( Irred )
            , _characteristic(  (Residu_t) _bF.characteristic() )
            , _extension_order( (Residu_t) _pD.degree(Irred).value() )
            , _exponent(        (Residu_t)( _extension_order * (Residu_t)Exponent_Trait(_bF)) )
                , _cardinality(     (Integer) pow( Integer(_bF.cardinality()) , (uint64_t)_extension_order ) )
            , zero(             (Element)(_pD.zero))
            , one (             (Element)(_pD.one))
            , mOne (             (Element)(_pD.mOne))
            {
                if (polydomain.isOne(_irred)) {
                    if (_cardinality < (1<<20) )
                        _pD.creux_random_irreducible( _irred,  (int64_t) _extension_order);
                    else
                        _pD.random_irreducible( _irred,  (int64_t) _extension_order);
                }
            }

        Extension ( const Self_t& eF) :
                _bF( eF._bF ), _pD( eF._pD ), _irred( eF._irred )
            , _characteristic( eF._characteristic )
            , _extension_order( eF._extension_order )
            , _exponent( eF._exponent ), _cardinality( eF._cardinality )
            , zero (_pD.zero), one (_pD.one), mOne (_pD.mOne)
            { }

	Self_t & operator=(const Self_t& eF)
	{
		if (this != &eF) {
			_bF = eF._bF;
			_pD = eF._pD;
			_irred = eF._irred;
			_characteristic = eF._characteristic;
			_exponent = eF._exponent;
			_extension_order = eF._extension_order;
			_cardinality = eF._cardinality;
			zero = eF.zero;
			one = eF.one;
			mOne = eF.mOne;
		}
		return *this;
	}

	PolElement& init( PolElement& e) const
            {
		return _pD.init(e) ;
            }

	template<class XXX>
	PolElement& init( PolElement& e, const XXX& i) const
            {
		return _pD.modin( _pD.init(e, i), _irred) ;
            }

	PolElement& assign( PolElement& e, const BFElement& a) const
	{
		return _pD.assign(e, a) ;
	}

	PolElement& assign( PolElement& e, const PolElement& a) const
            {
		return _pD.assign(e, a) ;
            }

        Integer& convert(Integer& i, const PolElement& e) const {
            
            return (Poly1PadicDom<BaseField_t>(_pD)).eval(i, e);
        }
        PolElement& init(PolElement& e, const Integer& i) const {
            return (Poly1PadicDom<BaseField_t>(_pD)).radix(e, i, _extension_order);
        }
        
                
                
	template<class XXX>
	XXX& convert( XXX& i, const PolElement& e) const
            {
		return _pD.convert( i, e) ;
            }

	PolElement& add (PolElement& r, const PolElement& a, const PolElement& b) const
            {
		return _pD.add( r, a, b);
            }

	PolElement& sub (PolElement& r, const PolElement& a, const PolElement& b) const
            {
		return _pD.sub( r, a, b);
            }

	PolElement& neg (PolElement& r, const PolElement& a) const
            {
		return _pD.neg( r, a );
            }

	PolElement& mul (PolElement& r, const PolElement& a, const PolElement& b) const
            {
		return _pD.modin( _pD.mul( r, a, b), _irred );
            }

	PolElement& inv (PolElement& r, const PolElement& a) const
            {
                    //          _pD.write(_pD.write(_pD.write( std::cerr << "(", _pD.invmod( r, a, _irred)) << ") * (", a) << ")   == 1 + V * (", _irred) << std::endl;
		return  _pD.invmod( r, a, _irred);
            }

	PolElement& div (PolElement& r, const PolElement& a, const PolElement& b) const
            {
		return _pD.modin( _pD.mulin( inv(r, b), a), _irred );
            }

	PolElement& axpy (PolElement& r, const PolElement& a, const PolElement& b, const PolElement& c) const
            {
                    //         return _pD.modin( _pD.addin(_pD.mul( r, a, b), c), _irred );
                    //          return _pD.modin( _pD.axpy(r, a, b, c), _irred );
		return addin(mul(r,a,b),c);
            }

            // -- maxpy: r <- c - a * b mod p
	PolElement& maxpy (PolElement& r, const PolElement a, const PolElement b, const PolElement c) const
            {
		return _pD.modin( _pD.maxpy( r, a, b, c), _irred );
            }

            // -- maxpyin: r <- r - a * b mod p
	PolElement& maxpyin(PolElement& r, const PolElement a, const PolElement b) const
            {
		return _pD.modin( _pD.maxpyin( r, a, b), _irred );
            }

            // -- axmy: r <- a * x - y mod p
	PolElement& axmy  (PolElement& r, const PolElement a, const PolElement b, const PolElement c) const
            {
		return subin(mul(r,a,b),c);
            }

            // -- axmyin: r <- a * x - r mod p
	PolElement& axmyin(PolElement& r, const PolElement a, const PolElement b) const
            {
		maxpyin(r,a,b);
		return negin(r);
            }

	PolElement& addin(PolElement& r, const PolElement& b) const
            {
		return _pD.addin( r, b);
            }

	PolElement& subin(PolElement& r, const PolElement& b) const
            {
		return _pD.subin( r, b);
            }

	PolElement& negin(PolElement& r) const
            {
		return _pD.negin( r );
            }

	PolElement& mulin(PolElement& r, const PolElement& b) const
            {
		return _pD.modin( _pD.mulin( r, b), _irred );
            }

	PolElement& invin(PolElement& r) const
            {
		PolElement a(r);
		return _pD.invmod( r, a, _irred);
            }

	PolElement& divin(PolElement& r, const PolElement& b) const
            {
		PolElement tmp;
		inv(tmp,b);
		return _pD.modin( _pD.mulin( r, tmp), _irred );
            }

	PolElement& axpyin(PolElement& r, const PolElement& b, const PolElement& c) const
            {
		PolElement tmp; _pD.mul(tmp,b,c);
		return _pD.modin( _pD.addin( r, tmp), _irred );
            }

	bool areEqual (const PolElement& b, const PolElement& c) const
            {
		return _pD.areEqual( b, c) ;
            }

	bool isZero (const PolElement& b) const
            {
		return _pD.isZero(b) ;
            }

	bool isOne (const PolElement& b) const
            {
		return _pD.isOne(b) ;
            }
	bool isUnit (const PolElement& b) const
            {
		return _pD.isUnit(b) ;
            }
	bool isMOne (const PolElement& b) const
            {
		return _pD.isMOne(b) ;
            }


	template<class RandIter> Element& random(RandIter& g, Element& r) const
            {
	       	return _pD.random(g,r,Degree((int64_t)_extension_order-1));
            }
	template<class RandIter> Element& random(RandIter& g, Element& r, int64_t s) const
            {
	       	return _pD.random(g,r,(s>=_extension_order?_extension_order-1:s));
            }
	template<class RandIter> Element& random(RandIter& g, Element& r, const Element& b) const
            {
	      	return _pD.random(g,r,b.size());
            }
	template<class RandIter> Element& nonzerorandom(RandIter& g, Element& r) const
            {
	       	return _pD.nonzerorandom(g,r,Degree((int64_t)_extension_order-1));
            }
	template<class RandIter> Element& nonzerorandom(RandIter& g, Element& r, int64_t s) const
            {
	       	return _pD.nonzerorandom(g,r,(s>=_extension_order?_extension_order-1:s));
            }
	template<class RandIter> Element& nonzerorandom(RandIter& g, Element& r, const Element& b) const
            {
	      	return _pD.nonzerorandom(g,r,b.size());
            }

	typedef GIV_ExtensionrandIter< Self_t, Integer >  RandIter;



        // ----- Access to the modulus
    Residu_t residu() const { return _bF.residu(); }

	Integer &cardinality (Integer &c) const
            {
	       	return c=_cardinality;
            }

	Residu_t cardinality() const
            {
		return _cardinality ;
            }

	Integer &characteristic (Integer &c) const
            {
	       	return c=_characteristic;
            }

	Residu_t characteristic() const
            {
		return _characteristic;
            }

	int64_t & characteristic(int64_t & c) const
            {
		return c = (int64_t) _characteristic;
            }

	Residu_t exponent() const
            {
		return _exponent;
            }

	Residu_t order() const
            {
		return _extension_order;
            }

    	PolElement& irreducible(PolElement& P) const
            {
                return _pD.assign(P, _irred);
            }

    	const PolElement& irreducible() const {
            return _irred;
        }

	const BaseField_t& base_field() const
            {
		return _bF;
            }


	const Pol_t&  polynomial_domain() const
            {
		return _pD;
            }



	std::ostream&  write( std::ostream& o ) const
            {
		return _pD.write( _pD.write(o) << "/(", _irred) << ")";
            }


	std::istream& read ( std::istream& s, PolElement& a ) const
            {
		_pD.read( s, a);
		_pD.modin( a, _irred);
		return s;
            }

	std::ostream& write( std::ostream& o, const PolElement& R) const
            {
		return _pD.write( o, R );
            }


	std::istream&  read( std::istream& o ) const
            {
		std::cerr << "READ Extension, NOT YET IMPLEMENTED" << std::endl;
		return o;
            }

    };

//! Extension rand iters
    template <class ExtensionField, class Type>
    class GIV_ExtensionrandIter {

    public:

            /** @name Common Object Interface.
             * These methods are required of all LinBox random field Element generators.
             */
            //@{

            /** Field Element type.
             * The field Element must contain a default constructor,
             * a copy constructor, a destructor, and an assignment operator.
             */
	typedef typename ExtensionField::PolElement Element;
	typedef typename ExtensionField::Residu_t Residu_t;

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
	GIV_ExtensionrandIter(const  ExtensionField& F,
			      const Type& size = 0,
			      const Type& seed = 0) :
	       	_size(size), _givrand( GivRandom(seed) ), _field(F)
            {
		Type charact    = Type( F.characteristic() );
		if ((_size > charact) || (_size == 0) )
                    _size = charact;
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
	GIV_ExtensionrandIter(const GIV_ExtensionrandIter& R) :
	       	_size(R._size), _givrand(R._givrand) , _field(R._field)
            {}

            /** Destructor.
             * This destructs the random field Element generator object.
             * In this implementation, this destroys the generator by deleting
             * the random generator object to which _randIter_ptr points.
             */
	~GIV_ExtensionrandIter(void) {}

//             /** Assignment operator.
//              * Assigns ALP_randIter object R to generator.
//              * In this implementation, this means copying the generator to
//              * which R._randIter_ptr points.
//              * @param  R ALP_randIter object.
//              */
// 	GIV_ExtensionrandIter<ExtensionField,Type>& operator= ( const GIV_ExtensionrandIter< ExtensionField, Type >& R )
//             {
// 		if (this != &R) // guard against self-assignment
// 		{
//                     _size = R._size;
//                     _givrand = R._givrand;
//                     _field = R._field;
// 		}
// 		return *this;
//             }

            /** Random field Element creator with assignement.
             * This returns a random field Element from the information supplied
             * at the creation of the generator.
             * @return random field Element
             */
	Element& random(Element& elt) const
            {
                    // Create new random Elements
		elt.resize( (uint64_t)(_field.order()));
		for(typename Element::iterator it = elt.begin(); it != elt.end() ; ++ it) {
                    int64_t tmp = static_cast<int64_t>((double (_givrand()) / double(_GIVRAN_MODULO_)) * double(_size));
                    (_field.base_field()).init(*it , tmp);
			//(_field.base_field()) . random (*it);
		}
		return elt;
            } // Element& random(Element& )

            /** Random field Element creator with assignement.
             * This returns a random field Element from the information supplied
             * at the creation of the generator.
             * @return random field Element
             */
	Element& operator()(Element& elt) const
            {
		return this->random(elt);
            }

            /** Random field Element creator.
             * This returns a random field Element from the information supplied
             * at the creation of the generator.
             * @return random field Element
             */
	Element& operator() (void)
            {
		Element* x=new Element;
		return this->random(*x);

            } // Element& operator() (void)



            //@} Common Object Iterface

        const ExtensionField& ring() { return _field; }                

    private:

            /// Sampling size
	Type _size;

            /// Random generator
	GivRandom _givrand;

            /// ExtensionField
	const ExtensionField& _field;

    }; //  class GIV_ExtensionrandIter

} // namespace Givaro

#endif //__GIVARO_extension_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
