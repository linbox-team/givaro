#include <gmp.h>
#include <givaro/givgfq.h>
#include <givaro/givpoly1.h>
#include <givaro/givpoly1factor.h>
#include "givaro/givtablelimits.h"


template<class Rt> Rt FF_EXPONENT_MAX(const Rt p, const Rt e = 1) {
    Rt f = 0;
    for(Rt i = p; (i < FF_TABLE_MAX) && (f < e); ++f, i*=p) 
        ;
    return f;
}

#define NEED_POLYNOMIAL_REPRESENTATION(p,e) ((e) > FF_EXPONENT_MAX((p),(e)))

#define EXTENSION(q,expo) ( NEED_POLYNOMIAL_REPRESENTATION(q,expo) ? Extension<>(q, expo) : GFqDom<long>(q, expo) )

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
    const Integer _cardinality;

public:

    bool extension_type () const { return true; }

    typedef Polelement Element;
    typedef Polelement element;


    Extension ( const Residu_t p, const Residu_t e = 1)
            : _bF(p, FF_EXPONENT_MAX(p,e) ), _pD( _bF, "ý"  ), _characteristic( p ), _exponent ( e ) , _cardinality( pow(Integer(p),(unsigned long)(e)) ) {
/*     cerr << "Pol Cstor" << endl; */
        unsigned long basedegree = FF_EXPONENT_MAX(p,e) ;
        if (basedegree >= e) 
            cerr << "Try a direct extension field GFDom instead of a polynomial extension" << endl;
        else
            _pD.creux_random_irreducible( _irred, e-basedegree+1 );
    }

    Extension ( const BaseField_t& bF, const Residu_t e = 1)
            : _bF( bF ), _pD( _bF, "ý"  ), _characteristic( bF.characteristic() ), _exponent( e + bF.exponent() ), _cardinality( pow( Integer(bF.cardinality()), (unsigned long)(e) ) ) {
        _pD.creux_random_irreducible( _irred, e);
    }

    Extension ( const Self_t& eF)
            : _bF( eF._bF ), _pD( eF._pD ), _irred( eF._irred ), _characteristic( eF._characteristic ), _exponent( eF._exponent ) , _cardinality( eF._cardinality ) { }

    Self_t & operator=(const Self_t& eF) {
        if (this != &eF) {
            _bF = eF._bF;
            _pD = eF._pD;
            _irred = eF._irred;
            _characteristic = eF._characteristic;
            _exponent = eF._exponent;
            _cardinality = eF._cardinality;
        }
        return *this;
    }
    
    Polelement& init( Polelement& e, const Integer i) const { 
        return _pD.modin( _pD.init(e, i), _irred) ; 
    }

    Polelement& assign( Polelement& e, const Polelement& a) const { 
        return _pD.assign(e, a) ; 
    }

    Integer& convert( Integer& i, const Polelement& e) const { 
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
         Polelement g, d, v;
         _pD.gcd( g, r, v, a, _irred);
         return  r;
     }
  
    Polelement& div (Polelement& r, const Polelement& a, const Polelement& b) const {
        return _pD.modin( _pD.mul( inv(r, b), a), _irred );
    }
    
    Polelement& axpy (Polelement& r, const Polelement& a, const Polelement& b, const Polelement& c) const {
        return _pD.modin( _pD.axpy( r, a, b, c), _irred );
    }
  
    Polelement& addin(Polelement& r, const Polelement& b) const {
        return _pD.addin( r, b);
    }
  
    Polelement& sub (Polelement& r, const Polelement& b) const {
        return _pD.subin( r, b);
    }
  
    Polelement& negin(Polelement& r) const {
        return _pD.negin( r );
    }
  
    Polelement& mulin(Polelement& r, const Polelement& b) const {
        return _pD.modin( _pD.mulin( r, b), _irred );
    }
  
    Polelement& invin(Polelement& r) const {
         Polelement g, d, v, a(r);
         _pD.gcd( g, r, v, a, _irred);
         return  r;
     }
  
    Polelement& divin(Polelement& r, const Polelement& b) const {
        Polelement tmp;
        inv(tmp,b);
        return _pD.modin( _pD.mulin( r, tmp), _irred );
    }
    
    Polelement& axpyin(Polelement& r, const Polelement& b, const Polelement& c) const {
        return _pD.modin( _pD.axpyin( r, b, c), _irred );
    }
  
    bool areEqual (const Polelement& b, const Polelement& c) {
        return _pD.areEqual( b, c) ;
    }
            
    bool isZero (const Polelement& b) {
        return _pD.isZero(b) ;
    }
            
    bool isOne (const Polelement& b) {
        return _pD.isOne(b) ;
    }
            
    
    Integer &cardinality (Integer &c) const 
        { return c=_cardinality; }

    Integer &characteristic (Integer &c) const
        { return c=_characteristic; }

   

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
// read(istream)

    
};

