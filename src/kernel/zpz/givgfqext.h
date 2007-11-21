#ifndef _GIVARO_GFQ_EXTENSION_H_
#define _GIVARO_GFQ_EXTENSION_H_

// ==========================================================================
// file: givgfqext.h 
// Time-stamp: <21 Nov 07 10:27:14 Jean-Guillaume.Dumas@imag.fr>
// (c) Givaro Team
// date: 2007
// version: 
// author: Jean-Guillaume.Dumas
// Description:
//   Arithmetic on GF(p^k), with p a prime number less than 2^15
//   WARNING : k strictly greater than 1
//   Specialized for fast conversions to floating point numbers
//   See [JG Dumas, Q-adic Transform Revisited, 2007]
//   Main difference in interface is init/convert:
//       - In GFq
// ==========================================================================

#include "givaro/givgfq.h"
#include "givaro/givpower.h"
#include <limits>
#include <vector>
#include <deque>



template<class TT> class GFqExt : public GFqDom<TT> {
protected:
    typedef typename Signed_Trait<TT>::unsigned_type UTT;
    typedef TT Rep;
    typedef GFqDom<TT> Father_t;

    UTT _BITS; 	// the q-adic transform will be with q=2^_BITS
    UTT _BASE;	// this is 2^_BITS
    UTT _MASK;	// 2^_BITS - 1
    UTT _maxn;	// Worst case Maximal number of multiplications
        	// without reduction
    UTT _degree;// exponent-1
    UTT _pceil;	// smallest such that characteristic<2^_pceil, 
                // used for fast for table indexing
    UTT _MODOUT;// Largest accepted double for init
    
        // Conversion tables from exponent to double z-adic representation
    std::vector<double> _log2dbl;	// Exponent to double
    std::vector<UTT> 	_high2log;	// Half double to exponent
    std::vector<UTT> 	_low2log;	// Other half in Time-Memory Trade-Off
   

public:
    typedef GFqExt<TT> Self_t;
    
    typedef Rep Element;
    typedef UTT Residu_t;

    typedef Rep* Array;
    typedef const Rep* constArray;

    typedef GIV_randIter< GFqExt<TT> , Rep> RandIter; 

    GFqExt(): Father_t(), balanced(false) {}

        // Extension MUST be a parameter of the constructor
    GFqExt( const UTT P, const UTT e) : Father_t(P,e),
                                        balanced(false),
        _BITS( std::numeric_limits< double >::digits/( (e<<1)-1) ), 
        _BASE(1 << _BITS),
        _MASK( _BASE - 1),
        _maxn( _BASE/(P-1)/(P-1)/e),
        _degree( e-1 )
    {

        GIVARO_ASSERT(_maxn>0 , "[GFqExt]: field too large");
        builddoubletables();
        
    }

    Self_t operator=( const Self_t& F)
      {
	this->zero = F.zero;
	this->one = F.one;
	this->_characteristic = F._characteristic;
	this->_dcharacteristic = F._dcharacteristic;
	this->_inversecharacteristic = F._inversecharacteristic;
	this->_exponent = F._exponent;
	this->_q = F._q;
	this->_qm1 = F._qm1;
	this->_qm1o2 = F._qm1o2;
	this->_log2pol = F._log2pol;
	this->_pol2log = F._pol2log;
	this->_plus1 = F._plus1;
        this->_BITS = F._BITS;
        this->_BASE = F._BASE;
        this->_MASK = F._MASK;
        this->_maxn = F._maxn;
        this->_degree = F._degree;
        this->_log2dbl = F._log2dbl;
        this->_low2log = F._low2log;
        this->_high2log = F._high2log;
	return *this;
      }

    GFqExt( const GFqDom<TT>& F) : Father_t(F),
        _BITS( F._BITS ), _BASE( F._BASE ),_MASK( F._MASK ),
        _maxn( F._maxn ),_degree( F._degree ),
        _log2dbl ( F._log2dbl ), _low2log( F._low2log ),
        _high2log (F._high2log ), balanced(false) {
    }
    
        // Accesses 
    
    UTT bits() const { return _BITS;}
    UTT base() const { return _BASE;}
    UTT mask() const { return _MASK;}
    UTT maxdot() const { return _maxn; }
    UTT& characteristic(UTT& a) const { return a=this->_characteristic; }
    const bool balanced;
            
    Rep& init( Rep& r, const unsigned long l) const {
        return Father_t::init(r,l);
    }
    
           
    using Father_t::init;
    
            
    virtual double& convert(double& d, const Rep a) const {
        return d=_log2dbl[a];
    }        

    virtual float& convert(float& d, const Rep a) const {
        return d=(float)_log2dbl[a];
    }

    virtual Rep& init(Rep& pad, const double d) const {
//             // Precondition : d < _MODOUT
//         unsigned __GIVARO_INT64 rll(d); 
//         unsigned __GIVARO_INT64 tll(this->_inversecharacteristic*d);
            // JGD 16.11.2007 : Much slower but defensive
            // What shall we do ?
        unsigned __GIVARO_INT64 rll(d); 
        rll &= _MODOUT;
        unsigned __GIVARO_INT64 tll(this->_inversecharacteristic*rll);

        UTT prec(0); 
        UTT padl = (UTT)(rll - tll*this->_characteristic);
        if (padl == this->_characteristic) {
            padl -= this->_characteristic;
            tll += 1;
        }
        for(size_t j = 0;j<_degree;++j) {
            rll >>= _BITS;
            tll >>= _BITS;
            prec = (UTT)(rll-tll*this->_characteristic);

            padl <<= _pceil;
            padl += prec;
        }

        pad = prec;
        for(size_t j = 0;j<_degree;++j) {
            rll >>= _BITS;
            tll >>= _BITS;
            prec = (UTT)(rll-tll*this->_characteristic);

            pad <<= _pceil;
            pad += prec;
        }
        
        padl = this->_low2log[padl];
        pad = this->_high2log[pad];
        return this->addin(pad,padl);
    }
    
    virtual Rep& init(Rep& pad, const float d) const {
        return init(pad, (double)d);
    }
    
           


protected:

    void builddoubletables() {
        _log2dbl.resize(this->_log2pol.size());
        _pceil = 1;
        for(ulong ppow = 2; ppow < this->_characteristic; ppow <<= 1,++_pceil) {}
        ulong powersize = 1<<(_pceil * this->_exponent);
        _MODOUT = powersize - 1;
        _high2log.resize(powersize);
        _low2log.resize(powersize);

       

        typedef typename Father_t::Element ZElem;
        Father_t Zp(this->_characteristic,1);
        ZElem q,mq; Zp.init(q,2UL);
        dom_power(q,q,_BITS,Zp);
        Zp.neg(mq,q);

        Element xkmu;
            // This is X
        xkmu = this->_pol2log[this->_characteristic];
            // This is X^{e-1}
        dom_power(xkmu,xkmu,this->_exponent-1,*this);
        
        

        typedef Poly1FactorDom< Father_t, Dense > PolDom;
        PolDom Pdom( Zp );
        
        typedef Poly1PadicDom< Father_t, Dense > PadicDom;
        PadicDom PAD(Pdom);
         
        Father_t Z2B(2,_BITS); 
        PolDom P2dom( Z2B );
        PadicDom P2AD( P2dom );        

        std::vector<double>::iterator dblit = _log2dbl.begin();
        typename std::vector<UTT>::const_iterator polit = this->_log2pol.begin();
        

        for( ; polit != this->_log2pol.end(); ++polit, ++dblit) {
            
            std::vector<double> vect;
            std::deque<ZElem> low_ui;

            P2AD.evaldirect( *dblit, 
                             PAD.radixdirect(
                                 vect,
                                 (double)(*polit),
                                 this->_exponent)
                             );

            ulong binpolit = vect[0];
            for(size_t i =1; i<this->_exponent; ++i) {
                binpolit <<= _pceil;
                binpolit += vect[i];
            }

            ZElem tmp, prec, cour; 
            Zp.init(prec, vect[0]);
            for(size_t i = 1; i<this->_exponent; ++i) {
                Zp.init(cour, vect[i]);
                Zp.axpy(tmp, mq, cour, prec);    
                low_ui.push_back(tmp);
                prec = cour;
            }
            
            PAD.eval(tmp , low_ui );
            _low2log[binpolit] = this->_pol2log[tmp];

            low_ui.push_back(cour);
            PAD.eval( tmp, low_ui);
            Father_t::mul((Element&)_high2log[binpolit],(Element)this->_pol2log[tmp], xkmu);
        }
    }
};


#endif
