// ==========================================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// file: gfqext.h
// Time-stamp: <26 Jan 24 15:58:15 Jean-Guillaume.Dumas@imag.fr>
// date: 2007
// version:
// author: Jean-Guillaume.Dumas

/*! @file gfqext.h
 * @ingroup zpz
 * @brief   Arithmetic on GF(p^k), with p a prime number less than 2^15.
 *   Specialized for fast conversions to floating point numbers.
 *  Main difference in interface is init/convert.
 * @bib
 * - JG Dumas, <i>Q-adic Transform Revisited</i>, ISSAC 2008
 *  @warning k strictly greater than 1
 */

#ifndef __GIVARO_gfq_extension_H
#define __GIVARO_gfq_extension_H

#include "givaro/gfq.h"
#include "givaro/givpower.h"

#include <limits>
#include <deque>

namespace Givaro {

	// init with preconditions, bad entry could segfault
	template<class TT> class GFqExtFast;

	// init defensive, bad entry are transformed, to the cost of slowdown
	template<class TT> class GFqExt;


	//! GFq Ext
	template<class TT> class GFqExtFast : public GFqDom<TT> {
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
		typedef GFqExtFast<TT> Self_t;

		typedef Rep Element;
		typedef Element* Element_ptr ;
		typedef const Element* ConstElement_ptr;


		typedef UTT Residu_t;

		typedef Rep* Array;
		typedef const Rep* constArray;

		GFqExtFast():
			Father_t()
			// , balanced(false)
		{
		}

		// Extension MUST be a parameter of the constructor
		GFqExtFast( const UTT P, const UTT e) :
			Father_t(P,e),
			_BITS( std::numeric_limits< double >::digits/( (e<<1)-1) ),
			_BASE(1 << _BITS),
			_MASK( _BASE - 1),
			_maxn( _BASE/(P-1)/(P-1)/e),
			_degree( e-1 )
			// , balanced(false)
		{

			GIVARO_ASSERT(_maxn>0 , "[GFqExtFast]: field too large");
			builddoubletables();

		}

		// Extension MUST be a parameter of the constructor
		template<typename Vector>
		GFqExtFast( const UTT P, const UTT e, const Vector& modPoly) :
			Father_t(P,e, modPoly),
			_BITS( std::numeric_limits< double >::digits/( (e<<1)-1) ),
			_BASE(1 << _BITS),
			_MASK( _BASE - 1),
			_maxn( _BASE/(P-1)/(P-1)/e),
			_degree( e-1 )
			// , balanced(false)
		{
			GIVARO_ASSERT(_maxn>0 , "[GFqExtFast]: field too large");
			builddoubletables();

		}

		virtual ~GFqExtFast() {};

		Self_t operator=( const Self_t& F)
		{
			this->zero = F.zero;
			this->one = F.one;
			this->mOne = F.mOne;
			this->_characteristic = F._characteristic;
			this->_dcharacteristic = F._dcharacteristic;
			this->_exponent = F._exponent;
			this->_q = F._q;
			this->_qm1 = F._qm1;
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

		GFqExtFast( const GFqDom<TT>& F) :
			Father_t(F),
			_BITS( F._BITS ), _BASE( F._BASE ),_MASK( F._MASK ),
			_maxn( F._maxn ),_degree( F._degree ),
			_log2dbl ( F._log2dbl ), _low2log( F._low2log ),
			_high2log (F._high2log )
			// , balanced(false)
		{
		}

		// Accesses

		UTT bits() const
		{ return _BITS;}

		UTT base() const
		{ return _BASE;}

		UTT mask() const
		{ return _MASK;}

		UTT maxdot() const
		{ return _maxn; }

		UTT& characteristic(UTT& a) const
		{ return a=this->_characteristic; }

		UTT characteristic() const
		{ return this->_characteristic; }

		// const bool balanced;

		// Init/Convert

		using Father_t::init;

		Rep& init( Rep& r, const unsigned long l) const
		{
			return Father_t::init(r,l);
		}

		virtual Rep& init(Rep& pad, const double d) const
		{
			GIVARO_ASSERT(d>=0.0 , "[GFqExtFast]: init from a negative number");
			GIVARO_ASSERT(d<_MODOUT, "[GFqExtFast]: init from a too large number");
			// WARNING WARNING WARNING WARNING
			// Precondition : 0 <= d < _MODOUT
			// Can segfault if d is too large
			// WARNING WARNING WARNING WARNING
			uint64_t rll( static_cast<uint64_t>(d) );
			uint64_t tll( static_cast<uint64_t>(d/this->_dcharacteristic) );
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
				padl ^= prec;
			}

			pad = (Rep)prec;
			for(size_t j = 0;j<_degree;++j) {
				rll >>= _BITS;
				tll >>= _BITS;
				prec = (UTT)(rll-tll*this->_characteristic);

				pad <<= _pceil;
				pad ^= prec;
			}

			padl = this->_low2log[(size_t)padl];
			pad = (Rep)this->_high2log[(size_t)pad];
			return this->addin(pad,(Rep)padl);
		}

		virtual Rep& init(Rep& pad, const float d) const
		{
			return init(pad, (double)d);
		}

		using Father_t::convert;

		virtual double& convert(double& d, const Rep a) const
		{
			return d=_log2dbl[(size_t)a];
		}

		virtual float& convert(float& d, const Rep a) const
		{
			return d=(float)_log2dbl[(size_t)a];
		}

		template<class RandIter> Rep& random(RandIter& g, Rep& r) const
		{
			return init(r, static_cast<double>( (UTT)g() % _MODOUT));
		}

	protected:

		void builddoubletables()
		{
			_log2dbl.resize(this->_log2pol.size());
			_pceil = 1;
			for(unsigned long ppow = 2; ppow < this->_characteristic; ppow <<= 1,++_pceil) {}
			unsigned long powersize = 1<<(_pceil * this->_exponent);
			_MODOUT = UTT(powersize - 1);
			_high2log.resize(powersize);
			_low2log.resize(powersize);



			typedef typename Father_t::Element ZElem;
			Father_t Zp(this->_characteristic,1);
			ZElem q,mq; Zp.init(q,2);
			dom_power(q,q,_BITS,Zp);
			Zp.neg(mq,q);

			Element xkmu;
			// This is X
			xkmu = (Element)this->_pol2log[(size_t)this->_characteristic];
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

				unsigned long binpolit = static_cast<unsigned long>(vect[0]);
				for(size_t i =1; i<this->_exponent; ++i) {
					binpolit <<= _pceil;
					binpolit += static_cast<unsigned long>(vect[i]);
				}

				ZElem tmp, prec, cour; Zp.init(cour);
				Zp.init(prec, vect[0]);
				for(size_t i = 1; i<this->_exponent; ++i) {
					Zp.init(cour, vect[i]);
					Zp.axpy(tmp, mq, cour, prec);
					low_ui.push_back(tmp);
					prec = cour;
				}

				PAD.eval(tmp , low_ui );
				_low2log[binpolit] = this->_pol2log[(size_t)tmp];

				low_ui.push_back(cour);
				PAD.eval( tmp, low_ui);
				Father_t::mul((Element&)_high2log[(size_t)binpolit],(Element)this->_pol2log[(size_t)tmp], xkmu);
			}
		}
	};



	//! GFq Ext (other)
	template<class TT> class GFqExt : public GFqExtFast<TT> {
	protected:
		typedef typename Signed_Trait<TT>::unsigned_type UTT;
		typedef TT Rep;
		typedef GFqDom<TT> Father_t;
		typedef GFqExtFast<TT> DirectFather_t;

		double _fMODOUT;

	public:
		typedef GFqExt<TT> Self_t;

		typedef Rep Element;
		typedef UTT Residu_t;

		typedef Rep* Array;
		typedef const Rep* constArray;

		GFqExt(): DirectFather_t(),
		_fMODOUT(static_cast<double>(this->_MODOUT)) {}

		GFqExt( const UTT P, const UTT e) :
			DirectFather_t(P,e),
			_fMODOUT(static_cast<double>(this->_MODOUT)) {}

		GFqExt( const GFqDom<TT>& F) :
			DirectFather_t(F),
			_fMODOUT(static_cast<double>(this->_MODOUT)) {}

		~GFqExt() {}

		using Father_t::init;

		virtual Rep& init(Rep& pad, const double d) const
		{
			// Defensive init
			const double tmp(fmod(d,this->_fMODOUT));
			return DirectFather_t::init(pad, (tmp>0.0)?tmp:(tmp+_fMODOUT) );
		}
		virtual Rep& init(Rep& pad, const float d) const
		{
			return init(pad, (double)d);
		}

	};

} // namespace Givaro

#endif // __GIVARO_gfq_extension_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
