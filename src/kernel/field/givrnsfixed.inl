// Copyright(c)'1994-2011 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <08 Mar 11 16:23:50 Jean-Guillaume.Dumas@imag.fr>
// ==========================================================================
// Description:
//  Chinese Remainder Algorithm.
#ifndef __GIVARO_arithmodu_fixedprimes_INL
#define __GIVARO_arithmodu_fixedprimes_INL

#include "givaro/givpower.h"

namespace Givaro {

	template<class Ints>
	template<class smallIntVector>
	Ints& RNSsystemFixed<Ints>::RnsToRingLeft( Ints& I, const smallIntVector& residues, const int level, const int col )
	{
		if (level) {
			int lowercol=col<<1, lowercolnext=lowercol+1;
			int lowerlevel=level-1;
			Ints u0;
			RnsToRingLeft(u0, residues, lowerlevel, lowercol);
			RnsToRingRight(I, residues, lowerlevel, lowercolnext);

			I -= u0; 				// u1-u0
			I *= _primes[(size_t)lowerlevel][(size_t)lowercolnext];	// (u1-u0) * (p0(p0^{-1} mod p1))
			I += u0;				// (u1-u0)*M01+u0
			return Integer::modin(I, _primes[(size_t)level][(size_t)col]);// (u1-u0)*M01 +u0 mod p0p1, between 0 and p0p1-1
		}
		else {
			return I=residues[(size_t)col];
		}
	}

	template<class Ints>
	template<class smallIntVector>
	Ints& RNSsystemFixed<Ints>::RnsToRingRight( Ints& I, const smallIntVector& residues, const int level, const int col )
	{
		if (level) {
			int lowercol=col<<1, lowercolnext=lowercol+1;
			int lowerlevel=level-1;
			Ints u0;
			RnsToRingLeft(u0, residues, lowerlevel, lowercol);
			RnsToRingRight(I, residues, lowerlevel, lowercolnext);

			I -= u0; 				// u1-u0
			I *= _primes[(size_t)lowerlevel][(size_t)lowercolnext];	// (u1-u0) * (p0(p0^{-1} mod p1))
			return I += u0;				// (u1-u0)*M01+u0
		} else {
			return I=residues[(size_t)col];
		}
	}



	// Convert to an Ints:
	template<class Ints>
	template<class smallIntVector>
	Ints& RNSsystemFixed<Ints>::RnsToRing( Ints& I, const smallIntVector& rns)
	{
		int ir = (int)_RNS.Primes().size();
		typename RNS_t::array Reds( (size_t)ir );
		for(int i = (int)_primes.size(); i-- ;) {
			int is = (int)_primes[(size_t)i].size();
			if (is & 1) {
				RnsToRingLeft(Reds[--ir], rns, i, --is);
			}
		}

		return _RNS.RnsToRing( I, Reds);
	}


	// -- free memory allocated in array !
	template<class Ints>
	RNSsystemFixed<Ints>::~RNSsystemFixed()
	{}

	template<class Ints>
	RNSsystemFixed<Ints>::RNSsystemFixed () :
		_primes(0)
	{}

	template<class Ints>
	RNSsystemFixed<Ints>::RNSsystemFixed (const Self_t& R) :
		_primes(R._primes, givWithCopy())
	{}


	// -- Array of primes are given
	template<class Ints>
	RNSsystemFixed<Ints>::RNSsystemFixed( const array& inprimes) : _primes(0)
	{
		GIVARO_ASSERT( inprimes.size()>0, "[RNSsystemFixed<Ints>::RNSsystemFixed] bad size of array");
		_primes.reserve( GIVINTLOG(inprimes.size()) );
		_primes.resize(1);

		for(typename array::const_iterator pit = inprimes.begin(); pit != inprimes.end(); ++pit) {
			_primes.front().push_back( *pit );
			for(int i = 1; i< (int)_primes.size(); ++i) {
				int s = (int)_primes[(size_t)i-1].size();
				if (s & 1) break;
				else {
					Ints& p0(_primes[(size_t)i-1][(size_t)s-2]), & p1(_primes[(size_t)i-1][(size_t)s-1]);
					Ints prod(p0*p1);
					inv(p1, p0, p1) *= p0;
					_primes[(size_t)i].push_back( prod );
				}
			}

			int lastp = int(_primes.size()-1);

			if (! (_primes.back().size() & 1)) {

				//            std::cout << '[';
				//            for(int j=0; j<_primes.back().size(); ++j)
				//                std::cout << ' ' << _primes.back()[j];
				//            std::cout << ']'  << std::endl;



				array newlevel;

				int s = (int)_primes.back().size();

				Ints& p0(_primes[(size_t)lastp][(size_t)s-2]), & p1(_primes[(size_t)lastp][(size_t)s-1]);
				newlevel.push_back( p0*p1  );

				inv(p1, p0, p1) *= p0;
				_primes.push_back( newlevel );
			}

			//        std::cout  << std::endl << "----------------" << std::endl;
			//        std::cout << '[';
			//        for(int i=0; i<_primes.size(); ++i) {
			//            std::cout << '[';
			//            for(int j=0; j<_primes[i].size(); ++j)
			//                std::cout << ' ' << _primes[i][j];
			//            std::cout << ']'  << std::endl;
			//        }
			//        std::cout << ']' << std::endl;
			//        std::cout << "----------------" << std::endl;



		}

		int numodd=0;
		for(int i = (int)_primes.size(); i-- ; )
			if (_primes[(size_t)i].size() & 1) ++numodd;


		typename RNS_t::domains Mods((size_t)numodd);
		for(int i = (int)_primes.size(); i-- ; ) {
			if (_primes[(size_t)i].size() & 1)
				Mods[--numodd] = _primes[(size_t)i].back();
		}
		_RNS.setPrimes( Mods );
	}


	template<class Ints>
	const typename RNSsystemFixed<Ints>::tree& RNSsystemFixed<Ints>::Primes() const
	{
		return _primes;
	}


	template<class Ints>
	const Ints RNSsystemFixed<Ints>::ith(const size_t i) const
	{
		return _primes.front()[i];
	}

} // namespace Givaro

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
