/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
// ==========================================================================
// Copyright(c)'1994-2016 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: Bastien Vialla <bastien.vialla@lirmm.fr>
// ==========================================================================

#ifndef __GIVARO_MODULAR_EXTENDED_INL
#define __GIVARO_MODULAR_EXTENDED_INL

namespace Givaro {

	// ----------------
	// ----- IO methods

	template<>
	inline
	std::ostream &ModularExtended<float>::write (std::ostream &os) const
	{
		return os << "ModularExtended<float> mod " << _p;
	}

	template<>
	inline
	std::ostream &ModularExtended<double>::write (std::ostream &os) const
	{
		return os << "ModularExtended<double> mod " << _p;
	}

	template<typename _Element>
	inline
	std::istream &ModularExtended<_Element>::read (std::istream &is)
	{
		is >> _p;
		return is;
	}

	template<typename _Element>
	inline
	std::ostream &ModularExtended<_Element>::write (std::ostream &os, const Element &x) const
	{
		return os << static_cast<int64_t>(x);
	}

	template<typename _Element>
	inline
	std::istream &ModularExtended<_Element>::read (std::istream &is, Element &x) const
	{
		int64_t tmp;
		is >> tmp;
		init(x,tmp);
		return is;
	}

	// --------------------
	// ----- Initialisation de Modular<double>

	template<>
	template<>
	inline
	typename ModularExtended<double>::Element&
	ModularExtended<double>::init<const int64_t> (Element& x, const int64_t y) const
	{
		x = static_cast<Element>(std::abs(y) % _lp);
		if (y < 0) negin(x);
		return x;
	}

	template<>
	template<>
	inline
	typename ModularExtended<double>::Element&
	ModularExtended<double>::init<const uint64_t> (Element& x, const uint64_t y) const
	{
		return x = static_cast<Element>(y % (uint64_t)(_lp));
	}

	template<>
	template<>
	inline
	typename ModularExtended<double>::Element&
	ModularExtended<double>::init<const Integer &> (Element& x, const Integer& y) const
	{
		x = static_cast<Element>(y % _lp);
		if (x < 0) x += _p;
		return x;
	}

	template<>
	template<>
	inline
	typename ModularExtended<double>::Element&
	ModularExtended<double>::init<const typename ModularExtended<double>::Element &> (Element& x, const Element& y) const
	{
		return x = y;
	}

	// --------------------
	// ----- Initialisation de ModularExtended<float>

	template<>
	template<>
	inline
	typename ModularExtended<float>::Element&
	ModularExtended<float>::init(typename ModularExtended<float>::Element& r, const double a) const
	{
		r = static_cast<float>(std::fmod(a, _p));
		if (r < 0.f) r += _p;
		return r;
	}

	template<>
	template<>
	inline
	typename ModularExtended<float>::Element&
	ModularExtended<float>::init(typename ModularExtended<float>::Element& r, const int32_t a) const
	{
		r = static_cast<Element>(std::abs(a) % _lp);
		if (a < 0) negin(r);
		return r;
	}

	template<>
	template<>
	inline
	typename ModularExtended<float>::Element&
	ModularExtended<float>::init(typename ModularExtended<float>::Element& r, const uint32_t a) const
	{
		return r = static_cast<Element>(a % uint32_t(_lp));
	}

	template<>
	template<>
	inline
	typename ModularExtended<float>::Element&
	ModularExtended<float>::init(typename ModularExtended<float>::Element& r, const int64_t a) const
	{
		r = static_cast<Element>(std::abs(a) % int64_t(_lp));
		if (a < 0) negin(r);
		return r;
	}

	template<>
	template<>
	inline
	typename ModularExtended<float>::Element&
	ModularExtended<float>::init(typename ModularExtended<float>::Element& r, const uint64_t a) const
	{
		return r = static_cast<Element>(a % uint64_t(_lp));
	}

	template<>
	template<>
	inline
	typename ModularExtended<float>::Element&
	ModularExtended<float>::init(typename ModularExtended<float>::Element& r, const Integer& a) const
	{
		r = static_cast<Element>(a % _lp);
		if (a < 0) negin(r);
		return r;
	}

	// --------------
	// Multiplication
	// --------------
	template<>
	inline
	typename ModularExtended<float>::Element&
	ModularExtended<float>::mul(typename ModularExtended<float>::Element& r,
								const typename ModularExtended<float>::Element& a,
								const typename ModularExtended<float>::Element& b) const {
#ifdef FP_FAST_FMAF
		Element abh, abl, pql, q;
		abh = a * b;
		abl = fma(a, b, -abh);
		q = std::floor(abh*_invp);
		pql = fma (-q, _p, abh);
		r = abl + pql;
#else
		Element abh, abl, pql, pqh, q;
		mult_dekker(a, b, abh, abl);
		q = std::floor(abh*_invp);
		mult_dekker(-q, _p, pqh, pql);
		r = (abh + pqh) + (abl + pql);
#endif
		if(r >= _p)
			r-= _p;
		else if(r < 0)
			r += _p;
#ifndef NDEBUG
		assert((r < _p) && (r >= 0));
#endif
		return r;
	}

	template<>
	inline
	typename ModularExtended<double>::Element&
	ModularExtended<double>::mul(typename ModularExtended<double>::Element& r,
								const typename ModularExtended<double>::Element& a,
								const typename ModularExtended<double>::Element& b) const {
#ifdef FP_FAST_FMA
		Element abh, abl, pql, q;
		abh = a * b;
		abl = fma(a, b, -abh);
		q = std::floor(abh*_invp);
		pql = fma (-q, _p, abh);
		r = abl + pql;
#else
		Element abh, abl, pql, pqh, q;
		mult_dekker(a, b, abh, abl);
		q = std::floor(abh*_invp);
		mult_dekker(-q, _p, pqh, pql);
		r = (abh + pqh) + (abl + pql);
#endif
		if(r >= _p)
			r-= _p;
		else if(r < 0)
			r += _p;
#ifndef NDEBUG
		assert((r < _p) && (r >= 0));
#endif
		return r;
	}

	// --------------
	//    Reduction
	// --------------
	template<>
	inline
	typename ModularExtended<float>::Element&
	ModularExtended<float>::reduce (typename ModularExtended<float>::Element& a) const{
#ifdef FP_FAST_FMAF
		Element pql, q;
		q = std::floor(a*_invp);
		pql = fma (-q, _p, a);
		a = pql;
#else
		Element pql, pqh, q;
		q = std::floor(a*_invp);
		mult_dekker(-q, _p, pqh, pql);
		a = (a + pqh) + pql;
#endif
		if(a >= _p)
			a-= _p;
		else if(a < 0)
			a += _p;
#ifndef NDEBUG
		assert((a < _p) && (a >= 0));
#endif
		return a;
	}

	template<>
	inline
	typename ModularExtended<double>::Element&
	ModularExtended<double>::reduce (typename ModularExtended<double>::Element& a) const{
#ifdef FP_FAST_FMA
		Element pql, q;
		q = std::floor(a*_invp);
		pql = fma (-q, _p, a);
		a = pql;
#else
		Element pql, pqh, q;
		q = std::floor(a*_invp);
		mult_dekker(-q, _p, pqh, pql);
		a = (a + pqh) + pql;
#endif
		if(a >= _p)
			a-= _p;
		else if(a < 0)
			a += _p;
#ifndef NDEBUG
		assert((a < _p) && (a >= 0));
#endif
		return a;
	}

	template<>
    inline bool ModularExtended<double>::isUnit(const Element& a) const{ 
        Element u,d; 
        invext(u,d,a,_p); 
        return isOne(d) || isMOne(d); 
    }

} // Givaro

#endif // __GIVARO_MODULAR_EXTENDED_INL
