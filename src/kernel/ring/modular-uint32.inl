// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz32uns.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givzpz32uns.inl,v 1.15 2011-02-02 17:16:43 bboyer Exp $
// ==========================================================================

#ifndef __GIVARO_zpz32_uns_INL
#define __GIVARO_zpz32_uns_INL

// Description:

// ---------
// -- normalized operations
// ---------


// r = a*b
#define __GIVARO_ZPZ32_Uns_MUL(r,p,a,b) ( r = (a*b) % p )
// r *= a
#define __GIVARO_ZPZ32_Uns_MULIN(r,p,a) ( r = (r*a) % p )

// r = a - b
#define __GIVARO_ZPZ32_Uns_SUB(r,p,a,b) ( r = (a>=b) ? a-b: (p-b)+a )

// r -= a
#define __GIVARO_ZPZ32_Uns_SUBIN(r,p,a) { if (r<a) r+=(p-a); else r-=a; }

// r = a+b
#define __GIVARO_ZPZ32_Uns_ADD(r,p,a,b) { r = (a+b); r= (r < p ? r : r-p); }
// r += a
#define __GIVARO_ZPZ32_Uns_ADDIN(r,p,a) { r += a;  r= (r < p ? r : r-p); }

// r <- a*b+c % p
#define __GIVARO_ZPZ32_Uns_MULADD(r,p,a,b,c) \
{ r = (a*b+c) % p;  }

#define __GIVARO_ZPZ32_Uns_MULADDIN(r,p,a,b) \
{ r += a*b; r= (r < p ? r : r % p);  }

#define __GIVARO_ZPZ32_Uns_NEG(r,p,a) ( r = (a == 0 ? 0 : p-a) )
#define __GIVARO_ZPZ32_Uns_NEGIN(r,p) ( r = (r == 0 ? 0 : p-r) )

// a*b-c
#define __GIVARO_ZPZ32_Uns_MULSUB(r,p,a,b,c) \
{ r = (a*b+p-c); r= (r<p ? r : r % p);  }
// a*b-c
#define __GIVARO_ZPZ32_Uns_SUBMULIN(r,p,a,b) \
{ r = (a*b+p-r); r= (r<p ? r : r % p); __GIVARO_ZPZ32_Uns_NEGIN(r,p); }

namespace Givaro {

	inline Modular<uint32_t>::Modular( ) :
		zero(0), one(1), mOne(0), _p(0), _dp(0.0)
	{}

	inline Modular<uint32_t>::Modular( Residu_t p ) :
		zero(0), one(1), mOne(Rep(p-1)),_p(p), _dp((double)p)

	{}

	inline Modular<uint32_t>::Residu_t Modular<uint32_t>::residu( ) const
	{
		return _p;
	}

	inline Modular<uint32_t>::Modular(const Modular<uint32_t>& F) :
		zero(F.zero), one(F.one), mOne(F.mOne), _p(F._p), _dp(F._dp)
	{ }

	inline Modular<uint32_t>::Rep& Modular<uint32_t>::mul (Rep& r, const Rep a, const Rep b) const
	{
		return __GIVARO_ZPZ32_Uns_MUL(r,_p,a,b);
	}

	inline Modular<uint32_t>::Rep& Modular<uint32_t>::sub (Rep& r, const Rep a, const Rep b) const
	{
		return __GIVARO_ZPZ32_Uns_SUB(r,_p,a,b);
	}

	inline Modular<uint32_t>::Rep& Modular<uint32_t>::add (Rep& r, const Rep a, const Rep b) const
	{
		__GIVARO_ZPZ32_Uns_ADD(r,_p,a,b);
		return r;
	}

	inline Modular<uint32_t>::Rep& Modular<uint32_t>::neg (Rep& r, const Rep a) const
	{
		return __GIVARO_ZPZ32_Uns_NEG(r,_p,a);
	}

	inline Modular<uint32_t>::Rep& Modular<uint32_t>::inv (Rep& r,
								 const Rep a) const
	{
		return Modular<uint32_t>::invext(r, a, _p);
	}

	inline Modular<uint32_t>::Rep& Modular<uint32_t>::div (Rep& r, const Rep a, const Rep b) const
	{
		/*
		   uint32_t tmp;
		   uint32_t ib;
		   inv(ib, b);
		   __GIVARO_ZPZ32_Uns_MUL(tmp,_p,a,ib);
		   return r = (Modular<uint32_t>::Rep)tmp;
		   */
		return mulin( inv(r, b), a );
	}

	// -- inline array operations between Modular<uint32_t>::Rep
	inline void Modular<uint32_t>::mul (const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			uint32_t tmp;
			__GIVARO_ZPZ32_Uns_MUL(tmp, _p,a[i], b[i]);
			r[i] = (Modular<uint32_t>::Rep)tmp;
		}
	}

	inline void Modular<uint32_t>::mul (const size_t sz, Array r, constArray a, Rep b) const
	{
		for ( size_t i=sz ; --i ; ) {
			uint32_t tmp;
			__GIVARO_ZPZ32_Uns_MUL(tmp, _p, a[i], b);
			r[i] = (Modular<uint32_t>::Rep)tmp;
		}
	}

	inline void Modular<uint32_t>::div (const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			div( r[i], a[i], b[i]);
		}
	}

	inline void Modular<uint32_t>::div (const size_t sz, Array r, constArray a, Rep b) const
	{
		Modular<uint32_t>::Rep ib;
		inv(ib, b);
		mul(sz, r, a, ib);
	}

	inline void Modular<uint32_t>::add (const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			uint32_t tmp;
			__GIVARO_ZPZ32_Uns_ADD(tmp, _p, a[i], b[i]);
			r[i] = (Modular<uint32_t>::Rep)tmp;
		}
	}

	inline void Modular<uint32_t>::add (const size_t sz, Array r, constArray a, Rep b) const
	{
		for ( size_t i=sz ; --i ; ) {
			uint32_t tmp;
			__GIVARO_ZPZ32_Uns_ADD(tmp, _p, a[i], b);
			r[i] = (Modular<uint32_t>::Rep)tmp;
		}
	}

	inline void Modular<uint32_t>::sub (const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			uint32_t tmp;
			__GIVARO_ZPZ32_Uns_SUB(tmp, _p, a[i], b[i]);
			r[i] = (Modular<uint32_t>::Rep)tmp;
		}
	}

	inline void Modular<uint32_t>::sub (const size_t sz, Array r, constArray a, Rep b) const
	{
		for ( size_t i=sz ; --i ; ) {
			uint32_t tmp;
			__GIVARO_ZPZ32_Uns_SUB(tmp, _p, a[i], b);
			r[i] = (Modular<uint32_t>::Rep)tmp;
		}
	}

	inline void Modular<uint32_t>::neg (const size_t sz, Array r, constArray a) const
	{
		for ( size_t i=sz ; --i ; ) {
			uint32_t tmp;
			__GIVARO_ZPZ32_Uns_NEG(tmp, _p, a[i]);
			r[i] = (Modular<uint32_t>::Rep)tmp;
		}
	}


	inline Modular<uint32_t>::Rep& Modular<uint32_t>::mulin (Rep& r, const Rep a) const
	{
		return __GIVARO_ZPZ32_Uns_MULIN(r,_p, a);
	}

	inline Modular<uint32_t>::Rep& Modular<uint32_t>::divin (Rep& r, const Rep a) const
	{
		Modular<uint32_t>::Rep ia;
		inv(ia, a);
		return mulin(r, ia);
	}

	inline Modular<uint32_t>::Rep& Modular<uint32_t>::addin (Rep& r, const Rep a) const
	{
		uint32_t tmp = r;
		__GIVARO_ZPZ32_Uns_ADDIN(tmp,_p, a);
		return r = (Modular<uint32_t>::Rep)tmp;
	}

	inline Modular<uint32_t>::Rep& Modular<uint32_t>::subin (Rep& r, const Rep a) const
	{
		uint32_t tmp = r;
		__GIVARO_ZPZ32_Uns_SUBIN(tmp,_p, a);
		return r = (Modular<uint32_t>::Rep)tmp;
	}


	inline Modular<uint32_t>::Rep& Modular<uint32_t>::negin (Rep& r) const
	{
		return __GIVARO_ZPZ32_Uns_NEGIN(r,_p);
	}

	inline Modular<uint32_t>::Rep& Modular<uint32_t>::invin (Rep& r) const
	{
		return Modular<uint32_t>::invext(r, r, _p);
	}


	inline Modular<uint32_t>::Rep& Modular<uint32_t>::axpy (Rep& r,
								  const Rep a, const Rep b, const Rep c) const
	{
		uint32_t tmp;
		__GIVARO_ZPZ32_Uns_MULADD(tmp, _p, a, b, c);
		return r = (Modular<uint32_t>::Rep)tmp;
	}

	inline Modular<uint32_t>::Rep&  Modular<uint32_t>::axpyin (Rep& r,
								     const Rep a, const Rep b) const
	{
		uint32_t tmp = r;
		__GIVARO_ZPZ32_Uns_MULADDIN(tmp, _p, a, b);
		return r = (Modular<uint32_t>::Rep)tmp;
	}


	inline void Modular<uint32_t>::axpy (const size_t sz, Array r,
					      constArray a, constArray x, constArray y) const
	{
		for ( size_t i=sz ; --i ; ) {
			uint32_t tmp;
			__GIVARO_ZPZ32_Uns_MULADD(tmp, _p, a[i], x[i], y[i]);
			r[i] = (Modular<uint32_t>::Rep)tmp;
		}
	}

	inline void Modular<uint32_t>::axpyin (const size_t sz, Array r,
						constArray a, constArray x) const
	{
		for ( size_t i=sz ; --i ; ) {
			uint32_t tmp = r[i];
			__GIVARO_ZPZ32_Uns_MULADDIN(tmp, _p, a[i], x[i]);
			r[i] = (Modular<uint32_t>::Rep)tmp;
		}
	}

	inline Modular<uint32_t>::Rep&  Modular<uint32_t>::axmy (Rep& r,
								   const Rep a, const Rep b, const Rep c) const
	{
		uint32_t tmp;
		__GIVARO_ZPZ32_Uns_MULSUB(tmp, _p, a, b, c);
		return r = (Modular<uint32_t>::Rep)tmp;
	}

	// r = c-a*b
	inline Modular<uint32_t>::Rep&  Modular<uint32_t>::maxpy (Rep& r,
								    const Rep a, const Rep b, const Rep c) const
	{
		uint32_t tmp=c;
		__GIVARO_ZPZ32_Uns_SUBMULIN(tmp, _p, a, b);
		return r = (Modular<uint32_t>::Rep)tmp;
	}

	// r -= a*b
	inline Modular<uint32_t>::Rep&  Modular<uint32_t>::maxpyin (Rep& r,
								      const Rep a, const Rep b) const
	{
		uint32_t tmp = r;
		__GIVARO_ZPZ32_Uns_SUBMULIN(tmp, _p, a, b );
		return r = (Modular<uint32_t>::Rep)tmp;
	}

	// r = a*b - r
	inline Modular<uint32_t>::Rep&  Modular<uint32_t>::axmyin (Rep& r,
								     const Rep a, const Rep b) const
	{
		maxpyin(r,a,b);
		return negin(r);
	}


	inline void Modular<uint32_t>::axmy (const size_t sz, Array r,
					      constArray a, constArray x, constArray y) const
	{
		for ( size_t i=sz ; --i ; ) {
			uint32_t tmp;
			__GIVARO_ZPZ32_Uns_MULSUB(tmp, _p, a[i], x[i], y[i]);
			r[i] = (Modular<uint32_t>::Rep)tmp;
		}
	}

	// r -= a*b
	inline void Modular<uint32_t>::maxpyin (const size_t sz, Array r,
						 constArray a, constArray x) const
	{
		for ( size_t i=sz ; --i ; ) {
			uint32_t tmp = r[i];
			__GIVARO_ZPZ32_Uns_SUBMULIN(tmp, _p, a[i], x[i]);
			r[i] = (Modular<uint32_t>::Rep)tmp;
		}
	}

	// ------------------------- Miscellaneous functions

	inline int Modular<uint32_t>::isZero(const Rep a) const
	{
		return a == Modular<uint32_t>::zero;
	}

	inline int Modular<uint32_t>::isOne(const Rep a) const
	{
		return a == Modular<uint32_t>::one;
	}

inline int Modular<uint32_t>::isMOne(const Rep a) const
	{
		return a == Modular<uint32_t>::mOne;
	}

	inline size_t Modular<uint32_t>::length(const Rep ) const
	{
		return Modular<uint32_t>::size_rep;
	}

	// ---------
	// -- misc operations
	// ---------


	inline  Modular<uint32_t>::Rep&  Modular<uint32_t>::init ( Rep& r, const double a ) const
	{
		int sign; double ua;
		if (a < 0.0) { sign =-1; ua = -a;}
		else { ua = a; sign =1; }
		if ( ua > Signed_Trait<uint32_t>::max()){
			//     ua -= (double)floor(ua * _invdp)*_dp;
			ua = fmod(ua,_dp);
			r = (Rep) ua;
		} else
			r = (ua >=_p) ? (uint32_t) ua % _p : (uint32_t) ua;
		if (r && (sign ==-1)) r = _p - r;
		return r;
	}

	inline  Modular<uint32_t>::Rep&  Modular<uint32_t>::init ( Rep& r, const float a ) const
	{
		return init(r, (double)a);
	}



	inline  Modular<uint32_t>::Rep&  Modular<uint32_t>::init ( Rep& r, const unsigned long a ) const
	{
		return r = (Rep)( a >= (unsigned long)_p ? a % (unsigned long)_p : a);
	}

	inline  Modular<uint32_t>::Rep&  Modular<uint32_t>::init ( Rep& r, const long a ) const
	{
		int sign;
		unsigned long ua;
		if (a <0) {
			sign =-1;
			ua = (unsigned long)-a;
		}
		else {
			ua = (unsigned long)a;
			sign =1;
		}
		r = Rep( (ua >=_p) ? ua % _p : ua );
		if (r && (sign ==-1))
			r = _p - r;
		return r;
	}

	inline Modular<uint32_t>::Rep&  Modular<uint32_t>::init ( Rep& r, const Integer& Residu ) const
	{
		long tr;
		if (Residu <0) {
			// -a = b [p]
			// a = p-b [p]
			if ( Residu <= (Integer)(-_p) )
				tr = long( (-Residu) % _p) ;
			else
				tr = long(-Residu);
			if (tr)
				return r = Rep( _p - (unsigned long)tr );
			else
				return r = zero;
		} else {
			if (Residu >= (Integer)_p )
				tr =   long(Residu % _p) ;
			else
				tr = long(Residu);
			return r = Rep(tr);
		}
	}





	inline Modular<uint32_t>::Rep& Modular<uint32_t>::init( Rep& a, const int i) const
	{
		return init(a,(long)i);
	}

	inline Modular<uint32_t>::Rep& Modular<uint32_t>::init( Rep& a, const unsigned int i) const
	{
		return init(a,(unsigned long)i);
	}


	inline void Modular<uint32_t>::assign
	( const size_t sz, Array r, constArray a ) const
	{
		for ( size_t i=sz ; --i ; ) {
			if (a[i] <Modular<uint32_t>::zero) {
				r[i] = a[i] + _p;
				if (r[i] <Modular<uint32_t>::zero) r[i] = r[i] % _p;
			}
			else if ((uint32_t)a[i] >_p) {
				r[i] = a[i] - _p;
				if ((uint32_t)r[i] >=_p) r[i] = r[i] % _p;
			}
			else r[i] = a[i];
		}
	}

	inline  Modular<uint32_t>::Rep&  Modular<uint32_t>::assign ( Rep& r, const long a ) const
	{
		int sign;
		uint32_t ua;
		if (a <0) {
			sign =-1;
			ua = (uint32_t) -a;
		}
		else {
			ua = uint32_t(a);
			sign = 1;
		}
		r = Rep( (ua >=_p) ? (ua % _p) : ua );
		if (sign ==-1)
			r = _p - r;
		return r;
	}

	inline  Modular<uint32_t>::Rep&  Modular<uint32_t>::assign ( Rep& r, const short a ) const
	{
		return Modular<uint32_t>::assign( r, (long)a);
	}

	inline  Modular<uint32_t>::Rep&  Modular<uint32_t>::assign ( Rep& r, const unsigned long a ) const
	{
		return r = Rep( (a >=_p) ? (a % _p) : a );
	}

	inline  Modular<uint32_t>::Rep&  Modular<uint32_t>::assign
	( Rep& r, const unsigned short a ) const
	{
		return r = Rep( (a >=_p) ? (a % _p) : a );
	}

	inline  Modular<uint32_t>::Rep&  Modular<uint32_t>::assign
	( Rep& r, const Rep a ) const
	{
		return assign(r, (long)a);
	}


	inline void Modular<uint32_t>::init
	( const size_t sz, Array r, constArray a ) const
	{
		for ( size_t i=sz ; --i ; )
			r[i] = a[i];
	}

	inline Modular<uint32_t>::Rep& Modular<uint32_t>::init ( Rep& r ) const
	{
		return r = zero;
	}

	inline Modular<uint32_t>::Rep&  Modular<uint32_t>::dotprod
	( Rep& r, const int bound, const size_t sz, constArray a, constArray b ) const
	{
		unsigned int stride = 1;
		if ((unsigned long)bound < GIVARO_MAXUINT32)
			//    stride = GIVARO_MAXULONG/((unsigned long)bound * (unsigned long)bound);
			stride = (unsigned int) ( GIVARO_MAXULONG/((unsigned long)bound) / ((unsigned long)bound) );
		unsigned long dot = zero;
		if ((sz <10) && (sz <stride)) {
			for(  size_t i= sz; i--; )
				dot += a[i] * b[i];
			if (dot > _p)
				return r = (Rep)(dot % _p);
			else
				return r = (Rep)dot;
		}
		size_t i_begin=0;
		stride &= (unsigned int)~0x1;
		if (stride ==0) {
			for(  size_t i= sz; --i; ) {
				dot += a[i] * b[i];
				if (dot>_p) dot %= _p;
			}
			return r = (Rep)dot;
		}
		do {
			size_t min_sz = ((sz-i_begin) < stride ? (sz-i_begin) : stride);
			if ((min_sz & 0x1) !=0) {
				min_sz--;
				i_begin++;
				dot += a++[min_sz] * b++[min_sz];
			}
			if (min_sz > 1)
				for(  size_t i= min_sz; i>0; --i, --i, ++a, ++a, ++b, ++b )
				{
					dot += a[0] * b[0];
					dot += a[1] * b[1];
				}
			if (dot>_p)
				dot %= _p;
			i_begin += min_sz;
		} while (i_begin <sz);
		return r = (Rep)dot;
	}

	template< class RandIter >
	inline  Modular<uint32_t>::Rep& Modular<uint32_t>::random(RandIter& g, Rep& a) const
	{
		return init(a, g());
	}

	template< class RandIter >
	inline  Modular<uint32_t>::Rep& Modular<uint32_t>::random(RandIter& g, Rep& a, const Rep& ) const
	{
		return init(a, g());
	}
	template< class RandIter >
	inline  Modular<uint32_t>::Rep& Modular<uint32_t>::random(RandIter& g, Rep& a, long b) const
	{
		return init(a, g() %(uint32_t) b);

	}

	template< class RandIter >
	inline  Modular<uint32_t>::Rep& Modular<uint32_t>::nonzerorandom(RandIter& g, Rep& a) const
	{
		while (isZero(init(a, g()))) {};
		return a;
	}

	template< class RandIter >
	inline  Modular<uint32_t>::Rep& Modular<uint32_t>::nonzerorandom(RandIter& g, Rep& a, const Rep& ) const
	{
		while (isZero(init(a, g()))) {};
		return a;
	}

	template< class RandIter >
	inline  Modular<uint32_t>::Rep& Modular<uint32_t>::nonzerorandom(RandIter& g, Rep& a, long b) const
	{
		while (isZero(init(a, g() %(uint32_t) b))) {};
		return a;
	}

	inline Modular<uint32_t>::Rep&  Modular<uint32_t>::dotprod
	( Rep& r, const size_t sz, constArray a, constArray b ) const
	{
		return Modular<uint32_t>::dotprod(r, (int)_p, sz, a, b);
	}


	//  a -> r: uint32_t to double
	inline void
	Modular<uint32_t>::i2d ( const size_t sz, double* r, constArray a ) const
	{
		for (size_t i=0; i<sz; ++i) r[i] = a[i];
	}

	//  a -> r: double to uint32_t
	inline void
	Modular<uint32_t>::d2i ( const size_t sz, Array r, const double* a ) const
	{
		union d_2_l {
			double d;
			uint32_t r[2];
		};
		//  static const double offset = 4503599627370496.0; // 2^52
		double offset = 4503599627370496.0; // 2^52
		for (size_t i=0; i<sz; ++i)
		{
			d_2_l tmp;
			// - normalization: put fractional part at the end of the representation
			tmp.d = a[i] + offset;
			r[i] = tmp.r[1];
			if (r[i] <_p) r[i] %= _p;
		}
		//    r[i] = (tmp.r[1] <_p ? tmp.r[1] : tmp.r[1]-_p);
		//    r[i] = (r[i] <_p ? r[i] : r[i]%_p);
		//    r[i] = (tmp.r[1] <_p ? tmp.r[1] : tmp.r[1]%_p);
	}



	// -- Input: (z, <_p>)
	inline std::istream& Modular<uint32_t>::read (std::istream& s)
	{
		char ch;
		s >> std::ws >> ch;
		if (ch != '(')
			//    GivError::throw_error( GivBadFormat("Modular<uint32_t>::read: syntax error: no '('"));
			std::cerr << "GivBadFormat(Modular<uint32_t>::read: syntax error: no '('))" << std::endl;

		s >> std::ws >> ch;
		if (ch != 'z')
			//    GivError::throw_error( GivBadFormat("Modular<uint32_t>::read: bad domain object"));
			std::cerr << "GivBadFormat(Modular<uint32_t>::read: bad domain object))" << std::endl;

		s >> std::ws >> ch;
		if (ch != ',')
			//    GivError::throw_error( GivBadFormat("Modular<uint32_t>::read: syntax error: no ','"));
			std::cerr << "GivBadFormat(Modular<uint32_t>::read: syntax error: no ',')) " << std::endl;

		s >> std::ws >> _p;

		s >> std::ws >> ch;
		if (ch != ')')
			//    GivError::throw_error( GivBadFormat("Modular<uint32_t>::read: syntax error: no ')'"));
			std::cerr << "GivBadFormat(Modular<uint32_t>::read: syntax error: no ')')) " << std::endl;

		return s;
	}

	inline std::ostream& Modular<uint32_t>::write (std::ostream& s ) const
	{
		return s << "Uns32 Givaro Z/pZ modulo " << residu();
	}

	inline std::istream& Modular<uint32_t>::read (std::istream& s, Rep& a) const
	{
		Integer tmp;
		s >> tmp;
		init(a, tmp);
		return s;
	}

	inline std::ostream& Modular<uint32_t>::write (std::ostream& s, const Rep a) const
	{
		return s << a;
	}

} // namespace Givaro

#endif // __GIVARO_zpz32_uns_INL

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
