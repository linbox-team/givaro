// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz16std.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givzpz16std.inl,v 1.18 2011-02-02 17:16:43 bboyer Exp $
// ==========================================================================
// Description:

// ---------
// -- normalized operations
// ---------
#ifndef __GIVARO_zpz16std_INL
#define __GIVARO_zpz16std_INL

// r = a - b
//#define __GIVARO_ZPZ16_N_SUB(r,p,a,b) { r = (a-b); r= (r < 0 ? r+p : r);}
#define __GIVARO_ZPZ16_N_SUB(r,p,a,b) ( r = Rep(a>=b? a-b: (p-b)+a) )

// r -= a
#define __GIVARO_ZPZ16_N_SUBIN(r,p,a) { r = Rep(r-a); r= Rep(r < 0 ? r+p : r);}

// r = a+b
#define __GIVARO_ZPZ16_N_ADD(r,p,a,b) { r = Rep(a+b); r= Rep(r < p ? r : r-p);}
// r += a
#define __GIVARO_ZPZ16_N_ADDIN(r,p,a) { r = Rep(r+a);  r= Rep(r < p ? r : r-p);}

#define __GIVARO_ZPZ16_N_NEG(r,p,a) ( r = Rep(a == 0 ? 0 : p-a) )
#define __GIVARO_ZPZ16_N_NEGIN(r,p) ( r = Rep(r == 0 ? 0 : p-r) )


namespace Givaro {


	inline ZpzDom<Std16>::Residu_t ZpzDom<Std16>::residu( ) const
	{
		return _p;
	}

	inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::mul (Rep& r, const Rep a, const Rep b) const
	{
		int32_t tmp;
		__GIVARO_ZPZ32_N_MUL(tmp,(int32_t)_p,(int32_t)a,(int32_t)b);
		return r = (ZpzDom<Std16>::Rep)tmp;
	}

	inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::neg (Rep& r, const Rep a) const
	{
		return __GIVARO_ZPZ16_N_NEG(r,_p,a);
	}

	inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::inv (Rep& r, const Rep a) const
	{
		int32_t u;
		ZpzDom<Std16>::invext(u, a, _p);
		//   if ((d != 1) && (d != -1)) std::cerr << "GivMathDivZero(Zpz::inv)" << std::endl;
		return r = (u<0)?(ZpzDom<Std16>::Rep)(u+_p):(ZpzDom<Std16>::Rep)u;
	}

	inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::div (Rep& r, const Rep a, const Rep b) const
	{
		int32_t tmp;
		Rep ib;
		inv(ib, b);
		__GIVARO_ZPZ32_N_MUL(tmp,(int32_t)_p,(int32_t)a,(int32_t)ib);
		return r = (ZpzDom<Std16>::Rep)tmp;
	}

	inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::sub (Rep& r, const Rep a, const Rep b) const
	{
		return __GIVARO_ZPZ16_N_SUB(r,_p,a,b);
	}

	inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::add (Rep& r, const Rep a, const Rep b) const
	{
		__GIVARO_ZPZ16_N_ADD(r,_p,a,b);
		return r;
	}


	// -- inline array operations between ZpzDom<Std16>::Rep
	inline void ZpzDom<Std16>::mul (const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp;
			__GIVARO_ZPZ32_N_MUL(tmp, (int32_t)_p,(int32_t)a[i], (int32_t)b[i]);
			r[i] = (ZpzDom<Std16>::Rep)tmp;
		}
	}

	inline void ZpzDom<Std16>::mul (const size_t sz, Array r, constArray a, Rep b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp;
			__GIVARO_ZPZ32_N_MUL(tmp, (int32_t)_p, (int32_t)a[i], (int32_t)b);
			r[i] = (ZpzDom<Std16>::Rep)tmp;
		}
	}

	inline void ZpzDom<Std16>::div (const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			div( r[i], a[i], b[i]);
		}
	}

	inline void ZpzDom<Std16>::div (const size_t sz, Array r, constArray a, Rep b) const
	{
		ZpzDom<Std16>::Rep ib;
		inv(ib, b);
		mul(sz, r, a, ib);
	}

	inline void ZpzDom<Std16>::add (const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; )
			__GIVARO_ZPZ16_N_ADD(r[i], _p, a[i], b[i]);
	}

	inline void ZpzDom<Std16>::add (const size_t sz, Array r, constArray a, Rep b) const
	{
		for ( size_t i=sz ; --i ; )
			__GIVARO_ZPZ16_N_ADD(r[i], _p, a[i], b);
	}

	inline void ZpzDom<Std16>::sub (const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; )
			__GIVARO_ZPZ16_N_SUB(r[i], _p, a[i], b[i]);
	}

	inline void ZpzDom<Std16>::sub (const size_t sz, Array r, constArray a, Rep b) const
	{
		for ( size_t i=sz ; --i ; )
			__GIVARO_ZPZ16_N_SUB(r[i], _p, a[i], b);
	}

	inline void ZpzDom<Std16>::neg (const size_t sz, Array r, constArray a) const
	{
		for ( size_t i=sz ; --i ; )
			__GIVARO_ZPZ16_N_NEG(r[i], _p, a[i]);
	}


	inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::mulin (Rep& r, const Rep a) const
	{
		int32_t tmp = (int32_t)r;
		__GIVARO_ZPZ32_N_MULIN(tmp,(int32_t)_p, (int32_t)a);
		return r = (ZpzDom<Std16>::Rep)tmp;
	}

	inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::divin (Rep& r, const Rep a) const
	{
		ZpzDom<Std16>::Rep ia;
		inv(ia, a);
		return mulin(r, ia);
	}

	inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::addin (Rep& r, const Rep a) const
	{
		__GIVARO_ZPZ16_N_ADDIN(r,_p, a);
		return r;
	}

	inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::subin (Rep& r, const Rep a) const
	{
		__GIVARO_ZPZ16_N_SUBIN(r,_p,a);
		return r;
	}


	inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::negin (Rep& r) const
	{
		return __GIVARO_ZPZ16_N_NEGIN(r,_p);
	}

	inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::invin (Rep& r) const
	{
		int32_t u;
		ZpzDom<Std16>::invext(u, r, _p);
		//   if ((d != 1) && (d != -1)) std::cerr << "GivMathDivZero(Zpz::inv)" << std::endl;
		return r = (u<0)?(ZpzDom<Std16>::Rep)(u+_p):(ZpzDom<Std16>::Rep)u;
	}


	inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::axpy
	(Rep& r, const Rep a, const Rep b, const Rep c) const
	{
		int32_t tmp;
		__GIVARO_ZPZ32_N_MULADD(tmp, (int32_t)_p, (int32_t)a, (int32_t)b, (int32_t)c);
		return r = (ZpzDom<Std16>::Rep)tmp;
	}

	inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::axpyin
	(Rep& r, const Rep a, const Rep b) const
	{
		int32_t tmp = (int32_t)r;
		__GIVARO_ZPZ32_N_MULADDIN(tmp, (int32_t)_p, (int32_t)a, (int32_t)b);
		return r = (ZpzDom<Std16>::Rep)tmp;
	}


	inline void ZpzDom<Std16>::axpy
	(const size_t sz, Array r, constArray a, constArray x, constArray y) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp;
			__GIVARO_ZPZ32_N_MULADD(tmp, (int32_t)_p, (int32_t)a[i], (int32_t)x[i], (int32_t)y[i]);
			r[i] = (ZpzDom<Std16>::Rep)tmp;
		}
	}

	inline void ZpzDom<Std16>::axpyin
	(const size_t sz, Array r, constArray a, constArray x) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp = (int32_t)r[i];
			__GIVARO_ZPZ32_N_MULADDIN(tmp, (int32_t)_p, (int32_t)a[i], (int32_t)x[i]);
			r[i] = (ZpzDom<Std16>::Rep)tmp;
		}
	}

	inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::maxpy
	(Rep& r, const Rep a, const Rep b, const Rep c) const
	{
		int32_t tmp;
		__GIVARO_ZPZ32_N_MUL(tmp, (int32_t)_p, (int32_t)a, (int32_t)b);
		return __GIVARO_ZPZ16_N_SUB(r, _p, c, (Rep)tmp);
	}


	inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::axmy
	(Rep& r, const Rep a, const Rep b, const Rep c) const
	{
		int32_t tmp;
		__GIVARO_ZPZ32_N_MULSUB(tmp, (int32_t)_p, (int32_t)a, (int32_t)b, (int32_t)c);
		return r = (ZpzDom<Std16>::Rep)tmp;
	}

	// r -= a*b
	inline ZpzDom<Std16>::Rep&  ZpzDom<Std16>::maxpyin
	(Rep& r, const Rep a, const Rep b) const
	{
		int32_t tmp = (int32_t)r;
		__GIVARO_ZPZ32_N_SUBMULIN(tmp, (int32_t)_p, (int32_t)a, (int32_t)b );
		return r = (ZpzDom<Std16>::Rep)tmp;
		//    int32_t tmp = (int32_t)r;
		//   __GIVARO_ZPZ16_N_SUBMULIN(tmp, (int32_t)_p, (int32_t)a, (int32_t)b );
		//   return r = (ZpzDom<Std16>::Rep)tmp;
	}

	// r = a*b - r
	inline ZpzDom<Std16>::Rep&  ZpzDom<Std16>::axmyin (Rep& r,
							   const Rep a, const Rep b) const
	{
		int32_t tmp = (int32_t)r;
		__GIVARO_ZPZ32_N_MULSUB(tmp, (int32_t)_p, (int32_t)a, (int32_t)b , tmp);
		return r = (ZpzDom<Std16>::Rep)tmp;
	}


	inline void ZpzDom<Std16>::axmy (const size_t sz, Array r,
					 constArray a, constArray x, constArray y) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp;
			__GIVARO_ZPZ32_N_MULSUB(tmp, (int32_t)_p, (int32_t)a[i], (int32_t)x[i], (int32_t)y[i]);
			r[i] = (ZpzDom<Std16>::Rep)tmp;
		}
	}

	// r -= a*b
	inline void ZpzDom<Std16>::maxpyin (const size_t sz, Array r,
					    constArray a, constArray x) const
	{
		for ( size_t i=sz ; --i ; ) {
			int32_t tmp = (int32_t)r[i];
			__GIVARO_ZPZ32_N_SUBMULIN(tmp, (int32_t)_p, (int32_t)a[i], (int32_t)x[i]);
			r[i] = (ZpzDom<Std16>::Rep)tmp;
		}
	}

	// ------------------------- Miscellaneous functions

	inline int ZpzDom<Std16>::areEqual(const Rep a, const Rep b) const
	{
		return a == b;
	}

	inline int ZpzDom<Std16>::areNEqual(const Rep a, const Rep b) const
	{
		return a != b;
	}

	inline int ZpzDom<Std16>::isZero(const Rep a) const
	{
		return a == ZpzDom<Std16>::zero;
	}

	inline int ZpzDom<Std16>::isnzero(const Rep a) const
	{
		return a != ZpzDom<Std16>::zero;
	}

	inline int ZpzDom<Std16>::isOne(const Rep a) const
	{
		return a == ZpzDom<Std16>::one;
	}



	inline size_t ZpzDom<Std16>::length(const Rep ) const
	{
		return ZpzDom<Std16>::size_rep;
	}

	// ---------
	// -- misc operations
	// ---------
	inline  ZpzDom<Std16>::Rep&  ZpzDom<Std16>::init ( Rep& r, const unsigned long a ) const
	{
		return r = (Rep)( a >= (unsigned long)_p ? a % (unsigned long)_p : a);
	}

	inline  ZpzDom<Std16>::Rep&  ZpzDom<Std16>::init ( Rep& r, const long a ) const
	{
		int sign; long ua;
		if (a <0) {
			sign =-1;
			ua = -a;
		}
		else {
			ua = a;
			sign =1;
		}
		r = Rep( (ua >=_p) ? ua % (uint16_t)_p : ua );
		if (r && (sign ==-1))
			r = (Rep)(_p - r);
		return r;
	}

	inline ZpzDom<Std16>::Rep&  ZpzDom<Std16>::init ( Rep& r, const Integer& Residu ) const
	{
		Rep tr;
		if (Residu <0) {
			// -a = b [p]
			// a = p-b [p]
			if ( Residu <= (Integer)(-_p) ) tr = Rep( (-Residu) % (uint16_t)_p) ;
			else tr = Rep(-Residu);
			if (tr)
				return r = Rep((uint16_t)_p - (uint16_t)tr);
			else
				return r = zero;
		}
		else {
			if (Residu >= (Integer)_p ) tr =   Rep(Residu % _p) ;
			else tr = Rep(Residu);
			return r = (Rep)tr;
		}
	}




	inline  ZpzDom<Std16>::Rep& ZpzDom<Std16>::init( Rep& a, const int i) const
	{
		return init(a,(long)i);
	}
	inline  ZpzDom<Std16>::Rep& ZpzDom<Std16>::init( Rep& a, const double i) const
	{
		return init(a,(long)i);
	}
	inline  ZpzDom<Std16>::Rep& ZpzDom<Std16>::init( Rep& a, const float i) const
	{
		return init(a,(double)i);
	}
	inline  ZpzDom<Std16>::Rep& ZpzDom<Std16>::init( Rep& a, const unsigned int i) const
	{
		return init(a,(unsigned long)i);
	}


	inline void ZpzDom<Std16>::assign
	( const size_t sz, Array r, constArray a ) const
	{
		for ( size_t i=sz ; --i ; )
			r[i] = a[i];
	}

	inline  ZpzDom<Std16>::Rep&  ZpzDom<Std16>::assign ( Rep& r, const Rep a ) const
	{  return r=a;
	}

	template< class RandIter >
	inline  ZpzDom<Std16>::Rep& ZpzDom<Std16>::random(RandIter& g, Rep& a) const
	{
		return init(a, g());
	}

	template< class RandIter >
	inline  ZpzDom<Std16>::Rep& ZpzDom<Std16>::random(RandIter& g, Rep& a, const Rep& b) const
	{
		return init(a, g());
	}

	template< class RandIter >
	inline  ZpzDom<Std16>::Rep& ZpzDom<Std16>::random(RandIter& g, Rep& a, long b) const
	{
		return init(a, g() %(uint16_t) b);
	}

	template< class RandIter >
	inline  ZpzDom<Std16>::Rep& ZpzDom<Std16>::nonzerorandom(RandIter& g, Rep& a) const
	{
		while (isZero(init(a, g()))) {};
		return a;
	}

	template< class RandIter >
	inline  ZpzDom<Std16>::Rep& ZpzDom<Std16>::nonzerorandom(RandIter& g, Rep& a, const Rep& b) const
	{
		while (isZero(init(a, g()))) {};
		return a;
	}

	template< class RandIter >
	inline  ZpzDom<Std16>::Rep& ZpzDom<Std16>::nonzerorandom(RandIter& g, Rep& a, long b) const
	{
		while (isZero(init(a, g() %(uint16_t) b))) {};
		return a;
	}

	inline void ZpzDom<Std16>::init
	( const size_t sz, Array r, constArray a ) const
	{
		for ( size_t i=sz ; --i ; )
			r[i] = a[i];
	}

	inline ZpzDom<Std16>::Rep& ZpzDom<Std16>::init ( Rep& r ) const
	{
		return r = zero;
	}

	inline void ZpzDom<Std16>::dotprod
	( Rep& r, const int bound, const size_t sz, constArray a, constArray b ) const
	{
		unsigned int stride = 1;
		if ((unsigned long)bound < GIVARO_MAXUINT16)
			stride = (unsigned int) (GIVARO_MAXULONG/((unsigned long)bound * (unsigned long)bound));
		unsigned long dot = (unsigned long)zero;
		if ((sz <10) && (sz <stride)) {
			for(  size_t i= sz; i--; )
				dot += (unsigned long)a[i] * (unsigned long)b[i];
			if (dot > _p) r = (Rep)(dot % (uint16_t)_p);
			else r = (Rep)dot;
			return;
		}
		unsigned int i_begin=0;
		stride &= (unsigned int)~0x1;
		if (stride ==0) {
			for(  size_t i = sz; i--; ) {
				dot += (unsigned long)a[i] * (unsigned long)b[i];
				if (dot>_p) dot %= _p;
			}
			r = (Rep)dot;
			return;
		}
		do {
			size_t min_sz = ((sz-i_begin) < stride ? (sz-i_begin) : stride);
			if ((min_sz & 0x1) !=0) {
				min_sz--;
				i_begin++;
				dot += (unsigned long)a++[min_sz] * (unsigned long)b++[min_sz];
			}
			if (min_sz > 1)
				for(  size_t i= min_sz; i>0; --i, --i, ++a, ++a, ++b, ++b )
				{
					dot += (unsigned long)a[0] * (unsigned long)b[0];
					dot += (unsigned long)a[1] * (unsigned long)b[1];
				}
			if (dot>_p) dot %= (uint16_t)_p;
			i_begin += (unsigned int) min_sz;
		} while (i_begin <sz);
		r = (Rep)dot;
	}

	inline void ZpzDom<Std16>::dotprod
	( Rep& r, const size_t sz, constArray a, constArray b ) const
	{
		return ZpzDom<Std16>::dotprod(r, _p, sz, a, b);
	}


	//  a -> r: int16_t to double
	inline void
	ZpzDom<Std16>::i2d ( const size_t sz, double* r, constArray a ) const
	{
		for (size_t i=0; i<sz; ++i) r[i] = a[i];
	}

	//  a -> r: double to int16_t
	inline void
	ZpzDom<Std16>::d2i ( const size_t sz, Array r, const double* a ) const
	{
		union d_2_l {
			double d;
			int32_t r[2];
		};
		//  static const double offset = 4503599627370496.0; // 2^52
		double offset = 4503599627370496.0; // 2^52
		for (size_t i=0; i<sz; ++i)
		{
			d_2_l tmp;
			// - normalization: put fractional part at the end of the representation
			tmp.d = a[i] + offset;
			r[i] = (Rep) tmp.r[1];
			if (r[i] <_p)
				r[i] = Rep(r[i]%_p);
		}
		//    r[i] = (tmp.r[1] <_p ? tmp.r[1] : tmp.r[1]-_p);
		//    r[i] = (r[i] <_p ? r[i] : r[i]%_p);
		//    r[i] = (tmp.r[1] <_p ? tmp.r[1] : tmp.r[1]%_p);
	}



	// -- Input: (z, <_p>)
	inline std::istream& ZpzDom<Std16>::read (std::istream& s)
	{
		char ch;
		s >> std::ws >> ch;
		if (ch != '(')
			//    GivError::throw_error( GivBadFormat("ZpzDom<Std16>::read: syntax error: no '('"));
			std::cerr << "GivBadFormat(ZpzDom<Std16>::read: syntax error: no '('))" << std::endl;

		s >> std::ws >> ch;
		if (ch != 'z')
			//    GivError::throw_error( GivBadFormat("ZpzDom<Std16>::read: bad domain object"));
			std::cerr << "GivBadFormat(ZpzDom<Std16>::read: bad domain object))" << std::endl;

		s >> std::ws >> ch;
		if (ch != ',')
			//    GivError::throw_error( GivBadFormat("ZpzDom<Std16>::read: syntax error: no ','"));
			std::cerr << "GivBadFormat(ZpzDom<Std16>::read: syntax error: no ',')) " << std::endl;

		s >> std::ws >> _p;

		s >> std::ws >> ch;
		if (ch != ')')
			//    GivError::throw_error( GivBadFormat("ZpzDom<Std16>::read: syntax error: no ')'"));
			std::cerr << "GivBadFormat(ZpzDom<Std16>::read: syntax error: no ')')) " << std::endl;

		return s;
	}

	inline std::ostream& ZpzDom<Std16>::write (std::ostream& s ) const
	{
		return s << "Std16 Givaro Z/pZ modulo " << residu();
	}

	inline std::istream& ZpzDom<Std16>::read (std::istream& s, Rep& a) const
	{
        	Integer tmp;
		s >> tmp;
		init(a, tmp);
		return s;
	}

	inline std::ostream& ZpzDom<Std16>::write (std::ostream& s, const Rep a) const
	{
		return s << a;
	}

} // namespace Givaro

#endif // __GIVARO_zpz16std_INL
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
