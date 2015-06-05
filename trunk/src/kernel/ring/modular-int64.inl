// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz64std.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givzpz64std.inl,v 1.19 2011-02-02 17:16:43 bboyer Exp $
// ==========================================================================

#ifndef __GIVARO_zpz64std_INL
#define __GIVARO_zpz64std_INL

#include <givaro/modular-defines.h>

namespace Givaro {

	// -------------
	// ----- Modular

    template<>
    inline Modular<int64_t, uint64_t>::Residu_t
	Modular<int64_t, uint64_t>::getMaxModulus() { return 2147483647u; } // 2^31 - 1

    template<>
    inline Modular<int64_t, int64_t>::Residu_t
	Modular<int64_t, int64_t>::getMaxModulus() { return 2147483647u; }

	// ------------------------
	// ----- Classic arithmetic
	
	template<typename COMP>
    inline typename Modular<int64_t, COMP>::Element& Modular<int64_t, COMP>::mul
		(Element& r, const Element& a, const Element& b) const
	{
		return  __GIVARO_MODULAR_INTEGER_MUL(r,_p,a,b);
	}

	template<typename COMP>
    inline typename Modular<int64_t, COMP>::Element& Modular<int64_t, COMP>::sub
		(Element& r, const Element& a, const Element& b) const
	{
		return __GIVARO_MODULAR_INTEGER_SUB(r,_p,a,b);
	}

	template<typename COMP>
    inline typename Modular<int64_t, COMP>::Element& Modular<int64_t, COMP>::add
		(Element& r, const Element& a, const Element& b) const
	{
		__GIVARO_MODULAR_INTEGER_ADD(r,_p,a,b);
		return r;
	}

	template<typename COMP>
    inline typename Modular<int64_t, COMP>::Element& Modular<int64_t, COMP>::neg
		(Element& r, const Element& a) const
	{
		return __GIVARO_MODULAR_INTEGER_NEG(r,_p,a);
	}

	template<typename COMP>
    inline typename Modular<int64_t, COMP>::Element& Modular<int64_t, COMP>::inv
		(Element& r, const Element& a) const
	{
		invext(r, a, int64_t(_p));
		return (r < 0)? r += (int64_t)_p : r;
	}

	template<typename COMP>
    inline typename Modular<int64_t, COMP>::Element& Modular<int64_t, COMP>::div
		(Element& r, const Element& a, const Element& b) const
	{
		return mulin( inv(r,b), a );
	}
	template<typename COMP>
    inline typename Modular<int64_t, COMP>::Element& Modular<int64_t, COMP>::mulin
		(Element& r, const Element& a) const
	{
		return __GIVARO_MODULAR_INTEGER_MULIN(r,_p, a);
	}

	template<typename COMP>
    inline typename Modular<int64_t, COMP>::Element& Modular<int64_t, COMP>::divin
		(Element& r, const Element& a) const
	{
		typename Modular<int64_t, COMP>::Element ia;
		inv(ia, a);
		return mulin(r, ia);
	}

	template<typename COMP>
    inline typename Modular<int64_t, COMP>::Element& Modular<int64_t, COMP>::addin
		(Element& r, const Element& a) const
	{
		int64_t tmp = (int64_t)r;
		__GIVARO_MODULAR_INTEGER_ADDIN(tmp ,_p, a);
		return r = (typename Modular<int64_t, COMP>::Element)tmp;
	}

	template<typename COMP>
    inline typename Modular<int64_t, COMP>::Element& Modular<int64_t, COMP>::subin
		(Element& r, const Element& a) const
	{
		int64_t tmp = (int64_t)r;
		__GIVARO_MODULAR_INTEGER_SUBIN(tmp, _p, a);
		return r = (typename Modular<int64_t, COMP>::Element)tmp;
	}

	template<typename COMP>
    inline typename Modular<int64_t, COMP>::Element& Modular<int64_t, COMP>::negin
		(Element& r) const
	{
		return __GIVARO_MODULAR_INTEGER_NEGIN(r,_p);
	}

	template<typename COMP>
    inline typename Modular<int64_t, COMP>::Element& Modular<int64_t, COMP>::invin
		(Element& r) const
	{
		return inv(r, r);
	}
	
	template<typename COMP>
    inline typename Modular<int64_t, COMP>::Element& Modular<int64_t, COMP>::axpy
		(Element& r, const Element& a, const Element& b, const Element& c) const
	{
		return __GIVARO_MODULAR_INTEGER_MULADD(r, _p, a, b, c);
	}

	template<typename COMP>
    inline typename Modular<int64_t, COMP>::Element&  Modular<int64_t, COMP>::axpyin
		(Element& r, const Element& a, const Element& b) const
	{
		return __GIVARO_MODULAR_INTEGER_MULADDIN(r, _p, a, b);
	}
	
	template<typename COMP>
    inline typename Modular<int64_t, COMP>::Element& Modular<int64_t, COMP>::maxpy
		(Element& r, const Element& a, const Element& b, const Element& c) const
	{
		int64_t tmp;
		__GIVARO_MODULAR_INTEGER_MUL(tmp, _p, a, b);
		__GIVARO_MODULAR_INTEGER_SUB(r, _p, c, tmp);
		return r;
	}

	template<typename COMP>
    inline typename Modular<int64_t, COMP>::Element&  Modular<int64_t, COMP>::axmy
		(Element& r, const Element& a, const Element& b, const Element& c) const
	{
		int64_t tmp;
		__GIVARO_MODULAR_INTEGER_MULSUB(tmp, _p, a, b, c);
		return r = (typename Modular<int64_t, COMP>::Element)tmp;
	}

	template<typename COMP>
    inline typename Modular<int64_t, COMP>::Element&  Modular<int64_t, COMP>::maxpyin
		(Element& r, const Element& a, const Element& b) const
	{
		__GIVARO_MODULAR_INTEGER_SUBMULIN(r, _p, a, b );
		return r;
	}

	template<typename COMP>
    inline typename Modular<int64_t, COMP>::Element&  Modular<int64_t, COMP>::axmyin
		(Element& r, const Element& a, const Element& b) const
	{
		return __GIVARO_MODULAR_INTEGER_MULSUB(r, _p, a, b, r );
	}
	
	// ----------------------------------
	// ----- Classic arithmetic on arrays

	template<typename COMP>
    inline void Modular<int64_t, COMP>::mul
		(const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int64_t tmp;
			__GIVARO_MODULAR_INTEGER_MUL(tmp, _p,a[i], b[i]);
			r[i] = (typename Modular<int64_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int64_t, COMP>::mul
		(const size_t sz, Array r, constArray a, Element b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int64_t tmp;
			__GIVARO_MODULAR_INTEGER_MUL(tmp, _p, a[i], b);
			r[i] = (typename Modular<int64_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int64_t, COMP>::div
		(const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			div( r[i], a[i], b[i]);
		}
	}

	template<typename COMP>
    inline void Modular<int64_t, COMP>::div
		(const size_t sz, Array r, constArray a, Element b) const
	{
		typename Modular<int64_t, COMP>::Element ib;
		inv(ib, b);
		mul(sz, r, a, ib);
	}

	template<typename COMP>
    inline void Modular<int64_t, COMP>::add
		(const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int64_t tmp;
			__GIVARO_MODULAR_INTEGER_ADD(tmp, _p, a[i], b[i]);
			r[i] = (typename Modular<int64_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int64_t, COMP>::add
		(const size_t sz, Array r, constArray a, Element b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int64_t tmp;
			__GIVARO_MODULAR_INTEGER_ADD(tmp,_p, a[i], b);
			r[i] = (typename Modular<int64_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int64_t, COMP>::sub
		(const size_t sz, Array r, constArray a, constArray b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int64_t tmp;
			__GIVARO_MODULAR_INTEGER_SUB(tmp, _p, a[i], b[i]);
			r[i] = (typename Modular<int64_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int64_t, COMP>::sub
		(const size_t sz, Array r, constArray a, Element b) const
	{
		for ( size_t i=sz ; --i ; ) {
			int64_t tmp;
			__GIVARO_MODULAR_INTEGER_SUB(tmp, _p, a[i], b);
			r[i] = (typename Modular<int64_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int64_t, COMP>::neg
		(const size_t sz, Array r, constArray a) const
	{
		for ( size_t i=sz ; --i ; ) {
			int64_t tmp;
			__GIVARO_MODULAR_INTEGER_NEG(tmp, _p, a[i]);
			r[i] = (typename Modular<int64_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int64_t, COMP>::axpy
		(const size_t sz, Array r, constArray a, constArray x, constArray y) const
	{
		for ( size_t i=sz ; --i ; ) {
			int64_t tmp;
			__GIVARO_MODULAR_INTEGER_MULADD(tmp, _p, a[i], (int64_t)x[i], (int64_t)y[i]);
			r[i] = (typename Modular<int64_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int64_t, COMP>::axpyin
		(const size_t sz, Array r, constArray a, constArray x) const
	{
		for ( size_t i=sz ; --i ; ) {
			int64_t tmp = (int64_t)r[i];
			__GIVARO_MODULAR_INTEGER_MULADDIN(tmp, _p, a[i], (int64_t)x[i]);
			r[i] = (typename Modular<int64_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int64_t, COMP>::axmy
		(const size_t sz, Array r, constArray a, constArray x, constArray y) const
	{
		for ( size_t i=sz; i--; ) {
			int64_t tmp;
			__GIVARO_MODULAR_INTEGER_MULSUB(tmp, _p, a[i], (int64_t)x[i], (int64_t)y[i]);
			r[i] = (typename Modular<int64_t, COMP>::Element)tmp;
		}
	}

	template<typename COMP>
    inline void Modular<int64_t, COMP>::maxpyin
		(const size_t sz, Array r, constArray a, constArray x) const
	{
		for ( size_t i=sz ; --i ; ) {
			int64_t tmp = (int64_t)r[i];
			__GIVARO_MODULAR_INTEGER_SUBMULIN(tmp, _p, a[i], (int64_t)x[i]);
			r[i] = (typename Modular<int64_t, COMP>::Element)tmp;
		}
	}
	
	// --------------------
	// ----- Initialisation
	
	template<typename COMP>
    inline  typename Modular<int64_t, COMP>::Element&  Modular<int64_t, COMP>::init ( Element& r, const unsigned long a ) const
	{
		return r = (Element)( a >= (uint64_t)_p ? a % (uint64_t)_p : a);
	}

	template<typename COMP>
    inline  typename Modular<int64_t, COMP>::Element&  Modular<int64_t, COMP>::init ( Element& r, const long a ) const
	{
		int64_t sign; uint64_t ua;
		if (a <0) { sign =-1; ua = (uint64_t)-a;}
		else { ua = (uint64_t)a; sign =1; }
		r = (Element)((ua >=_p) ? ua % (uint64_t)_p : ua);
                if (r && (sign ==-1)) r = (Element)_p - r;
		return r;
	}


	template<typename COMP>
    inline typename Modular<int64_t, COMP>::Element&  Modular<int64_t, COMP>::init ( Element& r, const Integer& Residu ) const
	{
		if (Residu <0) {
			int64_t tr;
			// -a = b [p]
			// a = p-b [p]
			if ( (-Residu) >= (Integer)(_p) ) tr = int64_t( (-Residu) % (Integer)_p) ;
			else tr = int64_t(-Residu);
			if (tr) return r = (Element)( (uint64_t)_p - (uint64_t)tr ) ;
			else return r = zero;
		} else {
			Integer ip(_p);
			if (Residu >= ip ) return r =   int64_t(Residu % ip) ;
			else return r = int64_t(Residu);
		}
	}

	template<typename COMP>
    inline  typename Modular<int64_t, COMP>::Element& Modular<int64_t, COMP>::init( Element& a, const int i) const
	{
		return init(a,(long)i);
	}
	
	template<typename COMP>
    inline  typename Modular<int64_t, COMP>::Element& Modular<int64_t, COMP>::init( Element& a, const unsigned int i) const
	{
		return init(a,(unsigned long)i);
	}


	template<typename COMP>
    inline  typename Modular<int64_t, COMP>::Element&  Modular<int64_t, COMP>::init ( Element& r, const unsigned long long a ) const
	{
		return r = (Element)( a >= _p ? a % _p : a);
	}

	template<typename COMP>
    inline  typename Modular<int64_t, COMP>::Element&  Modular<int64_t, COMP>::init ( Element& r, const double a ) const
	{
		return init(r, (int64_t)a);
	}

	template<typename COMP>
    inline  typename Modular<int64_t, COMP>::Element&  Modular<int64_t, COMP>::init ( Element& r, const float a ) const
	{
		return init(r, (double)a);
	}

	template<typename COMP>
    inline  typename Modular<int64_t, COMP>::Element&  Modular<int64_t, COMP>::init ( Element& r, const long long a ) const
	{
		int sign; uint64_t ua;
		if (a <0) { sign =-1; ua = (unsigned int)-a;}
		else { ua = (unsigned int)a; sign =1; }
		r = (Element) ( (ua >=_p) ? ua % (uint64_t)_p : ua) ;
		if (r && (sign ==-1)) r = (Element)_p - r;
                return r;
	}

	template<typename COMP>
    inline void Modular<int64_t, COMP>::init ( const size_t sz, Array r, constArray a ) const
	{
		for ( size_t i = sz; --i; ) {
			r[i] = a[i];
		}
	}

	template<typename COMP>
    inline typename Modular<int64_t, COMP>::Element& Modular<int64_t, COMP>::init ( Element& r ) const
	{
		return r = zero;
	}
	
	// ---------
	// -- Assign

	template<typename COMP>
    inline void Modular<int64_t, COMP>::assign( const size_t sz, Array r, constArray a ) const
	{
		for ( size_t i = sz; --i; )
		    r[i] = a[i];
	}

	template<typename COMP>
    inline typename Modular<int64_t, COMP>::Element&  Modular<int64_t, COMP>::assign
	  ( Element& r, const Element a ) const
	{
		return r = a;
	}


    template<typename COMP>
        inline typename Modular<int64_t, COMP>::Element& Modular<int64_t, COMP>::dotprod
      ( Element& r, const int bound, const size_t sz, constArray a, constArray b ) const
    {
      unsigned int stride = 1;
      if (bound < Signed_Trait<Element>::max() )
       stride = (unsigned int) ( GIVARO_MAXULONG/((unsigned long)bound * (unsigned long)bound) );
      unsigned long dot = (unsigned long) zero; // this is intented !
      if ((sz <10) && (sz <stride)) {
        for(  size_t i= sz; i--; )
    #ifdef __x86_64__
          dot += (unsigned long)a[i] * (unsigned long)b[i];
    #else
          dot = (unsigned long) (dot + a[i] * b[i]);
    #endif
        if (dot > _p) r = (Element)(dot % (uint64_t)_p);
        else r = (Element)dot;
        return r;
      }
      unsigned int i_begin=0;
      stride &= (unsigned int)~0x1;
      if (stride ==0) {
        for(  size_t i= sz; --i; ) {
    #ifdef __x86_64__
          dot += (unsigned long)a[i] * (unsigned long)b[i];
          if (dot>_p) dot %= _p;
    #else
          dot = (unsigned long) (dot + a[i] * b[i]);
          if (dot>_p) dot = (unsigned long) (dot % _p);
    #endif
        }
        r = (Element)dot;
        return r;
      }
      do {
        size_t min_sz = ((sz-i_begin) < stride ? (sz-i_begin) : stride);
        if ((min_sz & 0x1) !=0)
          {
	          --min_sz;
	          ++i_begin;
    #ifdef __x86_64__
	          dot += (unsigned long)a++[min_sz] * (unsigned long)b++[min_sz];
    #else
	          dot = (unsigned long) (dot + a++[min_sz] * b++[min_sz]);
    #endif
          }
        if (min_sz > 1)
          for(  size_t i= min_sz; i>0; --i, --i, ++a, ++a, ++b, ++b ) //!@todo o_O
          {
    #ifdef __x86_64__
            dot += (unsigned long)a[0] * (unsigned long)b[0];
            dot += (unsigned long)a[1] * (unsigned long)b[1];
    #else
	    dot = (unsigned long) (dot +  a[0] * b[0] );
	    dot = (unsigned long) (dot +  a[1] * b[1] );
    #endif

          }
    #ifdef __x86_64__
        if (dot>(uint64_t)_p) dot %= (uint64_t)_p;
    #else
        if (dot>_p) dot = (unsigned long) (dot % _p);
    #endif
        i_begin += (unsigned int) min_sz;
      } while (i_begin <sz);
      r = (Element)dot;
      return r;
    }

    template<typename COMP>
    inline typename Modular<int64_t, COMP>::Element& Modular<int64_t, COMP>::dotprod
      ( Element& r, const size_t sz, constArray a, constArray b ) const
    {
	    return Modular<int64_t, COMP>::dotprod(r, int(_p), sz, a, b);
    }


      //  a -> r: int64_t to double
    template<typename COMP>
    inline void Modular<int64_t, COMP>::i2d ( const size_t sz, double* r, constArray a ) const
    {
      for (size_t i=0; i<sz; ++i)  {
	      r[i] = (double) a[i];
      }
    }

      //  a -> r: double to int64_t
    template<typename COMP>
    inline void Modular<int64_t, COMP>::d2i ( const size_t sz, Array r, const double* a ) const
    {
      union d_2_l {
        double d;
        int64_t r[2];
      };
    //  static const double offset = 4503599627370496.0; // 2^52
      double offset = 4503599627370496.0; // 2^52
      for (size_t i=0; i<sz; ++i)
      {
           d_2_l tmp;
          // - normalization: put fractional part at the end of the representation
          tmp.d = a[i] + offset;
          r[i] = tmp.r[1];
          if (Compute_t(r[i]) < Compute_t(_p))
          	r[i] %= Compute_t(_p);
      }
      //    r[i] = (tmp.r[1] <_p ? tmp.r[1] : tmp.r[1]-_p);
      //    r[i] = (r[i] <_p ? r[i] : r[i]%_p);
      //    r[i] = (tmp.r[1] <_p ? tmp.r[1] : tmp.r[1]%_p);
    }


     // -- Input: (z, <_p>)
    template<typename COMP>
    inline std::istream& Modular<int64_t, COMP>::read (std::istream& s)
    {
      char ch;
      s >> std::ws >> ch;
      if (ch != '(')
        std::cerr << "GivBadFormat(Modular<int64_t, COMP>::read: syntax error: no '('))" << std::endl;

      s >> std::ws >> ch;
      if (ch != 'z')
        std::cerr << "GivBadFormat(Modular<int64_t, COMP>::read: bad domain object))" << std::endl;

      s >> std::ws >> ch;
      if (ch != ',')
        std::cerr << "GivBadFormat(Modular<int64_t, COMP>::read: syntax error: no ',')) " << std::endl;

      s >> std::ws >> _p;


      s >> std::ws >> ch;
      if (ch != ')')
        std::cerr << "GivBadFormat(Modular<int64_t, COMP>::read: syntax error: no ')')) " << std::endl;

      return s;
    }

    template<>
    inline std::ostream& Modular<int64_t, int64_t>::write (std::ostream& s ) const
    {
      	return s << "Modular<int64_t, uint64_t> modulo " << residu();
    }

    template<>
    inline std::ostream& Modular<int64_t, uint64_t>::write (std::ostream& s ) const
    {
      	return s << "Modular<int64_t, uint64_t> modulo " << residu();
    }

    template<typename COMP>
    inline std::istream& Modular<int64_t, COMP>::read (std::istream& s, Element& a) const
    {
	    Integer tmp;
	    s >> tmp;
	    init(a, tmp);
	    return s;
    }

    template<typename COMP>
    inline std::ostream& Modular<int64_t, COMP>::write (std::ostream& s, const Element a) const
    {
      	return s << a;
    }
} // namespace Givaro

#endif // __GIVARO_zpz64std_INL
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
