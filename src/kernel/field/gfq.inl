// ==========================================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// file: gfq.inl
// Description:
//   Arithmetic on GF(q)
// Bugs:
// Authors : JG Dumas
//           Modified 20 Mar 03 by Clement Pernet
// Time-stamp: <09 Jul 08 08:47:17 Jean-Guillaume.Dumas@imag.fr>
// ==========================================================================

#ifndef __GIVARO_gfq_INL
#define __GIVARO_gfq_INL

#include <cmath>
#include <givaro/givinteger.h>
#include <givaro/givintnumtheo.h>
#include <givaro/givpower.h>
#include <givaro/givpoly1padic.h>


// Warning : valid iff b != c
#ifndef __GIVARO_COUNT__

#define _GIVARO_GFQ_ADD(c, a, b, mun, plun) { if ((b)==0) (c)=(a); else if ((a)==0) (c)=(b); else { \
	(c) = (a)-(b); \
	(c) = ((c)>0)?(c):(c)+ (TT)(mun); \
	(c) = (plun)[(UT)(c)]; \
	if (c) { \
		(c) = (c)+(b); \
		(c) = ((c)>0)?(c):(c)+(TT)(mun); \
	} } }

#define _GIVARO_GFQ_NEG(res, a, mo, mun) { if ( (a)==0 ) (res)=0; else\
	{ (res) = (Rep) ( (a) - (Rep) (mo) ) ; (res) = (Rep) ( ((res)>0)?(res):(res)+(Rep)(mun) ) ; } }

// Warning : valid iff a != c
// if not use AUTOSUB ...
#define _GIVARO_GFQ_SUB(c, a, b, mo, mun, plun) { if ((a)==0) {_GIVARO_GFQ_NEG(c,b,mo,mun);} else if ((b)==0) (c)=(a); else { \
	(c) = (b)-(a)-(TT)(mo); \
	(c) = ((c)>0)?(c):(c)+(TT)(mun); \
	(c) = ((c)>0)?(c):(c)+ (TT)(mun); \
	(c) = (plun)[(UT)(c)]; \
	if (c) { \
		(c) = (c)+(a); \
		(c) = ((c)>0)?(c):(c)+(TT)(mun); \
	} } }
#define _GIVARO_GFQ_AUTOSUB(c, b, mo, mun, plun) { if ((c)==0) {_GIVARO_GFQ_NEG(c,b,mo,mun);} else if ((b)!=0) { \
	(c) = (c)-(b)-(TT)(mo); \
	(c) = ((c)>0)?(c):(c)+(TT)(mun); \
	(c) = ((c)>0)?(c):(c)+ (TT)(mun); \
	(c) = (plun)[(UT)(c)]; \
	if (c) { \
		(c) = (c)+(b); \
		(c) = ((c)>0)?(c)-(TT)(mo):(c)+(TT)(mo); \
		(c) = ((c)>0)?(c):(c)+(TT)(mun); \
	} } }



#define _GIVARO_GFQ_MUL(res, a, b, mun) { if ( ((a)==0) || ((b)==0) ) { (res) =0; } else { (res) = (((res) = (a)+(b) )>(TT)(mun))?(res)-(TT)(mun):(res); } }

// JGD 02.04.1998 :  if a==1, a /= a used to be --> 0 !!!
#define _GIVARO_GFQ_INV(res, a, mun)    { (res) = (Rep)( (Rep)(mun)-(a) ); (res)= (Rep) ( (res)?(res):(Rep)(mun) ); }

#define _GIVARO_GFQ_DIV(res, a, b, mun) {  \
	if ( (a)==0 ) { (res)=0; } else { (res) = (((res)=(a)-(b))>0)?(res):(res)+(TT)(mun); } }



#define _GIVARO_GFQ_SQ(res, a, mun) { if ( (a)==0) (res) = 0; else \
	{ (res) = ( (a) << 1) - (mun); \
		(res) = ((res)>0)?(res):(res)+ (mun); } }

// plun -> 1+^c - (q-1) !!!
// Warning : valid iff b != c
#define _GIVARO_GFQ_SQADD(c,a,b,mun,plun) { \
	if ((a)==0) { (c)=(b); \
	} else if ((b)==0) { \
		(c) = ((  (c)=((a) << 1) - (mun)        )>0)?(c):(c) + (mun); \
	} else { \
		(c) = ((    (c) = ((a) << 1)-(b)-(mun)             )<0)?(c)+(mun):(c); \
		if (  (c) = (plun)[(UT)(((c)>0)?(c):(c)+(mun))]     ) { \
			(c) = ((    (c) = (c)+(b)         )>0)?(c):(c)+(mun); } \
	}\
}

// Warning : valid iff b != c
#define _GIVARO_GFQ_MULADD(c,a1,a2,b,mun,plun) { \
	if (((a1)==0) || ((a2)==0)) { (c)=(b); \
	} else if ((b)==0) { \
		(c) = ((    (c)=(a1)+(a2) - (TT)(mun)       )>0)?(c):(c) + (TT)(mun); \
	} else { \
		(c) = ((    (c) = (a1)+(a2)-(b)-(TT)(mun)        )<0)?(c)+(TT)(mun):(c); \
		if (( (c) = (plun)[(UT)( ((c)>0)?(c):(c)+(TT)(mun)   )])  ) { \
			(c) = ((    (c) = (c)+(b)        )>0)?(c):(c)+(TT)(mun); }\
	}\
}

// Warning : valid iff b != c
#define _GIVARO_GFQ_MULSUB(c,a1,a2,b,mo,mun,plun) { \
	if (((a1)==0) || ((a2)==0)) { (c)=(b); \
	} else if ((b)==0) { \
		(c) = ((    (c)=(a1)+(a2) - (mo) -(mun)       )>0)?(c):(c) + (mun); \
		(c) = (c)>0?(c):(c) + (mun); \
	} else { \
		(c) = ((    (c) = (a1)+(a2)-(b)-(mun) - (mo)       )<0)?(c)+(mun):(c); \
		(c) = (c)<0?(c)+(mun):(c); \
		if ( (c) = (plun)[(UT)( ((c)>0)?(c):(c)+(mun)   )]  ) { \
			(c) = ((    (c) = (c)+(b)        )>0)?(c):(c)+(mun); }\
	}\
}



#else


// Warning : valid iff b != c

#define _GIVARO_GFQ_ADD(c, a, b, mun, plun) { ++_add_call; if ((b)==0) (c)=(a); else if ((a)==0) (c)=(b); else { \
	(c) = (a)-(b); \
	(c) = ((c)>0)?(c):(c)+ (mun); \
	(c) = (plun)[(UT)(c)]; \
	if (c) { \
		(c) = (c)+(b); \
		(c) = ((c)>0)?(c):(c)+(mun); \
	} ++_add_count; } }

#define _GIVARO_GFQ_NEG(res, a, mo, mun) { ++_neg_call; if ( (a)==0 ) (res)=0; else\
	{ (res) = (Rep) ((a) - (mo)) ; (res) = (Rep) ( ((res)>0)?(res):(res)+(mun) ); ++_neg_count; } }

// Warning : valid iff a != c
// if not use AUTOSUB ...
#define _GIVARO_GFQ_SUB(c, a, b, mo, mun, plun) { ++_sub_call; if ((a)==0) {_GIVARO_GFQ_NEG(c,b,mo,mun);} else if ((b)==0) (c)=(a); else { \
	(c) = (b)-(a)-(mo); \
	(c) = ((c)>0)?(c):(c)+(mun); \
	(c) = ((c)>0)?(c):(c)+ (mun); \
	(c) = (plun)[(UT)(c)]; \
	if (c) { \
		(c) = (c)+(a); \
		(c) = ((c)>0)?(c):(c)+(mun); \
	} ++_sub_count; } }
#define _GIVARO_GFQ_AUTOSUB(c, b, mo, mun, plun) { ++_sub_call; if ((c)==0) {_GIVARO_GFQ_NEG(c,b,mo,mun);} else if ((b)!=0) { \
	(c) = (c)-(b)-(mo); \
	(c) = ((c)>0)?(c):(c)+(mun); \
	(c) = ((c)>0)?(c):(c)+ (mun); \
	(c) = (plun)[(UT)(c)]; \
	if (c) { \
		(c) = (c)+(b); \
		(c) = ((c)>0)?(c)-(mo):(c)+(mo); \
		(c) = ((c)>0)?(c):(c)+(mun); \
	} ++_sub_count; } }


#define _GIVARO_GFQ_MUL(res, a, b, mun) { ++_mul_call; if ( ((a)==0) || ((b)==0) ) { (res) =0; } else { (res) = (((res) = (a)+(b)-(mun) )>0)?(res):(res)+ (mun); ++_mul_count; } }

// JGD 02.04.1998 :  if a==1, a /= a used to be --> 0 !!!
#define _GIVARO_GFQ_INV(res, a, mun)    { ++_inv_call; (res) = (mun)-(a); (res)=(res)?(res):(mun); ++_inv_count; }

#define _GIVARO_GFQ_DIV(res, a, b, mun) {  ++_div_call; \
	if ( (a)==0 ) { (res)=0; } else { (res) = (((res)=(a)-(b))>0)?(res):(res)+(mun); ++_div_count; } }



#define _GIVARO_GFQ_SQ(res, a, mun) { ++_mul_call; if ( (a)==0) (res) = 0; else \
	{ (res) = ( (a) << 1) - (mun); \
		(res) = ((res)>0)?(res):(res)+ (mun); ++_mul_count; } }

// plun -> 1+^c - (q-1) !!!
// Warning : valid iff b != c
#define _GIVARO_GFQ_SQADD(c,a,b,mun,plun) { ++_mul_call; ++_add_call; \
	if ((a)==0) { (c)=(b); \
	} else if ((b)==0) { \
		(c) = ((  (c)=((a) << 1) - (mun)        )>0)?(c):(c) + (mun); \
		++_mul_count; } else { \
			(c) = ((    (c) = ((a) << 1)-(b)-(mun)             )<0)?(c)+(mun):(c); \
			if (  (c) = (plun)[(UT)(((c)>0)?(c):(c)+(mun))]     ) { \
				(c) = ((    (c) = (c)+(b)         )>0)?(c):(c)+(mun); } \
			++_mul_count; ++_add_count;      }\
}

// Warning : valid iff b != c
#define _GIVARO_GFQ_MULADD(c,a1,a2,b,mun,plun) { ++_mul_call; ++_add_call; \
	if (((a1)==0) || ((a2)==0)) { (c)=(b); \
	} else if ((b)==0) { \
		(c) = ((    (c)=(a1)+(a2) - (mun)       )>0)?(c):(c) + (mun); \
		++_mul_count; } else { \
			(c) = ((    (c) = (a1)+(a2)-(b)-(mun)        )<0)?(c)+(mun):(c); \
			if ( (c) = (plun)[(UT)( ((c)>0)?(c):(c)+(mun)   )]  ) { \
				(c) = ((    (c) = (c)+(b)        )>0)?(c):(c)+(mun); }\
			++_mul_count; ++_add_count;      }\
}

// Warning : valid iff b != c
#define _GIVARO_GFQ_MULSUB(c,a1,a2,b,mo,mun,plun) { ++_mul_call; ++_sub_call; \
	if (((a1)==0) || ((a2)==0)) { (c)=(b); \
	} else if ((b)==0) { \
		(c) = ((    (c)=(a1)+(a2) - (mo) -(mun)       )>0)?(c):(c) + (mun); \
		(c) = (c)>0?(c):(c) + (mun); \
		++_mul_count; ++_neg_count;      } else { \
			(c) = ((    (c) = (a1)+(a2)-(b)-(mun) - (mo)       )<0)?(c)+(mun):(c); \
			(c) = (c)<0?(c)+(mun):(c); \
			if ( (c) = (plun)[(UT)( ((c)>0)?(c):(c)+(mun)   )]  ) { \
				(c) = ((    (c) = (c)+(b)        )>0)?(c):(c)+(mun); }\
			++_mul_count; ++_sub_count;      }\
}


#endif

namespace Givaro {

    template<typename TT>
    inline typename GFqDom<TT>::Residu_t GFqDom<TT>::residu() const
    { return _q; }

    template<typename TT> inline typename GFqDom<TT>::Residu_t GFqDom<TT>::cardinality() const
    { return _q; }
    template<typename TT> inline typename GFqDom<TT>::Residu_t GFqDom<TT>::characteristic() const
    { return _characteristic; }

    template<typename TT>
    inline typename GFqDom<TT>::Residu_t GFqDom<TT>::generator() const
    { return _log2pol[1]; }

    template<typename TT>
    inline typename GFqDom<TT>::Rep& GFqDom<TT>::generator(Rep& g) const
    { return g=1; }

    template<typename TT>
    inline typename GFqDom<TT>::Rep& GFqDom<TT>::indeterminate(Rep& X) const
    {
        if (exponent()>1) {
            return X=(Rep)_pol2log[(size_t)_characteristic];
        } else {
            return X=one;
        }
    }

    template<typename TT>
    inline typename GFqDom<TT>::Rep GFqDom<TT>::indeterminate() const
    {
        Rep X; return indeterminate(X);
    }

    template<typename TT>
    inline typename GFqDom<TT>::Rep GFqDom<TT>::sage_generator() const
    {
        return indeterminate();
    }

    template<typename TT>
    inline typename GFqDom<TT>::Residu_t GFqDom<TT>::irreducible() const
    { return _irred; }

    template<typename TT> inline typename GFqDom<TT>::Residu_t GFqDom<TT>::exponent() const
    { return _exponent; }

    template<typename TT> inline typename GFqDom<TT>::Residu_t GFqDom<TT>::size() const
    { return _q; }



	// ------------------------- Miscellaneous functions

    template<typename TT>
    inline bool GFqDom<TT>::areEqual(const Rep a, const Rep b) const
    { return a == b ; }

    template<typename TT>
    inline bool GFqDom<TT>::areNEqual(const Rep a, const Rep b) const
    { return a != b ; }

    template<typename TT>
    inline bool GFqDom<TT>::isZero(const Rep a) const
    { return a == GFqDom<TT>::zero ; }

    template<typename TT>
    inline bool GFqDom<TT>::isnzero(const Rep a) const
    { return a != GFqDom<TT>::zero ; }

    template<typename TT>
    inline bool GFqDom<TT>::isOne(const Rep a) const
    { return a == GFqDom<TT>::one ; }

	template<typename TT>
    inline bool GFqDom<TT>::isMOne(const Rep a) const
    { return a == GFqDom<TT>::mOne ; }

    template<typename TT>
    inline bool GFqDom<TT>::isUnit(const Rep a) const
    {
            // Fermat : x^(p-1) = 1 whenever x is a unit
        return (a!=0) && (( ( a * (_characteristic-1) ) % _qm1 ) == 0);
    }

    template<typename TT>
    inline size_t GFqDom<TT>::length(const Rep ) const
    { return sizeof(TT) ;}

	// ----------- Usefull method :
    template<typename TT>
    inline typename GFqDom<TT>::Rep&  GFqDom<TT>::mul
    (Rep& r, const Rep a, const Rep b) const
    { _GIVARO_GFQ_MUL(r,a,b, GFqDom<TT>::_qm1) ; return r; }

    template<typename TT>
    inline typename GFqDom<TT>::Rep&  GFqDom<TT>::mulin
    (Rep& r, const Rep a) const
    { _GIVARO_GFQ_MUL(r,r,a, GFqDom<TT>::_qm1) ; return r; }

    template<typename TT>
    inline typename GFqDom<TT>::Rep&  GFqDom<TT>::div
    (Rep& r, const Rep a, const Rep b) const
    {
	   	_GIVARO_GFQ_DIV(r, a, b, GFqDom<TT>::_qm1) ;
		return r;
	}

    template<typename TT>
    inline typename GFqDom<TT>::Rep&  GFqDom<TT>::divin
    (Rep& r, const Rep a) const
    { _GIVARO_GFQ_DIV(r, r, a, GFqDom<TT>::_qm1) ; return r; }

    template<typename TT>
    inline typename GFqDom<TT>::Rep&  GFqDom<TT>::add
    (Rep& r, const Rep a, const Rep b) const
    { _GIVARO_GFQ_ADD(r, a, b, GFqDom<TT>::_qm1, GFqDom<TT>::_plus1) ; return r; }

    template<typename TT>
    inline typename GFqDom<TT>::Rep&  GFqDom<TT>::addin
    (Rep& r, const Rep a) const
    { _GIVARO_GFQ_ADD(r, r, a, GFqDom<TT>::_qm1, GFqDom<TT>::_plus1) ; return r; }

    template<typename TT>
    inline typename GFqDom<TT>::Rep&  GFqDom<TT>::sub
    (Rep& r, const Rep a, const Rep b) const
    { _GIVARO_GFQ_SUB(r, a, b, GFqDom<TT>::mOne, GFqDom<TT>::_qm1, GFqDom<TT>::_plus1) ; return r; }

    template<typename TT>
    inline typename GFqDom<TT>::Rep&  GFqDom<TT>::subin
    (Rep& r, const Rep a) const
    { _GIVARO_GFQ_AUTOSUB(r, a, GFqDom<TT>::mOne, GFqDom<TT>::_qm1, GFqDom<TT>::_plus1) ; return r; }

    template<typename TT>
    inline typename GFqDom<TT>::Rep&  GFqDom<TT>::neg
    (Rep& r, const Rep a) const
    { _GIVARO_GFQ_NEG(r, a, GFqDom<TT>::mOne, GFqDom<TT>::_qm1) ; return r; }

    template<typename TT>
    inline typename GFqDom<TT>::Rep&  GFqDom<TT>::negin
    (Rep& r) const
    { _GIVARO_GFQ_NEG(r, r, GFqDom<TT>::mOne, GFqDom<TT>::_qm1) ; return r; }

    template<typename TT>
    inline typename GFqDom<TT>::Rep&  GFqDom<TT>::inv
    (Rep& r, const Rep a) const
    { _GIVARO_GFQ_INV(r, a, GFqDom<TT>::_qm1) ; return r; }

    template<typename TT>
    inline typename GFqDom<TT>::Rep&  GFqDom<TT>::invin
    (Rep& r) const
    { _GIVARO_GFQ_INV(r, r, GFqDom<TT>::_qm1) ; return r; }

    template<typename TT>
    inline typename GFqDom<TT>::Rep&  GFqDom<TT>::axpy
    (Rep& r, const Rep a, const Rep b, const Rep c)
	const
    { _GIVARO_GFQ_MULADD(r,a,b,c, GFqDom<TT>::_qm1, GFqDom<TT>::_plus1) ; return r; }

    template<typename TT>
    inline typename GFqDom<TT>::Rep&  GFqDom<TT>::axpyin
    (Rep& r, const Rep a, const Rep b) const
    {
        Rep tmp = r;
        _GIVARO_GFQ_MULADD((r),a,b,tmp, (GFqDom<TT>::_qm1), (GFqDom<TT>::_plus1)) ;
        return r;
	}

        // r <- r-a*b
    template<typename TT>
    inline typename GFqDom<TT>::Rep&  GFqDom<TT>::maxpyin (Rep& r,
                                                           const Rep a,
                                                           const Rep b) const
    {
            //   Rep tmp = r;
            //   _GIVARO_GFQ_MULSUB(r,a,b,tmp, mOne, _qm1, _plus1) ;
        Rep tmp; _GIVARO_GFQ_MUL(tmp,a,b, _qm1) ;
        _GIVARO_GFQ_AUTOSUB(r,tmp, mOne, _qm1, _plus1) ;
        return r;
    }

    template<typename TT>
    inline typename GFqDom<TT>::Rep&  GFqDom<TT>::axmyin (Rep& r,
                                                          const Rep a,
                                                          const Rep b) const
    {
        this->maxpyin(r,a,b);
        return this->negin(r);
    }

        // r <- a*b-c
    template<typename TT>
    inline typename GFqDom<TT>::Rep&  GFqDom<TT>::axmy
    (Rep& r, const Rep a, const Rep b, const Rep c)
        const
    {
        _GIVARO_GFQ_MUL(r,a,b, GFqDom<TT>::_qm1) ;
        _GIVARO_GFQ_AUTOSUB(r,c, GFqDom<TT>::mOne, GFqDom<TT>::_qm1, GFqDom<TT>::_plus1) ;
        return r; }

    template<typename TT>
    inline typename GFqDom<TT>::Rep&  GFqDom<TT>::maxpy
    (Rep& r, const Rep a, const Rep b, const Rep c)
        const
    {
        _GIVARO_GFQ_MUL(r,a,b, GFqDom<TT>::_qm1) ;
        _GIVARO_GFQ_SUB(r,c,r, GFqDom<TT>::mOne, GFqDom<TT>::_qm1, GFqDom<TT>::_plus1) ;
        return r; }


        // -- inline array operations between Reps
    template<typename TT>
    inline void GFqDom<TT>::mul
    (const size_t sz, Array r, constArray a, constArray b) const
    {
        for ( size_t i=sz ; --i ; ) {
            _GIVARO_GFQ_MUL(r[i],a[i],b[i], GFqDom<TT>::_qm1) ;
        }
    }

    template<typename TT>
    inline void GFqDom<TT>::mul
    (const size_t sz, Array r, constArray a, const Rep b) const
    {
        for ( size_t i=sz ; --i ; ) {
            _GIVARO_GFQ_MUL(r[i],a[i],b, GFqDom<TT>::_qm1) ;
        }
    }

    template<typename TT>
    inline void GFqDom<TT>::div
    (const size_t sz, Array r, constArray a, constArray b) const
    {
        for ( size_t i=sz ; --i ; ) {
            _GIVARO_GFQ_DIV(r[i],a[i],b[i], GFqDom<TT>::_qm1) ;
        }
    }

    template<typename TT>
    inline void GFqDom<TT>::div
    (const size_t sz, Array r, constArray a, const Rep b) const
    {
        for ( size_t i=sz ; --i ; ) {
            _GIVARO_GFQ_DIV(r[i],a[i],b, GFqDom<TT>::_qm1) ;
        }
    }

    template<typename TT>
    inline void GFqDom<TT>::add
    (const size_t sz, Array r, constArray a, constArray b) const
    {
        for ( size_t i=sz ; --i ; ) {
            _GIVARO_GFQ_ADD(r[i], a[i], b[i], GFqDom<TT>::_qm1, GFqDom<TT>::_plus1) ;
        }
    }

    template<typename TT>
    inline void GFqDom<TT>::add
    (const size_t sz, Array r, constArray a, const Rep b) const
    {
        for ( size_t i=sz ; --i ; ) {
            _GIVARO_GFQ_ADD(r[i], a[i], b, GFqDom<TT>::_qm1, GFqDom<TT>::_plus1) ;
        }
    }

    template<typename TT>
    inline void GFqDom<TT>::sub
    (const size_t sz, Array r, constArray a, constArray b) const
    {
        for ( size_t i=sz ; --i ; ) {
            _GIVARO_GFQ_SUB(r[i], a[i], b[i], GFqDom<TT>::mOne, GFqDom<TT>::_qm1, GFqDom<TT>::_plus1) ;
        }
    }

    template<typename TT>
    inline void GFqDom<TT>::sub
    (const size_t sz, Array r, constArray a, const Rep b) const
    {
        for ( size_t i=sz ; --i ; ) {
            _GIVARO_GFQ_SUB(r[i], a[i], b, GFqDom<TT>::mOne, GFqDom<TT>::_qm1, GFqDom<TT>::_plus1) ;
        }
    }

    template<typename TT>
    inline void GFqDom<TT>::neg
    (const size_t sz, Array r, constArray a) const
    {
        for ( size_t i=sz ; --i ; ) {
            _GIVARO_GFQ_NEG(r[i], a[i], GFqDom<TT>::mOne, GFqDom<TT>::_qm1) ;
        }
    }

    template<typename TT>
    inline void GFqDom<TT>::inv
    (const size_t sz, Array r, constArray a) const
    {
        for ( size_t i=sz ; --i ; ) {
            _GIVARO_GFQ_INV(r[i], a[i], GFqDom<TT>::_qm1) ;
        }
    }

    template<typename TT>
    inline void GFqDom<TT>::axpy
    (const size_t sz, Array r, const Rep a, constArray x, constArray y) const
    {
        for ( size_t i=sz ; --i ; ) {
            _GIVARO_GFQ_MULADD(r[i], a, x[i], y[i], GFqDom<TT>::_qm1, GFqDom<TT>::_plus1) ;
        }
    }

    template<typename TT>
    inline void GFqDom<TT>::axpyin
    (const size_t sz, Array r, const Rep a, constArray x) const
    {
        Rep tmp;
        for ( size_t i=sz ; --i ; ) {
            tmp = r[i];
            _GIVARO_GFQ_MULADD(r[i], a, x[i], tmp, GFqDom<TT>::_qm1, GFqDom<TT>::_plus1) ;
        }
    }

    template<typename TT>
    inline void GFqDom<TT>::axpy
    (const size_t sz, Array r, const Rep a, constArray x, const Rep y) const
    {
        for ( size_t i=sz ; --i ; ) {
            _GIVARO_GFQ_MULADD(r[i], a, x[i], y, GFqDom<TT>::_qm1, GFqDom<TT>::_plus1) ;
        }
    }

    template<typename TT>
    inline void GFqDom<TT>::axmy
    (const size_t sz, Array r, const Rep a, constArray x, constArray y) const
    {
        for ( size_t i=sz ; --i ; ) {
            _GIVARO_GFQ_MUL(r[i], a, x[i], GFqDom<TT>::_qm1) ;
            _GIVARO_GFQ_AUTOSUB(r[i], y[i], GFqDom<TT>::mOne, GFqDom<TT>::_qm1, GFqDom<TT>::_plus1) ;
        }
    }

    template<typename TT>
    inline void GFqDom<TT>::axmy
    (const size_t sz, Array r, const Rep a, constArray x, const Rep y) const
    {
        for ( size_t i=sz ; --i ; ) {
            _GIVARO_GFQ_MUL(r[i], a, x[i], GFqDom<TT>::_qm1) ;
            _GIVARO_GFQ_AUTOSUB(r[i], y, GFqDom<TT>::mOne, GFqDom<TT>::_qm1, GFqDom<TT>::_plus1) ;
        }
    }

    template<typename TT>
    inline void GFqDom<TT>::maxpyin (const size_t sz, Array r,
                                     const Rep a, constArray x) const
    {
        Rep tmp;
        for ( size_t i=sz ; --i ; ) {
            _GIVARO_GFQ_MUL(tmp, a, x[i], GFqDom<TT>::_qm1) ;
            _GIVARO_GFQ_AUTOSUB(r[i], tmp, GFqDom<TT>::mOne, GFqDom<TT>::_qm1, GFqDom<TT>::_plus1) ;
        }
    }

        // ------------------------------------
        // Input - Output  of the Domain
        //
    template<typename TT>
    inline std::istream& GFqDom<TT>::read (std::istream& s) {
        char ch;
        s >> std::ws >> ch;
        if (ch != '(')
            std::cerr << "GFqDom::read: syntax error: no '('" << std::endl;
        UTT p;
        s >> p;
        s >> std::ws >> ch;
        if (ch == ')')
            *this = GFqDom<TT>(p,UTT(1));
        else {
            if (ch != '^')
                std::cerr << "GFqDom::read: syntax error: no '^'" << std::endl;
            UTT k;
            s >> std::ws >> k;
            if (ch != ')')
                std::cerr << "GFqDom::read: syntax error: no ')'" << std::endl;

            *this = GFqDom<TT>(p,k); // Seems like a useless copy... better have a reinit function.
        }
        return s;
    }

    template<typename TT>
    inline std::ostream& GFqDom<TT>::write (std::ostream& o) const
    {
        return o << "GFqDom<> (" <<  _characteristic << '^' << _exponent << ")";
    }
    template<typename TT>
    inline std::ostream& GFqDom<TT>::write (std::ostream& o, const std::string& s) const
    {
		return this->write(o) << s;
    }

        // ------------------------------------
        // Input - Output  of the Elements
        //
    template<typename TT>
    inline std::istream& GFqDom<TT>::read (std::istream& i, Rep& a) const
    {
        TT t;
        i >> t;
        init(a,t);
        return i;
    }

    template<typename TT>
    inline typename GFqDom<TT>::Rep& GFqDom<TT>::reduce( Rep& r, const Rep e) const {
            return r = e;
    }
    template<typename TT>
    inline typename GFqDom<TT>::Rep& GFqDom<TT>::reduce( Rep& r ) const {
        return r;
    }
       

 //    template<typename TT>
//     template<typename XXX>
//     inline typename GFqDom<TT>::Element& GFqDom<TT>::init(Element& r, const XXX& value) const {
//         return r = (Rep)_pol2log[ (UT)(value) ];
//     }
        

    template<typename TT>
    inline typename GFqDom<TT>::Rep& GFqDom<TT>::init( Rep& r, const double Residu ) const
    {
        double tr = Residu ;
        if (tr <0) {
                // -a = b [p]  <==>  a = p-b [p]
            tr = -tr;
            if (tr > Signed_Trait<UTT>::max() )
                tr = fmod(tr,(double)_q);
                //tr -= (double)floor(tr * _inversecharacteristic)*_dcharacteristic;
            else{
                if (tr >= (TT)_q )
                    tr = double((UTT)tr % _q) ;
            }

            if ((bool)tr)
                return r = (Rep)_pol2log[ UT(_q - (UTT)tr) ];
            else
                return r = zero;
        } else {
            if (tr > Signed_Trait<UTT>::max() )
                tr = fmod(tr, (double)_q);
                //tr -= (double)floor(tr * _inversecharacteristic)*_dcharacteristic;
            else{
                if (tr >= (TT)_q )
                    tr = double((UTT)tr % _q) ;
            }
            return r = (Rep)_pol2log[ (UT)tr ];
        }
    }

    template<typename TT>
    inline typename GFqDom<TT>::Rep& GFqDom<TT>::init( Rep& r, const float Residu ) const
    {
        return init(r, static_cast<double>(Residu));
    }



     template<typename TT>
     inline typename GFqDom<TT>::Rep& GFqDom<TT>::init( Rep& r, const int32_t Residu ) const
     {
         int32_t tr = Residu ;
         if (tr <0) {
                 // -a = b [p]
                 // a = p-b [p]
             tr = -tr;
             if (tr >= (int32_t)_q )
                 tr =(int32_t)( (UT)tr % _q ) ;
             if (tr)
                 return r = (Rep) _pol2log[(UT) _q - (UT)tr ];
             else
                 return r = zero;
         }
         else {
             if (tr >= (int32_t)_q )
                 tr = int32_t((uint32_t)tr % _q ) ;
             return r = (Rep)_pol2log[ (UT)tr ];
         }
     }
     template<typename TT>
     inline typename GFqDom<TT>::Rep& GFqDom<TT>::init( Rep& r, const int64_t Residu ) const
     {
         int64_t tr = Residu ;
         if (tr <0) {
                 // -a = b [p]
                 // a = p-b [p]
             tr = -tr;
             if (tr >= (int64_t)_q )
 				tr = tr % (int64_t)_q ;
             if (tr)
                 return r = (typename GFqDom<TT>::Rep) _pol2log[ (size_t)_q - (size_t)tr ];
             else
                 return r = zero;
         } else {
             if (tr >= (int64_t)_q )
 				tr = tr % (int64_t)_q ;
             return r = (Rep)_pol2log[ (size_t)tr ];
         }
     }

    template<typename TT>
    inline typename GFqDom<TT>::Rep& GFqDom<TT>::init( Rep& r, const Integer Residu ) const
    {
        UTT tr;
        if (Residu <0) {
                // -a = b [p]
                // a = p-b [p]
            if ( Residu <= (Integer)(-_q) )
				tr = (UTT) (  (-Residu) % (UTT)_q );
            else
                tr = UTT(-Residu);
            if (tr)
                return r = (Rep)_pol2log[ _q - (UTT)tr ];
            else
                return r = zero;
        }
		else { /* Residu >=0 */
            if (Residu >= (Integer)_q )
				tr =  (UTT)(Residu % (UTT)_q );
            else
				tr = UTT(Residu);
            return r = (Rep)_pol2log[ (size_t)tr ];
        }
    }

     template<typename TT>
     inline typename GFqDom<TT>::Rep& GFqDom<TT>::init( Rep& r, const uint64_t Residu ) const
     {
         uint64_t tr = Residu ;
         if (tr >= _q )
 			tr =tr %  (uint64_t) _q ;
         return r = (Rep)_pol2log[ (size_t)tr ];
     }
     template<typename TT>
     inline typename GFqDom<TT>::Rep& GFqDom<TT>::init( Rep& r, const uint32_t Residu ) const
     {
         uint64_t tr = static_cast<uint64_t>(Residu) ;
         if (tr >= _q ) tr = tr % (uint64_t) _q ;
         return r = (Rep)_pol2log[ (size_t)tr ];
     }


    template<typename TT>
    inline double& GFqDom<TT>::convert (double& r, const Rep a) const
    {
        return r = (double)_log2pol[ (UT)a] ;
    }

    template<typename TT>
    inline float& GFqDom<TT>::convert (float& r, const Rep a) const
    {
        return r = (float)_log2pol[ (UT)a] ;
    }

    template<typename TT>
    inline std::ostream& GFqDom<TT>::write (std::ostream& o, const Rep a) const
    {
        return o << _log2pol[ (UT)a] ;
    }



    template<typename TT>
    inline int64_t& GFqDom<TT>::convert (int64_t& r, const Rep a) const
    {
        return r = (int64_t)_log2pol[ (uint64_t)a] ;
    }

    template<typename TT>
    inline uint64_t& GFqDom<TT>::convert (uint64_t& r, const Rep a) const
    {
        return r = (uint64_t)_log2pol[ (uint64_t)a] ;
    }

    template<typename TT>
    inline int32_t& GFqDom<TT>::convert (int32_t& r, const Rep a) const
    {
        return r = (int32_t)_log2pol[ (UT)a] ;
    }

    template<typename TT>
    inline uint32_t& GFqDom<TT>::convert (uint32_t& r, const Rep a) const
    {
        return r = (uint32_t)_log2pol[ (UT)a] ;
    }

    template<typename TT>
    inline TT GFqDom<TT>::convert (const Rep a) const
    {
        return (TT)_log2pol[ (UT)a] ;
    }

    template<typename TT>
    inline Integer& GFqDom<TT>::convert (Integer& r, const Rep a) const
    {
        return r = (Integer)_log2pol[ (UT)a] ;
    }


        // ---------
        // -- Initialization operations
        // ---------
    template<typename TT>
    template<typename val_t, template<class, class> class Vector, template <class> class Alloc>
    inline typename GFqDom<TT>::Rep& GFqDom<TT>::init( Rep& r, const Vector<val_t, Alloc<val_t> >& P) const {
        static Self_t PrimeField(this->_characteristic);
        typedef Poly1Dom< Self_t, Dense > PolDom;
        static PolDom Pdom( PrimeField );
        typedef Poly1PadicDom< GFqDom<TT>, Dense > PadicDom;
        static PadicDom PAD(Pdom);
        Degree d;  Pdom.degree(d, P);
        if (d >= (int64_t)this->_exponent) {
            static typename PadicDom::Element tmp;
            static typename PadicDom::Element Irreducible = PAD.radix(tmp, this->_irred);
                // All this was to get the irreducible polynomial
                // Now we can mod it out
            typename PolDom::Element modP; Pdom.mod(modP, P, Irreducible);
            TT tr;
            PAD.eval(tr, modP);
            return r = (Rep) this->_pol2log[(size_t) tr ];
        } else {
            TT tr;
            PAD.eval(tr, P);
            return r = (Rep) this->_pol2log[ (size_t)tr ];
        }
    }


    template<typename TT>
    inline typename GFqDom<TT>::Rep& GFqDom<TT>::assign( Rep& r, const Integer a) const
    { return init (r, a); }

    template<typename TT>
    inline typename GFqDom<TT>::Rep& GFqDom<TT>::assign( Rep& r, const Rep a) const
    { return r = a; }


    template<typename TT> inline void GFqDom<TT>::assign( const size_t sz, Array r, constArray a ) const
    {
        TT tr;
            //    for ( size_t i=sz ; --i ; )
        for ( size_t i=sz; i--;) {
            tr = a[i] ;
            if (tr <0) {
                    // -a = b [p]
                    // a = p-b [p]
                tr = -tr;
                if (tr >=_characteristic ) tr = tr % _characteristic ;
                if (tr)
                    r[i] = _pol2log[ _characteristic - tr ];
                else
                    r[i] = 0;
            } else {
                if (tr >=_characteristic ) tr = tr % _characteristic ;
                r[i] = _pol2log[ tr ];
            }
        }
    }

    template<typename TT>
    inline typename  GFqDom<TT>::Rep& GFqDom<TT>::dotprod
    ( Rep& r, const size_t sz, constArray a, constArray b ) const
    {
        if (sz) {
            _GIVARO_GFQ_MUL(r,a[0],b[0],_qm1);
            Rep tmp;
            for(  int i= (int)sz; --i; ) {
                _GIVARO_GFQ_MUL(tmp,a[i],b[i],_qm1);
                _GIVARO_GFQ_ADD(r,r,tmp,_qm1,_plus1);
            }
            return r;
        } else
            return r = zero;
    }


        // ----- random generators
    template<typename TT> template<typename randIter>
    inline typename GFqDom<TT>::Rep& GFqDom<TT>::nonzerorandom(randIter& g, Rep& a, const Residu_t& s) const
    {
        a = Rep( ((UTT)(g()) % (s-1)) + 1);
        return a = (a<0?a+(Rep)_q:a); //@fixme can it really be <0?
    }

    template<typename TT> template<typename randIter>
    inline typename GFqDom<TT>::Rep& GFqDom<TT>::random(randIter& g, Rep& a, const Residu_t& s) const
    {
        a = Rep( (UTT)(g()) % s);
        return a = (a<0?a+(Rep)_q:a); //@fixme can it really be <0?
    }

    template<typename TT> template<typename randIter>
    inline typename GFqDom<TT>::Rep& GFqDom<TT>::random(randIter& g, Rep& r) const
    {
        return random(g,r,_q);
    }

    template<typename TT> template<typename randIter>
    inline typename GFqDom<TT>::Rep& GFqDom<TT>::nonzerorandom(randIter& g, Rep& r) const
    {
        return nonzerorandom(g,r,_q);
    }

    template<typename TT>
    inline GFqDom<TT>::GFqDom(const UTT P, const UTT e)
            // Precondition P prime
            :  zero(0)
        , one ( (TT)power(P,e) - 1  )
        , mOne(  (P==2)?  (one)  :  (one >> 1) )   // 1 == -1 in GF(2^k)
        , _characteristic(P)
        , _exponent(e)
        , _q( (UTT) one + 1 )
        , _qm1 ( (UTT) one )
        , _log2pol((UT) _q )
        , _pol2log( (UT)_q )
        , _plus1( (UT)_q )
        , _dcharacteristic( (double)P )
    {

            // 1 is represented by q-1, zero by 0
        _log2pol[0] = (UTT) zero;

        if (e <= 1) {
            IntNumTheoDom<> NTD;
            IntNumTheoDom<>::Rep IP(P), pr;
                //         UTT seed = (UTT) ( NTD.Integer2long( NTD.lowest_prim_root(pr, IP) ) );
            UTT seed;
            NTD.convert(seed, NTD.lowest_prim_root(pr, IP) );
            UTT accu = 1;
            for(UTT i=1; i<P; i++) {
                accu = (accu * seed) % P;
                _log2pol[(UT)i] = accu;
            }
        } else {
                // Fisrt compute an irreductible polynomial F over Z/pZ of degree e
                // Then a primitive root G (i.e. a generator of GF(q))
            GFqDom<TT> Zp(P,1);
                //         typedef CyclotomicTable<  GFqDom<TT>, Dense > PolDom;
                //         PolDom Pdom( Zp, e );
            typedef Poly1FactorDom< Self_t, Dense > PolDom;
            PolDom Pdom( Zp );
            typename PolDom::Element F, G, H;

                // F is irreducible of degree e over Zp
                // G is a primitive polynomial for F
                //         Pdom.random_prim_root(F,G, Degree(e));

                // F is an irreducible factor of the
                // (p^e-1) th cyclotomic polynomial
                // G is a primitive polynomial for F : X
                //         Pdom.getcyclo(F);
                //         Pdom.init(G, Degree(1), Zp.one);

                // F is irreducible of degree e over Zp
                // with X as a primitive polynomial
#ifndef GIVARO_RANDOM_IRREDUCTIBLE_PRIMITIVE_ROOT
            Pdom.ixe_irreducible(F, Degree((long)e));
                //         Pdom.init(G, Degree(1), Zp.one);
                //         Pdom.assign(G, Degree(1), Zp.one);
            Pdom.init(G, Degree(1));
#else
            Pdom.random_irreducible(F, Degree((int64_t)e));
            Pdom.give_random_prim_root(G,F);
#endif

            Pdom.assign(H, G);

            typedef Poly1PadicDom< GFqDom<TT>, Dense > PadicDom;
            PadicDom PAD(Pdom);

            PAD.eval(_log2pol[1], H);
            PAD.eval(_irred, F);

            for (UTT i = 2; i < _qm1; ++i) {
                Pdom.mulin(H, G);
                Pdom.modin(H, F);
                PAD.eval(_log2pol[(UT)i], H);
            }

            _log2pol[(UT)_qm1] = 1;

        }

        _log2pol[0] = 0;

            // pol2log[ j ] = i such that log2pol[i] = j
        for (UTT i = 0; i < _q; ++i)
            _pol2log[ (UT)_log2pol[(UT)i] ] = i;

            // plus1[i] = k such that G^i + 1 = G^k
            // WARNING : in the plus1 table, we now pre-substract (_q - 1)
        _plus1[0] = 0;

        UTT a,b,r;
        for (UTT i = 1; i < _q; ++i) {
            a = _log2pol[(UT)i];
            r = a % P;
            if (r == (P - 1))
                b = a - r;
            else
                b = a + 1;
                // WARNING : in the plus1 table we pre-substract (_q - 1)
            _plus1[(UT)i] = (TT)(_pol2log[(UT)b] - _qm1);
        }
            // -1 + 1 == 0
        _plus1[(UT)mOne] = 0;
    }

        // Dan Roche 6-15-04, adapted my/ported back to Givaro
        // by Martin Albrecht 10-06-06
        // This constructor takes a vector of ints that represent the polynomial
        // to use (for modular arithmetic on the extension field).
    template<typename TT>
    template<typename Vector>
    inline GFqDom<TT>::GFqDom(const UTT P, const UTT e, const Vector& modPoly):
            zero(0)
        , one ((TT) power(P,e) - 1  )
        , mOne(  (P==2)?  (one)  :  ( one >> 1) )   // 1 == -1 in GF(2^k)
        , _characteristic(P)
        , _exponent(e)
        , _q( (UTT) one + 1 )
        , _qm1 ( (UTT)one )
        , _log2pol( (UT)_q )
        , _pol2log( (UT)_q )
        , _plus1( (UT)_q )
        , _dcharacteristic( (double)P )
    {

            // 1 is represented by q-1, zero by 0
        _log2pol[0] = (UTT)zero;

        GFqDom<TT> Zp(P,1);
        typedef Poly1FactorDom< GFqDom<TT>, Dense > PolDom;
        PolDom Pdom( Zp );
        typename PolDom::Element Ft, F(e+1), G, H;

        for( size_t i = 0; i < F.size(); ++i )
            Zp.init( F[i], modPoly[i]);

        Pdom.give_prim_root(G,F);
        Pdom.assign(H,G);

        typedef Poly1PadicDom< GFqDom<TT>, Dense > PadicDom;
        PadicDom PAD(Pdom);

        PAD.eval(_log2pol[1],H);
        PAD.eval(_irred, F);

        for (UTT i = 2; i < _qm1; ++i) {
            Pdom.mulin(H, G);
            Pdom.modin(H, F);
            PAD.eval(_log2pol[i], H);
        }
        _log2pol[_qm1] = 1;

        _log2pol[0] = 0;

        for (UTT i = 0; i < _q; ++i)
            _pol2log[ _log2pol[i] ] = i;

        _plus1[0] = 0;

        UTT a,b,r;
        for (UTT i = 1; i < _q; ++i) {
            a = _log2pol[i];
            r = a % P;
            if (r == (P - 1))
                b = a - r;
            else
                b = a + 1;
            _plus1[i] = (TT)_pol2log[b] - (TT)_qm1;
        }

        _plus1[size_t(mOne)] = 0;
    }

         // Construction with prescribed irreducible polynomial
        //   and with prescribed generator polynomial
        //   coefficients of the vector should be integers-like
        //   there will be a call to this->init to build the
        //   representation of both polynomials
    template<typename TT>
    template<typename Vector>
    inline GFqDom<TT>::GFqDom(const UTT P, const UTT e, 
                              const Vector& modPoly, const Vector& genPoly):
            zero(0)
        , one ((TT) power(P,e) - 1  )
        , mOne(  (P==2)?  (one)  :  ( one >> 1) )   // 1 == -1 in GF(2^k)
        , _characteristic(P)
        , _exponent(e)
        , _q( (UTT) one + 1 )
        , _qm1 ( (UTT)one )
        , _log2pol( (UT)_q )
        , _pol2log( (UT)_q )
        , _plus1( (UT)_q )
        , _dcharacteristic( (double)P )
    {

            // 1 is represented by q-1, zero by 0
        _log2pol[0] = (UTT)zero;

        GFqDom<TT> Zp(P,1);
        typedef Poly1FactorDom< GFqDom<TT>, Dense > PolDom;
        PolDom Pdom( Zp );
        typename PolDom::Element Ft, F(e+1), G(genPoly.size()), H;

        for( size_t i = 0; i < F.size(); ++i )
            Zp.init( F[i], modPoly[i]);

        for( size_t i = 0; i < G.size(); ++i )
            Zp.init( G[i], genPoly[i]);

        Pdom.assign(H,G);

        typedef Poly1PadicDom< GFqDom<TT>, Dense > PadicDom;
        PadicDom PAD(Pdom);

        PAD.eval(_log2pol[1],H);
        PAD.eval(_irred, F);

        for (UTT i = 2; i < _qm1; ++i) {
            Pdom.mulin(H, G);
            Pdom.modin(H, F);
            PAD.eval(_log2pol[i], H);
        }
        _log2pol[_qm1] = 1;

        _log2pol[0] = 0;

        for (UTT i = 0; i < _q; ++i)
            _pol2log[ _log2pol[i] ] = i;

        _plus1[0] = 0;

        UTT a,b,r;
        for (UTT i = 1; i < _q; ++i) {
            a = _log2pol[i];
            r = a % P;
            if (r == (P - 1))
                b = a - r;
            else
                b = a + 1;
            _plus1[i] = (TT)_pol2log[b] - (TT)_qm1;
        }

        _plus1[size_t(mOne)] = 0;
    }

    template<typename TT> inline void GFqDom<TT>::Init() {}

    template<typename TT> inline void GFqDom<TT>::End() {}


#ifdef __GIVARO_COUNT__
    template<typename TT> int64_t GFqDom<TT>::_mul_count = 0;
    template<typename TT> int64_t GFqDom<TT>::_add_count = 0;
    template<typename TT> int64_t GFqDom<TT>::_div_count = 0;
    template<typename TT> int64_t GFqDom<TT>::_sub_count = 0;
    template<typename TT> int64_t GFqDom<TT>::_neg_count = 0;
    template<typename TT> int64_t GFqDom<TT>::_inv_count = 0;
    template<typename TT> int64_t GFqDom<TT>::_mul_call = 0;
    template<typename TT> int64_t GFqDom<TT>::_add_call = 0;
    template<typename TT> int64_t GFqDom<TT>::_div_call = 0;
    template<typename TT> int64_t GFqDom<TT>::_sub_call = 0;
    template<typename TT> int64_t GFqDom<TT>::_neg_call = 0;
    template<typename TT> int64_t GFqDom<TT>::_inv_call = 0;
#endif

} // namespace Givaro

#endif // __GIVARO_gfq_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
