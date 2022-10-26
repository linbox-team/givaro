// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/system/givconfig.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givconfig.h,v 1.22 2011-02-04 14:50:07 briceboyer Exp $
// ==========================================================================
/** @file givconfig.h
 * @ingroup system
 * @brief configuration file for Givaro
 */
#ifndef _GIVARO_INTERNAL_CONFIG_H_
#define _GIVARO_INTERNAL_CONFIG_H_ 1


// * GIVARO_HAVE_ANSI_NAMESPACE:
//   value: defined/undefined
//   purpose: define it in order to use name space in place of prefix
// * GIVARO_HAVE_ANSI_SPECIALIZED:
//   value: defined/undefined
//   purpose: define it in order to use ANSI declaration for specialization
//   of template functions.
// * GIVARO_HAVE_ANSI_EXCEPTION:
//   value: defined/undefined
//   purpose: define it in order to use ANSI exception classes hierarchy
//   as based classes for exceptions.
// * GIVARO_HAVE_ANSI_LIBRARY:
//   value: defined/undefined
//   purpose: define it in order to use ANSI library utilities (except for
//   container) as allocation without throwing exception, new str stream,
//   std header files...
// * GIVARO_HAVE_LBLAS:
//   value: defined/undefined
//   purpose: define it in order to use blas numerical library for some operations
//   as floating point matrice operations.
// * GIVARO_HAVE_TYPENAME:
//   value: defined/undefined
//   purpose: define it if the C++ compiler have typename keyword.
// * GIVARO_HAVE_LONG_LONG:
//   value: defined/undefined
//   purpose: define it if the C++ compiler have long long definition.
// * GIVARO_ASSERT_MACRO:
//   value: defined/undefined
//   purpose: expand additional code for assertion.


// -- Currently my machine & compiler:
#define GIVARO_HAVE_TYPENAME 1
#define GIVARO_USE_SIXTYFOUR 1


// ==========================================================================
// -- Version of the library :
// value: xxyyzz, where
// - xx: major version number
// - yy: minor version number
// - zz: revision number
#define GIVARO_MAJOR_VERSION    4
#define GIVARO_MINOR_VERSION    2
#define GIVARO_REVISION_VERSION 1
#define GIVARO_VERSION          40201

// ==========================================================================
// -- Defines the basic integer arithmetics available on this machine
#include <givaro-config.h>
#ifdef __GIVARO_HAVE_STDINT_H
#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif
#include <stdint.h>
#if !defined(INT64_MAX)
#if !defined(__CYGWIN__) && !defined(__PATHCC__)
#pragma message "#warning somebody nasty previously included <stdint.h> without __STDC_LIMIT_MACROS :)"
#endif
#include <limits>
#define INT64_MAX std::numeric_limits<int64_t>::max()
#endif

#if !defined(UINT64_MAX)
#if !defined(__CYGWIN__) && !defined(__PATHCC__)
#pragma message "#warning somebody nasty previously included <stdint.h> without __STDC_LIMIT_MACROS :)"
#endif
#include <limits>
#define UINT64_MAX std::numeric_limits<uint64_t>::max()
#endif

#if !defined(INT32_MAX)
#if !defined(__CYGWIN__) && !defined(__PATHCC__)
#pragma message "#warning somebody nasty previously included <stdint.h> without __STDC_LIMIT_MACROS :)"
#endif
#include <limits>
#define INT32_MAX std::numeric_limits<int32_t>::max()
#endif

#if !defined(UINT32_MAX)
#if !defined(__CYGWIN__) && !defined(__PATHCC__)
#pragma message "#warning somebody nasty previously included <stdint.h> without __STDC_LIMIT_MACROS :)"
#endif
#include <limits>
#define UINT32_MAX std::numeric_limits<uint32_t>::max()
#endif

#if !defined(INT16_MAX)
#if !defined(__CYGWIN__) && !defined(__PATHCC__)
#pragma message "#warning somebody nasty previously included <stdint.h> without __STDC_LIMIT_MACROS :)"
#endif
#include <limits>
#define INT16_MAX std::numeric_limits<int16_t>::max()
#endif

#if !defined(UINT16_MAX)
#if !defined(__CYGWIN__) && !defined(__PATHCC__)
#pragma message "#warning somebody nasty previously included <stdint.h> without __STDC_LIMIT_MACROS :)"
#endif
#include <limits>
#define UINT16_MAX std::numeric_limits<uint16_t>::max()
#endif

#if !defined(INT8_MAX)
#if !defined(__CYGWIN__) && !defined(__PATHCC__)
#pragma message "#warning somebody nasty previously included <stdint.h> without __STDC_LIMIT_MACROS :)"
#endif
#include <limits>
#define INT8_MAX std::numeric_limits<int8_t>::max()
#endif

#if !defined(UINT8_MAX)
#if !defined(__CYGWIN__) && !defined(__PATHCC__)
#pragma message "#warning somebody nasty previously included <stdint.h> without __STDC_LIMIT_MACROS :)"
#endif
#include <limits>
#define UINT8_MAX std::numeric_limits<uint8_t>::max()
#endif


#endif


# define GIVARO_BITS_PER_LONGINT 	SIZEOF_LONG
# define GIVARO_BITS_PER_INT		SIZEOF_INT
# define GIVARO_BITS_PER_SHORTINT	SIZEOF_SHORT
# define GIVARO_BITS_PER_CHAR		SIZEOF_CHAR

#ifndef __GIVARO_HAVE_STDINT_H
# define GIVARO_MAXUINT8		255U 		    // 2^8-1
# define GIVARO_MAXUINT16		65535U 		    // 2^16-1
# define GIVARO_MAXUINT32		4294967295U 	// 2^32-1
# define GIVARO_MAXULONG		4294967295U 	// 2^32-1
#else
# define GIVARO_MAXUINT8		UINT8_MAX 	    // 2^8-1
# define GIVARO_MAXUINT16		UINT16_MAX 		// 2^16-1
# define GIVARO_MAXUINT32		UINT32_MAX 	    // 2^32-1
# define GIVARO_MAXULONG		4294967295U 	//! @bug 2^32-1
#endif


// ==========================================================================
// -- Code expansion depending on the previous defined macros


// -- specialized template

#ifdef GIVARO_HAVE_ANSI_SPECIALIZED
#  define GIVARO_SPECIALIZED template<>
#else
#  define GIVARO_SPECIALIZED
#endif


// -- typename for Givaro

#ifdef GIVARO_HAVE_TYPENAME
#else
#define typename
#endif


// -- Macros for debug

#define GIV_XVALTOSTR(msg)   #msg
#define GIV_VALTOSTR(msg)   GIV_XVALTOSTR(msg)
// #define GIV_ERROR(msg) GivError( ##msg " File:" __FILE__ ", line:" GIV_VALTOSTR(__LINE__))

#ifdef __GIVARO_DEBUG
#  ifdef GIVARO_HAVE_ANSI_LIBRARY  // here is ANSI C++ header definition !!!
#    include <sstream>
#    include <givaro/giverror.h>
#    define GIVARO_ASSERT(cond, msg) { \
    if (!(cond)) {\
        std::ostringstream ostr;\
        ostr << msg << "\nFile:"##__FILE__##", Line:" << __LINE__;\
        GivError::throw_error( GivError(ostr.str().c_str()) );\
    }}
#    define GIVARO_ASSERT2(cond, msg1, msg2) { \
    if (!(cond)) {\
        std::ostringstream ostr;\
        ostr << msg1 << msg2 << "\nFile:"##__FILE__##", Line:" << __LINE__;\
        GivError::throw_error( GivError(ostr.str().c_str()) );\
    }}
#    define GIVARO_ASSERT_MATHDIV0(cond, msg) { \
    if (!(cond)) {\
        std::ostringstream ostr;\
        ostr << msg << "\nFile:"##__FILE__##", Line:" << __LINE__;\
        throw  GivMathDivZero(ostr.str().c_str()) ;\
    }}

#  else
#    include <sstream>
#    include <givaro/giverror.h>
#    define GIVARO_ASSERT(cond, msg) { \
    if (!(cond)) {\
        std::ostringstream ostr;\
        ostr << msg << "\nFile:" << __FILE__ << ", Line:" << __LINE__;\
        GivError::throw_error( GivError(ostr.str().c_str()) );\
    }}
#    define GIVARO_ASSERT2(cond, msg1, msg2) { \
    if (!(cond)) {\
        std::ostringstream ostr;\
        ostr << msg1 << msg2 << "\nFile:" << __FILE__ << ", Line:" << __LINE__;\
        GivError::throw_error( GivError(ostr.str().c_str()) );\
    }}
#    define GIVARO_ASSERT_MATHDIV0(cond, msg) { \
    if (!(cond)) {\
        std::ostringstream ostr;\
        ostr << msg << "\nFile:"<<__FILE__<<", Line:" << __LINE__;    \
        throw GivMathDivZero(ostr.str().c_str());\
    }}
#  endif

#define GIVARO_REQUIRE(bool_expression, msg) GIVARO_ASSERT2(bool_expression, "Givaro required: ", msg)
#define GIVARO_ENSURE(bool_expression, msg) GIVARO_ASSERT2(bool_expression, "Givaro ensured: ", msg)
#define GIVARO_STATE(expression)  expression
#else
#define GIVARO_ASSERT(cond, msg)
#define GIVARO_ASSERT2(cond, msg1, msg2)
#define GIVARO_REQUIRE(ignore, msg) ((void) 0)
#define GIVARO_ENSURE(ignore, msg) ((void) 0)
#define GIVARO_ASSERT_MATHDIV0(cond, msg)
// #define GIVARO_STATE(expression) do { (expression) } while (0)
#define GIVARO_STATE(ignore)  ((void) 0)

#endif

// ==========================================================================
// -- System features
#define _SYS_UNDEF 0
#define _SYS_MACOS 1

// ==========================================================================
// -- System variable
#ifndef GIVARO_SYS
#  ifdef __APPLE__
#    define GIVARO_SYS _SYS_MACOS
#  else
#    define GIVARO_SYS _SYS_UNDEF
#  endif
#endif

// ==========================================================================
// -- Misc features. Should be deleted

// -- Define this macro to store a log of memory address
// allocated during computation
//#define GIVARO_MAPMEM

// -- Define this variable to compute statistics about memory usage
#ifndef GIVARO_STATMEM
#define GIVARO_STATMEM  1
#endif



// -- Signed Traits (JGD 15.12.1999)
// JGD 11.06.03
#if !defined( __GNUC__ ) || ( __GNUC__ != 2 )
// JGD 21.03.03
#include <limits>
template<class XXX>
struct GIVARO_numeric_limits {
    typedef XXX self_type;
    static XXX max() { return  std::numeric_limits<XXX>::max(); }
};
#else
#include <limits.h>
template<class XXX>
struct GIVARO_numeric_limits {
    typedef XXX self_type;
    static XXX max() { return (XXX)(pow(2,8*sizeof(XXX))-1); }
};

/* Declared in <float.h> on ANSI C systems.  */
#ifndef DBL_MIN
#define DBL_MIN 1e-37
#endif
#ifndef DBL_MAX
#define DBL_MAX 1e+37
#endif
#ifndef FLT_MIN
#define FLT_MIN 1e-37
#endif
#ifndef FLT_MAX
#define FLT_MAX 1e+37
#endif

#ifndef __GIVARO__DONOTUSE_longlong__
#  ifndef LLONG_MAX
#    ifdef LONG_LONG_MAX
#	define LLONG_MAX	LONG_LONG_MAX
#	define LLONG_MIN	LONG_LONG_MIN
#	define ULLONG_MAX	ULONG_LONG_MAX
#    else
#	define LLONG_MAX	9223372036854775807LL
#	define LLONG_MIN	(-LLONG_MAX - 1LL)
#	define ULLONG_MAX	18446744073709551615ULL
#    endif
#  endif
#endif

//! @bug why not use numeric limits from stl ?
template<> inline float GIVARO_numeric_limits<float>::max() { return FLT_MAX; }
template<> inline double GIVARO_numeric_limits<double>::max() { return DBL_MAX; }
template<> inline short GIVARO_numeric_limits<short>::max() { return SHRT_MAX; }
template<> inline unsigned short GIVARO_numeric_limits<unsigned short>::max() { return USHRT_MAX; }
template<> inline unsigned char GIVARO_numeric_limits<unsigned char>::max() { return CHAR_MAX; }
template<> inline signed char GIVARO_numeric_limits<signed char>::max() { return UCHAR_MAX; }
template<> inline int GIVARO_numeric_limits<int>::max() { return INT_MAX; }
template<> inline unsigned int GIVARO_numeric_limits<unsigned int>::max() { return UINT_MAX; }
template<> inline long GIVARO_numeric_limits<long>::max() { return LONG_MAX; }
template<> inline unsigned long GIVARO_numeric_limits<unsigned long>::max() { return ULONG_MAX; }
#ifndef __GIVARO__DONOTUSE_longlong__
template<> inline long long GIVARO_numeric_limits<long long>::max() { return LLONG_MAX; }
template<> inline unsigned long long GIVARO_numeric_limits<unsigned long long>::max() { return ULLONG_MAX; }
#endif

#endif

template<class XXX> struct Signed_Trait : public GIVARO_numeric_limits<XXX> {
    typedef XXX signed_type;
    typedef XXX unsigned_type;
};

template<> struct Signed_Trait<float>  : public GIVARO_numeric_limits<float> {
    typedef float signed_type;
    typedef unsigned long unsigned_type;
};

template<> struct Signed_Trait<double>  : public GIVARO_numeric_limits<double> {
    typedef double signed_type;
    typedef double unsigned_type;
    //    typedef unsigned long unsigned_type;
};


template<> struct Signed_Trait<unsigned short>  : public GIVARO_numeric_limits<unsigned short> {
    typedef short signed_type;
    typedef unsigned short unsigned_type;
};

template<> struct Signed_Trait<short>  : public GIVARO_numeric_limits<short> {
    typedef short signed_type;
    typedef unsigned short unsigned_type;
};

template<> struct Signed_Trait<unsigned char>  : public GIVARO_numeric_limits<unsigned char> {
    typedef short signed_type;
    typedef unsigned char unsigned_type;
};

template<> struct Signed_Trait<signed char>  : public GIVARO_numeric_limits<signed char> {
    typedef signed char signed_type;
    typedef unsigned char unsigned_type;
};

template<> struct Signed_Trait<int>  : public GIVARO_numeric_limits<int> {
    typedef int signed_type;
    typedef unsigned int unsigned_type;
};

template<> struct Signed_Trait<unsigned int>  : public GIVARO_numeric_limits<unsigned int> {
    typedef int signed_type;
    typedef unsigned int unsigned_type;
};

template<> struct Signed_Trait<long>  : public GIVARO_numeric_limits<long> {
    typedef long signed_type;
    typedef unsigned long unsigned_type;
};

template<> struct Signed_Trait<unsigned long>  : public GIVARO_numeric_limits<unsigned long> {
    typedef long signed_type;
    typedef unsigned long unsigned_type;
};


#ifndef __GIVARO__DONOTUSE_longlong__
template<> struct Signed_Trait<long long>  : public GIVARO_numeric_limits<long long> {
    typedef long long signed_type;
    typedef unsigned long long unsigned_type;
};


template<> struct Signed_Trait<unsigned long long>  : public GIVARO_numeric_limits<unsigned long long> {
    typedef long long signed_type;
    typedef unsigned long long unsigned_type;
};
#endif


#if defined(_OPENMP) || defined(OMP_H) || defined(__OMP_H) || defined(__pmp_omp_h)
#  ifndef __GIVARO_USE_OPENMP
#    define __GIVARO_USE_OPENMP 1
#  endif
//#else
//#  undef __GIVARO_USE_OPENMP
#endif

#ifdef __GIVARO_HAVE_INT128
/* Define int128 type */
#define int128_t __int128_t

/* Define uint128 type */
#define uint128_t __uint128_t
#endif

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
