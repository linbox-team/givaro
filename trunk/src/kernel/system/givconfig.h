#ifndef _GIVARO_CONFIG_H_
#define _GIVARO_CONFIG_H_
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/system/givconfig.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givconfig.h,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description: configuration file for Givaro


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
// * GIVARO_DEBUG:
//   value: defined/undefined
//   purpose: expand additional code for verification and debuging.
// * GIVARO_DEBUG_LEVEL:
//   value: integer > 0
//   purpose: <to be defined>


// -- Currently my machine & compiler:
#define GIVARO_HAVE_TYPENAME


// ==========================================================================
// -- Version of the library :
// value: xxyyzz, where
// - xx: major version number
// - yy: minor version number
// - zz: revision number
#define GIVARO_VERSION  0x020000

// -- Define this value both to compile the library of user program
// value: integer that defines debug level trace information (not well defined)
//#define GIVARO_DEBUG 1


// ==========================================================================
// -- Defined the basic integer arithmetics available on this machine
#include "givaro-config.h"

# define GIVARO_BITS_PER_LONGINT 	SIZEOF_LONG
# define GIVARO_BITS_PER_INT		SIZEOF_INT
# define GIVARO_BITS_PER_SHORTINT	SIZEOF_SHORT
# define GIVARO_BITS_PER_CHAR		SIZEOF_CHAR

typedef signed    __GIVARO_INT8       int8;
typedef signed    __GIVARO_INT16      int16;
typedef signed 	  __GIVARO_INT32      int32;
typedef unsigned  __GIVARO_INT8      uint8;
typedef unsigned  __GIVARO_INT16     uint16;
typedef unsigned  __GIVARO_INT32     uint32;

# define GIVARO_MAXUINT8		255U 		// 2^8-1
# define GIVARO_MAXUINT16		65535U 		// 2^16-1
# define GIVARO_MAXUINT32		4294967295U 	// 2^32-1
# define GIVARO_MAXULONG		4294967295U 	// 2^32-1



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

#ifdef GIVARO_DEBUG
#  ifdef GIVARO_HAVE_ANSI_LIBRARY  // here is ANSI C++ header definition !!!
#    include <sstream>
#    define GIVARO_ASSERT(cond, msg) { \
      if (!(cond)) {\
        ostringstream ostr;\
        ostr << msg << "\nFile:"##__FILE__##", Line:" << __LINE__;\
        GivError::throw_error( GivError(ostr.str().c_str()) );\
      }}
#    define GIVARO_ASSERT2(cond, msg1, msg2) { \
      if (!(cond)) {\
        ostringstream ostr;\
        ostr << msg1 << msg2 << "\nFile:"##__FILE__##", Line:" << __LINE__;\
        GivError::throw_error( GivError(ostr.str().c_str()) );\
      }}
#  else
#    include <strstream>
#    define GIVARO_ASSERT(cond, msg) { \
      if (!(cond)) {\
        ostrstream ostr;\
        ostr << msg << "\nFile:" << __FILE__ << ", Line:" << __LINE__;\
        GivError::throw_error( GivError(ostr.str()) );\
      }}
#    define GIVARO_ASSERT2(cond, msg1, msg2) { \
      if (!(cond)) {\
        ostrstream ostr;\
        ostr << msg1 << msg2 << "\nFile:" << __FILE__ << ", Line:" << __LINE__;\
        GivError::throw_error( GivError(ostr.str()) );\
      }}
#  endif
#else
#define GIVARO_ASSERT(cond, msg)
#define GIVARO_ASSERT2(cond, msg1, msg2)
#endif

// ==========================================================================
// -- System features
#define _SYS_UNDEF 0
#define _SYS_MACOS 1



// ==========================================================================
// -- System variable
#define GIVARO_SYS _SYS_UNDEF



// ==========================================================================
// -- Misc features. Should be deleted

// -- Define this macro to store a log of memory address
// allocated during computation
//#define GIVARO_MAPMEM 

// -- Define this variable to compute statistics about memory usage
#define GIVARO_STATMEM  1



// -- Signed Traits (JGD 15.12.1999)
// JGD 11.06.03
#if !defined(__GNUC__) || (__GNUC__ != 2)
// JGD 21.03.03
#include <limits>
template<class XXX> struct GIVARO_numeric_limits {
    	typedef XXX self_type;
	static XXX max() { return  std::numeric_limits<XXX>::max(); }
};
#else
#include <limits.h>
template<class XXX> struct GIVARO_numeric_limits {
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
  #ifndef __DONOTUSE_longlong__
template<> inline long long GIVARO_numeric_limits<long long>::max() { return LLONG_MAX; }
template<> inline unsigned long long GIVARO_numeric_limits<unsigned long long>::max() { return ULLONG_MAX; }
  #endif

#endif

template<class XXX> struct Signed_Trait : public GIVARO_numeric_limits<XXX> {
    typedef XXX signed_type;
};

template<> struct Signed_Trait<float>  : public GIVARO_numeric_limits<float> {
    typedef float signed_type;
    typedef unsigned long unsigned_type;
};

template<> struct Signed_Trait<double>  : public GIVARO_numeric_limits<double> {
    typedef double signed_type;
    typedef unsigned long unsigned_type;
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


  #ifndef __DONOTUSE_longlong__
template<> struct Signed_Trait<long long>  : public GIVARO_numeric_limits<long long> {
    typedef long long signed_type;
    typedef unsigned long long unsigned_type;
};


template<> struct Signed_Trait<unsigned long long>  : public GIVARO_numeric_limits<unsigned long long> {
    typedef long long signed_type;
    typedef unsigned long long unsigned_type;
};
  #endif 

#endif
