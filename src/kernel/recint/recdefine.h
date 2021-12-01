/* consig/define.h - Debug and useful defines

   Copyright Université Grenoble Alpes
Contributors :
Alexis BREUST (alexis.breust@gmail.com 2014)

This software is a computer program whose purpose is to provide a
fixed precision arithmetic library.

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/


#ifndef __RECINTCONFIG_DEFINE_H
#define __RECINTCONFIG_DEFINE_H

// --------------------------------------------------------------
// --------------------- Various defines ------------------------

#include <cstdlib> /* For size_t */
#include <cstdint> /* For uint64_t and so on */

// FIXME Get info at configure-time - A.B.
// Check for anonymous unions + anonymous structs + __uint128_t
// + constructor variables in aggregate
// NOTE : __GIVARO_HAVE_INT128 is now available to test for __uint128_t
// #define __RECINT_USE_FAST_128

namespace RecInt
{
    /* Size of a minimal recint */
    typedef uint64_t limb;
}

/* Set here the threshold above which Karatsuba method of multiplication is used.
   Changing this value may be interesting for some machines. (Default: 10) */
#if not defined(__RECINT_THRESHOLD_KARA)
#define __RECINT_THRESHOLD_KARA 10
#endif

/* All computations assume this is a 64 bits machine.
   However, the code will work fine (but slower) on a 32 bits machine. */
#define __RECINT_LIMB_BITS 64
#define __RECINT_LIMB_SIZE 6

/* Some useful defines. */
#define __RECINT_MINUSONE limb(0xffffffffffffffff)
#define __RECINT_MAXPOWTWO limb(0x8000000000000000)
#define __RECINT_TYPENOTMAXPOWTWO(T) ~(T(1) << (8 * sizeof(T) - 1))
#define __RECINT_THIRTYONEPOINTFIVE 3037000499u

// --------------------------------------------------------------
// ----------------- Variables for quick access -----------------

namespace RecInt
{
    /* NBLIMB<K>::value gives the number of limbs of a ruint<K> */
    template <size_t K> struct NBLIMB   { static const uint32_t value = NBLIMB<K-1>::value << 1; };
    template<> struct NBLIMB<__RECINT_LIMB_SIZE> { static const uint32_t value = 1; };

    /* NBBITS<K>::value gives the total size of a ruint<K> */
    template <size_t K> struct NBBITS   { static const uint32_t value = NBBITS<K-1>::value << 1; };
    template<> struct NBBITS<__RECINT_LIMB_SIZE> { static const uint32_t value = __RECINT_LIMB_BITS; };
}


// --------------------------------------------------------------
// ------------------- Template compatibility -------------------

#include <type_traits>

/* If typename T is an arithmetic type,
   then template enable and return value is RET */
#define __RECINT_IS_ARITH(T, ...)    typename std::enable_if<std::is_arithmetic<T>::value, __VA_ARGS__>::type

/* If typename T is an unsigned type,
   then template enable and return value is RET */
#define __RECINT_IS_UNSIGNED(T, ...) typename std::enable_if<std::is_unsigned<T>::value, __VA_ARGS__>::type

/* If typename T is a signed type,
   then template enable and return value is RET */
#define __RECINT_IS_SIGNED(T, ...)   typename std::enable_if<std::is_signed<T>::value, __VA_ARGS__>::type

/* If typename T is not a fundamental type,
   then template enable and return value is RET */
#define __RECINT_IS_NOT_FUNDAMENTAL(T, ...)   typename std::enable_if<!std::is_fundamental<T>::value, __VA_ARGS__>::type


// --------------------------------------------------------------
// --------------------- Debug stuff ----------------------------

/* Speed print in hexadecimal and line */
#define __RECINT_DEBUG_LINE()  std::cout << "------" << std::endl
#define __RECINT_DEBUG_SHOW(X) std::cout << std::hex << (X) << std::dec << std::endl


// --------------------------------------------------------------
// ----------------- Long long from GMP library -----------------

namespace RecInt
{
    /* See extern/longlong.h */
    typedef uint64_t UWtype;
    typedef uint64_t UHWtype;
    typedef uint64_t UDWtype;

    typedef uint64_t UDItype;
    typedef uint32_t USItype;
    typedef int64_t  DItype;
    typedef int32_t  SItype;

#ifdef W_TYPE_SIZE
#undef W_TYPE_SIZE
#endif

#ifndef NO_ASM
#define NO_ASM
#endif

#define W_TYPE_SIZE 64
#include "reclonglong.h"
}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
