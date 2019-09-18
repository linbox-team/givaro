/* longlong.h -- definitions for mixed size 32/64 bit arithmetic.

   Copyright 1991, 1992, 1993, 1994, 1996, 1997, 1999, 2000, 2001, 2002, 2003,
   2004, 2005, 2007, 2008, 2009, 2011, 2012 Free Software Foundation, Inc.

   This file is free software; you can redistribute it and/or modify it under the
   terms of the GNU Lesser General Public License as published by the Free
   Software Foundation; either version 3 of the License, or (at your option) any
   later version.

   This file is distributed in the hope that it will be useful, but WITHOUT ANY
   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
   PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
   details.

   You should have received a copy of the GNU Lesser General Public License
   along with this file.  If not, see http://www.gnu.org/licenses/.  */

/* You have to define the following before including this file:

   UWtype -- An unsigned type, default type for operations (typically a "word")
   UHWtype -- An unsigned type, at least half the size of UWtype
   UDWtype -- An unsigned type, at least twice as large a UWtype
   W_TYPE_SIZE -- size in bits of UWtype

   SItype, USItype -- Signed and unsigned 32 bit types
   DItype, UDItype -- Signed and unsigned 64 bit types

   On a 32 bit machine UWtype should typically be USItype;
   on a 64 bit machine, UWtype should typically be UDItype.

   Optionally, define:

   LONGLONG_STANDALONE -- Avoid code that needs machine-dependent support files
   NO_ASM -- Disable inline asm


   CAUTION!  Using this version of longlong.h outside of GMP is not safe.  You
   need to include gmp.h and gmp-impl.h, or certain things might not work as
   expected.
*/

/* longlong.h may already be included from another library (e.g. flint) */

#define __BITS4 (W_TYPE_SIZE / 4)
#define __ll_B ((UWtype) 1 << (W_TYPE_SIZE / 2))
#define __ll_lowpart(t) ((UWtype) (t) & (__ll_B - 1))
#define __ll_highpart(t) ((UWtype) (t) >> (W_TYPE_SIZE / 2))

/* This is used to make sure no undesirable sharing between different libraries
   that use this file takes place.  */
#ifndef __MPN
#define __MPN(x) __##x
#endif

/* Define auxiliary asm macros.

   1) recint_umul_ppmm(high_prod, low_prod, multiplier, multiplicand) multiplies two
   UWtype integers MULTIPLIER and MULTIPLICAND, and generates a two UWtype
   word product in HIGH_PROD and LOW_PROD.

   2) __umulsidi3(a,b) multiplies two UWtype integers A and B, and returns a
   UDWtype product.  This is just a variant of recint_umul_ppmm.

   3) recint_udiv_qrnnd(quotient, remainder, high_numerator, low_numerator,
   denominator) divides a UDWtype, composed by the UWtype integers
   HIGH_NUMERATOR and LOW_NUMERATOR, by DENOMINATOR and places the quotient
   in QUOTIENT and the remainder in REMAINDER.  HIGH_NUMERATOR must be less
   than DENOMINATOR for correct operation.  If, in addition, the most
   significant bit of DENOMINATOR must be 1, then the pre-processor symbol
   UDIV_NEEDS_NORMALIZATION is defined to 1.

   4) sdiv_qrnnd(quotient, remainder, high_numerator, low_numerator,
   denominator).  Like recint_udiv_qrnnd but the numbers are signed.  The quotient
   is rounded towards 0.

   5) recint_count_leading_zeros(count, x) counts the number of zero-bits from the
   msb to the first non-zero bit in the UWtype X.  This is the number of
   steps X needs to be shifted left to set the msb.  Undefined for X == 0,
   unless the symbol RECINT_COUNT_LEADING_ZEROS_0 is defined to some value.

   6) recint_count_trailing_zeros(count, x) like recint_count_leading_zeros, but counts
   from the least significant end.

   7) recint_add_ssaaaa(high_sum, low_sum, high_addend_1, low_addend_1,
   high_addend_2, low_addend_2) adds two UWtype integers, composed by
   HIGH_ADDEND_1 and LOW_ADDEND_1, and HIGH_ADDEND_2 and LOW_ADDEND_2
   respectively.  The result is placed in HIGH_SUM and LOW_SUM.  Overflow
   (i.e. carry out) is not stored anywhere, and is lost.

   8) recint_sub_ddmmss(high_difference, low_difference, high_minuend, low_minuend,
   high_subtrahend, low_subtrahend) subtracts two two-word UWtype integers,
   composed by HIGH_MINUEND_1 and LOW_MINUEND_1, and HIGH_SUBTRAHEND_2 and
   LOW_SUBTRAHEND_2 respectively.  The result is placed in HIGH_DIFFERENCE
   and LOW_DIFFERENCE.  Overflow (i.e. carry out) is not stored anywhere,
   and is lost.

   If any of these macros are left undefined for a particular CPU,
   C macros are used.


   Notes:

   For recint_add_ssaaaa the two high and two low addends can both commute, but
   unfortunately gcc only supports one "%" commutative in each asm block.
   This has always been so but is only documented in recent versions
   (eg. pre-release 3.3).  Having two or more "%"s can cause an internal
   compiler error in certain rare circumstances.

   Apparently it was only the last "%" that was ever actually respected, so
   the code has been updated to leave just that.  Clearly there's a free
   choice whether high or low should get it, if there's a reason to favour
   one over the other.  Also obviously when the constraints on the two
   operands are identical there's no benefit to the reloader in any "%" at
   all.

*/

/* The CPUs come in alphabetical order below.

   Please add support for more CPUs here, or improve the current support
   for the CPUs below!  */


/* recint_count_leading_zeros_gcc_clz is recint_count_leading_zeros implemented with gcc
   3.4 __builtin_clzl or __builtin_clzll, according to our limb size.
   Similarly recint_count_trailing_zeros_gcc_ctz using __builtin_ctzl or
   __builtin_ctzll.

   These builtins are only used when we check what code comes out, on some
   chips they're merely libgcc calls, where we will instead want an inline
   in that case (either asm or generic C).

   These builtins are better than an asm block of the same insn, since an
   asm block doesn't give gcc any information about scheduling or resource
   usage.  We keep an asm block for use on prior versions of gcc though.

   For reference, __builtin_ffs existed in gcc prior to __builtin_clz, but
   it's not used (for recint_count_leading_zeros) because it generally gives extra
   code to ensure the result is 0 when the input is 0, which we don't need
   or want.  */

#ifdef _LONG_LONG_LIMB
#define recint_count_leading_zeros_gcc_clz(count,x)    \
    do {                                        \
        (count) = __builtin_clzll (x);          \
    } while (0)
#else
#define recint_count_leading_zeros_gcc_clz(count,x)    \
    do {                                        \
        (count) = __builtin_clzl (x);           \
    } while (0)
#endif

#ifdef _LONG_LONG_LIMB
#define recint_count_trailing_zeros_gcc_ctz(count,x)   \
    do {                                        \
        (count) = __builtin_ctzll (x);          \
    } while (0)
#else
#define recint_count_trailing_zeros_gcc_ctz(count,x)   \
    do {                                        \
        (count) = __builtin_ctzl (x);           \
    } while (0)
#endif



    /* FIXME: "sidi" here is highly doubtful, should sometimes be "diti".  */
#if !defined (recint_umul_ppmm) && defined (__umulsidi3)
#define recint_umul_ppmm(ph, pl, m0, m1)               \
    {                                           \
        UDWtype __ll = __umulsidi3 (m0, m1);    \
        ph = (UWtype) (__ll >> W_TYPE_SIZE);    \
        pl = (UWtype) __ll;                     \
    }
#endif

#if !defined (__umulsidi3)
#define __umulsidi3(u, v)                               \
    ({UWtype __hi, __lo;                                \
        recint_umul_ppmm (__hi, __lo, u, v);                   \
        ((UDWtype) __hi << W_TYPE_SIZE) | __lo; })
#endif


    /* Use mpn_recint_umul_ppmm or mpn_recint_udiv_qrnnd functions, if they exist.  The "_r"
       forms have "reversed" arguments, meaning the pointer is last, which
       sometimes allows better parameter passing, in particular on 64-bit
       hppa. */

#define mpn_recint_umul_ppmm  __MPN(recint_umul_ppmm)
extern UWtype mpn_recint_umul_ppmm (UWtype *, UWtype, UWtype);

#if ! defined (recint_umul_ppmm) && HAVE_NATIVE_mpn_recint_umul_ppmm  \
    && ! defined (LONGLONG_STANDALONE)
#define recint_umul_ppmm(wh, wl, u, v)                                         \
    do {                                                                \
        UWtype __recint_umul_ppmm__p0;                                         \
        (wh) = mpn_recint_umul_ppmm (&__recint_umul_ppmm__p0, (UWtype) (u), (UWtype) (v)); \
        (wl) = __recint_umul_ppmm__p0;                                         \
    } while (0)
#endif

#define mpn_recint_umul_ppmm_r  __MPN(recint_umul_ppmm_r)
extern UWtype mpn_recint_umul_ppmm_r (UWtype, UWtype, UWtype *);

#if ! defined (recint_umul_ppmm) && HAVE_NATIVE_mpn_recint_umul_ppmm_r	\
    && ! defined (LONGLONG_STANDALONE)
#define recint_umul_ppmm(wh, wl, u, v)                                         \
    do {                                                                \
        UWtype __recint_umul_ppmm__p0;                                         \
        (wh) = mpn_recint_umul_ppmm_r ((UWtype) (u), (UWtype) (v), &__recint_umul_ppmm__p0); \
        (wl) = __recint_umul_ppmm__p0;                                         \
    } while (0)
#endif

#define mpn_recint_udiv_qrnnd  __MPN(recint_udiv_qrnnd)
extern UWtype mpn_recint_udiv_qrnnd (UWtype *, UWtype, UWtype, UWtype);

#if ! defined (recint_udiv_qrnnd) && HAVE_NATIVE_mpn_recint_udiv_qrnnd	\
    && ! defined (LONGLONG_STANDALONE)
#define recint_udiv_qrnnd(q, r, n1, n0, d)					\
    do {                                                                \
        UWtype __recint_udiv_qrnnd__r;						\
        (q) = mpn_recint_udiv_qrnnd (&__recint_udiv_qrnnd__r,				\
                              (UWtype) (n1), (UWtype) (n0), (UWtype) d); \
        (r) = __recint_udiv_qrnnd__r;						\
    } while (0)
#endif

#define mpn_recint_udiv_qrnnd_r  __MPN(recint_udiv_qrnnd_r)
extern UWtype mpn_recint_udiv_qrnnd_r (UWtype, UWtype, UWtype, UWtype *);

#if ! defined (recint_udiv_qrnnd) && HAVE_NATIVE_mpn_recint_udiv_qrnnd_r	\
    && ! defined (LONGLONG_STANDALONE)
#define recint_udiv_qrnnd(q, r, n1, n0, d)					\
    do {                                                                \
        UWtype __recint_udiv_qrnnd__r;						\
        (q) = mpn_recint_udiv_qrnnd_r ((UWtype) (n1), (UWtype) (n0), (UWtype) d, \
                                &__recint_udiv_qrnnd__r);                      \
        (r) = __recint_udiv_qrnnd__r;						\
    } while (0)
#endif


    /* If this machine has no inline assembler, use C macros.  */

#if !defined (recint_add_ssaaaa)
#define recint_add_ssaaaa(sh, sl, ah, al, bh, bl)      \
    do {                                        \
        UWtype __x;                             \
        __x = (al) + (bl);                      \
        (sh) = (ah) + (bh) + (__x < (al));      \
        (sl) = __x;                             \
    } while (0)
#endif

#if !defined (recint_sub_ddmmss)
#define recint_sub_ddmmss(sh, sl, ah, al, bh, bl)      \
    do {                                        \
        UWtype __x;                             \
        __x = (al) - (bl);                      \
        (sh) = (ah) - (bh) - ((al) < (bl));     \
        (sl) = __x;                             \
    } while (0)
#endif

    /* If we lack recint_umul_ppmm but have recint_smul_ppmm, define recint_umul_ppmm in terms of
       recint_smul_ppmm.  */
#if !defined (recint_umul_ppmm) && defined (recint_smul_ppmm)
#define recint_umul_ppmm(w1, w0, u, v)                                 \
    do {                                                        \
        UWtype __w1;                                            \
        UWtype __xm0 = (u), __xm1 = (v);                        \
        recint_smul_ppmm (__w1, w0, __xm0, __xm1);                     \
        (w1) = __w1 + (-(__xm0 >> (W_TYPE_SIZE - 1)) & __xm1)   \
            + (-(__xm1 >> (W_TYPE_SIZE - 1)) & __xm0);		\
    } while (0)
#endif

    /* If we still don't have recint_umul_ppmm, define it using plain C.

       For reference, when this code is used for squaring (ie. u and v identical
       expressions), gcc recognises __x1 and __x2 are the same and generates 3
       multiplies, not 4.  The subsequent additions could be optimized a bit,
       but the only place GMP currently uses such a square is mpn_sqr_basecase,
       and chips obliged to use this generic C umul will have plenty of worse
       performance problems than a couple of extra instructions on the diagonal
       of sqr_basecase.  */

#if !defined (recint_umul_ppmm)
#define recint_umul_ppmm(w1, w0, u, v)						\
    do {                                                                \
        UWtype __x0, __x1, __x2, __x3;					\
        UHWtype __ul, __vl, __uh, __vh;					\
        UWtype __u = (u), __v = (v);					\
                                                                        \
        __ul = __ll_lowpart (__u);                                      \
        __uh = __ll_highpart (__u);                                     \
        __vl = __ll_lowpart (__v);                                      \
        __vh = __ll_highpart (__v);                                     \
                                                                        \
        __x0 = (UWtype) __ul * __vl;					\
        __x1 = (UWtype) __ul * __vh;					\
        __x2 = (UWtype) __uh * __vl;					\
        __x3 = (UWtype) __uh * __vh;					\
                                                                        \
        __x1 += __ll_highpart (__x0);/* this can't give carry */        \
        __x1 += __x2;		/* but this indeed can */		\
        if (__x1 < __x2)		/* did we get it? */            \
            __x3 += __ll_B;		/* yes, add it in the proper pos. */ \
                                                                        \
        (w1) = __x3 + __ll_highpart (__x1);                             \
        (w0) = (__x1 << W_TYPE_SIZE/2) + __ll_lowpart (__x0);		\
    } while (0)
#endif

    /* If we don't have recint_smul_ppmm, define it using recint_umul_ppmm (which surely will
       exist in one form or another.  */
#if !defined (recint_smul_ppmm)
#define recint_smul_ppmm(w1, w0, u, v)                                 \
    do {                                                        \
        UWtype __w1;                                            \
        UWtype __xm0 = (u), __xm1 = (v);                        \
        recint_umul_ppmm (__w1, w0, __xm0, __xm1);                     \
        (w1) = __w1 - (-(__xm0 >> (W_TYPE_SIZE - 1)) & __xm1)   \
            - (-(__xm1 >> (W_TYPE_SIZE - 1)) & __xm0);		\
    } while (0)
#endif

    /* Define this unconditionally, so it can be used for debugging.  */
#define __recint_udiv_qrnnd_c(q, r, n1, n0, d)                                 \
    do {                                                                \
        UWtype __d1, __d0, __q1, __q0, __r1, __r0, __m;			\
                                                                        \
        __d1 = __ll_highpart (d);                                       \
        __d0 = __ll_lowpart (d);                                        \
                                                                        \
        __q1 = (n1) / __d1;                                             \
        __r1 = (n1) - __q1 * __d1;                                      \
        __m = __q1 * __d0;                                              \
        __r1 = __r1 * __ll_B | __ll_highpart (n0);                      \
        if (__r1 < __m)							\
        {                                                               \
            __q1--, __r1 += (d);                                        \
            if (__r1 >= (d)) /* i.e. we didn't get carry when adding to __r1 */ \
                if (__r1 < __m)						\
                    __q1--, __r1 += (d);                                \
        }                                                               \
        __r1 -= __m;							\
                                                                        \
        __q0 = __r1 / __d1;                                             \
        __r0 = __r1  - __q0 * __d1;                                     \
        __m = __q0 * __d0;                                              \
        __r0 = __r0 * __ll_B | __ll_lowpart (n0);                       \
        if (__r0 < __m)							\
        {                                                               \
            __q0--, __r0 += (d);                                        \
            if (__r0 >= (d))						\
                if (__r0 < __m)						\
                    __q0--, __r0 += (d);                                \
        }                                                               \
        __r0 -= __m;							\
                                                                        \
        (q) = __q1 * __ll_B | __q0;                                     \
        (r) = __r0;                                                     \
    } while (0)

    /* If the processor has no recint_udiv_qrnnd but sdiv_qrnnd, go through
       __udiv_w_sdiv (defined in libgcc or elsewhere).  */
/*
#if !defined (recint_udiv_qrnnd) && defined (sdiv_qrnnd)
   #define recint_udiv_qrnnd(q, r, nh, nl, d)              \
     do {                                                \
         UWtype __r;                                     \
         (q) = __MPN(udiv_w_sdiv) (&__r, nh, nl, d);     \
         (r) = __r;                                      \
     } while (0)
 UWtype __MPN(udiv_w_sdiv) (UWtype *, UWtype, UWtype, UWtype);
 #endif
*/
    /* If recint_udiv_qrnnd was not defined for this processor, use __recint_udiv_qrnnd_c.  */
#if !defined (recint_udiv_qrnnd)
#define UDIV_NEEDS_NORMALIZATION 1
#define recint_udiv_qrnnd __recint_udiv_qrnnd_c
#endif

#if !defined (recint_count_leading_zeros)
#define recint_count_leading_zeros(count, x)                                   \
    do {                                                                \
        UWtype __xr = (x);                                              \
        UWtype __a;                                                     \
                                                                        \
        if (W_TYPE_SIZE == 32)						\
        {                                                               \
            __a = __xr < ((UWtype) 1 << 2*__BITS4)                      \
                ? (__xr < ((UWtype) 1 << __BITS4) ? 1 : __BITS4 + 1)    \
                : (__xr < ((UWtype) 1 << 3*__BITS4) ? 2*__BITS4 + 1     \
                   : 3*__BITS4 + 1);                                    \
        }                                                               \
        else								\
        {                                                               \
            for (__a = W_TYPE_SIZE - 8; __a > 0; __a -= 8)              \
                if (((__xr >> __a) & 0xff) != 0)                        \
                    break;                                              \
            ++__a;                                                      \
        }                                                               \
                                                                        \
        (count) = W_TYPE_SIZE + 1 - __a - __clz_tab[__xr >> __a];       \
    } while (0)
    /* This version gives a well-defined value for zero. */
#define RECINT_COUNT_LEADING_ZEROS_0 (W_TYPE_SIZE - 1)
#define RECINT_COUNT_LEADING_ZEROS_NEED_CLZ_TAB
#define RECINT_COUNT_LEADING_ZEROS_SLOW
#endif

    /* clz_tab needed by mpn/x86/pentium/mod_1.asm in a fat binary */
#if HAVE_HOST_CPU_FAMILY_x86 && WANT_FAT_BINARY
#define RECINT_COUNT_LEADING_ZEROS_NEED_CLZ_TAB
#endif

#ifdef RECINT_COUNT_LEADING_ZEROS_NEED_CLZ_TAB
extern const unsigned char __clz_tab[129];
#endif

#if !defined (recint_count_trailing_zeros)
#if !defined (RECINT_COUNT_LEADING_ZEROS_SLOW)
    /* Define recint_count_trailing_zeros using an asm recint_count_leading_zeros.  */
#define recint_count_trailing_zeros(count, x)                          \
    do {                                                        \
        UWtype __ctz_x = (x);                                   \
        UWtype __ctz_c;                                         \
        recint_count_leading_zeros (__ctz_c, __ctz_x & -__ctz_x);      \
        (count) = W_TYPE_SIZE - 1 - __ctz_c;                    \
    } while (0)
#else
    /* Define recint_count_trailing_zeros in plain C, assuming small counts are common.
       We use clz_tab without ado, since the C recint_count_leading_zeros above will have
       pulled it in.  */
#define recint_count_trailing_zeros(count, x)					\
    do {                                                                \
        UWtype __ctz_x = (x);						\
        int __ctz_c;							\
                                                                        \
        if (LIKELY ((__ctz_x & 0xff) != 0))                             \
            (count) = __clz_tab[__ctz_x & -__ctz_x] - 2;                \
        else								\
        {                                                               \
            for (__ctz_c = 8 - 2; __ctz_c < W_TYPE_SIZE - 2; __ctz_c += 8) \
            {								\
                __ctz_x >>= 8;						\
                if (LIKELY ((__ctz_x & 0xff) != 0))                     \
                    break;                                              \
            }								\
                                                                        \
            (count) = __ctz_c + __clz_tab[__ctz_x & -__ctz_x];		\
        }                                                               \
    } while (0)
#endif
#endif

#ifndef UDIV_NEEDS_NORMALIZATION
#define UDIV_NEEDS_NORMALIZATION 0
#endif

    /* Whether recint_udiv_qrnnd is actually implemented with recint_udiv_qrnnd_preinv, and
       that hence the latter should always be used.  */
#ifndef UDIV_PREINV_ALWAYS
#define UDIV_PREINV_ALWAYS 0
#endif

    /* Give defaults for UMUL_TIME and UDIV_TIME.  */
#ifndef UMUL_TIME
#define UMUL_TIME 1
#endif

#ifndef UDIV_TIME
#define UDIV_TIME UMUL_TIME
#endif
    /* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
    // vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
