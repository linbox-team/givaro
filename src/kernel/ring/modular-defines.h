// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: A. Breust <Alexis.Breust@imag.fr>
//         Brice Boyer (briceboyer) <boyer.brice@gmail.com>
// ==========================================================================

#ifndef __GIVARO_modular_defines_H
#define __GIVARO_modular_defines_H

// r <- a^-1
#define __GIVARO_MODULAR_RECINT_NEG(r,p,a) { if (a == 0) RecInt::reset(r); else RecInt::sub(r, p, a); }
#define __GIVARO_MODULAR_INTEGER_NEG(r,p,a) ( \
    r = a == 0 ? static_cast<Element>(0) : \
    static_cast<Element>(static_cast<Compute_t>(p)-static_cast<Compute_t>(a)) )
#define __GIVARO_MODULAR_FLOATING_NEG(r,p,a) ( r = (a == 0 ? 0 : p-a) )

// r <- r^-1
#define __GIVARO_MODULAR_RECINT_NEGIN(r,p) { if (r == 0) RecInt::reset(r); else RecInt::sub(r, p, r); }
#define __GIVARO_MODULAR_INTEGER_NEGIN(r,p) ( \
    r = r == 0 ? static_cast<Element>(0) : \
    static_cast<Element>(static_cast<Compute_t>(p)-static_cast<Compute_t>(r)) )
#define __GIVARO_MODULAR_FLOATING_NEGIN(r,p) ( r = (r == 0 ? 0 : p-r) )

// r <- a * b
#define __GIVARO_MODULAR_RECINT_MUL(r,p,a,b) { RecInt::mod_n(RecInt::mul(r, a, b), p); }
#define __GIVARO_MODULAR_RECINT_LMUL(r,p,a,b) { Compute_t tmp;RecInt::lmul(tmp, a, b); RecInt::mod_n(r, tmp, p); }
#define __GIVARO_MODULAR_INTEGER_MUL(r,p,a,b) ( \
    r = static_cast<Element>(static_cast<Compute_t>(a)*static_cast<Compute_t>(b) % static_cast<Compute_t>(p)) )
#define __GIVARO_MODULAR_FLOATING_MUL(r,p,a,b) ( r = std::fmod(a*b, p) )

// r <- r * a
#define __GIVARO_MODULAR_RECINT_MULIN(r,p,a) { RecInt::mod_n(RecInt::mul(r, a), p); }
#define __GIVARO_MODULAR_RECINT_LMULIN(r,p,a) { Compute_t tmp;RecInt::lmul(tmp, r, a); RecInt::mod_n(r, tmp, p); }
#define __GIVARO_MODULAR_INTEGER_MULIN(r,p,a) ( r = static_cast<Element>(static_cast<Compute_t>(r)*static_cast<Compute_t>(a) % static_cast<Compute_t>(p)) )
#define __GIVARO_MODULAR_FLOATING_MULIN(r,p,a) ( r = std::fmod((r*a), p) )

// r <- a - b
#define __GIVARO_MODULAR_RECINT_SUB(r,p,a,b) { \
    if (a < b) { RecInt::sub(r, p, b); RecInt::add(r, a); } \
    else { RecInt::sub(r, a, b); } }
#define __GIVARO_MODULAR_INTEGER_SUB(r,p,a,b) ( \
    r = static_cast<Element>((static_cast<Compute_t>(a) >= static_cast<Compute_t>(b)) ? \
    static_cast<Compute_t>(a)-static_cast<Compute_t>(b) : \
    static_cast<Compute_t>(p)-static_cast<Compute_t>(b)+static_cast<Compute_t>(a)) )
#define __GIVARO_MODULAR_FLOATING_SUB(r,p,a,b) ( r = (a>=b) ? a-b: (p-b)+a )

// r <- r - a
#define __GIVARO_MODULAR_RECINT_SUBIN(r,p,a) { \
    if (r < a) { RecInt::add(r, p - a); } \
    else { RecInt::sub(r, a); } }
#define __GIVARO_MODULAR_INTEGER_SUBIN(r,p,a) { \
	if (static_cast<Compute_t>(r) < static_cast<Compute_t>(a)) \
	r = static_cast<Element>(static_cast<Compute_t>(r) + static_cast<Compute_t>(p) - static_cast<Compute_t>(a)); \
	else r = static_cast<Element>(static_cast<Compute_t>(r) - static_cast<Compute_t>(a)); }
#define __GIVARO_MODULAR_FLOATING_SUBIN(r,p,a) { if (r<a) r+=(p-a); else r-=a; }

// r <- a + b
#define __GIVARO_MODULAR_RECINT_ADD(r,p,a,b) { \
    RecInt::add(r, a, b); \
    if (r >= p) RecInt::sub(r, p); }
#define __GIVARO_MODULAR_INTEGER_ADD(r,p,a,b) { \
    r = static_cast<Element>(static_cast<Compute_t>(a) + static_cast<Compute_t>(b)); \
    r = static_cast<Element>((static_cast<Compute_t>(r) < static_cast<Compute_t>(p)) ? r : \
    static_cast<Compute_t>(r) - static_cast<Compute_t>(p)); }
#define __GIVARO_MODULAR_FLOATING_ADD(r,p,a,b) { r = (a+b); r= (r < p ? r : r-p); }

// r <- r + a
#define __GIVARO_MODULAR_RECINT_ADDIN(r,p,a) { \
    RecInt::add(r, a); \
    if (r >= p) RecInt::sub(r, p); }
#define __GIVARO_MODULAR_INTEGER_ADDIN(r,p,a) { \
    r = static_cast<Element>(r + a); \
    r = static_cast<Element>((static_cast<Compute_t>(r) < static_cast<Compute_t>(p)) ? r : \
    static_cast<Compute_t>(r) - static_cast<Compute_t>(p)); }
#define __GIVARO_MODULAR_FLOATING_ADDIN(r,p,a) { r += a;  r= (r < p ? r : r-p); }

// r <- c + a*b
#define __GIVARO_MODULAR_RECINT_MULADD(r,p,a,b,c) {	\
    RecInt::copy(r, c);					\
    RecInt::addmul(r, a, b);				\
    RecInt::mod_n(r, p); }
#define __GIVARO_MODULAR_RECINT_LMULADD(r,p,a,b,c) {	\
    __GIVARO_MODULAR_RECINT_LMUL(r,p,a,b); \
    __GIVARO_MODULAR_RECINT_ADDIN(r,p,c); }
#define __GIVARO_MODULAR_INTEGER_MULADD(r,p,a,b,c) ( \
    r = static_cast<Element>((static_cast<Compute_t>(a)*static_cast<Compute_t>(b) \
    + static_cast<Compute_t>(c)) % static_cast<Compute_t>(p)))
#define __GIVARO_MODULAR_FLOATING_MULADD(r,p,a,b,c) ( r = std::fmod((a*b+c), p) )

// r <- r + a*b
#define __GIVARO_MODULAR_RECINT_MULADDIN(r,p,a,b) { \
    RecInt::addmul(r, a, b); \
    RecInt::mod_n(r, p); }
#define __GIVARO_MODULAR_RECINT_LMULADDIN(r,p,a,b) { \
    Element tmp=r;				     \
    __GIVARO_MODULAR_RECINT_LMULADD(r,p,a,b,tmp);} 
#define __GIVARO_MODULAR_INTEGER_MULADDIN(r,p,a,b) ( \
    r = static_cast<Element>((static_cast<Compute_t>(a)*static_cast<Compute_t>(b) \
    + static_cast<Compute_t>(r)) % static_cast<Compute_t>(p)))
#define __GIVARO_MODULAR_FLOATING_MULADDIN(r,p,a,b) ( r = std::fmod((a*b+r), p) )

// r <- a*b - c
#define __GIVARO_MODULAR_RECINT_MULSUB(r,p,a,b,c) { \
    __GIVARO_MODULAR_RECINT_MUL(r,p,a,b); \
    __GIVARO_MODULAR_RECINT_SUBIN(r,p,c); }
#define __GIVARO_MODULAR_RECINT_LMULSUB(r,p,a,b,c) { \
    __GIVARO_MODULAR_RECINT_LMUL(r,p,a,b); \
    __GIVARO_MODULAR_RECINT_SUBIN(r,p,c); }
#define __GIVARO_MODULAR_INTEGER_MULSUB(r,p,a,b,c) (			\
    r = static_cast<Element>((static_cast<Compute_t>(a)*static_cast<Compute_t>(b) \
    + static_cast<Compute_t>(p)-static_cast<Compute_t>(c)) % static_cast<Compute_t>(p)))
#define __GIVARO_MODULAR_FLOATING_MULSUB(r,p,a,b,c) ( r = std::fmod((a*b+p-c), p) )

// r <- a*b - r
#define __GIVARO_MODULAR_RECINT_SUBMULIN(r,p,a,b) { \
    __GIVARO_MODULAR_RECINT_NEGIN(r,p); \
    RecInt::addmul(r, a, b); \
    RecInt::mod_n(r, p); }
#define __GIVARO_MODULAR_RECINT_LSUBMULIN(r,p,a,b) { \
    Element tmp=r;				     \
    __GIVARO_MODULAR_RECINT_LMULSUB(r,p,a,b,tmp);} 
#define __GIVARO_MODULAR_INTEGER_SUBMULIN(r,p,a,b) {			\
    r = static_cast<Element>((static_cast<Compute_t>(a)*static_cast<Compute_t>(b) \
    +static_cast<Compute_t>(p)-static_cast<Compute_t>(r)) % static_cast<Compute_t>(p)); \
    __GIVARO_MODULAR_INTEGER_NEGIN(r,p); }
#define __GIVARO_MODULAR_FLOATING_SUBMULIN(r,p,a,b) { \
    r = (a*b+p-r); \
    r= (r<p ? r : std::fmod(r, p)); \
    __GIVARO_MODULAR_FLOATING_NEGIN(r,p); }

#endif // __GIVARO_modular_defines_H

