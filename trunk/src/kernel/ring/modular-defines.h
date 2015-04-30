// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: A. Breust <Alexis.Breust@imag.fr>
//          BB <brice.boyer@lip6.fr>
// ==========================================================================

#ifndef __GIVARO_modular_defines_H
#define __GIVARO_modular_defines_H

// r <- a^-1
#define __GIVARO_MODULAR_INTEGER_NEG(r,p,a) ( r = Element(a == 0 ? 0 : Compute_t(p)-Compute_t(a)) )
#define __GIVARO_MODULAR_FLOATING_NEG(r,p,a) ( r = (a == 0 ? 0 : p-a) )
// r <- r^-1
#define __GIVARO_MODULAR_INTEGER_NEGIN(r,p) ( r = Element(r == 0 ? 0 : Compute_t(p)-Compute_t(r)) )
#define __GIVARO_MODULAR_FLOATING_NEGIN(r,p) ( r = (r == 0 ? 0 : p-r) )

// r <- a * b
#define __GIVARO_MODULAR_INTEGER_MUL(r,p,a,b) ( r = Element(Compute_t(a)*Compute_t(b) % Compute_t(p)) )
#define __GIVARO_MODULAR_FLOATING_MUL(r,p,a,b) ( r = std::fmod(a*b, p) )
// r <- r * a
#define __GIVARO_MODULAR_INTEGER_MULIN(r,p,a) ( r = Element(Compute_t(r)*Compute_t(a) % Compute_t(p)) )
#define __GIVARO_MODULAR_FLOATING_MULIN(r,p,a) ( r = std::fmod((r*a), p) )

// r <- a - b
#define __GIVARO_MODULAR_INTEGER_SUB(r,p,a,b) \
( r = Element((Compute_t(a) >= Compute_t(b)) ? \
Compute_t(a)-Compute_t(b) : Compute_t(p)-Compute_t(b)+Compute_t(a)) )
#define __GIVARO_MODULAR_FLOATING_SUB(r,p,a,b) ( r = (a>=b) ? a-b: (p-b)+a )
// r <- r - a
#define __GIVARO_MODULAR_INTEGER_SUBIN(r,p,a) { \
	if (Compute_t(r) < Compute_t(a)) \
	r += Element(Compute_t(p)-Compute_t(a)); \
	else r = Element(Compute_t(r) - Compute_t(a)); }
#define __GIVARO_MODULAR_FLOATING_SUBIN(r,p,a) { if (r<a) r+=(p-a); else r-=a; }

// r <- a + b
#define __GIVARO_MODULAR_INTEGER_ADD(r,p,a,b) \
{ r = Element(Compute_t(a) + Compute_t(b)); \
r = Element((Compute_t(r) < Compute_t(p)) ? r : Compute_t(r) - Compute_t(p)); }
#define __GIVARO_MODULAR_FLOATING_ADD(r,p,a,b) { r = (a+b); r= (r < p ? r : r-p); }
// r <- r + a
#define __GIVARO_MODULAR_INTEGER_ADDIN(r,p,a) \
{ r += Element(a); r = Element((Compute_t(r) < Compute_t(p)) ? r : Compute_t(r) - Compute_t(p)); }
#define __GIVARO_MODULAR_FLOATING_ADDIN(r,p,a) { r += a;  r= (r < p ? r : r-p); }

// r <- c + a*b
#define __GIVARO_MODULAR_INTEGER_MULADD(r,p,a,b,c) \
( r = Element((Compute_t(a)*Compute_t(b)+Compute_t(c)) % Compute_t(p)))
#define __GIVARO_MODULAR_FLOATING_MULADD(r,p,a,b,c) ( r = std::fmod((a*b+c), p) )
// r <- r + a*b
#define __GIVARO_MODULAR_INTEGER_MULADDIN(r,p,a,b) \
( r = Element((Compute_t(a)*Compute_t(b)+Compute_t(r)) % Compute_t(p)))
#define __GIVARO_MODULAR_FLOATING_MULADDIN(r,p,a,b) ( r = std::fmod((a*b+r), p) )

// r <- a*b - c
#define __GIVARO_MODULAR_INTEGER_MULSUB(r,p,a,b,c) \
( r = Element((Compute_t(a)*Compute_t(b)+Compute_t(p)-Compute_t(c)) % Compute_t(p)))
#define __GIVARO_MODULAR_FLOATING_MULSUB(r,p,a,b,c) ( r = std::fmod((a*b+p-c), p) )

// r <- r - a*b
#define __GIVARO_MODULAR_INTEGER_SUBMULIN(r,p,a,b) \
{ r = Element((Compute_t(a)*Compute_t(b)+Compute_t(p)-Compute_t(r)) % Compute_t(p)); \
__GIVARO_MODULAR_INTEGER_NEGIN(r,p); }
#define __GIVARO_MODULAR_FLOATING_SUBMULIN(r,p,a,b) \
{ r = (a*b+p-r); r= (r<p ? r : std::fmod(r, p)); __GIVARO_MODULAR_FLOATING_NEGIN(r,p); }

#endif // __GIVARO_modular_defines_H

