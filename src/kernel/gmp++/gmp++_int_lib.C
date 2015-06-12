// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_lib.C,v $
// Copyright(c)'1994-2011 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: B Boyer
// $Id: gmp++_int_lib.C,v 1.4 2009-09-17 14:28:22 jgdumas Exp $
// ==========================================================================
/** @file gmp++/gmp++_int_lib.C
 * libing stuff.
 */


#ifndef __GMPplusplus_LIB_C__
#define __GMPplusplus_LIB_C__
#include "gmp++/gmp++.h"
namespace Givaro
{
	//------------------------------------- predefined null and one
	const Integer Integer::zero(0U);
	const Integer Integer::one(1U);
	const Integer Integer::mOne(-1L);

} // Givaro
#endif

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen

