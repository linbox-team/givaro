// ==========================================================================
// Copyright(c)'2012 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// file: givcaster.h
// Time-stamp: <21 May 12 15:35:16 Jean-Guillaume.Dumas@imag.fr>
// ==========================================================================

/** @file givcaster.h
 * @ingroup system
 * @brief NO DOC
 */

#ifndef __GIVARO_givcaster_H
#define __GIVARO_givcaster_H

namespace Givaro {
    template <typename Target, typename Source>
    Target& Caster (Target& t, const Source& s) {
        return t = static_cast<Target>(s);
    }
}

#endif // __GIVARO_caster_H
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
