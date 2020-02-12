// ==========================================================================
// $Source
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id
// ==========================================================================
#include <givaro/gfq.h>

namespace Givaro {

    template<>
    inline void VectorDom<GFqDom<int64_t>,Dense>::dot
    ( Type_t& res, const Element& op1, const Element& op2) const
    {
        size_t sz = dim(op1);
        const GFqDom<int64_t>& domain = subdomain();
        domain.dotprod( res, sz, op1.baseptr(), op2.baseptr() );
    }


} // Givaro
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
