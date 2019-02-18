// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/system/givbasictype.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givbasictype.h,v 1.4 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================

/** @file givbasictype.h
 * @ingroup system
 * @brief NO DOC
 */

#ifndef __GIVARO_basictype_H
#define __GIVARO_basictype_H
#include "givaro/givconfig.h"

#include <stdlib.h> // for size_t

#include <sys/types.h> // needed on MacOS X 10.5 for uint type

namespace Givaro {

    /** Neutral type.
     * definition of zero and one
     */
    class Neutral {
    public:
        static Neutral zero;
        static Neutral one;
        inline operator int() const { return _val; }
        inline int operator==( const Neutral& n) const { return _val==n._val; }
        inline int operator!=( const Neutral& n) const { return _val!=n._val; }
    private:
        Neutral( int val ) : _val(val) {};
        int _val;
    };

    //! Used to build no initialized object as static object
    class givNoInit {};
    //! Used to call cstor without copy
    class givNoCopy {};
    //! Used to call cstor with copy
    class givWithCopy {};

} // namespace Givaro

#endif // __GIVARO_basictype_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
