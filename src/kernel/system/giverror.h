// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/system/giverror.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: giverror.h,v 1.4 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================

/*! @file giverror.h
 * @ingroup system
 * @brief error exception.
 */

#ifndef __GIVARO_error_H
#define __GIVARO_error_H

#include <iostream>
namespace Givaro {
    // ------------------------------- GivError
    //! Base class for exeception handling in Givaro
    class GivError {
    public:
        GivError(const char* msg =0 )
        : strg(msg) {};

        virtual ~GivError() ;
        // -- virtual print of the error message
        virtual std::ostream& print( std::ostream& o )  const;

        // -- non virtual output operator
        friend std::ostream& operator<< (std::ostream& o, const GivError& E) ;

        // - useful to setup a break point on it
        static void throw_error( const GivError& err );

    protected:
        const char* strg;
    };

    //! Math error.
    class GivMathError : public GivError {
    public:
        virtual ~GivMathError() ;

        GivMathError(const char* msg = 0) : GivError(msg) { }
    };

    //! Exception thrown in input of data structure
    class GivBadFormat : public GivError {
    public:
        virtual ~GivBadFormat();
        GivBadFormat(const char* msg = 0) : GivError(msg) { }
    };

    //! Div by 0.
    class GivMathDivZero : public GivError {
    public:
        virtual ~GivMathDivZero();
        GivMathDivZero(const char* msg = 0) : GivError(msg) { }
    };

} // namespace Givaro

#endif // __GIVARO_error_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
