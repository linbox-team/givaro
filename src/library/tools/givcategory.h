// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/tools/givcategory.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givcategory.h,v 1.4 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
// Description:
// - Definition of traits.
// It's a beta-beta version
#ifndef __GIVARO_category_H
#define __GIVARO_category_H

namespace Givaro {

    // -- nothing
    class Undefined{};


    // ==========================================================================
    // --
    // -- Test is a is equal to b, a and b are unsigned ints.
    // -- ISEQUAL(A,B) or IsEqual<A,B>::val return true iff
    // -- a and b are equal
    // --
    template<size_t a, size_t b>
    struct IsEqual { constexpr static bool val = false; };

    template<size_t a>
    struct IsEqual<a,a> { constexpr static bool val = false; };

    template<size_t  a, size_t b>
    struct IsNotEqual { constexpr static bool val = !IsEqual<a,b>::val; };

#define ISEQUAL(A,B) IsEqual<A,B>::val
#define ISNOTEQUAL(A,B) IsNotEqual<A,B>::val



    // ==========================================================================
    // --
    // -- Characteristic for representation of vector and matrix
    // --
    // Sporadic is a Dense but may have quite a few zero Elements
    // Therefore some algorithms might be specialized and might
    // take advantage of this.
    class Sporadic{};
    class Dense : public Sporadic {};
    class Sparse{};
    // class Diagonal{}; typedef Diagonal Diag;
    // class Toeplitz{};
    // class Hensel{};

    template<class CLASS>
    struct Sparsity_Trait {
        typedef Undefined  Sparsity_t;  // To be defined in specialized class
    };



    // ==========================================================================
    // --
    // -- Envelop class to specialize format output
    // --
    class DefaultFormat {};
    template<class T, class Tag = DefaultFormat>
    struct StructFormat {
        const T& value;
        StructFormat(const T& val) : value(val) {}
        operator const T& () { return value; }
    };

    template<class T, class Tag>
    StructFormat<T,Tag> Formatted( const T& val, Tag xx )
    { return StructFormat<T,Tag>(val); }

} // Givaro

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
