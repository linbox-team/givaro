// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/bstruct/givstack.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Author: T. Gautier
// $Id: givstack.h,v 1.3 2011-02-02 16:23:55 briceboyer Exp $
// ==========================================================================
/** @file givstack.h
 * @ingroup bstruct
 * @brief no doc.
 */
#ifndef __GIVARO_stack_H
#define __GIVARO_stack_H

namespace Givaro {

    //! Stack
    template <class THING>
    class Stack {
    public :
        inline Stack() ;
        inline ~Stack() ;

        inline void push(const THING&) ;
        inline void pop() ;

        inline int isEmpty() const { return (ThePointer == NULL) ; }
        inline THING top() const ;

    private :
        struct inner_stack
        {
            THING thething ;
            inner_stack * next ;
        } *ThePointer ;
    } ;

} // Givaro

#include "givaro/givstack.inl"

#endif // __GIVARO_stack_H


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
