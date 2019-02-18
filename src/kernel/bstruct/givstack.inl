// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/bstruct/givstack.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Author: T. Gautier
// $Id: givstack.inl,v 1.3 2011-02-02 16:23:55 briceboyer Exp $
// ==========================================================================

#ifndef __GIVARO_stack_INL
#define __GIVARO_stack_INL

namespace Givaro {

    template <class THING>
    Stack<THING>::~Stack()
    { }

    template <class THING>
    Stack<THING>::Stack()
    {
        ThePointer = NULL ;
    }

    template <class THING>
    void Stack<THING>::push(const THING& T)
    {
        struct inner_stack * Newpt ;
        Newpt = new struct inner_stack ;
        Newpt->thething = T ;
        Newpt->next = ThePointer ;
        ThePointer = Newpt ;
    }

    template <class THING>
    void Stack<THING>::pop()
    {
        if (ThePointer == NULL)
        {
            cerr << "*** Error: Empty Stack" << endl ;
        }
        else {
            inner_stack* tmp = ThePointer ;
            ThePointer = ThePointer->next ;
            delete tmp ;
        }
    }

    template <class THING>
    THING Stack<THING>::top() const
    {
        return ThePointer->thething ;
    }

} // namespace Givaro


#endif // __GIVARO_stack_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
