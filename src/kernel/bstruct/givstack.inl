// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/bstruct/givstack.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.
// Author: T. Gautier
// $Id: givstack.inl,v 1.2 2009-09-17 14:28:22 jgdumas Exp $
// ==========================================================================

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

