// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/bstruct/givstack.inl,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Author: T. Gautier
// $Id: givstack.inl,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
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

