// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/bstruct/givstack.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.
// Author: T. Gautier
// $Id: givstack.h,v 1.2 2009-09-17 14:28:22 jgdumas Exp $
// ==========================================================================
#ifndef _STACK_H_
#define _STACK_H_

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

#include "givaro/givstack.inl"

#endif 

	
