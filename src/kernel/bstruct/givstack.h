#ifndef _STACK_H_
#define _STACK_H_
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/bstruct/givstack.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Author: T. Gautier
// $Id: givstack.h,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================

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

	
