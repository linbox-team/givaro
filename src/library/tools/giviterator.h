#ifndef _GIV_ITERATOR_H_
#define _GIV_ITERATOR_H_
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/tools/giviterator.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// Description:
// - Definition of traits for iterators.
// - The purpose of this trait class is for specialization of some algorithms
//   depending on the iteration mechanism provides by container.
// $Id: giviterator.h,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// It's a beta-beta version.
//

// -- isUndefined trait iterator 
class isUndefinedIterator{};

// --
// -- Forward iterator, should implement:
// - operator++(), operator++(int), operator+=(int): increment
// - operator*(): deference
// --
class isForwardIterator{}; 

// --
// -- BiDirectional iterator, should implement:
// - operator-- operator--int), operator-=(int): decrement
// --
class isBidirectionalIterator: public isForwardIterator{}; 

// --
// -- Random iterator, should implement: 
// - operator()(int), operator[](int): random access
// --
class isRandomIterator: public isBidirectionalIterator{}; 


// --
// -- Iterator trait: each iterator typename should provide its 
// -- categrory using this trait class.  
// --
template<class Iterator>
struct IteratorTraits {
  typedef isUndefinedIterator Category_t;
};


// --
// -- A container of name OO should provide one of the previous iterators:
// - * OO::Iterator: the name of the default iterator associated to OO. 
// - * OO::Iterator OO::begin(): 
// -- Depending on the trait associated to OO::Iterator, OO class should
// -- also provides the following interface:
// - * [ForwardIteratorTrait, BidirectionalIteratorTrait]:
// -   - OO::Iterator OO::end(): 
// - * [RandomIteratorTrait]:
// -   - size_t OO::bound(): return the number of elements in the sequence
// --
template<class Container>
struct IteratorInterface {
  typedef isUndefinedIterator 			Category_t;	// - category of iterator
  typedef typename Container::Iterator_t 	Iterator_t;	// - type of iterator
  typedef typename Container::constIterator_t	constIterator_t;// - type of constiterator
  typedef typename Container::Type_t 		Type_t; 	// - type of element 
  typedef typename Container::Indice_t 		Indice_t;	// - type of indice for RndIter

  // -- other operations that should be defined in specialized trait classes:
  // -- [depending of the category implemented by the container class].
  // -- If [ForwardIterator, BidirectionalIterator, RandomIterator]
  // * static Iterator_t begin(Container& cc);
  // * static const Iterator_t begin(const Container& cc);
  // -- If [ForwardIteratorTrait, BidirectionalIteratorTrait]
  // * static Iterator_t end(Container& cc);
  // * static const Iterator_t end(const Container& cc);
  // -- If [RandomIteratorTrait]
  // * static size_t bound(const Container& cc);
};


#endif
