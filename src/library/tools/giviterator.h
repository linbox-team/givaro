// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/tools/giviterator.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// Description:
// - Definition of traits for iterators.
// - The purpose of this trait class is for specialization of some algorithms
//   depending on the iteration mechanism provides by container.
// $Id: giviterator.h,v 1.3 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
// It's a beta-beta version.
#ifndef __GIVARO_iterator_H
#define __GIVARO_iterator_H

namespace Givaro {
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
    // -   - size_t OO::bound(): return the number of Elements in the sequence
    // --
    template<class Container>
    struct IteratorInterface {
        typedef isUndefinedIterator 			Category_t;	// - category of iterator
        typedef typename Container::Iterator_t 	Iterator_t;	// - type of iterator
        typedef typename Container::constIterator_t	constIterator_t;// - type of constiterator
        typedef typename Container::Type_t 		Type_t; 	// - type of Element
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

} // Givaro

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
