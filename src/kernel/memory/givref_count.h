// ==================================================================== //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// author: Th. Gautier
// version : 2.7
// date: 1995
// ==================================================================== //
/** @file givref_count.h
 * @ingroup memory
 * @brief Definition of the Counter class, Counter.
 * This class definition objects to handle reference
 * counter for memory allocation (eg array0).
 */

#ifndef __GIVARO_ref_counter_H
#define __GIVARO_ref_counter_H
#include <stddef.h>

namespace Givaro {

    //! Ref counter.
    class RefCounter {
    public:
        // Cstor and Dstor
        inline RefCounter( long l = 0) : counter(l) {}
        //inline RefCounter( const RefCounter& ) : counter(C.counter) {}
        inline ~RefCounter() {}

        //  Return the value
        inline long  getvalue() const { return counter ; }
        inline long  val() const { return counter ; }
        // Return a ref to the counter
        inline long& refvalue() { return counter ; }
        // Increments the counter and returns the new value
        inline long  incr() { return ++counter ; }
        // Decrements the value and returns the new value
        inline long  decr() { return --counter ; }

    protected:
        long counter ;
    } ;

} // namespace Givaro
#endif // __GIVARO_ref_counter_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
