// ==========================================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// Time-stamp: <08 Feb 02 16:33:39 Jean-Guillaume.Dumas@imag.fr>
// ==========================================================================
/** @file givintrns.h
 * @ingroup CRA
 * @brief arithmetic for RNS representations.
 *  Modular arithmetic for GIVARO. Here is defined arithmetic functions
 * on rns representation with Givaro Integers.
 */

#ifndef __GIVARO__arithmodu_intrns_H
#define __GIVARO__arithmodu_intrns_H

#include "givaro/givconfig.h"
#include "givaro/givinteger.h"

// ---------------------------------------------  class RNSsystem
// Structure which manages list of primes in order to do

namespace Givaro {
    // #ifndef __ECC

    //! RNS system class. No doc.
    template< template<class, class> class Container, template <class> class Alloc>
    class IntRNSsystem : public IntegerDom {
    public:
        //     typedef Element    Ring;
        //     typedef Element   Modulo;
        typedef Element   external;
        typedef Container< Element, Alloc<Element> > array;


        // Default Cstor, Dstor/Cstor of recopy:
        // -- free memory allocated in array !
        IntRNSsystem() : _primes(0), _prod(one), _ck(0) {}
        ~IntRNSsystem(){}
        IntRNSsystem(const IntRNSsystem& R) : _primes(R._primes), _prod(R._prod), _ck(R._primes) {}

        // -- Cstor with given primes
        IntRNSsystem( const array& primes );

        template<class TT>
        IntRNSsystem( const Container< TT, Alloc<TT> > & primes );

        // -- Computation of a mixed-radix representation of the residus.
        //     void RnsToMixedRadix(array&  mixrad, const array&  residu) const;
        template<class TT>
        void RnsToMixedRadix(array&  mixrad, const Container<TT, Alloc<TT> >&  residu) ;

        // -- Convert a mixed radix representation to an external
        void MixedRadixToRing( external& res,  const array& mixrad ) const;

        // -- Convert an Ring Element to a its representation
        // with the "this" rns system.
        void RingToRns( array& residu, const external& a ) ;

        // -- Fast conversion: requires pre-computation (first time it was called)
        void fastRingToRns( array& residu, const external& a ) const;

        // -- Convert a representation to an external Element
        template<class TT>
        void RnsToRing( external& a, const Container<TT, Alloc<TT> >& residu ) ;

        // -- Fast conversion: requires pre-computation (first time it was called)
        void fastRnsToRing( external& a, const array& residu ) const;

        // ------------- Access methods

        // -- Returns the number of primes of this ctxt
        int NumOfPrimes() const { return _primes.size(); }

        // -- Returns a array to the begin of the array of primes
        const array& Primes() const;
        // -- Returns the ith primes of the rns system
        const Element ith(const size_t i) const;

        // -- Returns a array of the reciprocal ck = (\prod_{j=0..k-1)p_j)^(-1) [pk]
        const array& Reciprocals() const;
        const Element reciprocal(const size_t i) const;
        const Element product() const;

    protected:
        // -- Compute some fields of the structure :
        void ComputeCk();

        // -- Compute product of primes
        void ComputeProd();

        // -- Compute the Qk for Ring -> RNS, allocate U
        void ComputeQk();

        array  _primes; 	// - array of the relatively primes numbers
        Element _prod;      // - product of primes
        array  _ck;     	// - reciprocals, _ck[0] = 1, same size as _primes

        // -- for fast conversion
        size_t _sizek;
        size_t _log2k;
        array  _qk;	// - cf algo Aho, Hopcroft & Ullman
        array  _u;	// - cf algo Aho, Hopcroft & Ullman
    };
    // #endif

#if 0 // defined(__ECC)
    //#else
    /* template<class Container> */
    /* class IntRNSsystem : public IntegerDom { */
    /* public: */
    /* //     typedef Element    Ring; */
    /* //     typedef Element   Modulo; */
    /*     typedef Element   external; */
    /*     typedef Container array; */


    /*         // Default Cstor, Dstor/Cstor of recopy:  */
    /* 	// -- free memory allocated in array ! */
    /*     IntRNSsystem() : _primes(0), _prod(one), _ck(0) {} */
    /*     ~IntRNSsystem(){}  */
    /*     IntRNSsystem(const IntRNSsystem& R) : _primes(R._primes), _prod(R._prod), _ck(R._primes) {} */

    /*         // -- Cstor with given primes  */
    /*     IntRNSsystem( const array& primes ); */

    /*     template<class ContTT> */
    /*     IntRNSsystem( const ContTT& primes ); */

    /*         // -- Computation of a mixed-radix representation of the residus. */
    /* //     void RnsToMixedRadix(array&  mixrad, const array&  residu) const;  */
    /*     template<class ContTT> */
    /*     void RnsToMixedRadix(array&  mixrad, const ContTT&  residu) const;  */

    /*         // -- Convert a mixed radix representation to an external */
    /*     void MixedRadixToRing( external& res,  const array& mixrad ) const; */

    /*         // -- Convert an Ring Element to a its representation */
    /*         // with the "this" rns system. */
    /*     void RingToRns( array& residu, const external& a ) const; */

    /*         // -- Fast conversion: requires pre-computation (first time it was called) */
    /*     void fastRingToRns( array& residu, const external& a ) const; */

    /*         // -- Convert a representation to an external Element */
    /*     template<class ContTT> */
    /*     void RnsToRing( external& a, const ContTT& residu ) const; */

    /*         // -- Fast conversion: requires pre-computation (first time it was called) */
    /*     void fastRnsToRing( external& a, const array& residu ) const; */

    /*         // ------------- Access methods */

    /*         // -- Returns the number of primes of this ctxt */
    /*     int NumOfPrimes() const { return _primes.size(); }  */

    /*         // -- Returns a array to the begin of the array of primes */
    /*     const array& Primes() const; */
    /*         // -- Returns the ith primes of the rns system */
    /*     const Element ith(const size_t i) const; */

    /*         // -- Returns a array of the reciprocal ck = (\prod_{j=0..k-1)p_j)^(-1) [pk] */
    /*     const array& Reciprocals() const; */
    /*     const Element reciprocal(const size_t i) const; */
    /*     const Element product() const; */

    /* protected: */
    /*         // -- Compute some fields of the structure : */
    /*     void ComputeCk(); */

    /*         // -- Compute product of primes */
    /*     void ComputeProd(); */

    /*         // -- Compute the Qk for Ring -> RNS, allocate U */
    /*     void ComputeQk(); */

    /*     array  _primes; 	// - array of the relatively primes numbers */
    /*     Element _prod;      // - product of primes */
    /*     array  _ck;     	// - reciprocals, _ck[0] = 1, same size as _primes  */

    /*         // -- for fast conversion */
    /*     size_t _sizek; */
    /*     size_t _log2k; */
    /*     array  _qk;	// - cf algo Aho, Hopcroft & Ullman */
    /*     array  _u;	// - cf algo Aho, Hopcroft & Ullman */
    /* }; */
    /* #endif */
#endif

    } // Givaro

#include "givaro/givintrns_cstor.inl"
#include "givaro/givintrns_convert.inl"

#endif // __GIVARO__arithmodu_intrns_H

    /* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
    // vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
