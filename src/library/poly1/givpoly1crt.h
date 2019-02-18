// ==========================================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: J-G Dumas
// Time-stamp: <06 May 10 13:47:28 Jean-Guillaume.Dumas@imag.fr>
// ==========================================================================

/** @file givpoly1crt.h
 * @ingroup poly1
 * @brief Polynomial Chinese Remaindering of degree 1
 */

#ifndef __GIVARO_poly1_crt_H
#define __GIVARO_poly1_crt_H
#include <givaro/givpoly1.h>
#include <givaro/givindeter.h>
#include <vector>

namespace Givaro {

    //! Poly1 CRT
    template<class Field>
    class Poly1CRT  {
        typedef Poly1CRT<Field> 			Self_t;
    public:
        typedef Field				Field_t;
        typedef Poly1Dom<Field, Dense>		Ring_t;
        typedef typename Field_t::Element 		Type_t;
        typedef typename Ring_t::Element 		Element;
        typedef Element 				Rep;
        typedef std::vector<Type_t> 			array_T;
        typedef std::vector<Element> 			array_E;

        // Default Cstor, Dstor/Cstor of recopy:
        Poly1CRT() ;
        ~Poly1CRT();
        Poly1CRT( const Self_t& R);

        // -- Cstor with given residues so that irreds are (X-primes[i])
        Poly1CRT( const Field& F, const array_T& primes, const Indeter& X = Indeter() );

        // -- Convert a Ring Element to a its RNS representation
        // with the "this" rns system.
        array_T& RingToRns( array_T& rns, const Element& a ) const;

        // -- Convert a RNS representation to a RING Element
        Element& RnsToRing( Element& a, const array_T& rns );

        // ------------- Access methods

        // -- Returns the number of primes of this ctxt
        int size() const { return _primes.size(); }

        // -- Returns a array to the begin of the array of primes
        const array_T& Primes() const;
        // -- Returns the ith primes of the rns system
        const Type_t& ith(const size_t i) const;

        // -- Returns an array of the reciprocal ck = (\prod_{j=0..k-1)p_j)^(-1) [pk]
        const array_E& Reciprocals() const;
        const Element& reciprocal(const size_t i) const;

        const Field_t& getdomain() { return _F; }
        const Ring_t&  getpolydom() { return _PolRing; }

        std::ostream& write( std::ostream& o ) const {
            return _PolRing.write(o << "CRT(") << ')';
        }

        std::istream& read ( std::istream& i, Element& n) const {
            return _PolRing.read(i,n);
        }

        std::ostream& write( std::ostream& o, const Element& n) const {
            return _PolRing.write(o,n);
        }

        std::istream& read ( std::istream& i, Type_t& n) const {
            return _F.read(i,n);
        }

        std::ostream& write( std::ostream& o, const Type_t& n) const {
            return _F.write(o,n);
        }



    protected:
        // -------------- Compute some fields of the structure :
        void ComputeCk();

        const Indeter  _XIndet;
        const Field_t& _F;
        const Ring_t   _PolRing;
        array_T   _primes; 	// - array of the primes
        array_E       _ck;  // - Radix list reciprocals
    };
    } // Givaro

#include "givaro/givpoly1crtcstor.inl"
#include "givaro/givpoly1crtconvert.inl"

#endif
    /* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
    // vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
