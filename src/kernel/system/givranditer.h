// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: Giorgi Pascal <pascal.giorgi@ens-lyon.fr>
// ==========================================================================

/** @file givranditer.h
 * @ingroup zpz
 * @brief NO DOC
 * Givaro ring Elements generator
 */

#ifndef __GIVARO_randiter_H
#define __GIVARO_randiter_H

#include "givaro/givconfig.h"
#include "givaro/givrandom.h"
#include "givaro/givtimer.h"

// For ModularBalancedRandIter
#include <sys/time.h>
#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE
#endif

#include <stdlib.h>
#include <limits>

namespace Givaro {

    /** Random ring Element generator.
     *   This class defines a ring Element generator for all givaro ring (Gfq and Zpz)
     *   throught a template argument as a ring.
     *   The random generator used is the givrandom.
     */
    template <class Ring , class Type>
    class GIV_randIter
    {
    public:

        /** @name Common Object Interface.
         * These methods are required of all LinBox random ring Element generators.
         */
        //@{

        /** Ring Element type.
         * The ring Element must contain a default constructor,
         * a copy constructor, a destructor, and an assignment operator.
         */
        typedef typename Ring::Element Element;
        typedef typename Ring::Residu_t Residu_t;

        /** Constructor from ring, sampling size, and seed.
         * The random ring Element iterator works in the ring F, is seeded
         * by seed, and it returns any one Element with probability no more
         * than 1/min(size, F.cardinality()).
         * A sampling size of zero means to sample from the entire ring.
         * A seed of zero means to use some arbitrary seed for the generator.
         * This implementation sets the sampling size to be no more than the
         * cardinality of the ring.
         * @param F LinBox ring archetype object in which to do arithmetic
         * @param size constant integer reference of sample size from which to
         *             sample (default = F.cardinality())
         * @param seed constant integer reference from which to seed random number
         *             generator (default = 0)
         */
        GIV_randIter(const Ring& F,
                     const uint64_t seed = 0,
                     const Residu_t size = 0)
        : _ring(F), _size(size?size:std::max(F.cardinality(),Residu_t(1))), _givrand(seed)
        {}

        /** Copy constructor.
         * Constructs ALP_randIter object by copying the random ring
         * Element generator.
         * This is required to allow generator objects to be passed by value
         * into functions.
         * In this implementation, this means copying the random ring Element
         * generator to which R._randIter_ptr points.
         * @param  R ALP_randIter object.
         */
        GIV_randIter(const GIV_randIter& R)
        : _ring(R._ring), _size(R._size), _givrand(R._givrand) {}

        /** Destructor.
         * This destructs the random ring Element generator object.
         * In this implementation, this destroys the generator by deleting
         * the random generator object to which _randIter_ptr points.
         */
        ~GIV_randIter(void) {}

        /** Assignment operator.
         * Assigns ALP_randIter object R to generator.
         * In this implementation, this means copying the generator to
         * which R._randIter_ptr points.
         * @param  R ALP_randIter object.
         */
        GIV_randIter<Ring,Type>& operator=(const GIV_randIter<Ring,Type>& R)
        {
            if (this != &R) // guard against self-assignment
            {
                _givrand = R._givrand;
                const_cast<Ring&>(_ring) = R._ring;
            }

            return *this;
        }

        /** Random ring Element creator with assignement.
         * This returns a random ring Element from the information supplied
         * at the creation of the generator.
         * @return random ring Element
         */
        Element& operator()(Element& elt) const
        {
            return ring().random (_givrand, elt, _size);
        }
        Element& random(Element& elt) const
        {
            return this->operator()(elt);
        }
        Element operator()() const
        {
            Element tmp;
            return this->operator()(tmp);
        }
        Element random() const
        {
            return this->operator()();
        }

        const Ring& ring() const { return _ring; }



        //@} Common Object Iterface

        /** @name Implementation-Specific Methods.
        */
        //@{

        //     /// Default constructor
        //   GIV_randIter(void) : _size(0), _givrand(), _ring() {}

        //@}

    private:

        /// Ring
        const Ring& _ring;

        /// Random generator
        const Residu_t _size;
        GivRandom _givrand;


    }; //  class GIV_randIter

    /** Random ring Element generator.
     *   This class defines a ring Element generator for all givaro modular rings (Gfq and Modular)
     *   throught a template argument as a ring.
     *   The random generator used is the givrandom.
     */
    template <class Ring>
    class ModularRandIter
    {
    public:

        /** @name Common Object Interface.
         * These methods are required of all LinBox random ring Element generators.
         */
        //@{

        /** Ring Element type.
         * The ring Element must contain a default constructor,
         * a copy constructor, a destructor, and an assignment operator.
         */
        typedef typename Ring::Element Element;
        typedef typename Ring::Residu_t Residu_t;

        /** Constructor from ring, sampling size, and seed.
         * The random ring Element iterator works in the ring F, is seeded
         * by seed, and it returns any one Element with probability no more
         * than 1/min(F.cardinality()).
         * size has no meaning in ModularRandIter.
         * A seed of zero means to use some arbitrary seed for the generator.
         * @param F LinBox ring archetype object in which to do arithmetic
         * @param seed constant integer reference from which to seed random number
         *             generator (default = 0)
         */
        ModularRandIter(const Ring& F, const uint64_t seed = 0, const Residu_t size = 0)
        : _givrand(seed), _size(size?size:F.cardinality()), _ring(F) {}

        /** Copy constructor.
         * Constructs ALP_randIter object by copying the random ring
         * Element generator.
         * This is required to allow generator objects to be passed by value
         * into functions.
         * In this implementation, this means copying the random ring Element
         * generator to which R._randIter_ptr points.
         * @param  R ALP_randIter object.
         */
        ModularRandIter(const ModularRandIter& R)
        : _givrand(R._givrand) , _size(R._size), _ring(R._ring) {}

        /** Destructor.
         * This destructs the random ring Element generator object.
         * In this implementation, this destroys the generator by deleting
         * the random generator object to which _randIter_ptr points.
         */
        ~ModularRandIter(void) {}

        /** Assignment operator.
         * Assigns ALP_randIter object R to generator.
         * In this implementation, this means copying the generator to
         * which R._randIter_ptr points.
         * @param  R ALP_randIter object.
         */
        ModularRandIter<Ring>& operator=(const ModularRandIter<Ring>& R)
        {
            // guard against self-assignment
            if (this != &R)
            {
                _givrand = R._givrand;
                const_cast<Ring&>(_ring) = R._ring;
            }
            return *this;
        }

        /** Random ring Element creator with assignement.
         * This returns a random ring Element from the information supplied
         * at the creation of the generator.
         * @return random ring Element
         */
        Element& operator()(Element& elt) const
        {
            // Create new random Elements
            return ring().random(_givrand, elt);
        }

        Element& random(Element& elt) const
        {
            return this->operator()(elt);

        }

        Element operator()() const
        {
            Element tmp;
            return this->operator()(tmp);
        }

        Element random() const
        {
            return this->operator()();
        }

        //@} Common Object Iterface

        /** @name Implementation-Specific Methods.
        */
        //@{

        //       /// Default constructor
        //     ModularRandIter(void) : _givrand(), _ring() {}

        //@}

        const Ring& ring() const { return _ring; }

    private:

        /// Random generator
        GivRandom _givrand;
        const Residu_t _size;

        /// Ring
        const Ring& _ring;

    }; //  class ModularRandIter

    /** UnparametricRandIter
     *
     * Imported from FFLAS-FFPACK UnparametricRandIter - AB 2015-01-12
     **/

    template <class Ring>
    class GeneralRingRandIter {
    public:
        typedef typename Ring::Element Element;
        typedef typename Ring::Residu_t Residu_t;

        GeneralRingRandIter(const Ring &F, uint64_t seed = 0, const Residu_t size = 0) : _F(F), _size(size?size:F.cardinality()), _givrand(seed)
        {}
        GeneralRingRandIter(const GeneralRingRandIter<Ring> &R) : _F(R._F), _size(R._size), _givrand(R._givrand) {}
        ~GeneralRingRandIter() {}

        Element& operator() (Element& a) const
        {
                // If no size given and cardinality is 0
            return ring().init(a, _size? _givrand() % static_cast<GivRandom::random_t>(_size) :_givrand());
        }
        Element& random (Element& a) const
        {
            return this->operator()(a);
        }
        Element operator() () const
        {
            Element a; return this->operator()(a);
        }
        Element random () const
        {
            return this->operator()();
        }

        const Ring& ring() const { return _F; }

    private:
        const Ring& _F;
        const Residu_t _size;
        /// Random generator
        GivRandom _givrand;

    };


    /** Random iterator for nonzero random numbers
     *
     * Wraps around an existing random iterator and ensures that the output
     * is entirely nonzero numbers.
     * Imported from FFLAS-FFPACK NonzeroRandIter - AB 2015-01-12
     **/
    template <class Ring, class RandIter = typename Ring::RandIter>
    class GeneralRingNonZeroRandIter
    {
    public:
        typedef typename Ring::Element Element;
        typedef typename RandIter::Residu_t Residu_t;

        GeneralRingNonZeroRandIter(RandIter& r) : _r(r) {}
        GeneralRingNonZeroRandIter(const GeneralRingNonZeroRandIter& R) : _r(R._r) {}
        ~GeneralRingNonZeroRandIter() {}

        Element& operator()(Element &a)  const
        {
            do _r.random(a); while ( ring().isZero(a));
            return a;
        }
        Element& random(Element &a) const
        {
            return this->operator()(a);
        }

        Element operator()()
        {
            Element a; return this->operator()(a);
        }
        Element random() const
        {
            return this->operator()();
        }

        const Ring& ring() const { return _r.ring(); }

    private:
        RandIter& _r;
    };

} // namespace Givaro

#endif // __GIVARO_randiter_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
