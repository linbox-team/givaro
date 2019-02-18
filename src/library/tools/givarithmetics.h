// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/tools/givarithmetics.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givarithmetics.h,v 1.2 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
// Description:
// - definition of standard arithmetics over Vector & Matrix
#ifndef __GIVARO_arithmethics_H
#define __GIVARO_arithmethics_H
//
#include "givaro/givconfig.h"

namespace Givaro {

    template<class Domain>
    struct KSpace {

        // -- exported types
        typedef 	   Domain		Domain_t;
        typedef typename Domain::Rep		Entry_t;
        typedef typename Domain::Type_t	Scalar_t;

        typedef 	   KSpace<Domain>	Self_t;

        // -- Cstor
        KSpace() : _domain(), _entry() {};
        KSpace(const Domain_t& D) : _domain(D), _entry() {};
        KSpace(const Domain_t& D, size_t dim ) : _domain(D), _entry()
        {
            _domain.init(_entry, dim);
        };
        KSpace(const Self_t& S) : _domain(S._domain)
        {
            _domain.assign(_entry,S._entry);
        };

        // -- Destor: implicit

        // -- assignment:
        Self_t& operator= (const Self_t& e ) {
            e._domain.assign(_entry, e._entry);
            return *this;
        }

        // -- arithmetic operators:
        Self_t& operator+= ( const Self_t& a ) {
            _domain.addin( _entry, a._entry );
            return *this;
        }
        Self_t& operator-= ( const Self_t& a ) {
            _domain.subin( _entry, a._entry );
            return *this;
        }
        Self_t& operator*= ( const Scalar_t& v ) {
            typedef typename Domain::Domain_t SubDomain_t;      // domain of the scalar
            Curried2<MulOp<SubDomain_t> > opcode( _domain.subdomain(), v);
            _domain.map( _entry, opcode, _entry );
            return *this;
        }
        Self_t operator+ ( const Self_t& a ) const {
            KSpace<Domain> res(_domain);
            _domain.init( res._entry, _domain.dim(_entry) );
            _domain.add( res._entry, _entry, a._entry );
            return res;
        }
        Self_t operator- ( const Self_t& a ) const {
            KSpace<Domain> res(_domain);
            _domain.init( res._entry, _domain.dim(_entry) );
            _domain.sub( res._entry, _entry, a._entry );
            return res;
        }
        Self_t operator* ( const Scalar_t& v ) const {
            typedef typename Domain::Domain_t SubDomain_t;      // domain of the scalar
            KSpace<Domain> res(_domain);
            _domain.init( res._entry, _domain.dim(_entry) );
            Curried2<MulOp<SubDomain_t> > opcode( _domain.subdomain(), (Scalar_t&)v);
            _domain.map( res._entry, opcode, _entry );
            return res;
        }

        Scalar_t& operator[] (int i) { return _entry[i]; }
        const Scalar_t& operator[] (int i) const { return _entry[i]; }


        // -- representation
        const Domain_t _domain;
        Entry_t	 _entry;
    };

    // -- friend operator
    template<class Domain>
    KSpace<Domain> operator*
    ( const typename KSpace<Domain>::Scalar_t& v, const KSpace<Domain>& U )
    {
        typedef typename KSpace<Domain>::Scalar_t Scalar_t;
        typedef typename Domain::Domain_t SubDomain_t;      // domain of the scalar
        KSpace<Domain> res(U._domain);
        U._domain.init( res._entry, U._domain.dim(U._entry) );
        Curried1<MulOp<SubDomain_t> > opcode( U._domain.subdomain(), (Scalar_t&)v);
        U._domain.map( res._entry, opcode, U._entry );
        return res;
    }

    template<class Domain>
    inline void dotprod(
                        typename Domain::Scalar_t& dot,
                        const KSpace<Domain>& U,
                        const KSpace<Domain>& V)
    {
        U._domain.dot(dot, U._entry, V._entry);
    }

    template<class Domain>
    inline ostream& operator<<(ostream& sout,  const KSpace<Domain>& U )
    {
        return U._domain.write(sout, U._entry);
    }

    template<class Domain>
    inline istream& operator>>(istream& sin,  KSpace<Domain>& U )
    {
        return U._domain.read(sin, U._entry);
    }

} // Givaro

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
