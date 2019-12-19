// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givpoly1io.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givpoly1io.inl,v 1.7 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================
// Description:
#ifndef __GIVARO_poly1_io_INL
#define __GIVARO_poly1_io_INL

#include <iostream>

namespace Givaro {
    // --
    template<class Domain>
    std::istream& Poly1Dom<Domain,Dense>::read ( std::istream& sin )
    {
        char ch;
        sin >> std::ws >> ch;
#ifdef __GIVARO_DEBUG
        if (ch != '(')
            GivError::throw_error(
                                  GivBadFormat("Poly1Dom<Domain,Dense>::read: syntax error no '('"));
#endif

        _domain.read(sin);

        sin >> std::ws >> ch;
#ifdef __GIVARO_DEBUG
        if (ch != ',')
            GivError::throw_error(
                                  GivBadFormat("Poly1Dom<Domain,Dense>::read: syntax error no ','"));
#endif

        sin >> _x;

        sin >> std::ws >> ch;
#ifdef __GIVARO_DEBUG
        if (ch != ')')
            GivError::throw_error(
                                  GivBadFormat("Poly1Dom<Domain,Dense>::read: syntax error no ')'"));
#endif
        return sin;
    }

    template<class Domain>
    std::ostream& Poly1Dom<Domain,Dense>::write( std::ostream& o ) const
    {
        return _domain.write(o) << '[' << _x << ']';
    }



    template<class Domain>
    std::ostream& Poly1Dom<Domain,Dense>::write( std::ostream& o, const Rep& R) const
    {
        if (R.size()) {
            Rep P; assign(P, R);
            setdegree(P);
            if (P.size()) {
                if (! _domain.isZero(P[0])) {
                    if (_domain.isOne(P[0]))
                        _domain.write(o,P[0]);
                    else
                        _domain.write(o << "(",P[0]) << ")";
                }
                if (P.size() > 1) {
                    if (! _domain.isZero(P[0])) o << " + ";
                    if (! _domain.isZero(P[1])) {
                        if (! _domain.isOne(P[1])) {
                            _domain.write(o << "(",P[1]) << ")*";
                        }
                        o << _x;
                    }
                    for(unsigned long l=2;l<P.size();++l) {
                        if (! _domain.isZero(P[l-1])) o << " + ";
                        if (! _domain.isZero(P[l])) {
                            if (! _domain.isOne(P[l])) {
                                _domain.write(o << "(",P[l]) << ")*";
                            }
                            o << _x << "^" << l;
                        }
                    }
                }
                return o;
            }
        }
        return o << "0";
    }

    template<class Domain>
    std::istream& Poly1Dom<Domain,Dense>::read ( std::istream& i, Rep& P) const
    {
        long deg;
        i >> deg;
        init(P,Degree(deg));
        // JGD 18.09.2002
        for(;deg>=0;--deg)
            _domain.read( i, P[(size_t)deg]);
        // i >> P[deg];
        return i;
    }
} // Givaro
#endif // __GIVARO_poly1_io_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
