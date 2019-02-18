// ==========================================================================
// Copyright(c)'2017 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: J. G. Dumas
// ==========================================================================
/** @file givprint.h
 * helper print for containers
 */
#ifndef __GIVARO_print_H
#define __GIVARO_print_H

#include <iostream>
#include <iterator>
#include <vector>
#include <list>
#include <set>

namespace std
{

    /*! Prints a vector on output.
     * @param o output stream
     * @param v vector
     * @warning <<(ostream&,T&) exists !
     */
    template <class T, typename A=std::allocator<T> >
    std::ostream& operator<<(std::ostream& out, const std::vector<T, A>& V)
    {
        std::copy(V.begin(), V.end(), std::ostream_iterator<T>(out << '['," "));
        return out << ']';
    }

    /*! Prints a pair.
     * @param o output stream
     * @param C a pair
     * @warning <<(ostream&,T&) exists !
     */
    template<class S, class T>
    std::ostream& operator<<(std::ostream& out, const std::pair<S, T> & C)
    {
        return out << '<' << C.first << ", " << C.second << '>';
    }


    /*! Prints a list.
     * @param o output stream
     * @param C a pair
     * @warning <<(ostream&,T&) exists !
     */
    template<class T, typename A=std::allocator<T> >
    std::ostream& operator<< (std::ostream& out, const std::list<T, A> & L)
    {
        std::copy(L.begin(), L.end(), std::ostream_iterator<T>(out << '('," "));
        return out << ')' ;
    }


    /*! Prints a set.
     * @param o output stream
     * @param C a pair
     * @warning <<(ostream&,T&) exists !
     */
    template<class T, typename A=std::allocator<T> >
    std::ostream& operator<< (std::ostream& out, const std::set<T, A> & S)
    {
        std::copy(S.begin(), S.end(), std::ostream_iterator<T>(out << '{'," "));
        return out << '}' ;
    }

}
#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
