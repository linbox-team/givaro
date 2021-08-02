//==================================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// -----------------------------------------------------------------
// Time-stamp: <24 Nov 17 17:27:34 Jean-Guillaume.Dumas@imag.fr>
// -----------------------------------------------------------------
// author: Jean-Guillaume.Dumas
// date: 2004
//==================================================================

/*! @file StaticElement.h
 * @ingroup zpz
 * @brief NO DOC
 */

#ifndef __GIVARO_static_element_H
#define __GIVARO_static_element_H

#include <iostream>
#include "gmp++/gmp++.h"

namespace Givaro {

	//! Static Element
template <class DomainStyle>
struct StaticElement {
    typedef DomainStyle Domain;
protected:
    static Domain _domain;
    typedef typename Domain::Element Rep;
    Rep _elem;
public:
    static void setDomain(const Domain& D) {
        StaticElement::_domain = D;
    }

    static const Domain & getDomain() {
        return StaticElement::_domain;
    }

    StaticElement()  { _domain.init(_elem); }


    StaticElement(const Integer& i) { _domain.init( (_elem), i); }
    StaticElement(const double& i) { _domain.init( (_elem), i); }
    StaticElement(const int32_t& i) { _domain.init( (_elem), (int64_t)i); }
    StaticElement(const uint32_t& i) { _domain.init( (_elem), (uint64_t)i); }
    StaticElement(const int64_t& i) { _domain.init( (_elem), i); }
    StaticElement(const uint64_t& i) { _domain.init( (_elem), i); }


    operator short() const		{ int64_t tmp; return (short)_domain.convert(tmp,_elem); }
    operator unsigned short() const 	{ uint64_t tmp; return (unsigned short)_domain.convert(tmp,_elem); }
    operator unsigned char() const 	{ uint64_t tmp; return (unsigned char)_domain.convert(tmp,_elem); }
    operator uint32_t() const	{ uint64_t tmp; return (uint32_t)_domain.convert(tmp,_elem); }
    operator int() const		{ int64_t tmp; return (int)_domain.convert(tmp,_elem); }
    operator float() const		{ double tmp; return (float)_domain.convert(tmp,_elem); }
    operator uint64_t() const	{ uint64_t tmp; return _domain.convert(tmp,_elem); }
    operator int64_t() const		{ int64_t tmp; return _domain.convert(tmp,_elem); }
    operator double() const		{ double tmp; return (double)_domain.convert(tmp,_elem); }
    operator Integer() const		{ Integer tmp; return _domain.convert(tmp,_elem); }
    template<class INITCST>
    operator INITCST() const 		{ INITCST tmp; return _domain.convert(tmp,_elem); }


    template<class INITCST>
    StaticElement(const INITCST& i) { _domain.init( (_elem), i); }

    template<class INITCST>
    StaticElement(const INITCST& i, const Domain& D) {
        setDomain(D);
        _domain.init((_elem), i);
    }

    StaticElement& operator=( const StaticElement& e) {
        _domain.assign((_elem), e._elem); return  *this;
    }

    bool operator==(const StaticElement& e) const {
        return _domain.areEqual((_elem), e._elem);
    }

    bool operator!=(const StaticElement& e) const {
        return !_domain.areEqual((_elem), e._elem);
    }

    bool isZero() const {
        return StaticElement::_domain.isZero(this->_elem);
    }

    static bool isZero(const StaticElement& e) {
        return e.isZero();
    }

    bool isOne() const {
        return StaticElement::_domain.isOne(this->_elem);
    }

    static bool isOne(const StaticElement& e) {
        return e.isOne();
    }

    bool isMOne() const {
        return StaticElement::_domain.isMOne(this->_elem);
    }

    static bool isMOne(const StaticElement& e) {
        return e.isMOne();
    }


    inline const StaticElement operator* (const StaticElement& e) const {
        StaticElement tmp;
        StaticElement(_domain.mul(tmp._elem,(_elem), e._elem));
        return tmp;
    }
    inline const StaticElement operator/ (const StaticElement& e) const {
        StaticElement tmp;
        StaticElement(_domain.div(tmp._elem,(_elem), e._elem));
        return tmp;
    }

    inline const StaticElement operator+ (const StaticElement& e) const {
        StaticElement tmp;
        StaticElement(_domain.add(tmp._elem,(_elem), e._elem));
        return tmp;
    }

    inline const StaticElement operator- (const StaticElement& e) const {
        StaticElement tmp;
        StaticElement(_domain.sub(tmp._elem,(_elem), e._elem));
        return tmp;
    }

    inline const StaticElement operator- () const {
        StaticElement tmp;
        StaticElement(_domain.neg(tmp._elem,(_elem)));
        return tmp;
    }



    inline StaticElement& operator*= (const StaticElement& e) {
        _domain.mulin((_elem), e._elem); return *this;
    }

    inline StaticElement& operator/= (const StaticElement& e) {
        _domain.divin((_elem), e._elem); return *this;
    }

    inline StaticElement& operator+= (const StaticElement& e) {
        _domain.addin((_elem), e._elem); return *this;
    }

    inline StaticElement& operator-= (const StaticElement& e) {
        _domain.subin((_elem), e._elem); return *this;
    }

    inline StaticElement& operator-- () {
        _domain.subin(this->_elem, _domain.one); return *this;
    }

    inline StaticElement& operator++ () {
        _domain.addin(this->_elem, _domain.one); return *this;
    }

    inline StaticElement operator-- (int) {
        StaticElement tmp(*this);
        _domain.subin(this->_elem, _domain.one); return tmp;
    }

    inline StaticElement operator++ (int) {
        StaticElement tmp(*this);
        _domain.addin(this->_elem, _domain.one); return tmp;
    }

    friend inline std::istream& operator>> (std::istream& i, StaticElement& a) {
        return _domain.read(i,(a._elem));
    }


    friend inline std::ostream& operator<< (std::ostream& o, const StaticElement& a) {
        return _domain.write(o,(a._elem));
    }

};

} // namespace Givaro

#endif // __GIVARO_static_element_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
