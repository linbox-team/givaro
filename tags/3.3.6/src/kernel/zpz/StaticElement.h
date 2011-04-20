//==================================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// -----------------------------------------------------------------
// Time-stamp: <19 Dec 06 10:55:08 Jean-Guillaume.Dumas@imag.fr>
// -----------------------------------------------------------------
// file: StaticElement.h
// author: Jean-Guillaume.Dumas
// date: 2004
//==================================================================

#ifndef __GIVARO_static_element_H
#define __GIVARO_static_element_H

#include <iostream>
#include <gmp++/gmp++.h>

template <class Domain>
class StaticElement {
    static Domain _domain;
    typedef typename Domain::Element Rep;
    Rep _elem;
public:
    static void setDomain(const Domain& D) {
        StaticElement::_domain = D;
    }

    static Domain getDomain() {
        return StaticElement::_domain;
    }

    StaticElement()  { _domain.init(_elem); }


    StaticElement(const Integer& i) { _domain.init( (_elem), i); }
    StaticElement(const double& i) { _domain.init( (_elem), i); }
    StaticElement(const int& i) { _domain.init( (_elem), (long)i); }
    StaticElement(const unsigned int& i) { _domain.init( (_elem), (unsigned long)i); }
    StaticElement(const long& i) { _domain.init( (_elem), i); }
    StaticElement(const unsigned long& i) { _domain.init( (_elem), i); }


    operator short() const		{ long tmp; return (short)_domain.convert(tmp,_elem); }
    operator unsigned short() const 	{ unsigned long tmp; return (unsigned short)_domain.convert(tmp,_elem); }
    operator unsigned char() const 	{ unsigned long tmp; return (unsigned char)_domain.convert(tmp,_elem); }
    operator unsigned int() const	{ unsigned long tmp; return (unsigned int)_domain.convert(tmp,_elem); }
    operator int() const		{ long tmp; return (int)_domain.convert(tmp,_elem); }
    operator float() const		{ double tmp; return (float)_domain.convert(tmp,_elem); }
    operator unsigned long() const	{ unsigned long tmp; return _domain.convert(tmp,_elem); }
    operator long() const		{ long tmp; return _domain.convert(tmp,_elem); }
    operator double() const		{ double tmp; return (double)_domain.convert(tmp,_elem); }
    operator Integer() const		{ Integer tmp; return _domain.convert(tmp,_elem); }
#ifndef __GIVARO__DONOTUSE_longlong__
    operator unsigned long long() const	{ Integer tmp; return (unsigned long long)_domain.convert(tmp,_elem); }
    operator long long() const		{ Integer tmp; return (long long)_domain.convert(tmp,_elem); }
#endif
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

    bool operator==(const StaticElement& e) {
        return _domain.areEqual((_elem), e._elem);
    }

    bool operator!=(const StaticElement& e) {
        return !_domain.areEqual((_elem), e._elem);
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

    friend inline std::istream& operator>> (std::istream& i, StaticElement& a) {
        return _domain.read(i,(a._elem));
    }


    friend inline std::ostream& operator<< (std::ostream& o, const StaticElement a) {
        return _domain.write(o,(a._elem));
    }

};
#endif // __GIVARO_static_element_H
