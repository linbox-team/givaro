/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// ------------------------------------------------------------------------
// Time-stamp: <08 Oct 04 19:02:35 Jean-Guillaume.Dumas@imag.fr> 
// ------------------------------------------------------------------------

/* StaticElement.h
 * Copyright (C) 2004 Jean-Guillaume.Dumas
 *
 * Written by Jean-Guillaume.Dumas
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __STATIC_ELEMENT_H__
#define __STATIC_ELEMENT_H__

#include <iostream>
#include <gmp++/gmp++_int.h>

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
#ifdef __USE_GMPPLUSPLUS_64__
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
#endif // __STATIC_ELEMENT_H__
