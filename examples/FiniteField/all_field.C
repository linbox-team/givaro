// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*! @file examples/FiniteField/all_field.C
 * @ingroup examples
 * @ingroup finitefields
 * @example examples/FiniteField/all_field.C
 * @brief NO DOC
 */

#include <iostream>
#include <givaro/gfq.h>
#include <givaro/montgomery.h>
#include <givaro/modular.h>
#include <givaro/StaticElement.h>

using namespace Givaro;


namespace Givaro {
    // Domain kind
    typedef Modular<uint32_t>	Field1;	typedef StaticElement< Field1 > Element1; 	template<> Field1 Element1::_domain(2);
    typedef GFqDom<int64_t>		Field2;	typedef StaticElement< Field2 > Element2;	template<> Field2 Element2::_domain(2);
    typedef Montgomery<int32_t>	Field3;	typedef StaticElement< Field3 > Element3;	template<> Field3 Element3::_domain(2);
    typedef Modular<Integer>		Field4;	typedef StaticElement< Field4 > Element4;	template<> Field4 Element4::_domain(2);
    typedef Modular<int32_t>		Field5;	typedef StaticElement< Field5 > Element5; 	template<> Field5 Element5::_domain(2);
    typedef Modular<int16_t>		Field6;	typedef StaticElement< Field6 > Element6; 	template<> Field6 Element6::_domain(2);
    typedef Modular<Log16>		Field7;	typedef StaticElement< Field7 > Element7; 	template<> Field7 Element7::_domain(2);
#ifdef GIVARO_USE_SIXTYFOUR
    typedef Modular<int64_t>		Field8;	typedef StaticElement< Field8 > Element8;
    template<>
    Field8 Element8::_domain(2);
#endif
}



template<class Field, class Element>
void TestField()
{
    uint64_t P = 251;

    // Initialization of static member
    Element::setDomain( Field( typename Field::Residu_t(P)) );

    // Initialisations of elements
    Element a(2),b(-29.8),c(33),d(Integer("123456789012345678901234567890")),e(0);



    e += (a = b);  	std::cout << a << " = " << b << " mod " << P << ";" << std::endl;
    e += (a = b + c);  	std::cout << a << " = " << b << " + " << c << " mod " << P << ";" << std::endl;
    e += (a = b - c);  	std::cout << a << " = " << b << " - " << c << " mod " << P << ";" << std::endl;
    e += (a = b * c);  	std::cout << a << " = " << b << " * " << c << " mod " << P << ";" << std::endl;
    e += (a = b / c);  	std::cout << a << " = " << b << " / " << c << " mod " << P << ";" << std::endl;

    std::cout << d << " + " << a << " mod " << P << " = "; e += (d += a); std::cout << d << ";" << std::endl;
    std::cout << d << " - " << a << " mod " << P << " = "; e += (d -= a); std::cout << d << ";" << std::endl;
    std::cout << d << " * " << a << " mod " << P << " = "; e += (d *= a); std::cout << d << ";" << std::endl;
    std::cout << d << " / " << a << " mod " << P << " = "; e += (d /= a); std::cout << d << ";" << std::endl;

    std::cout << a << " is non zero ? " << (a != Element(0) ) << std::endl;
    a = 0; std::cout << a << " is zero ? " << (a == Element(0) ) << std::endl;

    double dd(0.0); dd += (double)(e); dd += (float)e; dd += (uint32_t)e;
    Element::getDomain().write( std::cerr << "Test: " << dd << " within ") << std::endl;
}



int main(int argc, char ** argv) {

    TestField<Field1, Element1>();
    TestField<Field2, Element2>();
    TestField<Field3, Element3>();
    TestField<Field4, Element4>();
    TestField<Field5, Element5>();
    TestField<Field6, Element6>();
    TestField<Field7, Element7>();
#ifdef GIVARO_USE_SIXTYFOUR
    TestField<Field8, Element8>();
#endif

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
