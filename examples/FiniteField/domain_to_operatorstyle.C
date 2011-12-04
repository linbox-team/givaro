// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

/*! @file examples/FiniteField/domain_to_operatorstyle.C
 * @ingroup examples
 * @ingroup finitefields
 * @example examples/FiniteField/domain_to_operatorstyle.C
 * @brief NO DOC
 */

#include <iostream>
#include <givaro/givgfq.h>
#include <givaro/StaticElement.h>

using namespace Givaro;



// Finite Field with Domain style
typedef GFqDom<long> Field;

// Wrapper to give an operator style to the elements
typedef StaticElement< Field > Element;

// Mandatory declaration (because of static template)
// and an actual constructed field is mandatory (the "(2)") for g++ 3.4
namespace Givaro
{
template<>
Field Element::_domain(2);
}

int main(int argc, char ** argv) {
    unsigned long P = (argc>1 ? (unsigned long)atoi(argv[1]) : 5009UL);

        // Initialization of static member
    Element::setDomain( Field(P) );

        // Initialisations of elements
    Element a(2),b(-29.8),c(33),d(Integer("123456789012345678901234567890"));


        // Computations
    a = b;  	std::cerr << a << " = " << b << " mod " << P << ";" << std::endl;
    a = b + c;  std::cerr << a << " = " << b << " + " << c << " mod " << P << ";" << std::endl;
    a = b - c;  std::cerr << a << " = " << b << " - " << c << " mod " << P << ";" << std::endl;
    a = b * c;  std::cerr << a << " = " << b << " * " << c << " mod " << P << ";" << std::endl;
    a = b / c;  std::cerr << a << " = " << b << " / " << c << " mod " << P << ";" << std::endl;

        // Computations in place
    std::cerr << d << " + " << a << " mod " << P << " = ";
    d += a; std::cerr << d << ";" << std::endl;

    std::cerr << d << " - " << a << " mod " << P << " = ";
    d -= a; std::cerr << d << ";" << std::endl;

    std::cerr << d << " * " << a << " mod " << P << " = ";
    d *= a; std::cerr << d << ";" << std::endl;

    std::cerr << d << " / " << a << " mod " << P << " = ";
    d /= a; std::cerr << d << ";" << std::endl;

        // Tests
    std::cerr << a << " is non zero is " << (a != Element(0) ) << std::endl;

    a = 0; std::cerr << a << " is zero is " << (a == Element(0) ) << std::endl;


        // Access to Field object
    Field F = Element::getDomain();
    F.write( std::cerr << "Test: within ") << std::endl;

    return 0;
}
