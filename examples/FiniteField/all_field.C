#include <iostream>
#include <givaro/givzpz32uns.h>
#include <givaro/givgfq.h>
#include <givaro/givmontg32.h>
#include <givaro/givzpzInt.h>
#include "StaticElement.h"

// Domain kind
typedef ZpzDom<Unsigned32>	Field1;	typedef StaticElement< Field1 > Element1; 	Field1 Element1::_domain;
typedef GFqDom<long>		Field2;	typedef StaticElement< Field2 > Element2;	Field2 Element2::_domain;
typedef Montgomery<Std32>	Field3;	typedef StaticElement< Field3 > Element3;	Field3 Element3::_domain;
typedef ZpzDom<Integer>		Field4;	typedef StaticElement< Field4 > Element4;	Field4 Element4::_domain;
typedef ZpzDom<Std32>		Field5;	typedef StaticElement< Field5 > Element5; 	Field5 Element5::_domain;
typedef ZpzDom<Std16>		Field6;	typedef StaticElement< Field6 > Element6; 	Field6 Element6::_domain;
typedef ZpzDom<Log16>		Field7;	typedef StaticElement< Field7 > Element7; 	Field7 Element7::_domain;
#ifdef __USE_Givaro_64__
typedef ZpzDom<Std64>		Field8;	typedef StaticElement< Field8 > Element8; 	Field8 Element8::_domain;
#endif



template<class Field, class Element>
void TestField() {
    unsigned long P = 251;

        // Initialization of static member
    Element::setDomain( Field(P) );

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

    double dd(0.0); dd += (double)(e); dd += (float)e; dd += (unsigned int)e;    
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
#ifdef __USE_Givaro_64__
    TestField<Field8, Element8>();  
#endif

    return 0;
}
