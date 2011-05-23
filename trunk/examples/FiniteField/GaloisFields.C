// Copyright(c)'1994-2011 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

#include <iostream>
#include <givaro/givgfq.h>
#include <givaro/givextension.h>
#include <givaro/givzpz64std.h>

using namespace Givaro;

int main (int argc, char * * argv) {
    unsigned long q = (argc>1?atoi(argv[1]):13);
    unsigned long expo = (argc>2?atoi(argv[2]):8);
    
    {
            // This is the field with 11^2=121 elements
            // Using a generator representation and tables
        GFqDom<long> F1(11,2);
        F1.write(std::cout << "This is the field with 121 elements: ") << std::endl;
    }
    
    {
            // This is the field with q^expo elements using the best 
            // possible base field

        std::cout << "Exponent max for zech logs " << q << '^' << expo << " : " << FF_EXPONENT_MAX(q,expo) << std::endl;
        std::cout << "NEED polynomial representation : " << NEED_POLYNOMIAL_REPRESENTATION(q,expo) << std::endl;

        if ( NEED_POLYNOMIAL_REPRESENTATION(q,expo) ) {
                // The template parameter if the type of the base field
                // By default GFqDom<long> is used.
            Extension<> Fqe(q, expo);
            Fqe.write(std::cout << "This is the field with " << q << '^' << expo << " elements: ") << std::endl;
        } else {
            GFqDom<long> Fqe(q, expo);
            Fqe.write(std::cout << "This is the field with " << q << '^' << expo << " elements: ") << std::endl;
        }
    }
    

        // This is the field with 11^2=121 elements
        // Using a polynomial representation
    {
        GFqDom<long> F11(11);
        Extension< GFqDom<long> > F121(F11, 2);
        F121.write(std::cout << "This is the field with 121 elements: ") << std::endl;
    }

        // This is the field with 11^2=121 elements
        // Using a polynomial representation
        // And an alternative field, works only with Givaro >= 3.4.1
    {
        
        ZpzDom<Std64> F11(11);
        Extension< ZpzDom<Std64>  > F121(F11, 2);
        F121.write(std::cout << "This is the field with 121 elements: ") << std::endl;
    }


        // This is the field with 2^8 elements
    {
        GFqDom<long> F256(2,8);
        F256.write(std::cout << "This is the field with 256 elements: ") << ", using: " << F256.irreducible() << " as irreducible polynomial" << std::endl;
    }

        // This is the field with 2^8 elements
        // Using 1 + x +x^3 +x^4 +x^8 as irreducible polynomial
    {
        std::vector< GFqDom<long>::Residu_t > Irred(9);
        Irred[0] = 1; Irred[1] = 1; Irred[2] = 0; Irred[3] = 1;
        Irred[4] = 1; Irred[5] = 0; Irred[6] = 0; Irred[7] = 0;
        Irred[8] = 1;
        GFqDom<long> F256(2,8, Irred);
        F256.write(std::cout << "This is the field with 256 elements: ") << ", using: " << F256.irreducible() << " as irreducible polynomial" << std::endl;
    }


    return 0;
}
