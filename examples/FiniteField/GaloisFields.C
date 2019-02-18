// Copyright(c)'1994-2011 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

#include <iostream>
#include <givaro/gfq.h>
#include <givaro/extension.h>
#include <givaro/modular-integral.h>
#include <givaro/givprint.h>

using namespace Givaro;

int main (int argc, char * * argv) {
    uint64_t q = (argc>1?(uint64_t)atoi(argv[1]):13);
    uint64_t expo = (argc>2?(uint64_t)atoi(argv[2]):8);

    {
        // This is the field with 11^2=121 elements
        // Using a generator representation and tables
        GFqDom<int64_t> F1(11,2);
        F1.write(std::cout << "This is the field with 121 elements: ") << std::endl;
        std::cout << " using: " << F1.irreducible() << " as irreducible polynomial" << std::endl;
        GFqDom<int64_t>::Element primroot; F1.generator(primroot);
        F1.write(std::cout << " represented as indexes with respect to the generator: ", primroot) << std::endl;
        std::cout << " represented as indexes with respect to the generator: " << F1.generator() << std::endl;

    }

    {
        // This is the field with q^expo elements using the best
        // possible base field

        std::cerr << "Exponent max for zech logs with characteristic " << q << " : " << FF_EXPONENT_MAX(q,expo) << std::endl;
        std::cerr << "Sub-Exponent max for zech logs " << q << "^" << expo << " : " << FF_SUBEXPONENT_MAX(q,expo) << std::endl;
        std::cout << "NEED polynomial representation : " << NEED_POLYNOMIAL_REPRESENTATION(q,expo) << std::endl;

        if ( NEED_POLYNOMIAL_REPRESENTATION(q,expo) ) {
            // The template parameter if the type of the base field
            // By default GFqDom<int64_t> is used.
            Extension<> Fqe(q, expo);
            Fqe.write(std::cout << "This is the field with " << q << '^' << expo << " elements: ") << std::endl;
        } else {
            GFqDom<int64_t> Fqe(q, expo);
            Fqe.write(std::cout << "This is the field with " << q << '^' << expo << " elements: ") << std::endl;
        }
    }


    // This is the field with 11^2=121 elements
    // Using a polynomial representation
    {
        GFqDom<int64_t> F11(11);
        Extension< GFqDom<int64_t> > F121(F11, 2);
        F121.write(std::cout << "This is the field with 121 elements: ") << std::endl;
    }

    // This is the field with 11^2=121 elements
    // Using a polynomial representation
    // And an alternative field, works only with Givaro >= 3.4.1
    {

        Modular<int64_t> F11(11);
        Extension< Modular<int64_t> > F121(F11, 2);
        F121.write(std::cout << "This is the field with 121 elements: ") << ", using: " << F121.irreducible() << " as irreducible polynomial" << std::endl;


    }

    // This is the field with 11^2=121 elements
    // Using a polynomial representation
    // And an alternative field, works only with Givaro >= 3.4.1
    {

        Modular<int64_t> F11(11);
        Poly1Dom< Modular<int64_t>, Dense > PolF11(F11,"Z");
        Poly1Dom< Modular<int64_t>, Dense >::Element Irred;
        PolF11.init(Irred, Degree(2));
        F11.assign(Irred[0],F11.one);
        F11.assign(Irred[1],F11.one); // Irred is Y^2+Y+1
        Extension< Modular<int64_t> > F121(PolF11, Irred);
        F121.write(std::cout << "This is the field with 121 elements: ") << ", using: " << F121.irreducible() << " as irreducible polynomial" << std::endl;


    }


    // This is the field with 2^8 elements
    {
        GFqDom<int64_t> F256(2,8);
        F256.write(std::cout << "This is the field with 256 elements: ") << ", using: " << F256.irreducible() << " as irreducible polynomial" << std::endl;
    }

    // This is the field with 2^8 elements
    // Using 1 + x +x^3 +x^4 +x^8 as irreducible polynomial
    {
        std::vector< GFqDom<int64_t>::Residu_t > Irred(9);
        Irred[0] = 1; Irred[1] = 1; Irred[2] = 0; Irred[3] = 1;
        Irred[4] = 1; Irred[5] = 0; Irred[6] = 0; Irred[7] = 0;
        Irred[8] = 1;
        GFqDom<int64_t> F256(2,8, Irred);
        F256.write(std::cout << "This is the field with 256 elements: ") << ", using: " << F256.irreducible() << " as irreducible polynomial" << std::endl;

        std::cout << " this field is generated (in 2-adic) by: " << F256.generator() << std::endl;

        std::cout << " in this field the indeterminate is represented by: " << F256.sage_generator() << " (same as " << F256.indeterminate() << ')' << std::endl;


        Givaro::GFqDom<int64_t> F2(2);
        typedef Givaro::Poly1Dom< Givaro::GFqDom<int64_t>, Givaro::Dense> PolDomain;
        typedef Givaro::Poly1PadicDom< Givaro::GFqDom<int64_t>, Givaro::Dense> PadicDomain;


        PolDomain Pol2(F2,"Z");
        PadicDomain Padic2(Pol2);
        PolDomain::Element polGen, polIrred;
        Padic2.radix(polGen, F256.generator() );
        Padic2.radix(polIrred, F256.irreducible() );

        Pol2.write(Pol2.write(
                              std::cout << " that is with this representation, (", polGen)
                   << ")^" << F256.indeterminate()
                   << " == Z mod (", polIrred) << ')' << std::endl;

        {
            GFqDom<int64_t>::Element a, b, c;
            givvector<int64_t> Av(8); // A := X^7+X^6+X^3+X
            Av[0]=0; Av[1]=1; Av[2]=0; Av[3]=1;
            Av[4]=0; Av[5]=0; Av[6]=1; Av[7]=1;
            givvector<int64_t> Bv(8); // B:=X^6+X^4+X+1;
            Bv[0]=1; Bv[1]=1; Bv[2]=0; Bv[3]=0;
            Bv[4]=1; Bv[5]=0; Bv[6]=1; Bv[7]=0;

            F256.init(a, Av); F256.init(b, Bv);

            F256.mul(c,a,b);
            F256.write( std::cout, a) << '*';
            F256.write( std::cout, b) << '=';
            F256.write( std::cout, c) << std::endl;


            F256.add(c,a,b);
            F256.write( F256.write( F256.write(
                                               std::cout, a) << '+', b) << '=' , c ) << std::endl;
        }

        {
            GFqDom<int64_t>::Element a, b, c, x;
            x = F256.sage_generator();

            F256.add(a,x,F256.one); 	// X+1
            F256.mulin(a,x);
            F256.mulin(a,x);
            F256.mulin(a,x); 		// (X+1)X^3
            F256.addin(a,F256.one); 	// (X+1)X^3+1
            F256.mulin(a,x);
            F256.mulin(a,x); 		// ((X+1)X^3+1)X^2
            F256.addin(a,F256.one); 	// ((X+1)X^3+1)X^2+1
            F256.mulin(a,x); 		// (((X+1)X^3+1)X^2+1)X


            F256.mul(b,x,x);
            F256.addin(b,F256.one); 	// X^2+1
            F256.mulin(b,x);
            F256.mulin(b,x);
            F256.mulin(b,x); 		// (X^2+1)X^3
            F256.addin(b,F256.one); 	// (X^2+1)X^3+1
            F256.mulin(b,x);		// ((X^2+1)X^3+1)X
            F256.addin(b,F256.one); 	// ((X^2+1)X^3+1)X+1

            F256.mul(c,a,b);
            F256.write( std::cout, a) << '*';
            F256.write( std::cout, b) << '=';
            F256.write( std::cout, c) << std::endl;


            F256.add(c,a,b);
            F256.write( F256.write( F256.write(
                                               std::cout, a) << '+', b) << '=' , c ) << std::endl;
        }

        {
            Modular<int16_t> GF2(2);
            Poly1PadicDom< Modular<int16_t> > P2(GF2,"X");

            givvector<int64_t> vect202; P2.radixdirect(vect202, 202, 8);
            givvector<int64_t> vect83; P2.radixdirect(vect83, 83, 8);

            givvector< uint64_t > Irred2; P2.radixdirect(Irred2, 283, 9);
            GFqDom<int64_t> F256_(2,8, Irred2);


            GFqDom<int64_t>::Element a, b, c;
            F256_.init(a, vect202); F256_.init(b, vect83);

            F256_.mul(c,a,b);
            F256_.write( std::cout, a) << '*';
            F256_.write( std::cout, b) << '=';
            F256_.write( std::cout, c) << std::endl;


            F256_.add(c,a,b);
            F256_.write( F256_.write( F256_.write(
                                                  std::cout, a) << '+', b) << '=' , c ) << std::endl;

        }

    }


    return 0;
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
