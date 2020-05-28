// ========================================================== //
// Copyright(c)'2020 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <28 May 20 15:26:49 Jean-Guillaume.Dumas@imag.fr>
// ========================================================== //

/*! @file examples/Polynomial/AES.C
 * @ingroup examples
 * @ingroup Polynomial
 * @example examples/Polynomial/AES.C
 * @brief NO DOC
 */

#include <givaro/gf2.h>
#include <givaro/gfq.h>
#include <givaro/givpoly1factor.h>
#include <givaro/givquotientdomain.h>
#include <givaro/givtimer.h>
#include <sstream>

using namespace Givaro;

int AESbox[256] =  
 {0x63 ,0x7c ,0x77 ,0x7b ,0xf2 ,0x6b ,0x6f ,0xc5 ,0x30 ,0x01 ,0x67 ,0x2b ,0xfe ,0xd7 ,0xab ,0x76
 ,0xca ,0x82 ,0xc9 ,0x7d ,0xfa ,0x59 ,0x47 ,0xf0 ,0xad ,0xd4 ,0xa2 ,0xaf ,0x9c ,0xa4 ,0x72 ,0xc0
 ,0xb7 ,0xfd ,0x93 ,0x26 ,0x36 ,0x3f ,0xf7 ,0xcc ,0x34 ,0xa5 ,0xe5 ,0xf1 ,0x71 ,0xd8 ,0x31 ,0x15
 ,0x04 ,0xc7 ,0x23 ,0xc3 ,0x18 ,0x96 ,0x05 ,0x9a ,0x07 ,0x12 ,0x80 ,0xe2 ,0xeb ,0x27 ,0xb2 ,0x75
 ,0x09 ,0x83 ,0x2c ,0x1a ,0x1b ,0x6e ,0x5a ,0xa0 ,0x52 ,0x3b ,0xd6 ,0xb3 ,0x29 ,0xe3 ,0x2f ,0x84
 ,0x53 ,0xd1 ,0x00 ,0xed ,0x20 ,0xfc ,0xb1 ,0x5b ,0x6a ,0xcb ,0xbe ,0x39 ,0x4a ,0x4c ,0x58 ,0xcf
 ,0xd0 ,0xef ,0xaa ,0xfb ,0x43 ,0x4d ,0x33 ,0x85 ,0x45 ,0xf9 ,0x02 ,0x7f ,0x50 ,0x3c ,0x9f ,0xa8
 ,0x51 ,0xa3 ,0x40 ,0x8f ,0x92 ,0x9d ,0x38 ,0xf5 ,0xbc ,0xb6 ,0xda ,0x21 ,0x10 ,0xff ,0xf3 ,0xd2
 ,0xcd ,0x0c ,0x13 ,0xec ,0x5f ,0x97 ,0x44 ,0x17 ,0xc4 ,0xa7 ,0x7e ,0x3d ,0x64 ,0x5d ,0x19 ,0x73
 ,0x60 ,0x81 ,0x4f ,0xdc ,0x22 ,0x2a ,0x90 ,0x88 ,0x46 ,0xee ,0xb8 ,0x14 ,0xde ,0x5e ,0x0b ,0xdb
 ,0xe0 ,0x32 ,0x3a ,0x0a ,0x49 ,0x06 ,0x24 ,0x5c ,0xc2 ,0xd3 ,0xac ,0x62 ,0x91 ,0x95 ,0xe4 ,0x79
 ,0xe7 ,0xc8 ,0x37 ,0x6d ,0x8d ,0xd5 ,0x4e ,0xa9 ,0x6c ,0x56 ,0xf4 ,0xea ,0x65 ,0x7a ,0xae ,0x08
 ,0xba ,0x78 ,0x25 ,0x2e ,0x1c ,0xa6 ,0xb4 ,0xc6 ,0xe8 ,0xdd ,0x74 ,0x1f ,0x4b ,0xbd ,0x8b ,0x8a
 ,0x70 ,0x3e ,0xb5 ,0x66 ,0x48 ,0x03 ,0xf6 ,0x0e ,0x61 ,0x35 ,0x57 ,0xb9 ,0x86 ,0xc1 ,0x1d ,0x9e
 ,0xe1 ,0xf8 ,0x98 ,0x11 ,0x69 ,0xd9 ,0x8e ,0x94 ,0x9b ,0x1e ,0x87 ,0xe9 ,0xce ,0x55 ,0x28 ,0xdf
 ,0x8c ,0xa1 ,0x89 ,0x0d ,0xbf ,0xe6 ,0x42 ,0x68 ,0x41 ,0x99 ,0x2d ,0x0f ,0xb0 ,0x54 ,0xbb ,0x16};

typedef GFqDom<int64_t> Field;
typedef Poly1Dom< GF2, Dense>::Element Polynomial;
typedef Field::Element Byte;

std::string dec2hex(const size_t decimal_value) {
    std::stringstream ss;
    ss<< '[' << std::hex << decimal_value << ']'; // int decimal_value
    return std::string( ss.str() );
}

std::string byte2hex(const Field& F, const Byte octet) {
    int64_t bval; F.convert(bval, octet);
    GIVARO_ASSERT( bval == F.zech2padic(octet), "field inconsistency");
    return dec2hex(bval);
}


int main(int argc, char** argv)
{
    GF2 F2;
    Poly1Dom< GF2, Dense> PD2(F2,'X');
    Poly1PadicDom< GF2, Dense > P2adic(PD2);

    
        // P8 = X^8 + X^4 + X^3 + X + 1
    Field GF256(2, 8, std::vector<int64_t>{1,1,0,1,1,0,0,0,1});
 
    Polynomial Q; P2adic.radix(Q, GF256.irreducible() );
    std::cout << "GF256 with irreducible (2-adic): " << GF256.irreducible() << '=' << dec2hex(GF256.irreducible()) << std::endl;
    PD2.write(std::cout << "  representing: ", Q) << std::endl;
    
        // Generator is 1+X
    Field::Element gen;
    GF256.generator(gen);
    GF256.write(std::cout <<
                "In this field, the generator used is (in 2-adic): ", gen)
                          << std::endl
                          << "  whose internal representation is "
                          << gen << std::endl;

        // A random byte
    GivRandom randiter;
    Byte octet; GF256.random(randiter, octet); 
    int64_t bval; GF256.convert(bval, octet);
    
    GF256.write(std::cout << "Byte " << dec2hex(bval) << '=' << bval << " is ", octet) << " represented as " << octet << '=' << byte2hex(GF256, octet) << std::endl;
    GF256.write(std::cout << "Checking that AES S-box of ", octet) << '=' << byte2hex(GF256,octet) << " is " << AESbox[bval] << '=' << dec2hex(AESbox[bval]) << std::endl;
    
        // Field inverse
    Byte inverse; 
    if (GF256.isZero(octet))
        GF256.assign(inverse, GF256.zero);
    else
        GF256.inv(inverse, octet);
    
    std::cout << "  " << byte2hex(GF256,inverse) << " is the inverse of " << byte2hex(GF256,octet) << std::endl;
    Polynomial binv; P2adic.radix(binv, GF256.zech2padic(inverse) );
    PD2.write(std::cout << "  and represents; ", binv) << std::endl;
    

        // Affine function
    Polynomial matrix, c3; 
    P2adic.radix(matrix, 1+(1<<1)+(1<<2)+(1<<3)+(1<<4));	// cyclic(1F)=1+X+X^2+X^3+X^4+X^5
    P2adic.radix(c3, 1+(1<<1)+(1<<5)+(1<<6)); 			// C3=1+X+X^5+X^6
    PD2.write(std::cout << "  Linear multiplicator: ", matrix)  << std::endl;
    PD2.write(std::cout << "  Affine constant: ", c3)  << std::endl;
    
    Polynomial deg8; P2adic.radix(deg8, 1+(1<<8));	// 1+X^8
    QuotientDom< Poly1Dom< GF2, Dense> > Q2D8(PD2, deg8);
    
    Polynomial tmp;
    Q2D8.mul( tmp, matrix, binv);
    PD2.write(std::cout << "  multiplied: ", tmp) << std::endl;
    Q2D8.addin( tmp, c3);
    PD2.write(std::cout << "  added: ", tmp) << std::endl;
    uint64_t res; P2adic.eval(res, tmp);
    std::cout << "  represented by: " << res << '=' << dec2hex(res) << std::endl;

//     Field F2(2);
//     Poly1Dom< Field, Dense> Pol2(F2);

   
//     GF128.init(b, 5U);
//     GF128.init(c, 3);
//     GF128.write(std::cout, b) << std::endl;
//     GF128.write(std::cout, c) << std::endl;

//     Field::Element f,g,h,j;
//     Field F2(2);
//     Poly1Dom< Field, Dense> Pol2(F2);
//     Poly1Dom< Field, Dense>::Element P, Q, R;
//     Pol2.init(P,Degree(1));
//     F2.init(P[0],1);
//     F2.init(P[1],1);
//     GF128.init(f, P);
//     GF128.write(std::cout << "2-adic representation of 1+X is: ", f)
//     << std::endl
//     << " ... while its internal representation is: "
//     << f << std::endl;

//     GF128.write(std::cout << "Indeed, we are in ") <<std::endl;

//     GF128.generator(gen);
//     GF128.write(std::cout <<
//                 "In this field, the generator used is (in 2-adic): ", gen)
//     << std::endl
//     << " whose internal representation is "
//     << gen << std::endl;


//     Poly1PadicDom< Field, Dense > Padic2(Pol2);
//     //

//     std::cout << "Irreducible (in 2-adic): "
//     << GF128.irreducible() << std::endl;

//     GF128.init(g, Padic2.radix( Q, Integer(5) ));
//     GF128.write(std::cout
//                 << "2-adic representation of 1+X^2 is: ", g)
//     << std::endl;

//     GF128.init(h);
//     GF128.add(h, g, f);
//     GF128.write(std::cout
//                 << "2-adic representation of X+X^2 is: ", h)
//     << std::endl;

//     GF128.mul(h, g, f);
//     GF128.write(std::cout
//                 << "2-adic representation of 1+X+X^2+X^3 is: ", h)
//     << std::endl;

//     GF128.div(h, g, f);
//     GF128.write(std::cout
//                 << "2-adic representation of 1+X is: ", h)
//     << std::endl;

//     GF128.init(j, Padic2.radix( Q, Integer(213) ));
//     GF128.write(std::cout
//                 << "2-adic representation of the moding out of X^7+X^6+X^4+X^2+1 by the irreducible is: ", j)
//     << std::endl;

    return 0;

}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
