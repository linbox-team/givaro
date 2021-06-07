// =================================================================== //
// Copyright(c)'2021 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <07 Jun 21 18:57:12 Jean-Guillaume.Dumas@imag.fr>
//
// Fibonacci numbers    /////////////////////////
// =================================================================== //

/*! @file examples/Integer/Fibonacci.C
 * @ingroup examples
 * @ingroup integers
 * @example examples/Integer/Fibonacci.C
 * @brief NO DOC
 */
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <givaro/givinteger.h>

#define GMPFIBTABLESIZE 93
#define GMPCNSTUINT(C) ((uint64_t) C##u)

const uint64_t givfibtable[GMPFIBTABLESIZE+2] {
  GMPCNSTUINT (0x1),  /* -1 */
  GMPCNSTUINT (0x0),  /* 0 */
  GMPCNSTUINT (0x1),  /* 1 */
  GMPCNSTUINT (0x1),  /* 2 */
  GMPCNSTUINT (0x2),  /* 3 */
  GMPCNSTUINT (0x3),  /* 4 */
  GMPCNSTUINT (0x5),  /* 5 */
  GMPCNSTUINT (0x8),  /* 6 */
  GMPCNSTUINT (0xd),  /* 7 */
  GMPCNSTUINT (0x15),  /* 8 */
  GMPCNSTUINT (0x22),  /* 9 */
  GMPCNSTUINT (0x37),  /* 10 */
  GMPCNSTUINT (0x59),  /* 11 */
  GMPCNSTUINT (0x90),  /* 12 */
  GMPCNSTUINT (0xe9),  /* 13 */
  GMPCNSTUINT (0x179),  /* 14 */
  GMPCNSTUINT (0x262),  /* 15 */
  GMPCNSTUINT (0x3db),  /* 16 */
  GMPCNSTUINT (0x63d),  /* 17 */
  GMPCNSTUINT (0xa18),  /* 18 */
  GMPCNSTUINT (0x1055),  /* 19 */
  GMPCNSTUINT (0x1a6d),  /* 20 */
  GMPCNSTUINT (0x2ac2),  /* 21 */
  GMPCNSTUINT (0x452f),  /* 22 */
  GMPCNSTUINT (0x6ff1),  /* 23 */
  GMPCNSTUINT (0xb520),  /* 24 */
  GMPCNSTUINT (0x12511),  /* 25 */
  GMPCNSTUINT (0x1da31),  /* 26 */
  GMPCNSTUINT (0x2ff42),  /* 27 */
  GMPCNSTUINT (0x4d973),  /* 28 */
  GMPCNSTUINT (0x7d8b5),  /* 29 */
  GMPCNSTUINT (0xcb228),  /* 30 */
  GMPCNSTUINT (0x148add),  /* 31 */
  GMPCNSTUINT (0x213d05),  /* 32 */
  GMPCNSTUINT (0x35c7e2),  /* 33 */
  GMPCNSTUINT (0x5704e7),  /* 34 */
  GMPCNSTUINT (0x8cccc9),  /* 35 */
  GMPCNSTUINT (0xe3d1b0),  /* 36 */
  GMPCNSTUINT (0x1709e79),  /* 37 */
  GMPCNSTUINT (0x2547029),  /* 38 */
  GMPCNSTUINT (0x3c50ea2),  /* 39 */
  GMPCNSTUINT (0x6197ecb),  /* 40 */
  GMPCNSTUINT (0x9de8d6d),  /* 41 */
  GMPCNSTUINT (0xff80c38),  /* 42 */
  GMPCNSTUINT (0x19d699a5),  /* 43 */
  GMPCNSTUINT (0x29cea5dd),  /* 44 */
  GMPCNSTUINT (0x43a53f82),  /* 45 */
  GMPCNSTUINT (0x6d73e55f),  /* 46 */
  GMPCNSTUINT (0xb11924e1),  /* 47 */
  GMPCNSTUINT (0x11e8d0a40),  /* 48 */
  GMPCNSTUINT (0x1cfa62f21),  /* 49 */
  GMPCNSTUINT (0x2ee333961),  /* 50 */
  GMPCNSTUINT (0x4bdd96882),  /* 51 */
  GMPCNSTUINT (0x7ac0ca1e3),  /* 52 */
  GMPCNSTUINT (0xc69e60a65),  /* 53 */
  GMPCNSTUINT (0x1415f2ac48),  /* 54 */
  GMPCNSTUINT (0x207fd8b6ad),  /* 55 */
  GMPCNSTUINT (0x3495cb62f5),  /* 56 */
  GMPCNSTUINT (0x5515a419a2),  /* 57 */
  GMPCNSTUINT (0x89ab6f7c97),  /* 58 */
  GMPCNSTUINT (0xdec1139639),  /* 59 */
  GMPCNSTUINT (0x1686c8312d0),  /* 60 */
  GMPCNSTUINT (0x2472d96a909),  /* 61 */
  GMPCNSTUINT (0x3af9a19bbd9),  /* 62 */
  GMPCNSTUINT (0x5f6c7b064e2),  /* 63 */
  GMPCNSTUINT (0x9a661ca20bb),  /* 64 */
  GMPCNSTUINT (0xf9d297a859d),  /* 65 */
  GMPCNSTUINT (0x19438b44a658),  /* 66 */
  GMPCNSTUINT (0x28e0b4bf2bf5),  /* 67 */
  GMPCNSTUINT (0x42244003d24d),  /* 68 */
  GMPCNSTUINT (0x6b04f4c2fe42),  /* 69 */
  GMPCNSTUINT (0xad2934c6d08f),  /* 70 */
  GMPCNSTUINT (0x1182e2989ced1),  /* 71 */
  GMPCNSTUINT (0x1c5575e509f60),  /* 72 */
  GMPCNSTUINT (0x2dd8587da6e31),  /* 73 */
  GMPCNSTUINT (0x4a2dce62b0d91),  /* 74 */
  GMPCNSTUINT (0x780626e057bc2),  /* 75 */
  GMPCNSTUINT (0xc233f54308953),  /* 76 */
  GMPCNSTUINT (0x13a3a1c2360515),  /* 77 */
  GMPCNSTUINT (0x1fc6e116668e68),  /* 78 */
  GMPCNSTUINT (0x336a82d89c937d),  /* 79 */
  GMPCNSTUINT (0x533163ef0321e5),  /* 80 */
  GMPCNSTUINT (0x869be6c79fb562),  /* 81 */
  GMPCNSTUINT (0xd9cd4ab6a2d747),  /* 82 */
  GMPCNSTUINT (0x16069317e428ca9),  /* 83 */
  GMPCNSTUINT (0x23a367c34e563f0),  /* 84 */
  GMPCNSTUINT (0x39a9fadb327f099),  /* 85 */
  GMPCNSTUINT (0x5d4d629e80d5489),  /* 86 */
  GMPCNSTUINT (0x96f75d79b354522),  /* 87 */
  GMPCNSTUINT (0xf444c01834299ab),  /* 88 */
  GMPCNSTUINT (0x18b3c1d91e77decd),  /* 89 */
  GMPCNSTUINT (0x27f80ddaa1ba7878),  /* 90 */
  GMPCNSTUINT (0x40abcfb3c0325745),  /* 91 */
  GMPCNSTUINT (0x68a3dd8e61eccfbd),  /* 92 */
  GMPCNSTUINT (0xa94fad42221f2702),  /* 93 */
};


namespace Givaro {

        // Returns r=Fibo(n-1) & t=Fibo(n)
        // b is temporary
        // Assumes n>=0
    Integer& Fibonacci(Integer& r, Integer& t,
                       const uint64_t n) {
        if (n <= 93) {
            uint64_t k(n);
            r=givfibtable[k];
            return t=givfibtable[++k];
        }

        uint64_t k(n); ++k; k >>=1 ;// n=2k-->k or n=2k+1--> k+1

        Fibonacci(t,r,k);		// n=2k-->k or n=2k+1--> k+1
        if (n & 1u) {
            Integer b(r);
            r <<= 1;
            r -= t;
            r *= t;					// r=t(2t-r)
            t <<= 1;
            t += b;
            t *= b;					// t=t(t-2t)
            t -= r;					// t=t^2+r^2
       } else {
            Integer b(t);
            t <<= 1;
            t += r;
            t *= r;					// t=r(r+2t)
            r <<= 1;
            Integer::sub(r,b,r);	// r = t - r;
            r *= b;					// r=t(t-2r)
            r += t;					// r=t^2+r^2
        }
        return t;
    }

        // Assumes n>=0
        // Avoids one of the last two (largest) multiplications
    Integer& Fibonacci(Integer& t, const uint64_t n) {
        Integer a;
        if (n & 1u) {
            uint64_t k(n); ++k; k >>=1 ;// n=2k-->k or n=2k+1--> k+1
               // F[2k+1]=F[k]*(F[k+1]+2F[k]) +/- 1;
            Fibonacci(t,a,k);
            Integer b(t);
            t <<= 1;
            t += a;
            t *= b;
            return ( (k & 1u)? ++t : --t );
        } else {
                // F[2k]=F[k]*(F[k]+2F[k-1]);
            Fibonacci(t,a, n>>1 );
            t <<= 1;
            t += a;
            return t *= a;
        }

    }
}

int main (int argc, char * * argv) {
    uint64_t n = argc > 1 ? (uint64_t)atoi(argv[1]) : 42u;

    Givaro::Integer t;

    Givaro::Timer tim; tim.start();
    Givaro::Fibonacci(t,n);
    tim.stop();

    std::cout << t  << std::endl;
    std::clog << tim << std::endl;

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
