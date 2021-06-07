// =================================================================== //
// Copyright(c)'2021 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <07 Jun 21 11:59:11 Jean-Guillaume.Dumas@imag.fr>
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

namespace Givaro {
        // Returns r=Fibo(n-1) & t=Fibo(n)
        // a, b are temporary (values at half n);
        // Assumes n>=1
    Integer& Fibonacci(Integer& r, Integer& t, Integer& a, Integer& b,
                   const uint64_t n) {
        if (n <= 1) {
            r=Integer::zero;
            return t=Integer::one;
        }
        Fibonacci(a,b,r,t,n>>1); // n/2
            // (a+bX)^2 mod X^2-X-1 = (a^2+b^2) + (2ab+b^2)X
        if (n & 1u) {
            r = a<<1;
            r += b;
            r *= b;		// r=b(b+2a)
            t = b<<1;
            t = a - t;
            t *= a;		// t=a(a-2b)
            t += r;		// t=a^2+b^2
            t += r;		// next term
        } else {
            t = a<<1;
            t += b;
            t *= b;		// t=b(b+2a)
            r = b<<1;
            r = a - r;
            r *= a;		// r=a(a-2b)
            r += t;		// r=a^2+b^2
        }
        return t;
    }

        // Assumes n>=0
    Integer& Fibonacci(Integer& t, const uint64_t n) {
        if (n) {
            Integer a,c,d;
            return Fibonacci(a,t,c,d,n);
        } else {
            return t=Integer::zero;
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
