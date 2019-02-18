#include <iostream>
#include <cstdlib>

#include <recint/ruint.h>
#include <recint/rmintmg.h>

using namespace RecInt;

int main(void)
{
    ruint<7> p, q;
    ruint<8> d, e, n, g, phin;
    rmint<8> m, c, f;

    // p and q are two primes
    fill_with_1(p); p -= 796;
    fill_with_1(q); q -= 712;

    // n = p*q;
    lmul(n, p, q);
    m.init_module(n);

    std::cout << "p = " << p << std::endl;
    std::cout << "q = " << q << std::endl;
    std::cout << "Module n = p*q = " << n << std::endl;

    //  phin = (p-1)*(q-1);
    lmul(phin, p-1, q-1);

    RecInt::srand(time(NULL));

    // Looking for relatively prime e with phin so that inv(e) exists
    do {
        rand(e);
        gcd(g, e, phin);
    } while (g != 1);

    // d = inv(e) mod phin
    inv_mod(d, e, phin);

    std::cout << std::endl << "Encryption key e = " << e << std::endl;
    std::cout << "Decryption key d = " << d << std::endl;

    rand(m);
    std::cout << std::endl << "Message m = " << m << std::endl;

    exp(c, m, e);
    std::cout << "Encryption: c = m^e mod n --> " << c << std::endl;

    exp(f, c, d);
    std::cout << "Decryption: m = c^d mod n --> " << f << std::endl;

    if (m == f) std::cout << std::endl << "Decryption OK" << std::endl;
    else std::cout << std::endl << "Decryption failed" << std::endl;

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
