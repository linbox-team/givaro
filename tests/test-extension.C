// Copyright(c) by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

// Part of this is adapted from issues #203, #204 reported by dragomang87
#include <iostream>
#include <givaro/givpower.h>
#include <givaro/modular.h>
#include <givaro/gf2.h>
#include <givaro/gfq.h>
#include <givaro/extension.h>

using namespace Givaro;

#ifndef GIV_PASSED_MSG
#  define GIV_PASSED_MSG "\033[1;32mPASSED.\033[0m"
#endif
#ifndef GIV_ERROR_MSG
#  define GIV_ERROR_MSG "\033[1;31m****** ERROR ******\033[0m "
#endif

template<typename Field>
bool TestIdentity(GivRandom& generator, const int MOD, const int expo) {
    const uint64_t seed(generator.seed());
    Field basefield(MOD);
    Extension< Field > field(basefield, expo);
    typename Extension< Field >::Element a, b;
    field.random(generator, a);
    field.random(generator, b);
    field.write(std::clog << "Field is: \n  ") << std::endl;
    field.write(std::clog << "Random elements: \n    a = ", a);
    field.write(std::clog << "\n    b = ", b) << std::endl;

    typename Extension< Field >::Element  power, product;
    dom_power(power,a,0,field);
    bool pass(field.areEqual(field.one, power));
    if (! pass) {
        std::cerr << GIV_ERROR_MSG << seed << std::endl;
        field.write(std::cerr << "    field.one is:", field.one) << std::endl;
        field.write(std::cerr << "    a^0       is:", power) << std::endl;
        std::cerr << "Is power==field.one? "
                  << pass
                  << std::endl;
    } else {
        std::clog << "[DomPow] : " << GIV_PASSED_MSG << std::endl;
    }

    field.mul(product,b,field.one);
    bool success(field.areEqual(product,b)); pass &= success;
    if (! success) {
        std::cerr << GIV_ERROR_MSG << seed << std::endl;
        field.write(std::cerr
                    << "product = field.mul(product,b,field.one) = "
                    ,  product)
                    << std::endl
                    << "Is product==b? "
                    << success
                    << std::endl;
    } else {
        std::clog << "[MulOne] : " << GIV_PASSED_MSG << std::endl;
    }

    field.mul(product,b,power);
    success = field.areEqual(product,b); pass &= success;
    if (! success) {
        std::cerr << GIV_ERROR_MSG << seed << std::endl;
        field.write(std::clog
                    << "product = field.mul(product,b,power    ) = "
                    ,  b)
                    << std::endl
                    << "Is product==b? "
                    << success
                    << std::endl;
    } else {
        std::clog << "[MulPow] : " << GIV_PASSED_MSG << std::endl;
    }

    return pass;
}

template<typename Field>
bool TestExponent(const int seed, const Field& F, const int expo) {
    using EFq = Extension< Field >;
    EFq Ef(F, expo, 'X');
    Extension< EFq > EEf( Ef, 5, 'Y');
    EEf.write(std::clog << "Extension is: \n  ") << std::endl;

    bool success( (int)EEf.exponent() == 5*expo*Exponent_Trait(F) );

    if (! success) {
        std::cerr << GIV_ERROR_MSG << seed << std::endl;
        std::cerr << "Extensions of order " << Exponent_Trait(F) << '*'
                  << expo << "*5 is "
                  << EEf.exponent() << std::endl;
    } else {
        std::clog << "[ExtExp] : " << GIV_PASSED_MSG << std::endl;
    }

    return success;
}

template<typename Field>
bool TestSignedness(GivRandom& generator, const int expo) {
    Field field(2,expo);
    typename Field::Element random;
    field.random(generator, random); // was segfaulting with unsigned
    field.write(std::clog << "[FrcSgn] : " << GIV_PASSED_MSG) << std::endl;
    return true;
}


int main(int argc, char ** argv)
{
    const int seed = int (argc>1?atoi(argv[1]):BaseTimer::seed());
#ifdef __GIVARO_DEBUG
    std::cerr << "seed: " << seed << std::endl;
#endif
    GivRandom generator(seed);

    const int MOD = int (argc>2?atoi(argv[2]):7);
    const int expo = int (argc>3?atoi(argv[3]):5);

    bool pass( TestExponent<GFq<>>(seed, GFq<>(MOD,3), expo) );
    pass &= TestExponent< GF2 >(seed, GF2(), expo<<1);
    pass &= TestExponent< GFq<> >(seed, GFq<>(2,7), expo<<1);

    pass &= TestIdentity<GFq<int64_t>>(generator, MOD, expo);
    pass &= TestIdentity<Modular<double>>(generator, MOD, expo);

    pass &=  TestSignedness<GFq<int64_t>>(generator,expo);
        // Uses a signed representation anyway
    pass &=  TestSignedness<GFq<uint64_t>>(generator,expo);

    return (! pass);
}
