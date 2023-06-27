// Copyright(c)'1994-2025 by The Givaro group
// This file is part of Givaro.
// BB: adapted,enhanced from examples/FiniteFields/ff-arith.C
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

#include <iostream>
#include <iomanip>

#include "test-fieldarith.h"
#include <givaro/modular.h>
#include <givaro/modular-balanced.h>
#include <givaro/modular-extended.h>
#include <givaro/montgomery.h>

#include <givaro/gf2.h>
#include <givaro/gfq.h>
#include <givaro/gfqext.h>
#include <givaro/extension.h>
#include <givaro/givintprime.h>

#include <recint/recint.h>

using namespace Givaro;


int main(int argc, char ** argv)
{
    int seed = int (argc>1?atoi(argv[1]):BaseTimer::seed());
#ifdef __GIVARO_DEBUG
    std::cerr << "seed: " << seed << std::endl;
#endif
    std::cerr<<std::setprecision(17);
    Integer::seeding((uint64_t)seed);
    RecInt::srand(seed);

    using ModularCUS = Modular<int8_t, uint16_t>;
    using ModularSUZ = Modular<int16_t, uint32_t>;
    using ModularZULL = Modular<int32_t, uint64_t>;
    using ModularUCUS = Modular<uint8_t, uint16_t>;
    using ModularUSUZ = Modular<uint16_t, uint32_t>;
    using ModularUZULL = Modular<uint32_t, uint64_t>;
    // using ModularFD = Modular<float, double>;

    //--------------------//
    //----- Modulo 2 -----//

    TEST_SPECIFIC(Modular<int8_t>, C2, 2);
    TEST_SPECIFIC(Modular<int16_t>, S2, 2);
    TEST_SPECIFIC(Modular<int32_t>, Z2, 2);
    TEST_SPECIFIC(Modular<int64_t>, LL2, 2);
    TEST_SPECIFIC(Modular<uint8_t>, UC2, 2);
    TEST_SPECIFIC(Modular<uint16_t>, US2, 2);
    TEST_SPECIFIC(Modular<uint32_t>, UZ2, 2);
    TEST_SPECIFIC(Modular<uint64_t>, ULL2, 2);
    TEST_SPECIFIC(ModularCUS, CUS2, 2);
    TEST_SPECIFIC(ModularSUZ, SUZ2, 2);
    TEST_SPECIFIC(ModularZULL, ZULL2, 2);
    TEST_SPECIFIC(ModularUCUS, UCUS2, 2);
    TEST_SPECIFIC(ModularUSUZ, USUZ2, 2);
    TEST_SPECIFIC(ModularUZULL, UZULL2, 2);
    TEST_SPECIFIC(Modular<Log16>, L2, 2);
    TEST_SPECIFIC(Modular<float>, F2, 2);
    TEST_SPECIFIC(Modular<double>, D2, 2);
    TEST_SPECIFIC(ModularExtended<float>, MEF2, 2);
    TEST_SPECIFIC(ModularExtended<double>, MED2, 2);
    //TEST_SPECIFIC(ModularFD, FD2, 2);
    TEST_SPECIFIC(Modular<Integer>, I2, 2);
    TEST_SPECIFIC(Modular<RecInt::rint128>, R2, 2);
    TEST_SPECIFIC(Modular<RecInt::ruint128>, RU2, 2);
    TEST_SPECIFIC(GF2, GF_2, 2);

    //--------------------//
    //----- Modulo 3 -----//

    TEST_SPECIFIC(ModularBalanced<int32_t>, BZ3, 3);
    TEST_SPECIFIC(ModularBalanced<int64_t>, BLL3, 3);
    TEST_SPECIFIC(ModularBalanced<float>, BF3, 3);
    TEST_SPECIFIC(ModularBalanced<double>, BD3, 3);

    TEST_SPECIFIC(Montgomery<int32_t>, MZ3, 3);
    TEST_SPECIFIC(Montgomery<RecInt::ruint128>, MRU3, 3);

    //---------------------//
    //----- Modulo 13 -----//

    TEST_SPECIFIC(Modular<int8_t>, C13, 13);
    TEST_SPECIFIC(Modular<int16_t>, S13, 13);
    TEST_SPECIFIC(Modular<int32_t>, Z13, 13);
    TEST_SPECIFIC(Modular<int64_t>, LL13, 13);
    TEST_SPECIFIC(Modular<uint8_t>, UC13, 13);
    TEST_SPECIFIC(Modular<uint16_t>, US13, 13);
    TEST_SPECIFIC(Modular<uint32_t>, UZ13, 13);
    TEST_SPECIFIC(Modular<uint64_t>, ULL13, 13);
    TEST_SPECIFIC(ModularCUS, CUS13, 13);
    TEST_SPECIFIC(ModularSUZ, SUZ13, 13);
    TEST_SPECIFIC(ModularZULL, ZULL13, 13);
    TEST_SPECIFIC(ModularUCUS, UCUS13, 13);
    TEST_SPECIFIC(ModularUSUZ, USUZ13, 13);
    TEST_SPECIFIC(ModularUZULL, UZULL13, 13);
    TEST_SPECIFIC(Modular<Log16>, L13, 13);
    TEST_SPECIFIC(Modular<float>, F13, 13);
    TEST_SPECIFIC(Modular<double>, D13, 13);
    TEST_SPECIFIC(ModularExtended<float>, MEF13, 13);
    TEST_SPECIFIC(ModularExtended<double>, MED13, 13);
    //TEST_SPECIFIC(ModularFD, FD13, 13);
    TEST_SPECIFIC(Modular<Integer>, I13, 13);
    TEST_SPECIFIC(Modular<RecInt::rint128>, R13, 13);
    TEST_SPECIFIC(Modular<RecInt::ruint128>, RU13, 13);
    using MU128U256 = Modular<RecInt::ruint128,RecInt::ruint256> ;
    using MU64U128 = Modular<RecInt::ruint64,RecInt::ruint128>;
    using M128256 = Modular<RecInt::rint128,RecInt::rint256> ;
    using M64128 = Modular<RecInt::rint64,RecInt::rint128>;
    TEST_SPECIFIC(MU64U128, R67U13, 13);
    TEST_SPECIFIC(M64128, R67R13, 13);
    TEST_SPECIFIC(MU128U256, RUU13, 13);
    TEST_SPECIFIC(M128256, R78R13, 13);

    TEST_SPECIFIC(ModularBalanced<int32_t>, BZ13, 13);
    TEST_SPECIFIC(ModularBalanced<int64_t>, BLL13, 13);
    TEST_SPECIFIC(ModularBalanced<float>, BF13, 13);
    TEST_SPECIFIC(ModularBalanced<double>, BD13, 13);

    TEST_SPECIFIC(Montgomery<int32_t>, MZ13, 13);
    TEST_SPECIFIC(Montgomery<RecInt::ruint128>, MRU13, 13);

    //--------------------------------//
    //----- Modulo maximal prime -----//

    TEST_LAST_PRIME(Modular<int8_t>, Cpmax);
    TEST_LAST_PRIME(Modular<int16_t>, Spmax);
    TEST_LAST_PRIME(Modular<int32_t>, Zpmax);
    TEST_LAST_PRIME(Modular<int64_t>, LLpmax);
    TEST_LAST_PRIME(Modular<uint8_t>, UCpmax);
    TEST_LAST_PRIME(Modular<uint16_t>, USpmax);
    TEST_LAST_PRIME(Modular<uint32_t>, UZpmax);
    TEST_LAST_PRIME(Modular<uint64_t>, ULLpmax);
    TEST_LAST_PRIME(ModularCUS, CUSpmax);
    TEST_LAST_PRIME(ModularSUZ, SUZpmax);
    TEST_LAST_PRIME(ModularZULL, ZULLpmax);
    TEST_LAST_PRIME(ModularUCUS, UCUSpmax);
    TEST_LAST_PRIME(ModularUSUZ, USUZpmax);
    TEST_LAST_PRIME(ModularUZULL, UZULLpmax);
    TEST_LAST_PRIME(Modular<Log16>, Lpmax);
    TEST_LAST_PRIME(Modular<float>, Fpmax);
    TEST_LAST_PRIME(Modular<double>, Dpmax);
    TEST_LAST_PRIME(ModularExtended<float>, MEFpmax);
    TEST_LAST_PRIME(ModularExtended<double>, MEDpmax);
    //TEST_LAST_PRIME(ModularFD, FDpmax);
    TEST_LAST_PRIME(Modular<RecInt::rint64>, R6pmax);
    TEST_LAST_PRIME(Modular<RecInt::ruint64>, R6Upmax);
    TEST_LAST_PRIME(Modular<RecInt::rint128>, R7pmax);
    TEST_LAST_PRIME(Modular<RecInt::ruint128>, R7Upmax);
    TEST_LAST_PRIME(Modular<RecInt::rint256>, R8pmax);
    TEST_LAST_PRIME(Modular<RecInt::ruint256>, R8Upmax);
    TEST_LAST_PRIME(MU64U128, RU67pmax);
    TEST_LAST_PRIME(MU128U256, RU78pmax);
    TEST_LAST_PRIME(M64128, R67pmax);
    TEST_LAST_PRIME(M128256, R78pmax);

    TEST_LAST_PRIME(ModularBalanced<int32_t>, BZpmax);
    TEST_LAST_PRIME(ModularBalanced<int64_t>, BLLpmax);
    TEST_LAST_PRIME(ModularBalanced<float>, BFpmax);
    TEST_LAST_PRIME(ModularBalanced<double>, BDpmax);

    TEST_LAST_PRIME(Montgomery<int32_t>, MZpmax);
    TEST_LAST_PRIME(Montgomery<RecInt::ruint128>, MRUpmax);


    //-------------------------//
    //----- Galois fields -----//

    TEST_SPECIFIC(GFqDom<int32_t>, GF13, 13);
    TEST_SPECIFIC(GFqDom<int32_t>, GFpmax, 65521U);
    TEST_SPECIFIC(GFqDom<int64_t>, GFLLplarge, 4194301);

    TEST_LAST_PRIME(GFqDom<int32_t>, GFpmmax);
        // int64_t maxCardinality would require 96GB
        // and about 960 CPU seconds ...

    // CP: disabled as it makes the testsuite fails on machine with <1.5GB.
    //     see https://github.com/linbox-team/givaro/issues/194
    // TEST_SPECIFIC(GFqDom<int64_t>, GFLLpXXL, 67108859); // already 1.5 GB

    // Zech log finite field with 256 elements
    // and prescribed 1 + x +x^3 +x^4 +x^8 irreducible polynomial
    //     std::vector< GFqDom<int64_t>::Residu_t > Irred(9);
    //     Irred[0] = 1; Irred[1] = 1; Irred[2] = 0; Irred[3] = 1;
    //     Irred[4] = 1; Irred[5] = 0; Irred[6] = 0; Irred[7] = 0;
    //     Irred[8] = 1;
    std::vector< int64_t > Irred {1,1,0,1,1,0,0,0,1};
    TEST_SPECIFIC(GFqDom<int64_t>, GF256, 2, 8, Irred);

    // Zech log finite field with 343 elements
    // and prescribed 3 +x^3  irreducible polynomial
    // and prescribed 5 +3x +4x^2 generator polynomial
    std::vector< int64_t > Irred3 {3,0,0,1}, Gene3 {5,3,4};
    TEST_SPECIFIC(GFqDom<int64_t>, GF343, 7, 3, Irred3, Gene3);

    TEST_SPECIFIC(GFqDom<int32_t>, GF625, 5, 4);
    TEST_SPECIFIC(GFqExt<int32_t>, GF81, 3, 4);

    // Zech log finite field with 2Mb tables
    TEST_SPECIFIC(GFqDom<int64_t>, GF2M, 2, 20);
    TEST_SPECIFIC(GFqDom<int64_t>, GF2M1, 2, 2);
    TEST_SPECIFIC(GFqDom<int64_t>, GF11E3, 11, 3);
    TEST_SPECIFIC(Extension<GFqDom<int64_t> >, GF11E9, GF11E3, 3);
    TEST_SPECIFIC(Extension<>, GF13E8, 13, 8);

    return 0;
}


/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
