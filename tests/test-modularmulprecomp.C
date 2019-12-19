// Copyright(c)'1994-2009,2010 by The Givaro group
// This file is part of Givaro.
// Copyright (C) 2016 Romain Lebreton, adapted from test-ringarith.C
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

#include <iostream>
#include <givaro/modular.h>
#include <givaro/givinteger.h>

using namespace Givaro;

#define TESTE_EG( a, b )										\
if (!F.areEqual((a),(b))) {									\
    F.write( F.write(std::cout,a) << "!=",b)					\
    << " failed (at line " <<  __LINE__ << ")" << std::endl;	\
    throw std::string( "Message d'erreur" );					\
    return -1;													\
}

#define NBITER 100

template<class Ring>
int TestOneMulPrecomp(const Ring& F, const typename Ring::Element& x, const typename Ring::Element& y)
{
#ifdef __GIVARO_DEBUG
    std::cerr << "Testing " ;
    F.write(std::cerr) << " : " << std::endl;
#endif

    typename Ring::Element a, b, c, d;
    typename Ring::Compute_t invp, invb, invb2;
    size_t bitsizep;

    F.assign(a, x);
    F.assign(b, y);
    if (a < 0){
        std::cout << "x : " << x << "\n a : " << a << std::endl;
    }

    F.mul(c,a,b);

    F.precomp_p(invp, bitsizep);
    F.mul_precomp_p(d,a,b,invp, bitsizep);
    TESTE_EG(c,d);

    F.precomp_b(invb, b);
    F.mul_precomp_b(d,a,b,invb);
    TESTE_EG(c,d);

    F.precomp_b(invb2, b, invp);
    TESTE_EG(invb,invb2);
    F.mul_precomp_b(d,a,b,invb);
    TESTE_EG(c,d);

    return 0;
}

#define JETESTE( F, x, y )						\
if (TestOneMulPrecomp(F,x,y)) {						\
    std::cout << #x << " " << #y << " failed !" << std::endl;	\
    return -1;							\
}


template<class Ring>
int TestMulPrecomp(const Ring& F, const uint64_t seed)
{

    typename Ring::Element x, y;
    typename Ring::RandIter g(F, 0_ui64, seed);

    F.init(x, 7U);
    F.init(y, -29.0);
    JETESTE(F,x,y);

    F.init(x, Ring::maxCardinality()-1U);
    F.init(y, Ring::maxCardinality()-1U);
    JETESTE(F,x,y);

    F.assign(x, F.maxElement());
    F.assign(y, F.maxElement());
    JETESTE(F,x,y);

    F.assign(x, F.minElement());
    F.assign(y, F.maxElement());
    JETESTE(F,x,y);

    F.assign(x, F.minElement());
    F.assign(y, F.minElement());
    JETESTE(F,x,y);

    for (size_t i = 0; i< NBITER; ++i) {
        g.random(x); g.random(y);
        JETESTE(F,x,y);
    }

    return 0;
}

#undef JETESTE

int main(int argc, char ** argv) {

    auto seed = static_cast<uint64_t>(argc>1?atoi(argv[1]):BaseTimer::seed());

#ifdef __GIVARO_DEBUG
    std::cerr << "seed: " << seed << std::endl;
#endif

    Integer::seeding(seed);
    RecInt::srand(seed);

    //--------------------------------//
    //----- Test mult with precomp ---//

    /**********************/
    /*  Modular<8, ...>   */
    /**********************/
    Modular<int8_t, int8_t> M88p (3); //2-bit prime
    TestMulPrecomp(M88p,seed);
    Modular<int8_t, uint8_t> M8u8p (3);
    TestMulPrecomp(M8u8p,seed);
    Modular<uint8_t, int8_t> Mu88p (3);
    TestMulPrecomp(Mu88p,seed);
    Modular<uint8_t, uint8_t> Mu8u8p (3);
    TestMulPrecomp(Mu8u8p,seed);

    Modular<int8_t, int16_t> M816p (61); //6-bit prime
    TestMulPrecomp(M816p,seed);
    Modular<int8_t, uint16_t> M8u16p (61);
    TestMulPrecomp(M8u16p,seed);
    Modular<uint8_t, int16_t> Mu816p (61);
    TestMulPrecomp(Mu816p,seed);
    Modular<uint8_t, uint16_t> Mu8u16p (61);
    TestMulPrecomp(Mu8u16p,seed);

    /**********************/
    /*  Modular<16, ...>  */
    /**********************/
    Modular<int16_t, int16_t> M1616p (61); //6-bit prime
    TestMulPrecomp(M1616p,seed);
    Modular<int16_t, uint16_t> M16u16p (61);
    TestMulPrecomp(M16u16p,seed);
    Modular<uint16_t, int16_t> Mu1616p (61);
    TestMulPrecomp(Mu1616p,seed);
    Modular<uint16_t, uint16_t> Mu16u16p (61);
    TestMulPrecomp(Mu16u16p,seed);

    Modular<int16_t, int32_t> M1632p (16381); //14-bit prime
    TestMulPrecomp(M1632p,seed);
    Modular<int16_t, uint32_t> M16u32p (16381);
    TestMulPrecomp(M16u32p,seed);
    Modular<uint16_t, int32_t> Mu1632p (16381);
    TestMulPrecomp(Mu1632p,seed);
    Modular<uint16_t, uint32_t> Mu16u32p (16381);
    TestMulPrecomp(Mu16u32p,seed);

    /**********************/
    /*  Modular<32, ...>  */
    /**********************/
    Modular<int32_t, int32_t> M3232p (16381); //14-bit prime
    TestMulPrecomp(M3232p,seed);
    Modular<int32_t, uint32_t> M32u32p (16381);
    TestMulPrecomp(M32u32p,seed);
    Modular<uint32_t, int32_t> Mu3232p (16381);
    TestMulPrecomp(Mu3232p,seed);
    Modular<uint32_t, uint32_t> Mu32u32p (16381);
    TestMulPrecomp(Mu32u32p,seed);


    Modular<int32_t, int64_t> M3264p (1073741789); //30-bit prime
    TestMulPrecomp(M3264p,seed);
    Modular<int32_t, uint64_t> M32u64p (1073741789);
    TestMulPrecomp(M32u64p,seed);
    Modular<uint32_t, int64_t> Mu3264p (1073741789);
    TestMulPrecomp(Mu3264p,seed);
    Modular<uint32_t, uint64_t> Mu32u64p (1073741789);
    TestMulPrecomp(Mu32u64p,seed);

    Modular<uint32_t, uint64_t> Mu32u64p2 (143513);
    TestMulPrecomp(Mu32u64p2,seed);


    /**********************/
    /*  Modular<64, ...>  */
    /**********************/
    Modular<int64_t, int64_t> M6464p (1073741789); //30-bit prime
    TestMulPrecomp(M6464p,seed);
    Modular<int64_t, uint64_t> M64u64p (1073741789);
    TestMulPrecomp(M64u64p,seed);
    Modular<uint64_t, int64_t> Mu6464p (1073741789);
    TestMulPrecomp(Mu6464p,seed);
    Modular<uint64_t, uint64_t> Mu64u64p (1073741789);
    TestMulPrecomp(Mu64u64p,seed);

#ifdef __GIVARO_HAVE_INT128
    Modular<int64_t, int128_t> M64128p (4611686018427387847ul); //62-bit prime
    TestMulPrecomp(M64128p,seed);
    Modular<int64_t, uint128_t> M64u128p (4611686018427387847ul);
    TestMulPrecomp(M64u128p,seed);
    Modular<uint64_t, int128_t> Mu64128p (4611686018427387847ul);
    TestMulPrecomp(Mu64128p,seed);
    Modular<uint64_t, uint128_t> Mu64u128p (4611686018427387847ul);
    TestMulPrecomp(Mu64u128p,seed);
#endif

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
