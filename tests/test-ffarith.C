// Copyright(c)'1994-2009,2010 by The Givaro group
// This file is part of Givaro.
// BB: adapted,enhanced from examples/FiniteFields/ff-arith.C
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

#include <iostream>

#include <givaro/modular.h>
#include <givaro/modular-balanced.h>
#include <givaro/unparametric.h>
#include <givaro/montgomery.h>

#include <givaro/gfq.h>
#include <givaro/gfqext.h>
#include <givaro/extension.h>
#include <givaro/givintprime.h>

#include <recint/recint.h>

using namespace Givaro;

#define TESTE_EG( a, b ) \
if (!F.areEqual((a),(b))) {\
	F.write(F.write(std::cout,a) << "!=",b) << " failed (at line " <<  __LINE__ << ")" << std::endl; \
	return(-1); \
}

#define JETESTE( a, s ) \
if (TestField( (a), int(s)) ) {\
	std::cout << #a << " failed !" << std::endl;\
	return -1 ; \
}

#define JEONETESTE( F, a, x ) \
if (TestOneField(F,(int)a,(float)x)) {\
	std::cout << #a << " failed !" << std::endl;\
	return -1 ; \
}

template<class Int1, class Int2>
long long locgcd ( const Int1 a, const Int2 b ) {
//     std::cerr << std::endl << '|' << a << '^' << b << '|';
    long long u3, v3; u3 = a; v3 = (long long)b;
    while (v3 != 0) {
        long long q, t3;
        q = u3 / v3;
        t3 = u3 - q * v3;
        u3 = v3; v3 = t3;
    }
//     std::cerr << u3 << std::endl;

    return u3;
}




template<class Field>
int TestOneField(const Field& F, const int FIRSTINT, const float FIRSTFLOAT)
{
#ifdef GIVARO_DEBUG
	std::cerr << "testing " ;
	F.write(std::cerr );
        std::cerr << " (" << FIRSTINT << ',' << FIRSTFLOAT << ')';
	std::cerr  << " : " << std::flush;

#endif



	typename Field::Element a, b, c, d,a_,b_,c_,d_,ma;
	typename Field::Element e,e_;

        F.init(a, 0UL);
        TESTE_EG(a, F.zero);
        F.init(a, 1UL);
//         F.write(std::cerr) << std::endl;
//         F.write(std::cerr << "a: ", a) << std::endl;
//         F.write(std::cerr << "1: ", F.one) << std::endl;
        TESTE_EG(a, F.one);

        F.init(ma,-1L);
//         F.write(std::cerr) << std::endl;
//         F.write(std::cerr << "a: ", a) << std::endl;
//         F.write(std::cerr << "ma: ", ma) << std::endl;
//         F.write(std::cerr << "1: ", F.one) << std::endl;
//         F.write(std::cerr << "-1: ", F.mOne) << std::endl;
        TESTE_EG(ma, F.mOne);

	F.init(a, FIRSTINT);
    double invertible=floor(FIRSTFLOAT<0?-FIRSTFLOAT:FIRSTFLOAT);
    do {
        invertible += 1;
        F.init(b, invertible);
    } while ( locgcd( (long long)(invertible),F.characteristic()) != 1 );



	F.init(c);            // empty constructor
	F.init(d);            // empty constructor

	F.add(c, a, b);       // c = a+b
	F.init(c_);           //! @warning F.init(c_,c); ne marche pas !
	F.assign(c_,c);       // c_ <- c

//         F.write(std::cerr) << std::endl;
//         F.write(std::cerr << "a:=", a) << ';' << std::endl;
//         F.write(std::cerr << "b:=", b) << ';' << std::endl;
//         F.write(std::cerr << "c:=", c) << ';' << std::endl;
//         F.write(std::cerr << "c_:=", c_) << ';' << std::endl;

	TESTE_EG(c,c_);
	F.subin(c_,a);
//         F.write(std::cerr) << std::endl;
//         F.write(std::cerr << "a:=", a) << ';' << std::endl;
//         F.write(std::cerr << "b:=", b) << ';' << std::endl;
//         F.write(std::cerr << "c:=", c) << ';' << std::endl;
//         F.write(std::cerr << "c_:=", c_) << ';' << std::endl;

	TESTE_EG(b,c_);

	F.mul(c, a, b);     // c = a*b
	F.assign(c_,c);       // c_ <- c
	F.divin(c_,b);      // c_ == a ?

//         F.write(std::cerr) << std::endl;
//         F.write(std::cerr << "a: ", a) << std::endl;
//         F.write(std::cerr << "b: ", b) << std::endl;
//         F.write(std::cerr << "c: ", c) << std::endl;
//         F.write(std::cerr << "c_: ", c_) << std::endl;
	TESTE_EG(a,c_);

	F.axpy(d, a, b, c); // d = a*b + c;
	F.init(d_);
	F.axmy(d_,a,b,c); // d_ = a*b - c
	F.addin(d_,c);
	F.subin(d,c);

//         F.write(std::cerr) << std::endl;
//         F.write(std::cerr << "a:=", a) << ';' << std::endl;
//         F.write(std::cerr << "b:=", b) << ';' << std::endl;
//         F.write(std::cerr << "c:=", c) << ';' << std::endl;
//         F.write(std::cerr << "d:=", d) << ';' << std::endl;
//         F.write(std::cerr << "d_:=", d_) << ';' << std::endl;
	TESTE_EG(d_,d);

	F.sub(d,a,b); // d = a -b
	F.add(c,a,b); // c = a+b
	F.init(e);
	F.init(e_);
	F.mul(e,d,c); // e = d*c;
	F.mul(a_,a,a); // a_ = a*a
	F.mul(b_,b,b); // b_ = b*b
	F.sub(e_,a_,b_); // e_ = a_ - b_

//         F.write(std::cerr) << std::endl;
//         F.write(std::cerr << "a:=", a) << ';' << std::endl;
//         F.write(std::cerr << "b:=", b) << ';' << std::endl;
//         F.write(std::cerr << "c:=", c) << ';' << std::endl;
//         F.write(std::cerr << "d:=", d) << ';' << std::endl;
//         F.write(std::cerr << "e:=", e) << ';' << std::endl;
//         F.write(std::cerr << "e_:=", e_) << ';' << std::endl;
//         F.write(std::cerr << "a_:=", a_) << ';' << std::endl;
//         F.write(std::cerr << "b_:=", b_) << ';' << std::endl;
	TESTE_EG(e,e_) // a^2 - b^2 = (a-b)(a+b) ;)

	// Four operations
	F.init(a_);
	F.assign(a_,a);
	F.addin(a, b) ;
	F.subin(a, b) ;
	F.mulin(a, b) ;
	F.divin(a, b) ;

	TESTE_EG(a_,a);


	F.maxpy(e, a, b, d); // e = d-a*b
	F.assign(e_,d);
	F.maxpyin(e_, a, b); // e = d - a*b;

//         F.write(std::cerr << "a:=", a) << ';' << std::endl;
//         F.write(std::cerr << "b:=", b) << ';' << std::endl;
//         F.write(std::cerr << "d:=", d) << ';' << std::endl;
//         F.write(std::cerr << "e:=", e) << ';' << std::endl;
//         F.write(std::cerr << "e_:=", e_) << ';' << std::endl;
	TESTE_EG(e,e_);


	F.axmy(e, a, b, d); // e = a*b -d;

	F.assign(e_,d);
	F.maxpyin(e_, a, b); // e = d - a*b;

	F.negin(e_);

	TESTE_EG(e,e_);


	if ((unsigned) F.characteristic() != (unsigned)2) {
		F.init(e,1);
		F.init(a,22993);
		F.inv(b,a);
		F.mul(c,b,a);

		TESTE_EG(e,c);

		F.init(a,22993);
		F.init(b,22993);
		F.invin(a);
		F.mulin(a,b);
	}

	TESTE_EG(e,a);

	F.init(a,37397);
	F.inv(b,a);
	F.mul(c,b,a);

	TESTE_EG(e,c);

	F.init(a,37397);
	F.init(b,37397);
	F.invin(a);
	F.mulin(a,b);

	TESTE_EG(e,a);

#ifdef GIVARO_DEBUG
	F.write(std::cerr );
	std::cerr  << " done." << std::endl;
#endif
	return 0 ;

}

#ifndef NBITER
#define NBITER 50
#endif

template<class Field>
int TestField(const Field& F, const int seed)
{
    long ch = (long) F.characteristic();
    JEONETESTE(F,7UL,-29.3);
    srand48(seed);
    for(size_t i=0; i< NBITER; ++i) {
        typename Field::Element x;
        float d;
	do {
		d = float((double)ch*drand48());
            F.init(x,d );
        } while(F.isZero(x));
        int a; do {
            F.init(x, a = (int)lrand48());
        } while(F.isZero(x));
        JEONETESTE(F,a,d);
    }
    return 0;
}

template<class Ints>
Ints previousprime(const Ints& a) {
    static IntPrimeDom IPD;
    Integer aI(a);
    IPD.prevprimein(aI);

	return (Ints)aI;
}


int main(int argc, char ** argv)
{
    int seed = int (argc>1?atoi(argv[1]):BaseTimer::seed());
#ifdef GIVARO_DEBUG
    std::cerr << "seed: " << seed << std::endl;
#endif
    Integer::seeding((unsigned long)seed);

#ifdef NDEBUG
    assert(0);
#endif

	//---------------------//
	//----- Modulo 13 -----//

	Modular<int8_t> S13(13);
	JETESTE(S13,seed);

	Modular<uint8_t> US13(13);
	JETESTE(US13,seed);

	Modular<int16_t> C13(13);
	JETESTE(C13,seed);

	Modular<uint16_t> UC13(13);
	JETESTE(UC13,seed);

	Modular<int32_t> Z13(13);
	JETESTE(Z13,seed);

	ModularBalanced<int32_t> BZ13(13);
	JETESTE(Z13,seed);

	Modular<uint32_t> U13(13);
	JETESTE(U13,seed);

	Modular<Log16> L13(13);
	JETESTE(L13,seed);

	Montgomery<int32_t> M13(13);
	JETESTE(M13,seed);

#ifdef __USE_Givaro_SIXTYFOUR__
	Modular<int64_t> LL13(13UL);
	JETESTE(LL13,seed);

	Modular<uint64_t> ULL13(13UL);
	JETESTE(ULL13,seed);

	ModularBalanced<int64_t> BULL13(13UL);
	JETESTE(BULL13,seed);
#endif

	Modular<Integer> I13(13);
	JETESTE(I13,seed);

	Modular<float> F13(13);
	JETESTE(F13,seed);

	ModularBalanced<float> BF13(13);
	JETESTE(BF13,seed);

	Modular<double> D13(13);
	JETESTE(D13,seed);

	ModularBalanced<double> BD13(13);
	JETESTE(BD13,seed);

	Modular<RecInt::ruint128> R13(13);
	JETESTE(R13,seed);

	//--------------------------------//
	//----- Modulo maximal prime -----//

	Modular<int8_t> Spmax(previousprime(Modular<int8_t>::getMaxModulus() ) );
	JETESTE(Spmax,seed);

	Modular<uint8_t> USpmax(previousprime(Modular<uint8_t>::getMaxModulus() ) );
	JETESTE(USpmax,seed);

	Modular<int16_t> CUpmax( previousprime(Modular<int16_t>::getMaxModulus()) );
	JETESTE(CUpmax,seed);

	Modular<uint16_t> UCUpmax( previousprime(Modular<int16_t>::getMaxModulus()) );
	JETESTE(UCUpmax,seed);

	Modular<Log16> Lpmax( previousprime(Modular<Log16>::getMaxModulus() ) );
	JETESTE(Lpmax,seed);

	Modular<int32_t> Zpmax( previousprime(Modular<int32_t>::getMaxModulus() ) );
	JETESTE(Zpmax,seed);

	ModularBalanced<int32_t> BZpmax( previousprime(ModularBalanced<int32_t>::getMaxModulus() ) );
	JETESTE(BZpmax,seed);

	Modular<uint32_t> Upmax(previousprime(Modular<uint32_t>::getMaxModulus() ) );
	JETESTE(Upmax,seed);

	Montgomery<int32_t> Mpmax(previousprime(Montgomery<int32_t>::getMaxModulus() ) );
	JETESTE(Mpmax,seed);

#ifdef __USE_Givaro_SIXTYFOUR__
	Modular<int64_t> LLpmax(previousprime(Modular<int64_t>::getMaxModulus() ) );
	JETESTE(LLpmax,seed);

	ModularBalanced<int64_t> BLLpmax(previousprime(ModularBalanced<int64_t>::getMaxModulus() ) );
	JETESTE(BLLpmax,seed);

	Modular<uint64_t> ULLpmax(previousprime(Modular<uint64_t>::getMaxModulus() ) );
	JETESTE(ULLpmax,seed);
#endif

	Modular<float> Fpmax(previousprime(Modular<float>::getMaxModulus() ) );
	JETESTE(Fpmax,seed);

	ModularBalanced<float> BFpmax(previousprime(ModularBalanced<float>::getMaxModulus() ) );
	JETESTE(BFpmax,seed);

	Modular<double> Dpmax(previousprime(Modular<double>::getMaxModulus() ) );
	JETESTE(Dpmax,seed);

	ModularBalanced<double> BDpmax(previousprime(ModularBalanced<double>::getMaxModulus() ) );
	JETESTE(BDpmax,seed);

	Modular<RecInt::ruint128> Rpmax(previousprime(Modular<RecInt::ruint128>::getMaxModulus()));
	JETESTE(Rpmax,seed);

	//--------------------------//
	//----- Modulo maximal -----//

	Modular<int8_t> Smax(Modular<int8_t>::getMaxModulus() );
	JETESTE(Smax,seed);

	Modular<uint8_t> USmax(Modular<uint8_t>::getMaxModulus() );
	JETESTE(USmax,seed);

	Modular<int16_t> CUmax(Modular<int16_t>::getMaxModulus() );
	JETESTE(CUmax,seed);

	Modular<uint16_t> UCUmax(Modular<uint16_t>::getMaxModulus() );
	JETESTE(UCUmax,seed);

	Modular<Log16> Lmax( Modular<Log16>::getMaxModulus()  );
	JETESTE(Lmax,seed);

	Modular<int32_t> Zmax(Modular<int32_t>::getMaxModulus());
	JETESTE(Zmax,seed);

	ModularBalanced<int32_t> BZmax(ModularBalanced<int32_t>::getMaxModulus());
	JETESTE(BZmax,seed);

	Modular<uint32_t> Umax(Modular<uint32_t>::getMaxModulus() );
	JETESTE(Umax,seed);

	Montgomery<int32_t> Mmax(Montgomery<int32_t>::getMaxModulus() );
	JETESTE(Mmax,seed);

#ifdef __USE_Givaro_SIXTYFOUR__
	Modular<int64_t> LLmax(Modular<int64_t>::getMaxModulus());
	JETESTE(LLmax,seed);

	ModularBalanced<int64_t> BLLmax(ModularBalanced<int64_t>::getMaxModulus());
	JETESTE(BLLmax,seed);

	Modular<uint64_t> ULLmax(Modular<uint64_t>::getMaxModulus());
	JETESTE(ULLmax,seed);
#endif

	Modular<float> Fmax(Modular<float>::getMaxModulus());
	JETESTE(Fmax,seed);

	ModularBalanced<float> BFmax(ModularBalanced<float>::getMaxModulus());
	JETESTE(BFmax,seed);

	Modular<double> Dmax(Modular<double>::getMaxModulus());
	JETESTE(Dmax,seed);

	ModularBalanced<double> BDmax(ModularBalanced<double>::getMaxModulus());
	JETESTE(BDmax,seed);

    // Not a prime
	// Modular<RecInt::ruint128> Rmax(Modular<RecInt::ruint128>::getMaxModulus());
	// JETESTE(Rmax,seed);

	//--------------------//
	//----- Modulo 2 -----//

	Modular<int8_t> S2(2);
	JETESTE(S2,seed);

	Modular<uint8_t> US2(2);
	JETESTE(US2,seed);

	Modular<int16_t> C2(2);
	JETESTE(C2,seed);

	Modular<uint16_t> UC2(2);
	JETESTE(UC2,seed);

	Modular<int32_t> Z2(2);
	JETESTE(Z2,seed);

	Modular<uint32_t> U2(2);
	JETESTE(U2,seed);

	Modular<Log16> L2(2);
	JETESTE(L2,seed);

	Modular<Log16> L2b( L2 );
	JETESTE(L2b,seed);

#ifdef __USE_Givaro_SIXTYFOUR__
	Modular<int64_t> LL2(2UL);
	JETESTE(LL2,seed);

	Modular<uint64_t> ULL2(2UL);
	JETESTE(ULL2,seed);
#endif

	Modular<Integer> I2(2);
	JETESTE(I2,seed);

	Modular<float> F2(2);
	JETESTE(F2,seed);

	Modular<double> D2(2);
	JETESTE(D2,seed);

	Modular<RecInt::ruint128> R2(2);
	JETESTE(R2,seed);

	//--------------------//
	//----- Modulo 3 -----//

	Montgomery<int32_t> M3(3);
	JETESTE(M3,seed);

	ModularBalanced<int32_t> BZ3(3);
	JETESTE(BZ3,seed);

#ifdef __USE_Givaro_SIXTYFOUR__
	ModularBalanced<int64_t> BLL3(3UL);
	JETESTE(BLL3,seed);
#endif

	ModularBalanced<float> BF3(3);
	JETESTE(BF3,seed);

	ModularBalanced<double> BD3(3);
	JETESTE(BD3,seed);

	Modular<RecInt::ruint128> R3(3);
	JETESTE(R3,seed);

	//---------------------------------//
	//----- Other Characteristics -----//

	Montgomery<int32_t> MR(39989);
	JETESTE(MR,seed);

	GFqDom<int> GF13( 13 );
	JETESTE(GF13,seed);

	// Zech log prime field with prime max
	GFqDom<int> GFpmax( 65521UL );
	JETESTE(GFpmax,seed);

#ifndef __GIVARO__DONOTUSE_longlong__
	// Zech log prime field with prime max (memory limited)
	GFqDom<long long> GFLLpmax( 4194301ULL );
	JETESTE(GFLLpmax,seed);
#endif

#ifdef __USE_Givaro_SIXTYFOUR__
	Modular<uint64_t> GenZ101(101);
	JETESTE(GenZ101,seed);
#endif

    // Zech log finite field with 256 elements
    // and prescribed 1 + x +x^3 +x^4 +x^8 irreducible polynomial
    std::vector< GFqDom<long>::Residu_t > Irred(9);
    Irred[0] = 1; Irred[1] = 1; Irred[2] = 0; Irred[3] = 1;
    Irred[4] = 1; Irred[5] = 0; Irred[6] = 0; Irred[7] = 0;
    Irred[8] = 1;
    GFqDom<long> GF256(2,8, Irred);
    JETESTE(GF256,seed);

	// Zech log finite field with 5^4 elements
	GFqDom<int> GF625( 5, 4 );
	JETESTE(GF625,seed);

	// Zech log finite field with 3^4 elements
	// Using the Q-adic Transform
	GFqExt<int> GF81( 3, 4 );
	JETESTE(GF81,seed);

	// Zech log finite field with 2Mb tables
#ifndef __GIVARO__DONOTUSE_longlong__
	GFqDom<long long> GF2M( 2, 20 );
	GFqDom<long long> GF2M1( 2, 2 );
    GFqDom<long long> GF11e3( 11, 3 );
    Extension<GFqDom<long long> > GF11e9(GF11e3,3);
#else
	GFqDom<long> GF2M( 2, 20) ;
	GFqDom<long> GF2M1( 2, 2) ;
    GFqDom<long> GF11e3( 11, 3) ;
    Extension<> GF11e9(GF11e3,3);
#endif

	JETESTE(GF2M,seed);
	JETESTE(GF2M1,seed);
	JETESTE(GF11e3,seed);
	JETESTE(GF11e9,seed);

    Extension<> GF13E8(13,8);
	JETESTE(GF13E8,seed);

#ifdef GIVARO_DEBUG
	std::cerr << std::endl ;
#endif

	return 0;
}

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
