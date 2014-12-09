// Copyright(c)'1994-2009,2010 by The Givaro group
// This file is part of Givaro.
// BB: adapted,enhanced from examples/FiniteFields/ff-arith.C
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

#include <iostream>
#include <givaro/modular.h>
#include <givaro/montgomery.h>
#include <givaro/givgfq.h>
#include <givaro/givintprime.h>
#include <givaro/givgfqext.h>
#include <givaro/givextension.h>

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
    double invertible=(FIRSTFLOAT<0?-FIRSTFLOAT:FIRSTFLOAT);
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

	// modulo 13 over 16 bits
	Modular<int16_t> C13(13);
	JETESTE(C13,seed);

	// modulo 13 over 32 bits
	Modular<int32_t> Z13(13);
	JETESTE(Z13,seed);

	// modulo 13 over unsigned 32 bits
	Modular<uint32_t> U13(13);
	JETESTE(U13,seed);

#ifdef __USE_Givaro_SIXTYFOUR__
	// modulo 13 over 64 bits
	Modular<int64_t> LL13(13UL);
	JETESTE(LL13,seed);
#endif

	// modulo 13 fully tabulated
	Modular<Log16> L13(13);
	JETESTE(L13,seed);

	// modulo 13 over 32 bits with Montgomery reduction
	Montgomery<int32_t> M13(13);
	JETESTE(M13,seed);

// Maximal prime values

	// prime modulo max over 15 bits
	Modular<int16_t> CUpmax( previousprime(Modular<int16_t>::getMaxModulus()) ); // 16381
	JETESTE(CUpmax,seed);

	// previous prime modulo max fully tabulated
	Modular<Log16> Lpmax( previousprime(Modular<Log16>::getMaxModulus() ) ); // 16369
	JETESTE(Lpmax,seed);

	// prime modulo max over 31 bits
	Modular<int32_t> Zpmax( previousprime(Modular<int32_t>::getMaxModulus() ) ); // 46337
	JETESTE(Zpmax,seed);

	// prime modulo max over 32 bits
	Modular<uint32_t> Upmax(previousprime(Modular<uint32_t>::getMaxModulus() ) ); // 65521
	JETESTE(Upmax,seed);

	// prime modulo max over 32 bits with Montgomery reduction
	Montgomery<int32_t> Mpmax(previousprime(Montgomery<int32_t>::getMaxModulus() ) ); // 40499
	JETESTE(Mpmax,seed);

#ifdef __USE_Givaro_SIXTYFOUR__
	// prime modulo max over 63 bits
	Modular<int64_t> LLpmax(previousprime(Modular<int64_t>::getMaxModulus() ) ); // 3037000493ULL
	JETESTE(LLpmax,seed);
#endif

	// Zech log prime field with prime max
	GFqDom<int> GFpmax( 65521UL );
	JETESTE(GFpmax,seed);

	// Zech log prime field with prime max (memory limited)
	GFqDom<long long> GFLLpmax( 4194301ULL );
	JETESTE(GFLLpmax,seed);


// Maximal values
	// prime modulo max over 15 bits
	Modular<int16_t> CUmax(Modular<int16_t>::getMaxModulus() );
	JETESTE(CUmax,seed);

	// prime modulo max fully tabulated
	Modular<Log16> Lmax( Modular<Log16>::getMaxModulus()  );
	JETESTE(Lmax,seed);

	// prime modulo max over 31 bits
	Modular<int32_t> Zmax(Modular<int32_t>::getMaxModulus());
	JETESTE(Zmax,seed);

	// modulo max over 32 bits
	Modular<uint32_t> Umax(Modular<uint32_t>::getMaxModulus() );
	JETESTE(Umax,seed);
	// prime modulo max over 32 bits with Montgomery reduction
	Montgomery<int32_t> Mmax(Montgomery<int32_t>::getMaxModulus() );
	JETESTE(Mmax,seed);

#ifdef __USE_Givaro_SIXTYFOUR__
	// prime modulo max over 63 bits
	Modular<int64_t> LLmax(Modular<int64_t>::getMaxModulus());
	JETESTE(LLmax,seed);
#endif



// Characteristic 2


	// modulo 2 over 16 bits
	Modular<int16_t> C2(2);
	JETESTE(C2,seed);

	// modulo 2 over 32 bits
	Modular<int32_t> Z2(2);
	JETESTE(Z2,seed);

	// modulo 2 over unsigned 32 bits
	Modular<uint32_t> U2(2);
	JETESTE(U2,seed);

#ifdef __USE_Givaro_SIXTYFOUR__
	// modulo 2 over 64 bits
	Modular<int64_t> LL2(2UL);
	JETESTE(LL2,seed);
#endif

	// modulo 2 fully tabulated
	Modular<Log16> L2(2);
	JETESTE(L2,seed);

	Modular<Log16> L2b( L2 );
	JETESTE(L2b,seed);

// Other Characteristics


	// modulo 3 over 32 bits with Montgomery reduction
	Montgomery<int32_t> M2(3);
	JETESTE(M2,seed);

	Montgomery<int32_t> M3(39989);
	JETESTE(M3,seed);

	// modulo 13 with primitive root representation
	GFqDom<int> GF13( 13 );
	JETESTE(GF13,seed);

	// modulo 13 over arbitrary size
	Modular<Integer> IntZ13(13);
	JETESTE(IntZ13,seed);

	// modulo 13 with generic implementation over signed integral type
	Modular<long long> GenZ13(13);
	JETESTE(GenZ13,seed);

	// modulo 101 with generic implementation over unsigned signed integral type
	Modular<unsigned long long> GenZ101(101);
	JETESTE(GenZ101,seed);

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
#else
	GFqDom<long> GF2M( 2, 20) ;
#endif
	JETESTE(GF2M,seed);

#ifndef __GIVARO__DONOTUSE_longlong__
	GFqDom<long long> GF2M1( 2, 2 );
#else
	GFqDom<long> GF2M1( 2, 2) ;
#endif
	JETESTE(GF2M1,seed);


        Extension<> GF13E8(13,8);
	JETESTE(GF13E8,seed);

#ifndef __GIVARO__DONOTUSE_longlong__
        GFqDom<long long> GF11e3( 11, 3 );
	JETESTE(GF11e3,seed);
        Extension<GFqDom<long long> > GF11e9(GF11e3,3);
	JETESTE(GF11e9,seed);
#else
        GFqDom<long> GF11e3( 11, 3) ;
	JETESTE(GF11e3,seed);
        Extension<> GF11e9(GF11e3,3);
	JETESTE(GF11e9,seed);
#endif

#ifdef GIVARO_DEBUG
	std::cerr << std::endl ;
#endif


	return 0;
}

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
