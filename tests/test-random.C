// Copyright(c)'2010 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// written by BB.

/** @file tests/test-random.C
 * @ingroup tests
 * @brief we test  bounds for random Integers
 * @test we test  bounds for random Integers
 */

#include <iostream>
#include <givaro/givinteger.h>

using namespace Givaro;
//#define GIVARO_DEBUG
using std::cout ; using std::endl;

//! tests <code>ret= .func(arg,arg);</code> ...
int test1()
{
	Integer toto  ;
	GivRandom gr;
	gr.seed();
#ifdef GIVARO_DEBUG
	cout << "this is a random() number : " << toto.random(gr) << endl;
#else
        toto.random(gr);
#endif

	gr.seed();
#ifdef GIVARO_DEBUG
	cout << "this is a random() number : " << toto.random(gr) << endl;
#else
        toto.random(gr);
#endif

	gr.seed();
#ifdef GIVARO_DEBUG
	cout << "this is a random() number : " << toto.random(gr) << endl;
#else
        toto.random(gr);
#endif

	Integer un(26);
	Integer autre(511);
	Integer large( Integer::one<<3000);
#ifdef GIVARO_DEBUG
	cout << "random...............OK" << endl;
#endif
	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer tata = toto.random_between(un,autre,gr);
		//cout << tata << endl;
		if (tata < un || tata >= autre) {
			cout << "random_between  failed" << endl;
			return -1  ;
		}
	}
#ifdef GIVARO_DEBUG
	cout << "random_between.......OK" << endl;
#endif

	unsigned long trois = 3 ;
	unsigned long petits = 6 ;
	unsigned long grands = 300 ;
	//std::vector<int> T(1<<petits) ;

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer tata = toto.random_between(trois,petits,gr);
		//cout << tata << endl;
		//T[tata] += 1 ;
		if (tata < (1<<trois) || tata >= (1<<petits) ) {
			cout << "random_between_exp  failed" << endl;
			return -1  ;
		}
	}
	//    for (size_t i = 0 ; i < 1<<petits ; ++i) cout << T[i] << " " ;
#ifdef GIVARO_DEBUG
	cout << "random_between_exp...OK" << endl;
#endif

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer tata = toto.random_exact(petits,gr);
		if ( tata.bitsize() != petits ) {
			//        cout << tata << endl;
			cout << "random_exact_exp  failed" << endl;
			return -1  ;
		}
	}
#ifdef GIVARO_DEBUG
	cout << "random_exact_exp.....OK" << endl;
#endif

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer tata = toto.random_exact(autre,gr);
		//        cout << tata << endl;
		if ( tata.bitsize() != autre.bitsize() ){
			cout << "random_exact  failed" << endl;
			return -1  ;
		}
	}
#ifdef GIVARO_DEBUG
	cout << "random_exact.........OK" << endl;
#endif

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer tata = toto.random_lessthan_2exp(grands,gr);
		if ( tata.bitsize() > grands ) {
			//        cout << tata << endl;
			cout << "random_lessthan_2exp  failed" << endl;
			return -1  ;
		}
	}
#ifdef GIVARO_DEBUG
	cout << "random_lessthan_2exp.....OK" << endl;
#endif

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer tata = toto.random_lessthan(large,gr);
		//        cout << tata << endl;
		if ( tata.bitsize() > large.bitsize() ){
			cout << "random_lessthan  failed" << endl;
			return -1  ;
		}
	}
#ifdef GIVARO_DEBUG
	cout << "random_lessthan.........OK" << endl;
#endif

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer tata = toto.nonzerorandom(petits,gr);
		if (tata == 0 || tata >= 1<<petits) {
			//        cout << tata << endl;
			cout << "nonzerorandom  failed" << endl;
			return -1  ;
		}
	}
#ifdef GIVARO_DEBUG
	cout << "nonzerorandom_exp....OK" << endl;
#endif

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer tata = toto.nonzerorandom(autre,gr);
		if (tata == 0 || tata >= autre) {
			       cout << tata << endl;
			cout << tata << "nonzerorandom  failed" << endl;
			return -1  ;
		}
	}
#ifdef GIVARO_DEBUG
	cout << "nonzerorandom........OK" << endl;
#endif

	return 0;

}

//! tests \verbatim ret= ::func(arg,arg); \endverbatim
int test2()
{
	GivRandom gr;
	Integer un(26);
	Integer autre(511);
	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer tata = Integer::random_between(un,autre,gr);
		//cout << tata << endl;
		if (tata < un || tata >= autre) {
			cout << "random_between  failed" << endl;
			return -1  ;
		}
	}
#ifdef GIVARO_DEBUG
	cout << "random_between.......OK" << endl;
#endif

	unsigned long trois = 3 ;
	unsigned long petits = 6 ;
	//std::vector<int> T(1<<petits) ;

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer tata = Integer::random_between(trois,petits,gr);
		//cout << tata << endl;
		//T[tata] += 1 ;
		if (tata < (1<<trois) || tata >= (1<<petits) ) {
			cout << "random_between_exp  failed" << endl;
			return -1  ;
		}
	}
#ifdef GIVARO_DEBUG
	//    for (size_t i = 0 ; i < 1<<petits ; ++i) cout << T[i] << " " ;
	cout << "random_between_exp...OK" << endl;
#endif

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer tata = Integer::random_exact(petits,gr);
		if ( tata.bitsize() != petits ) {
			//        cout << tata << endl;
			cout << "random_exact_exp  failed" << endl;
			return -1  ;
		}
	}
#ifdef GIVARO_DEBUG
	cout << "random_exact_exp.....OK" << endl;
#endif

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer tata = Integer::random_exact(autre,gr);
		//        cout << tata << endl;
		if ( tata.bitsize() != autre.bitsize() ){
			cout << "random_exact  failed" << endl;
			return -1  ;
		}
	}
#ifdef GIVARO_DEBUG
	cout << "random_exact.........OK" << endl;
#endif

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer tata = Integer::nonzerorandom(petits,gr);
		if (tata == 0 || tata >= 1<<petits) {
			//        cout << tata << endl;
			cout << "nonzerorandom  failed" << endl;
			return -1  ;
		}
	}
#ifdef GIVARO_DEBUG
	cout << "nonzerorandom_exp....OK" << endl;
#endif

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer tata = Integer::nonzerorandom(autre,gr);
		if (tata == 0 || tata >= autre) {
			//        cout << tata << endl;
			cout << "nonzerorandom  failed" << endl;
			return -1  ;
		}
	}
#ifdef GIVARO_DEBUG
	cout << "nonzerorandom........OK" << endl;
#endif

	return 0;

}

//! tests \verbatim ::func(ret,arg,arg); \endverbatim
int test3()
{
	GivRandom gr;
	Integer un(26);
	Integer autre(511);
	Integer tata ;
	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer::random_between(tata,un,autre,gr);
		//cout << tata << endl;
		if (tata < un || tata >= autre) {
			cout << "random_between  failed" << endl;
			return -1  ;
		}
	}
#ifdef GIVARO_DEBUG
	cout << "random_between.......OK" << endl;
#endif

	unsigned long trois = 3 ;
	unsigned long petits = 6 ;
	//std::vector<int> T(1<<petits) ;

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer::random_between(tata,trois,petits,gr);
		//cout << tata << endl;
		//T[tata] += 1 ;
		if (tata < (1<<trois) || tata >= (1<<petits) ) {
			cout << "random_between_exp  failed" << endl;
			return -1  ;
		}
	}
#ifdef GIVARO_DEBUG
	//    for (size_t i = 0 ; i < 1<<petits ; ++i) cout << T[i] << " " ;
	cout << "random_between_exp...OK" << endl;
#endif

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer::random_exact(tata,petits,gr);
		if ( tata.bitsize() != petits ) {
			//        cout << tata << endl;
			cout << "random_exact_exp  failed" << endl;
			return -1  ;
		}
	}
#ifdef GIVARO_DEBUG
	cout << "random_exact_exp.....OK" << endl;
#endif

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer::random_exact(tata,autre,gr);
		//        cout << tata << endl;
		if ( tata.bitsize() != autre.bitsize() ){
			cout << "random_exact  failed" << endl;
			return -1  ;
		}
	}
#ifdef GIVARO_DEBUG
	cout << "random_exact.........OK" << endl;
#endif
	const Integer tutu(autre);
	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer::random_exact(tata,tutu,gr);
		//        cout << tata << endl;
		if ( tata.bitsize() != tutu.bitsize() ){
			cout << "random_exact  failed" << endl;
			return -1  ;
		}
	}
#ifdef GIVARO_DEBUG
	cout << "random_exact.........OK" << endl;
#endif

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer::nonzerorandom(tata,petits,gr);
		if (tata == 0 || tata >= 1<<petits) {
			//        cout << tata << endl;
			cout << "nonzerorandom  failed" << endl;
			return -1  ;
		}
	}
#ifdef GIVARO_DEBUG
	cout << "nonzerorandom_exp....OK" << endl;
#endif

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer::nonzerorandom(tata,autre,gr);
		if (tata == 0 || tata >= autre) {
			//        cout << tata << endl;
			cout << "nonzerorandom  failed" << endl;
			return -1  ;
		}
	}
#ifdef GIVARO_DEBUG
	cout << "nonzerorandom........OK" << endl;
#endif

	return 0;

}

//! test possibly <0 random numbers
int test4()
{
	GivRandom gr;
	Integer un(26);
	Integer autre(511);
	Integer tata ;
	int count = 0 ;
	for (size_t i = 0 ; i < 5000 ; ++i) {
		if (Integer::RandBool(gr)) Integer::negin(un) ;
		if (Integer::RandBool(gr)) Integer::negin(autre) ;
		if (un>autre) std::swap(un,autre); // un < autre
		Integer::random_between(tata,un,autre,gr);
		//cout << tata << endl;
		if (tata < un || tata >= autre) {
			cout << "random_between  failed" << endl;
			return -1  ;
		}
		if (tata<0) ++count;
	}
	if ( !count ) {
		cout << "random_between  failed" << endl;
		return -1  ;
	}

#ifdef GIVARO_DEBUG
	cout << "random_between.......OK" << endl;
#endif

	unsigned long petits = 6 ;
	//std::vector<int> T(1<<petits) ;

	count = 0 ;
	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer::random_exact<false>(tata,petits,gr);
		if (tata<0) ++count ;
		if ( tata.bitsize() != petits) {
			//cout << tata << endl;
			cout << "random_exact_exp  failed" << endl;
			return -1  ;
		}
	}
	if ( !count ) {
		cout << "random_exact_exp  failed" << endl;
		return -1  ;
	}

#ifdef GIVARO_DEBUG
	cout << "random_exact_exp.....OK" << endl;
#endif

	count = 0 ;
	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer::random_exact<false>(tata,autre,gr);
		//        cout << tata << endl;
		if (tata<0) ++count ;
		if ( tata.bitsize() != autre.bitsize() ){
			cout << "random_exact  failed" << endl;
			return -1  ;
		}
	}
	if ( ! count ){
		cout << "random_exact  failed" << endl;
		return -1  ;
	}

#ifdef GIVARO_DEBUG
	cout << "random_exact.........OK" << endl;
#endif

	count = 0 ;
	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer::nonzerorandom<false>(tata,petits,gr);
		if (tata<0) count ++ ;
		if (tata == 0 || tata >= 1<<petits) {
//             cout << tata << endl;
			cout << "nonzerorandom  failed" << endl;
			return -1  ;
		}
	}
	if (!count){
		cout << "nonzerorandom  failed" << endl;
		return -1  ;
	}

#ifdef GIVARO_DEBUG
	cout << "nonzerorandom_exp....OK" << endl;
#endif
    if (autre<0) Integer::negin(autre);


	count = 0 ;
	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer::nonzerorandom<false>(tata,autre,gr);
		if (tata<0) ++ count ;
		if (tata == 0 || tata >= autre || (-tata) >= autre) {
			cout << "nonzerorandom  failed" << endl;
			return -1  ;
		}
	}
	if (!count)
	{
		cout << "nonzerorandom  failed" << endl;
		return -1  ;
	}
#ifdef GIVARO_DEBUG
	cout << "nonzerorandom........OK" << endl;
#endif

	return 0;

}

//! tests standard interface
int test5()
{
	GivRandom gr;
	Integer tata = 0 ;
	size_t count = 0 ;
#if (__GIVARO_SIZEOF_LONG==8) // BB : todo init to int64_t (either 'long' (64) or 'long long' (32))
	Integer toto(int64_t(15615486489765487));
#else
	Integer toto(15615486);
#endif
	unsigned long l = 5 ;
	for (size_t i = 0 ; i < 5000 ; ++i) {
		tata = Integer::random<false>(gr) ;
		if (tata<0) ++ count ;
		if (tata.bitsize() > 8*sizeof(mp_limb_t)) {
			//cout << tata << endl;
			cout << "random()  failed" << endl;
			return -1  ;
		}
	}
	if (!count)
	{
		cout << "random() failed" << endl;
		return -1  ;
	}
	for (size_t i = 0 ; i < 5000 ; ++i) {
		tata = Integer::random(gr) ;
		if (tata<0 || tata.bitsize() > 8*sizeof(mp_limb_t)) {
			//cout << tata << endl;
			cout << "random()  failed" << endl;
			return -1  ;
		}
	}
	if (!count)
	{
		cout << "random() failed" << endl;
		return -1  ;
	}
#ifdef GIVARO_DEBUG
	cout << "random().............OK" << endl;
#endif

	count = 0 ;
	for (size_t i = 0 ; i < 5000 ; ++i) {
		tata = Integer::random<false>(l,gr) ;
		if (tata<0) ++ count ;
		if (tata.bitsize() > l) {
			// cout << tata << endl;
			cout << "random  failed" << endl;
			return -1  ;
		}
	}
	if (!count)
	{
		cout << "random() failed" << endl;
		return -1  ;
	}
	for (size_t i = 0 ; i < 5000 ; ++i) {
		tata = Integer::random(l,gr) ;
		if (tata<0 || tata.bitsize() > l) {
			// cout << tata << endl;
			cout << "random  failed" << endl;
			return -1  ;
		}
	}
	if (!count)
	{
		cout << "random failed" << endl;
		return -1  ;
	}
#ifdef GIVARO_DEBUG
	cout << "random...............OK" << endl;
#endif

	count = 0 ;
	for (size_t i = 0 ; i < 5000 ; ++i) {
		tata = Integer::random<false>(toto,gr) ;
		if (tata<0) ++ count ;
		if (tata>toto || tata < -toto) {
			//cout << tata << endl;
			cout << "random  failed" << endl;
			return -1  ;
		}
	}
	if (!count)
	{
		cout << "random failed" << endl;
		return -1  ;
	}
	for (size_t i = 0 ; i < 5000 ; ++i) {
		tata = Integer::random(toto,gr) ;
		if (tata<0 || tata>toto) {
			//cout << tata << endl;
			cout << "random  failed" << endl;
			return -1  ;
		}
	}
	if (!count)
	{
		cout << "random failed" << endl;
		return -1  ;
	}
#ifdef GIVARO_DEBUG
	cout << "random...............OK" << endl;
#endif
	return 0 ;

}

#include <givaro/modular.h>
#include <givaro/montgomery.h>
#include <givaro/modular-ruint.h>
#include <givaro/modular-balanced.h>

template<typename Field>
int modtest(Field& F, uint64_t seed) {
    typedef typename Field::RandIter RI;
    RI gen(F);
    
    typename Field::Element a,b,c,d;
    F.init(a);F.init(b);F.init(c); F.init(d);
    
    gen.random(a);
//    gen(b);
    
    F.random(gen, c);

    return 0;
    
}



int main()
{
#ifdef GIVARO_DEBUG
	std::cout << "T1" << std::endl;
	if (test1()) return -1;
	std::cout << "T2" << std::endl;
	if (test2()) return -1;
	std::cout << "T3" << std::endl;
	if (test3()) return -1;
	std::cout << "T4" << std::endl;
	if (test4()) return -1;
	std::cout << "T5" << std::endl;
	if (test5()) return -1;
#else
	if (test1()) return -1;
	if (test2()) return -1;
	if (test3()) return -1;
	if (test4()) return -1;
	if (test5()) return -1;
#endif

	/* check seeding */
	GivRandom gr;
	Integer toto, tata ;
	gr.seed(3);
	toto = Integer::random(gr);
	gr.seed(3);
	tata = Integer::random(gr);
	if (toto != tata) return -1 ;

	// long unsigned int rangerand = (1<<20);
	// Integer seed;
	// seed = Integer::random_lessthan_2exp(rangerand); // <- l'erreur est dûe à cette ligne
	// cout << "seed : " << seed << endl;
	// Integer seed(1<<20);
	// cout << "seed : " << seed << endl;
        int pass=0;

        { Modular<double> F1(11); pass += modtest(F1,3); }
        { Modular<float> F1(11); pass += modtest(F1,3); }
        { Modular<uint64_t> F1(11); pass += modtest(F1,3); }
        { Modular<uint32_t> F1(11); pass += modtest(F1,3); }
        { Modular<int64_t> F1(11); pass += modtest(F1,3); }
        { Modular<int32_t> F1(11); pass += modtest(F1,3); }
        { Modular<uint16_t> F1(11); pass += modtest(F1,3); }
        { Modular<int16_t> F1(11); pass += modtest(F1,3); }
        { Modular<int8_t> F1(11); pass += modtest(F1,3); }
        { Modular<uint8_t> F1(11); pass += modtest(F1,3); }
        { Modular<Log16> F1(11); pass += modtest(F1,3); }
        { Montgomery<int32_t> F1(11); pass += modtest(F1,3); }
        { Montgomery<RecInt::ruint128> F1(11); pass += modtest(F1,3); }
        // { Modular<RecInt::rint128> F1(11); pass += modtest(F1,3); }
        { Modular<RecInt::ruint128> F1(11); pass += modtest(F1,3); }
        { Modular<Integer> F1(11); pass += modtest(F1,3); }
        { Modular<Integer> F1(11); pass += modtest(F1,3); }
        { ModularBalanced<double> F1(11); pass += modtest(F1,3); }
        { ModularBalanced<float> F1(11); pass += modtest(F1,3); }
        { ModularBalanced<int64_t> F1(11); pass += modtest(F1,3); }
        { ModularBalanced<int32_t> F1(11); pass += modtest(F1,3); }
        { Modular<Integer> F1(11); pass += modtest(F1,3); }
        { Modular<Integer> F1(11); pass += modtest(F1,3); }
        

	return pass ;

}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
