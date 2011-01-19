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

//#define GIVARO_DEBUG
using std::cout ; using std::endl;

//! tests <code>ret= .func(arg,arg);</code> ...
int test1()
{/*{{{*/
	Integer toto  ;
	toto.seeding((long unsigned int)0);
#ifdef GIVARO_DEBUG
	cout << "this is a random() number : " << toto.random() << endl;
#else
        toto.random();
#endif

	Integer un(26);
	Integer autre(511);
#ifdef GIVARO_DEBUG
	cout << "random...............OK" << endl;
#endif
	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer tata = toto.random_between(un,autre);
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
		Integer tata = toto.random_between(trois,petits);
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
		Integer tata = toto.random_exact(petits);
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
		Integer tata = toto.random_exact(autre);
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
		Integer tata = toto.nonzerorandom(petits);
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
		Integer tata = toto.nonzerorandom(autre);
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

}/*}}}*/

//! tests <code>ret= ::func(arg,arg);</code>
int test2()
{/*{{{*/
	Integer un(26);
	Integer autre(511);
	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer tata = Integer::random_between(un,autre);
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
		Integer tata = Integer::random_between(trois,petits);
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
		Integer tata = Integer::random_exact(petits);
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
		Integer tata = Integer::random_exact(autre);
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
		Integer tata = Integer::nonzerorandom(petits);
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
		Integer tata = Integer::nonzerorandom(autre);
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

}/*}}}*/

//! tests <code>::func(ret,arg,arg);</code>
int test3()
{/*{{{*/
	Integer un(26);
	Integer autre(511);
	Integer tata ;
	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer::random_between(tata,un,autre);
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
		Integer::random_between(tata,trois,petits);
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
		Integer::random_exact(tata,petits);
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
		Integer::random_exact(tata,autre);
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
		Integer::random_exact(tata,tutu);
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
		Integer::nonzerorandom(tata,petits);
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
		Integer::nonzerorandom(tata,autre);
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

}/*}}}*/

//! test possibly <0 random numbers
int test4()
{/*{{{*/
	Integer un(26);
	Integer autre(511);
	Integer tata ;
	int count = 0 ;
	for (size_t i = 0 ; i < 5000 ; ++i) {
		if (Integer::RandBool()) Integer::negin(un) ;
		if (Integer::RandBool()) Integer::negin(autre) ;
		if (un>autre) std::swap(un,autre); // un < autre
		Integer::random_between(tata,un,autre);
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
		Integer::random_exact<false>(tata,petits);
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
		Integer::random_exact<false>(tata,autre);
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
		Integer::nonzerorandom<false>(tata,petits);
		if (tata<0) count ++ ;
		if (tata == 0 || tata >= 1<<petits) {
			//        cout << tata << endl;
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

	count = 0 ;
	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer::nonzerorandom<false>(tata,autre);
		if (tata<0) ++ count ;
		if (tata == 0 || tata >= autre) {
			//        cout << tata << endl;
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

}/*}}}*/

//! tests standard interface
int test5()
{
	Integer tata = 0 ;
	size_t count = 0 ;
#if (__GIVARO_SIZEOF_LONG==8) // BB : todo init to int64 (either 'long' (64) or 'long long' (32))
	Integer toto(int64(15615486489765487));
#else
	Integer toto(15615486);
#endif
	unsigned long l = 5 ;
	for (size_t i = 0 ; i < 5000 ; ++i) {
		tata = Integer::random<false>() ;
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
		tata = Integer::random() ;
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
		tata = Integer::random<false>(l) ;
		if (tata<0) ++ count ;
		if (tata.bitsize() > l) {
			//cout << tata << endl;
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
		tata = Integer::random(l) ;
		if (tata<0 || tata.bitsize() > l) {
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

	count = 0 ;
	for (size_t i = 0 ; i < 5000 ; ++i) {
		tata = Integer::random<false>(toto) ;
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
		tata = Integer::random(toto) ;
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

int main()
{/*{{{*/
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
	Integer toto, tata ;
	Integer::seeding(3);
	toto = Integer::random();
	Integer::seeding(3);
	tata = Integer::random();
	if (toto != tata) return -1 ;

	return 0 ;

}/*}}}*/

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
