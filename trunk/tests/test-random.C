// Copyright(c)'2010 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.
// written by BB.

#include <iostream>
#include <givaro/givinteger.h>

using std::cout ; using std::endl;

int test1()
{/*{{{*/
	Integer toto  ;
	toto.seeding((long unsigned int)0);
	cout << "this is a random() number : " << toto.random() << endl;

	Integer un(26);
	Integer autre(511);
	cout << "random...............OK" << endl;
	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer tata = toto.random_between(un,autre);
		//cout << tata << endl;
		if (tata < un || tata >= autre) {
			cout << "random_between  failed" << endl;
			return -1  ;
		}
	}
	cout << "random_between.......OK" << endl;

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
	cout << "random_between_exp...OK" << endl;

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer tata = toto.random_exact(petits);
		if ( tata.bitsize() != petits ) {
			//        cout << tata << endl;
			cout << "random_exact_exp  failed" << endl;
			return -1  ;
		}
	}
	cout << "random_exact_exp.....OK" << endl;

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer tata = toto.random_exact(autre);
		//        cout << tata << endl;
		if ( tata.bitsize() != autre.bitsize() ){
			cout << "random_exact  failed" << endl;
			return -1  ;
		}
	}
	cout << "random_exact.........OK" << endl;

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer tata = toto.nonzerorandom(petits);
		if (tata == 0 || tata >= 1<<petits) {
			//        cout << tata << endl;
			cout << "nonzerorandom  failed" << endl;
			return -1  ;
		}
	}
	cout << "nonzerorandom_exp....OK" << endl;

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer tata = toto.nonzerorandom(autre);
		if (tata == 0 || tata >= autre) {
			//        cout << tata << endl;
			cout << "nonzerorandom  failed" << endl;
			return -1  ;
		}
	}
	cout << "nonzerorandom........OK" << endl;

	return 0;

}/*}}}*/

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
	cout << "random_between.......OK" << endl;

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
	//    for (size_t i = 0 ; i < 1<<petits ; ++i) cout << T[i] << " " ;
	cout << "random_between_exp...OK" << endl;

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer tata = Integer::random_exact(petits);
		if ( tata.bitsize() != petits ) {
			//        cout << tata << endl;
			cout << "random_exact_exp  failed" << endl;
			return -1  ;
		}
	}
	cout << "random_exact_exp.....OK" << endl;

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer tata = Integer::random_exact(autre);
		//        cout << tata << endl;
		if ( tata.bitsize() != autre.bitsize() ){
			cout << "random_exact  failed" << endl;
			return -1  ;
		}
	}
	cout << "random_exact.........OK" << endl;

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer tata = Integer::nonzerorandom(petits);
		if (tata == 0 || tata >= 1<<petits) {
			//        cout << tata << endl;
			cout << "nonzerorandom  failed" << endl;
			return -1  ;
		}
	}
	cout << "nonzerorandom_exp....OK" << endl;

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer tata = Integer::nonzerorandom(autre);
		if (tata == 0 || tata >= autre) {
			//        cout << tata << endl;
			cout << "nonzerorandom  failed" << endl;
			return -1  ;
		}
	}
	cout << "nonzerorandom........OK" << endl;

	return 0;

}/*}}}*/

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
	cout << "random_between.......OK" << endl;

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
	//    for (size_t i = 0 ; i < 1<<petits ; ++i) cout << T[i] << " " ;
	cout << "random_between_exp...OK" << endl;

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer::random_exact(tata,petits);
		if ( tata.bitsize() != petits ) {
			//        cout << tata << endl;
			cout << "random_exact_exp  failed" << endl;
			return -1  ;
		}
	}
	cout << "random_exact_exp.....OK" << endl;

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer::random_exact(tata,autre);
		//        cout << tata << endl;
		if ( tata.bitsize() != autre.bitsize() ){
			cout << "random_exact  failed" << endl;
			return -1  ;
		}
	}
	cout << "random_exact.........OK" << endl;

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer::nonzerorandom(tata,petits);
		if (tata == 0 || tata >= 1<<petits) {
			//        cout << tata << endl;
			cout << "nonzerorandom  failed" << endl;
			return -1  ;
		}
	}
	cout << "nonzerorandom_exp....OK" << endl;

	for (size_t i = 0 ; i < 5000 ; ++i) {
		Integer::nonzerorandom(tata,autre);
		if (tata == 0 || tata >= autre) {
			//        cout << tata << endl;
			cout << "nonzerorandom  failed" << endl;
			return -1  ;
		}
	}
	cout << "nonzerorandom........OK" << endl;

	return 0;

}/*}}}*/

int main() 
{/*{{{*/
	std::cout << "T1" << std::endl;
	if (test1()) return -1;
	std::cout << "T2" << std::endl;
	if (test2()) return -1;
	std::cout << "T3" << std::endl;
	if (test2()) return -1;

}/*}}}*/

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
