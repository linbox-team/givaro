// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz16table1.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: J.G. Dumas
// $Id: givzpz16table1.C,v 1.8 2011-02-04 14:11:46 jgdumas Exp $
// ==========================================================================
// Description:

#include <iostream>
#include "givaro/givzpz16table1.h"

namespace Givaro {

ZpzDom<Log16>::ZpzDom( Residu_t p ) :
	_p(p),_pmone(Residu_t(p-1)),zero(Rep(_pmone << 1)), one(0),mOne(Rep(_pmone>>1))
{
	int32_t i,j;

	// tab value: Domain -> Rep, or something very similar
	_tab_value2rep = new Power_t[_p];
	// tab power: Rep -> Domain
	_tab_rep2value = new Residu_t[_p];

	_tab_rep2value[0] = 1;
	_tab_value2rep[1] = 0;

	int32_t fourp = ((int32_t)p) << 2, fourpmone= ((int32_t)_pmone)<<2;

	int not_found = 1;
	Residu_t accu = 1;
	Residu_t seed =2;

	// -- Find a generator of the multiplicative group
	while (_p > 2 && not_found == 1)
	{
		for(i=1; i<_p; i++)
		{
			accu = Residu_t((accu * seed) % _p);
			//cout << i << " :: " << accu << endl;
			_tab_rep2value[i] = accu;
			if (accu == 1)
				break;
			_tab_value2rep[accu] = Rep(i);
		}
		if (accu != 1){
			std::cerr << "attempted to build Log16 field with non-prime base "<<_p<<", halting\n";
			return;
		}
		if (i ==_p-1) not_found = 0;
		else {
			do {
				seed = Residu_t(rand() % _p);
			} while ((seed ==0) && (seed !=1));
			//      if (seed < 0) seed += _p;
		}
	}
	// -- Set the zero at position 2 * _p - 2 in table
	_tab_value2rep[0] = zero;

	//cout << "Generateur: " << seed << endl;
	/*cout << "Table: Value -> Rep" << endl;
	  for(i=0; i<_p; i++) {
	  cout << i << " -> " << _tab_value2rep[i] << endl;
	  }
	  cout << "\nTable: Rep -> Value" << endl;
	  for(i=0; i<_p; i++) {
	  cout << i << " -> " << _tab_rep2value[i] << endl;
	  }
	  */
	// -- Table for multiplication
	_tab_mul = new Power_t[(size_t)fourp];
	for(j=0; j<_pmone; j++)
		_tab_mul[j] = (Rep)j;
	for(j=_pmone; j< (int32_t)zero; j++)
	       	_tab_mul[j] = Rep(j-_pmone);
	for(j=zero; j<= fourpmone; j++)
		_tab_mul[j] = zero;

	// -- Table for division and neg:
	_tab_div = &_tab_mul[_pmone];
	_tab_neg = &_tab_mul[_pmone/2];

	// -- Table for 1+value
	_tab_pone = new Power_t[(size_t)fourp];
	_tab_addone = &_tab_pone[(int32_t)(zero)];

	/* Pascal Giorgi 24/04/02
	   Error between _tab_rep2value and _tab_value2rep

	   for(j=0; j<_pmone; j++)
	   _tab_addone[j] = _tab_rep2value[ 1 + _tab_value2rep[j] ];
	   for(j=1-_pmone; j<0; j++)
	   _tab_addone[j] = _tab_rep2value[ 1 + _tab_value2rep[j + _pmone] ];

	   corrected by inversing the array
	   */
	for(j=0; j<_pmone; j++){
		if (_tab_rep2value[j] < _pmone)
			_tab_addone[j] = _tab_value2rep[ 1 + _tab_rep2value[j] ];
		else
			_tab_addone[j] = _tab_value2rep[0];
	}
	for(j=1-_pmone; j<0; j++){
		if (_tab_rep2value[j+_pmone] < _pmone)
			_tab_addone[j] = _tab_value2rep[ 1 + _tab_rep2value[j + _pmone] ];
		else
			_tab_addone[j] = _tab_value2rep[0];

	}
	for(j=_pmone; j<=(int32_t)zero; j++)
		_tab_addone[j] = 0;
	for(j=(int32_t)-zero; j<(int32_t)(1-_pmone); j++)
		_tab_addone[j] = (Rep)j;

	_tab_addone[_pmone / 2] = zero;
	_tab_addone[-_pmone / 2] = zero;


	// -- Table for 1-value
	_tab_mone = new Power_t[(size_t)fourp];
	_tab_subone = &_tab_mone[(int32_t)zero];

	for(j=_pmone; j<=(int32_t)zero; j++)
		_tab_subone[j] = 0;
	for(j=-zero; j<(int32_t)(1-3*_pmone/2); j++)
		_tab_subone[j] = Rep(j+_pmone/2);
	for(j=-3*_pmone/2; j<(1-_pmone); j++)
		_tab_subone[j] = Rep(j-_pmone/2);
	for(j=1-_pmone; j<(1-_pmone/2); j++)
		_tab_subone[j] = _tab_addone[j + _pmone/2 + _pmone];
	for(j=_pmone/2; j<_pmone; j++)
		_tab_subone[j] = _tab_addone[j - _pmone/2];

	for(j=-_pmone/2; j<_pmone/2; j++)
		_tab_subone[j] = _tab_addone[j+_pmone/2];

	numRefs = new int;
	(*numRefs) = 1;
#ifdef GIVARO_DEBUG
	std::cout << *(numRefs) << " Arefs, p="<<_p<<" \n";
#endif

	// -- temporary
}

ZpzDom<Log16>::ZpzDom(const ZpzDom<Log16>& F) :
	_p ( F._p),
	_pmone ( F._pmone),
	_tab_value2rep ( F._tab_value2rep),
	_tab_rep2value ( F._tab_rep2value),
	_tab_mul ( F._tab_mul),
	_tab_div ( F._tab_div),
	_tab_neg ( F._tab_neg),
	_tab_addone ( F._tab_addone),
	_tab_subone ( F._tab_subone),
	_tab_mone ( F._tab_mone),
	_tab_pone ( F._tab_pone),
	numRefs ( F.numRefs),

	zero(F.zero), one(F.one),mOne(F.mOne)
{
  (*numRefs)++;
#ifdef GIVARO_DEBUG
  std::cout << *(numRefs) << " Brefs, p="<<_p<<" \n";
#endif
}


ZpzDom<Log16>& ZpzDom<Log16>::operator=( const ZpzDom<Log16>& F)
{

	F.assign(const_cast<Element&>(one),F.one);
	F.assign(const_cast<Element&>(zero),F.zero);
	F.assign(const_cast<Element&>(mOne),F.mOne);


  if (this->numRefs) {
    (*(this->numRefs))--;
#ifdef GIVARO_DEBUG
    std::cout << *(this->numRefs) << " Crefs, p="<<this->_p<<" \n";
#endif
    if ((*(this->numRefs))==0) {
#ifdef GIVARO_DEBUG
      std::cout << "zero : " << zero << std::endl;
      std::cout << "Ddestroying, p="<<residu()<<"\n";
#endif
      delete [] _tab_value2rep;
      delete [] _tab_rep2value;
      delete [] _tab_mul;
      delete [] _tab_mone;
      delete [] _tab_pone;
//       delete [] (&_tab_addone[-zero]);
//       delete [] (&_tab_subone[-zero]);
      delete numRefs;
    }
  }

  this->_p = F.residu();
  this->_pmone = F._pmone;
  this->_tab_value2rep = F._tab_value2rep;
  this->_tab_rep2value = F._tab_rep2value;
  this->_tab_mul = F._tab_mul;
  this->_tab_div = F._tab_div;
  this->_tab_neg = F._tab_neg;
  this->_tab_mone = F._tab_mone;
  this->_tab_pone = F._tab_pone;
  this->_tab_addone = F._tab_addone;
  this->_tab_subone = F._tab_subone;
  this->numRefs = F.numRefs;
  (*(this->numRefs))++;
#ifdef GIVARO_DEBUG
  std::cout << *(this->numRefs) << " Erefs, p = "<<this->_p<<"\n";
#endif

  return *this;
}


ZpzDom<Log16>::~ZpzDom()
{
  (*numRefs)--;
  if (*numRefs == 0) {
#ifdef GIVARO_DEBUG
    std::cout << "zero : " << zero << std::endl;
    std::cout << "Fdestroying, p="<<residu()<<"\n";
#endif
    delete [] _tab_value2rep;
    delete [] _tab_rep2value;
    delete [] _tab_mul;
    delete [] _tab_mone;
    delete [] _tab_pone;
//     delete [] (&_tab_addone[-(int32_t)zero]);
//     delete [] (&_tab_subone[-(int32_t)zero]);
    delete numRefs;
  }
}

void ZpzDom<Log16>::Init()
{
}

void ZpzDom<Log16>::End()
{
}

} // namespace Givaro

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
