// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/system/givperf.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givperf.h,v 1.2 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
// Description:
// - performance analysis
#ifndef _GIV_PERF_H_
#define _GIV_PERF_H_

#ifdef GIVARO_PERF
#include <stddef>
#include <iostream>


struct __CoutCounter {
  void (*print)(ostream&);
  __CoutCounter( void (*prn)(ostream& ) ) : print(prn) 
  {
    cout << "New Counter" << endl;
  }
  ~__CoutCounter( ) 
  { 
    cout << "Destroy a counter" << endl;
    (*print)(cout << endl); cout << endl; 
  }
};


// --- Class that store a set of counters:
// - _count_cstor: #cstor calls except the recopy constructor calls
// - _count_cstor_recopy: #recopy cstor calls
// - _count_assign: #assignment calls
// - _count_dstor: #dstor calls

#define GIVARO_PERF_DEFCLASS(Name,Type)					\
template<class Type>							\
struct _Giv_perf##Name {						\
static size_t _count_cstor;						\
static size_t _count_cstor_recopy;					\
static size_t _count_assign;						\
static size_t _count_dstor;						\
static void print(ostream& o);						\
static __CoutCounter _coutcout;						\
  _Giv_perf##Name() { ++_count_cstor; }					\
  _Giv_perf##Name(const _Giv_perf##Name<Type>& S ) 			\
   { ++_count_cstor_recopy; }						\
  ~_Giv_perf##Name() { ++_count_dstor; }				\
  _Giv_perf##Name<Type>& operator=(const _Giv_perf##Name<Type>& S ) 	\
   { ++_count_assign; }							\
};									\
template<class Type>							\
size_t _Giv_perf##Name<Type>::_count_cstor =0;				\
template<class Type>							\
size_t _Giv_perf##Name<Type>::_count_cstor_recopy =0;			\
template<class Type>							\
size_t _Giv_perf##Name<Type>::_count_assign =0;				\
template<class Type>							\
size_t _Giv_perf##Name<Type>::_count_dstor =0;				\
template<class Type>							\
void _Giv_perf##Name<Type>::print(ostream& o ) {			\
  o << #Name ":  #cstor=" << _count_cstor 				\
    << ", #recopy=" << _count_cstor_recopy				\
    << ", #destor=" << _count_dstor					\
    << ", #assign=" << _count_assign;					\
}									\
template<class Type>							\
__CoutCounter _Giv_perf##Name<Type>::_coutcout= &_Giv_perf##Name<Type>::print;

#define GIVARO_PERF_CSTOR(Name,Type)   _Giv_perf##Name<Type>::_count_cstor++;
#define GIVARO_PERF_RECOPY(Name,Type)  _Giv_perf##Name<Type>::_count_cstor_recopy++;
#define GIVARO_PERF_DSTOR(Name,Type)   _Giv_perf##Name<Type>::_count_dstor++;
#define GIVARO_PERF_ASSIGN(Name,Type)  _Giv_perf##Name<Type>::_count_assign++;

#define GIVARO_PERF_INEHERIT(Name,Type)		\
: public _Giv_perf##Name<Type>

#define GIVARO_PERF_DISPLAY(Name,Type)  _Giv_perf##Name<Type>::print(cout), cout << endl;


#else // #ifdef GIVARO_PERF

#define GIVARO_PERF_DEFCLASS(N,T)	
#define GIVARO_PERF_INEHERIT(N,T)
#define GIVARO_PERF_CSTOR(Name,Type)  
#define GIVARO_PERF_RECOPY(Name,Type) 
#define GIVARO_PERF_DSTOR(Name,Type) 
#define GIVARO_PERF_ASSIGN(Name,Type)
#define GIVARO_PERF_DISPLAY(Name,Type)  


#endif

#endif
