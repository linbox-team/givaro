/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/system/givtimer.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givtimer.C,v 1.4 2011-02-02 14:15:45 jgdumas Exp $
// ==========================================================================
// Description:
// - various timer objects
// - to be rewritten to be more efficient

#include <cmath>
#include "givaro/givconfig.h"

#if (GIVARO_SYS == _SYS_MACOS)
#include <time.h> // for use of clock()
#define USER_TIME ((double) ( (double)clock()/ (double)CLOCKS_PER_SEC))
#else
extern "C" {
# include <sys/time.h>
# include <sys/resource.h>
	//  int getrusage (int, struct rusage*) ;
}
#endif

#include <iostream>
#include "givaro/givtimer.h"

namespace Givaro {

// Return a value to initialize random generator
int64_t BaseTimer::seed()
{
#if (GIVARO_SYS == _SYS_MACOS)
	return clock() ;
#else
	struct timeval tp;
	gettimeofday(&tp, 0) ;
	return static_cast<int64_t>(tp.tv_usec);
#endif
}

// Output the value of the timer :
std::ostream& BaseTimer::print( std::ostream& o ) const
{
	return o << _t ;
}

// Some arithmetic operator :
BaseTimer& BaseTimer::operator = (const BaseTimer & T)
{
	_start_t = T._start_t ;
	_t = T._t ;
	return *this ;
}

// Computes and returns interval of time
// beteween *this and T
const BaseTimer BaseTimer::operator - (const BaseTimer & T) const
{
	BaseTimer Tmp ;
	Tmp._t = _t - T._t ;
	return Tmp ;
}

const BaseTimer BaseTimer::operator - ()
{
	BaseTimer Tmp ;
	Tmp._t = -_t ;
	return Tmp ;
}

const BaseTimer BaseTimer::operator + (const BaseTimer & T)  const
{
	BaseTimer Tmp ;
	Tmp._t = _t + T._t ;
	return Tmp ;
}

// Start timer
void RealTimer::start()
{
#if (GIVARO_SYS == _SYS_MACOS)
	_t = _start_t = USER_TIME ;
#else
	struct timeval tmp2 ;
	gettimeofday (&tmp2, 0) ;

	// real time
	_t = _start_t = (double) tmp2.tv_sec +
	((double) tmp2.tv_usec)/ (double)BaseTimer::MSPSEC ;
#endif
}


// Stop timer
void RealTimer::stop()
{
#if (GIVARO_SYS == _SYS_MACOS)
	_t = USER_TIME - _start_t ;
#else
	struct timeval tmp2 ;
	gettimeofday (&tmp2, 0) ;

	// real time
	_t = (double) tmp2.tv_sec +
	((double) tmp2.tv_usec)/ (double)BaseTimer::MSPSEC - _start_t ;
#endif
}

// Start timer
void UserTimer::start()
{
#if (GIVARO_SYS == _SYS_MACOS)
	_t = _start_t = USER_TIME ;
#else
	struct rusage  tmp1 ;  // to getrusage (sys+user times)
	getrusage (RUSAGE_SELF, &tmp1) ;
	// user time
	_t = _start_t = (double) tmp1.ru_utime.tv_sec + ((double) tmp1.ru_utime.tv_usec)/ (double)MSPSEC ;

#endif
}

// Stop timer
void UserTimer::stop()
{
#if (GIVARO_SYS == _SYS_MACOS)
	_t = USER_TIME - _start_t;
#else
	struct rusage  tmp1 ;  // to getrusage (sys+user times)
	getrusage (RUSAGE_SELF, &tmp1) ;
	// user time
	_t = (double) tmp1.ru_utime.tv_sec + ((double) tmp1.ru_utime.tv_usec) / (double) MSPSEC - _start_t;
#endif
}


// Start timer
void SysTimer::start()
{
#if (GIVARO_SYS == _SYS_MACOS)
	_t = _start_t = USER_TIME ; //! @bug here.
#else
	struct rusage  tmp1 ;  // to getrusage (sys+user times)
	getrusage (RUSAGE_SELF, &tmp1) ;
	// user time
	_t = _start_t = (double) tmp1.ru_stime.tv_sec + ((double) tmp1.ru_stime.tv_usec)/ (double)MSPSEC ;
#endif
}


// Stop timer
void SysTimer::stop()
{
#if (GIVARO_SYS == _SYS_MACOS)
	_t = USER_TIME - _start_t ;
#else
	struct rusage  tmp1 ;  // to getrusage (sys+user times)
	getrusage (RUSAGE_SELF, &tmp1) ;
	// user time
	_t = (double) tmp1.ru_stime.tv_sec +
	((double) tmp1.ru_stime.tv_usec)/ (double)MSPSEC - _start_t ;
#endif
}



// Clear timer :
void Timer::clear()
{
	rt.clear();
	ut.clear();
	st.clear();
	_count = 0;
}

// Start timer
void Timer::start()
{
	rt.start() ;
	ut.start() ;
	st.start() ;
}

// Stop timer
void Timer::stop()
{
	rt.stop() ;
	ut.stop() ;
	st.stop() ;
	_count = 1;
}


std::ostream& Timer::print( std::ostream& o ) const
{
	o << "user time: " << usertime() << '\n' ;
	o << "sys. time: " << systime() << '\n' ;
	return o << "real time: " << realtime() << std::endl ;
}

// Some arithmetic operator :
Timer& Timer::operator = (const Timer & T)
{
	ut = T.ut ;
	st = T.st ;
	rt = T.rt ;
	_count = T._count;
	return *this ;
}

// Comput._tes and returns interval of time
// beteween *this and T
const Timer Timer::operator - (const Timer & T)  const
{
	Timer Tmp ;
	Tmp.ut = ut - T.ut ;
	Tmp.st = st - T.st ;
	Tmp.rt = rt - T.rt ;
	Tmp._count = _count - T._count;
	return Tmp ;
}

const Timer Timer::operator - ()
{
	Timer Tmp ;
	Tmp.ut = -ut ;
	Tmp.st = -st ;
	Tmp.rt = -rt ;
	Tmp._count = - _count;
	return Tmp ;
}

const Timer Timer::operator + (const Timer & T)  const
{
	Timer Tmp ;
	Tmp.ut = ut + T.ut ;
	Tmp.st = st + T.st ;
	Tmp.rt = rt + T.rt ;
	Tmp._count = _count + T._count;
	return Tmp ;
}
} // namespace Givaro
