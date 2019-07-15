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

extern "C" {
# include <sys/time.h>
# include <sys/resource.h>
    //  int getrusage (int, struct rusage*) ;
}

#include <iostream>
#include "givaro/givtimer.h"

namespace Givaro {

    // Return a value to initialize random generator
    int64_t BaseTimer::seed()
    {
        struct timespec ts;
#ifdef TIME_UTC
        timespec_get(&ts, TIME_UTC);
#else
# ifdef CLOCK_MONOTONIC
        clock_gettime(CLOCK_MONOTONIC, &ts);
# else
        clock_gettime(CLOCK_REALTIME, &ts);
# endif
#endif
        return static_cast<int64_t>(ts.tv_nsec);
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
    // between *this and T
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

    // Average timer
    const BaseTimer BaseTimer::operator / (const double nbiter) const
    {
        BaseTimer Tmp ;
        Tmp._t = _t / nbiter ;
        return Tmp ;
    }

    // Start timer
    void RealTimer::start()
    {
        struct timeval tmp2 ;
        gettimeofday (&tmp2, 0) ;

        // real time
        _t = _start_t = (double) tmp2.tv_sec +
        ((double) tmp2.tv_usec)/ (double)BaseTimer::MSPSEC ;
    }


    // Stop timer
    void RealTimer::stop()
    {
        struct timeval tmp2 ;
        gettimeofday (&tmp2, 0) ;

        // real time
        _t = (double) tmp2.tv_sec +
        ((double) tmp2.tv_usec)/ (double)BaseTimer::MSPSEC - _start_t ;
    }

    // Start timer
    void UserTimer::start()
    {
        struct rusage  tmp1 ;  // to getrusage (sys+user times)
        getrusage (RUSAGE_SELF, &tmp1) ;
        // user time
        _t = _start_t = (double) tmp1.ru_utime.tv_sec + ((double) tmp1.ru_utime.tv_usec)/ (double)MSPSEC ;
    }

    // Stop timer
    void UserTimer::stop()
    {
        struct rusage  tmp1 ;  // to getrusage (sys+user times)
        getrusage (RUSAGE_SELF, &tmp1) ;
        // user time
        _t = (double) tmp1.ru_utime.tv_sec + ((double) tmp1.ru_utime.tv_usec) / (double) MSPSEC - _start_t;
    }


    // Start timer
    void SysTimer::start()
    {
        struct rusage  tmp1 ;  // to getrusage (sys+user times)
        getrusage (RUSAGE_SELF, &tmp1) ;
        // user time
        _t = _start_t = (double) tmp1.ru_stime.tv_sec + ((double) tmp1.ru_stime.tv_usec)/ (double)MSPSEC ;
    }


    // Stop timer
    void SysTimer::stop()
    {
        struct rusage  tmp1 ;  // to getrusage (sys+user times)
        getrusage (RUSAGE_SELF, &tmp1) ;
        // user time
        _t = (double) tmp1.ru_stime.tv_sec +
        ((double) tmp1.ru_stime.tv_usec)/ (double)MSPSEC - _start_t ;
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
        ++_count;
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
    // between *this and T
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

    const Timer Timer::operator / (const double nbiter)  const
    {
        Timer Tmp ;
        Tmp.ut = ut / nbiter ;
        Tmp.st = st / nbiter ;
        Tmp.rt = rt / nbiter ;
        Tmp._count = 1;
        return Tmp ;
    }

} // namespace Givaro
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
