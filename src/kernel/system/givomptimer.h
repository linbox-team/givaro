/* kernel/system/givtimer.h
 * Copyright (C) 2014 Givaro Team
 *
 * Written by JG Dumas
 *
 * ------------------------------------
 * Modified by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Added _start_t member to BaseTimer, so that stop () does not clobber the
 * class' memory of its start time. This allows it to be called repeatedly to
 * get elapsed times.
 * ------------------------------------
 *
 * See COPYING for license information.
 *
 */

/** @file givtimer.h
 * @ingroup system
 * @brief timer
 */

#ifndef __GIVARO_OMP_timer_H
#define __GIVARO_OMP_timer_H

#include <iostream>
#include <givaro/givconfig.h>
#include <omp.h>

namespace Givaro {
    //! OMP timer
    struct OMPTimer {
        double _c;
        void start() { _c = omp_get_wtime(); }
        void stop() { _c = omp_get_wtime() - _c; }
        void clear() { _c = 0.0; }
        double realtime() { return _c; }
        double usertime() { return _c; }
        double time() const { return _c; }
        friend std::ostream& operator<<(std::ostream& o, const OMPTimer& t) {
            return o << t._c << 's';
        }

        OMPTimer& operator =(const OMPTimer& t) { _c = t._c; return *this; }
        OMPTimer& operator+=(const OMPTimer& t) { _c += t._c; return *this; }
        OMPTimer& operator-=(const OMPTimer& t) { _c -= t._c; return *this; }
        OMPTimer  operator +(const OMPTimer& t) const
        {
            OMPTimer r; r._c = _c + t._c; return r;
        }
        OMPTimer  operator -(const OMPTimer& t) const
        {
            OMPTimer r; r._c = _c - t._c; return r;
        }
        OMPTimer  operator -() { OMPTimer r; r._c = - _c; return r; }
    };
} // namespace Givaro

#endif // __GIVARO_OMP_timer_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
