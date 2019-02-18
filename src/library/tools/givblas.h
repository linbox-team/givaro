//=============================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Description: BLAS function prototypes
//=============================================================

#ifndef __GIVARO_blas_H
#define __GIVARO_blas_H

#ifdef BLAS__SUNPRO_CC
#define dcopy dcopy_
#define daxpy daxpy_
#define dgemv dgemv_
#define dtrmv dtrmv_
#define dtrsv dtrsv_
#define dgemm dgemm_
#define dtrsm dtrsm_
#define dsyrk dsyrk_
#define dger dger_
#define dpotrf dpotrf_
#endif

#ifdef __cplusplus
extern "C" {
#endif
    void dcopy( const int*, const double*, const int*, double*, const int*);
    void daxpy( const int*, const double*, const double*, const int*, double*, const int*);


    void dgemv( const char*, const int*, const int*,
                const double*, const double*, const int*,
                const double*, const int*, const double*,
                double*, const int*);
    void dtrmv(const char*, const char*, const char*, const int*, const double*, const int*, double*, const int*);
    void dtrsv(const char*, const char*, const char*, const int*, const double*, const int*, double*, const int*);


    void dtrsm( const char*, const char*, const char*, const char*, const int*, const int*, const double*, const double*, const int*, double*, const int*);
    void dgemm(const char*, const char*, const int*, const int*, const int*, const double*, const double*, const int*, const double*, const int*, const double*, const double*, const int*) ;void dsyrk( const char*, const char*, const int*, const int*, const double*, const double*, const int*, double*, double*, const int*);

    void dpotrf( const char*, const int*, const double*, const int*, const int*) ;

#ifdef __cplusplus
}
#endif



#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
