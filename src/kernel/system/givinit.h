// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/system/givinit.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givinit.h,v 1.3 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================
/** @file givinit.h
 * @ingroup system
 * @brief NO DOC
 */
#ifndef __GIVARO_init_H
#define __GIVARO_init_H
#include "givaro/givconfig.h"
#include <iostream>
#include <string>

//! Namespace in which the whole Givaro library resides.
namespace Givaro {

    // ==========================================================================
    /**  Initialisation of GIVARO .
     * - handler to manage signal
     * - init the memory manager
     * - init all other modules
     * .
     */
    class GivaroMain {
    public:
        //- Init of Givaro kernel :
        static void Init(int* argc, char*** argv) ;
        static void Init() ;

        //- End of Givaro kernel :
        static void End() ;

        //- Return the version of the library
        static const std::string Version()  ;

        // Display the prompt of Givaro
        static std::ostream& DisplayVersion( std::ostream& ) ;
        static void DisplayVersion();
    } ;


    // ==========================================================================
    // --
    //! Main application class
    //! Could be not used
    class GivaroAppli : public GivaroMain {
    public:
        //- Cstor, destor
        GivaroAppli() {};
        virtual ~GivaroAppli(){};

        //- main: must redefined by derived class
        virtual int main(int argc=0, char**argv=0) = 0;

        //- run: must be call by the user on its application object
        int run( int argc=0, char** argv=0);
    };

} // namespace Givaro
#endif // __GIVARO_init_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
