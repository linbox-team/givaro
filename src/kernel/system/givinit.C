// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/system/givinit.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givinit.C,v 1.2 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
// Description:
//  - Initialisation method

#include <iostream>
#include <string>
#include "givaro/givconfig.h"
#include "givaro/givinit.h"
#include "givaro/givmodule.h"
#include "givaro/givaromm.h"
#include "givaro/giverror.h"

#if (GIVAROSYS == MacOS)
#  if (__profile__ == 1)
#    pragma only_std_keywords off
#    include <profiler.h>
#    pragma only_std_keywords reset
#  endif
#endif

namespace Givaro {

    void GivaroMain::DisplayVersion()
    {
        GivaroMain::DisplayVersion( std::cout ) ;
    }

    std::ostream& GivaroMain::DisplayVersion(std::ostream& o)
    {

        o<<'\n' ;
        o<<"        /\\ \n" ;
        o<<"       /  \\    /\\      GIVARO : Parallel Algebraic Computing\n" ;
        o<<"      /\\__/\\  /  \\     by the Givaro Team\n" ;
        o<<"     /      \\/\\__/\\    All rights reserved, see copyright file.\n" ;
        o<<"    /                \\ " ;
        o<<"   /     Givaro-1.0   \\  Authors:\n" ;
        o<<"  /    (c) 1987-1998   \\    Th. Gautier, J.L. Roch, G.Villard\n" ;
        o<<" /       Givaro-4.0     \\  main co-Authors:\n" ;
        o<<"/      (c) 1998-2019     \\   J-G. Dumas, P. Giorgi, C. Pernet\n" ;
        o<<"--   -   -  -  -  --\n" ;
        return o<< "version: " << GivaroMain::Version() << std::endl;
    }

    const std::string GivaroMain::Version()
    {
        std::string Givaro_version("$Revision: ");
        Givaro_version += std::to_string(GIVARO_VERSION);
        Givaro_version.append(" GIVAROSYS",10);
        return Givaro_version;
    }

    void GivaroMain::Init(int* argc, char***argv)
    {
        GivModule::InitApp(argc, argv) ;
    }

    void GivaroMain::Init()
    {
        GivaroMain::Init(0,0) ;
    }

    // End of Givaro kernel :
    void GivaroMain::End()
    {
        GivModule::EndApp() ;
    }

    // -- run: must be call by the user on its application object
    int GivaroAppli::run( int argc, char** argv)
    {
        int res = 0;
#if (GIVAROSYS == MacOS)
#if ( __profile__ == 1)
        OSErr err = ProfilerInit(collectDetailed, bestTimeBase, 100,10);
        if (err != noErr) {
            cout << GivError("[Givaro::Init]: cannot initiliaze profile.");
            return 0;
        }
        ProfilerSetStatus(true);
#endif
#endif
        // -- call main function
        try {
            GivaroMain::Init(&argc, &argv);
            res = main(argc,argv);
            GivaroMain::End();
        }
        catch (const GivError& E) { std::cout << E << std::endl; }
        catch (...) { std::cout << "[GivaroAppli::run]: an error has occurred." << std::endl; }

#if (GIVAROSYS == MacOS)
#if (__profile__ == 1)
        ProfilerSetStatus(false);
        err = ProfilerDump("\pprofile.log");
        ProfilerTerm( );
        if (err != noErr) throw GivError("[Givaro::End]: cannot dump profile results.");
#endif
#endif

        return res;
    }

} // namespace Givaro
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
