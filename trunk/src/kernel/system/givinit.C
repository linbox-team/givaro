// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/system/givinit.C,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givinit.C,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:
//  - Initialisation method

#include <iostream>
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


void Givaro::DisplayVersion()
{
  Givaro::DisplayVersion( std::cout ) ;
}

void Givaro::DisplayVersion(std::ostream& o)
{
 
o<<'\n' ;
o<<"        /\\ \n" ;
o<<"       /  \\    /\\      GIVARO : Parallel Algebraic Computing\n" ;
o<<"      /\\__/\\  /  \\     by the Givaro Team\n" ;
o<<"     /      \\/\\__/\\    All rights reserved, see copyright file.\n" ;
o<<"    /                \\ " ;
o<<"   /     Givaro-1.0   \\  Authors:\n" ;
o<<"  /    (c) 1987-1998   \\    Th. Gautier, J.L. Roch, M.Samama, G.Villard\n" ;
o<<" /       Givaro-3.0     \\  co-Authors:\n" ;
o<<"/      (c) 1998-2002     \\    J-G. Dumas, P. Giorgi\n" ;
o<<"--   -   -  -  -  --\n" ;
o<< "version: " << Givaro::Version() << std::endl;
} 

const char* Givaro::Version()
{
  static const char* Givaro_version ="$Revision: 1.1.1.1 $ for ""GIVAROSYS"; 
  return Givaro_version;
}

void Givaro::Init(int* argc, char***argv) 
{
  GivModule::InitApp(argc, argv) ;
}

void Givaro::Init() 
{
  Givaro::Init(0,0) ;
}

  // End of Givaro kernel :
void Givaro::End()
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
    Givaro::Init(&argc, &argv);
    res = main(argc,argv);
    Givaro::End();
  }
  catch (GivError E) { std::cout << E << std::endl; }
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
