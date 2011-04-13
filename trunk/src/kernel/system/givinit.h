// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/system/givinit.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givinit.h,v 1.3 2011-02-02 16:23:56 bboyer Exp $
// ==========================================================================
#ifndef __GIVARO_init_H
#define __GIVARO_init_H
#include "givaro/givconfig.h"
#include <iostream>

namespace Givaro {

// ==========================================================================
// --
// -- Description:
//   - Initialisation of GIVARO :
//     * handler to manage signal
//     * init the memory manager
//     * init all other modules

class Givaro {
public:
  //- Init of Givaro kernel :
static void Init(int* argc, char*** argv) ;
static void Init() ;

  //- End of Givaro kernel :
static void End() ;

  //- Return the version of the library
static const char* Version()  ;

  // Display the prompt of Givaro
static void DisplayVersion( std::ostream& ) ;
static void DisplayVersion();
} ;


// ==========================================================================
// --
// -- Main application class
// -- Could be not used
class GivaroAppli : public Givaro {
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
