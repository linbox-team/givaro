#ifndef _GIV_INIT_H_
#define _GIV_INIT_H_
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/system/givinit.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givinit.h,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
#include "givaro/givconfig.h"
#include <iostream>


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

#endif
