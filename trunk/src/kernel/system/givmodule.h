#ifndef _MODULE_INIT_H_
#define _MODULE_INIT_H_
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/system/givmodule.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givmodule.h,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:

// -------------------------------------------------------- Fwd declaration

class GivModule;


// ----------------------------------------------------------- GivaroNoInit:
// Purpose: used to delay construction of object in the init function
// of a module definition.

class GivaroNoInit {};


// -------------------------------------------------------------- InitAfter:
// Purpose: define a precedence relation between two modules.

class InitAfter {
public:
  InitAfter( const GivModule& MI );
static InitAfter Default;
static InitAfter First;
static InitAfter Last;
private:
  InitAfter( int p );
  const GivModule* M;
  int priority;
  friend class GivModule;
  int operator < ( const InitAfter& M ) const;
};


// -------------------------------------------------------------- GivModule:
// Purpose: definition of module with precedence relation use to initialize
// them between different units compilation.

class GivModule {
public:
  typedef void (*ptFuncInit)(int* argc, char** *argv);
  typedef void (*ptFuncEnd)();
  enum {
    MaxPriority  = -100000 ,          // - maximum priority
    MinPriority  = -MaxPriority,      // - minimum priority
    DfltPriority = 0,                 // - default priority
    UndefPriority = MaxPriority-1     // - use to build depedences
  };
  // - Cstor of a module with a priority
  GivModule ( ptFuncInit init, ptFuncEnd end,
               const int p, const char* n=0 );

  // - Cstor of a module with precedence relation between an other module
  GivModule ( ptFuncInit init, ptFuncEnd end,
               const InitAfter& M, const char* n=0 );
  ~GivModule ();

private:
  // - Call by the Givaro::Init and Givaro::End functions
static void InitApp(int* argc, char***argv);
static void EndApp();
friend class Givaro;

private: 
  // - Internal data of a module
  int priority;
  InitAfter which;
  ptFuncInit f_init;
  ptFuncEnd f_end;
  const char* name;

friend class InitAfter;
static void SortGivModule();
};


// -------------------------------------------------------------- GivModule:
// Purpose: definition of object to be initialized after all modules
// initialization

class ObjectInit {
public:
  // -- when call: link in a global list, then ...
  ObjectInit();
  // -- ... call init during the initialization phase
  virtual void objinit() {};
private:
  ObjectInit* _next;
  friend class GivModule;
};

#endif
