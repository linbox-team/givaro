#ifndef _GIVERROR_H_
#define _GIVERROR_H_
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/system/giverror.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: giverror.h,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:
// - error exception 

#include <iostream>

// ------------------------------- GivError
// - Base class for execption handling in Givaro
class GivError {
public:
  GivError(const char* msg =0 ) 
  : strg(msg) {};

  // -- virtual print of the error message
  virtual std::ostream& print( std::ostream& o )  const;
  
  // -- non virtual output operator
  friend std::ostream& operator<< (std::ostream& o, const GivError& E) ;

  // - useful to setup a break point on it
  static void throw_error( const GivError& err );

protected:
  const char* strg;  
};

class GivMathError : public GivError {
public:
  GivMathError(const char* msg = 0) : GivError(msg) { }
};

// -- Exception thrown in input of data structure 
class GivBadFormat : public GivError {
public:
  GivBadFormat(const char* msg = 0) : GivError(msg) { }
};

class GivMathDivZero : public GivError {
public:
  GivMathDivZero(const char* msg = 0) : GivError(msg) { }
};



#endif
