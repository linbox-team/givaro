// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/system/giverror.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software. 
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: giverror.h,v 1.3 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
// Description:
// - error exception 
#ifndef _GIVERROR_H_
#define _GIVERROR_H_

#include <iostream>

// ------------------------------- GivError
// - Base class for execption handling in Givaro
class GivError {
public:
  GivError(const char* msg =0 ) 
  : strg(msg) {};

  virtual ~GivError() ;
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
  virtual ~GivMathError() ;

  GivMathError(const char* msg = 0) : GivError(msg) { }
};

// -- Exception thrown in input of data structure 
class GivBadFormat : public GivError {
public:
  virtual ~GivBadFormat();
  GivBadFormat(const char* msg = 0) : GivError(msg) { }
};

class GivMathDivZero : public GivError {
public:
  virtual ~GivMathDivZero();
  GivMathDivZero(const char* msg = 0) : GivError(msg) { }
};



#endif
