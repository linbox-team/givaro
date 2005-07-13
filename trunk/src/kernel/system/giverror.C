// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/system/giverror.C,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: giverror.C,v 1.2 2005-07-13 09:59:37 pernet Exp $
// ==========================================================================
// Description:
// - error exception 

#include "givaro/giverror.h"
#include <iostream>

std::ostream& GivError::print( std::ostream& o ) const
{ return o << strg ; }


GivError::~GivError(){}

GivMathError::~GivMathError(){}

GivBadFormat::~GivBadFormat(){}

GivMathDivZero::~GivMathDivZero(){}

void GivError::throw_error( const GivError& err ) 
{
  throw err;
}

std::ostream& operator<< (std::ostream& o, const GivError& E) 
{
   return E.print(o) ; 
}

