// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/system/giverror.C,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: giverror.C,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:
// - error exception 

#include "givaro/giverror.h"
#include <iostream>

std::ostream& GivError::print( std::ostream& o ) const
{ return o << strg ; }


void GivError::throw_error( const GivError& err ) 
{
  throw err;
}

std::ostream& operator<< (std::ostream& o, const GivError& E) 
{
   return E.print(o) ; 
}

