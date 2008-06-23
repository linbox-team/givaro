// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givdegree.C,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givdegree.C,v 1.2 2008-06-23 13:44:02 jgdumas Exp $
// ==========================================================================

#include "givaro/givdegree.h"

// -- Degree of zero polynomial
#ifndef __ECC
const long Degree::deginfty = DEGPOLYZERO;
#endif

