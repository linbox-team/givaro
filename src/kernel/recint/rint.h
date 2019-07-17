/* recint.h - RecInt library for fixed precision arithmetic

   Copyright Université Grenoble Alpes
Contributors :
Alexis BREUST (alexis.breust@gmail.com 2014)
Christophe CHABOT (christophechabotcc@gmail.com 2011)
Jean-Guillaume Dumas

Time-stamp: <20 Jun 12 10:31:24 Jean-Guillaume.Dumas@imag.fr>

This software is a computer program whose purpose is to provide an fixed precision arithmetic library.

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/


#ifndef RINT_H
#define RINT_H

#include "givaro-config.h"

/* Class definition */
#include "rrint.h"

/* Arithmetics in rint */
#include "radd.h"
#include "rsub.h"
#include "rmul.h"
#include "rdiv.h"
//#include "rexp.h"
//#include "rgcd.h"
//#include "rinvmod.h"
//#include "raddmul.h"

/* Comparisons of rint */
#include "rcmp.h"

/* Bit manipulation of rint */
#include "rfiddling.h"

/* Limb manipulation in rint */
//#include "rmanip.h"

/* Random in rint */
#include "rrandom.h"

/* Shift for rint */
//#include "rshift.h"

/* Display for rint */
#include "rdisplay.h"

/* Internal use for rint */
//#include "rtools.h"

/* GMP conversion system */
#include "rconvert.h"

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
