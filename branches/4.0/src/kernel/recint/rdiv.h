/* rint/arith.h - Arithmetic functions for rint

Copyright Université Joseph Fourier - Grenoble
Contributors :
    Alexis BREUST (alexis.breust@gmail.com 2014)
	Christophe CHABOT (christophechabotcc@gmail.com 2011)
    Jean-Guillaume DUMAS

Time-stamp: <20 Jun 12 10:28:29 Jean-Guillaume.Dumas@imag.fr>

This software is a computer program whose purpose is to provide an
fixed precision arithmetic library.

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


#ifndef RINT_ARITH_DIV_H
#define RINT_ARITH_DIV_H

#include <cassert>

#include "rrint.h"
#include "rfiddling.h"
#include "rudiv.h"
#include "rcmp.h"
#include <iostream>
// --------------------------------------------------------------
// ----------------------- DEFINTIONS ---------------------------

namespace RecInt
{
    template <size_t K> rint<K>& operator%=(rint<K>&, const rint<K>&);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, rint<K>&) operator%=(rint<K>&, const T&);

    template <size_t K> rint<K>& operator/=(rint<K>&, const rint<K>&);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, rint<K>&) operator/=(rint<K>&, const T&);

    template <size_t K> rint<K> operator%(const rint<K>&, const rint<K>&);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, rint<K>) operator%(const rint<K>&, const T&);

    template <size_t K> rint<K> operator/(const rint<K>&, const rint<K>&);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, rint<K>) operator/(const rint<K>&, const T&);

    // q = floor(a/b)
    template <size_t K> rint<K>& div_q(rint<K>& q, const rint<K>& a, const rint<K>& b);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, rint<K>&) div_q(rint<K>& q, const rint<K>& a, const T& b);
    
    // r = a mod b
    template <size_t K> rint<K>& div_r(rint<K>& r, const rint<K>& a, const rint<K>& b);
    template <size_t K, typename T> __RECINT_IS_ARITH(T, T&) div_r(T& r, const rint<K>& a, const T& b);
}


// --------------------------------------------------------------
// ------------------------ Operators ---------------------------

namespace RecInt
{
    // Operator %=
    template <size_t K>
    inline rint<K>& operator%=(rint<K>& a, const rint<K>& b) {
        return div_r(a, a, b);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, rint<K>&) operator%=(rint<K>& a, const T& b) {
        T aa;
        div_r(aa, a, b);
        return (a = aa);
    }

    // Operator /=
    template <size_t K>
    inline rint<K>& operator/=(rint<K>& a, const rint<K>& b) {
        return div_q(a, a, b);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, rint<K>&) operator/=(rint<K>& a, const T& b) {
        return div_q(a, a, b);
    }

    // Operator %
    template <size_t K>
    inline rint<K> operator%(const rint<K>& b, const rint<K>& c) {
        rint<K> a;
        div_r(a, b, c);
        return a;
    }

    // Operator /
    template <size_t K>
    inline rint<K> operator/(const rint<K>& b, const rint<K>& c) {
        rint<K> a;
        div_q(a, b, c);
        return a;
    }
}


// --------------------------------------------------------------
// ------------------------ Division ---------------------------

namespace RecInt
{
    // r = a mod b
    template <size_t K>
    inline rint<K>& div_q(rint<K>& q, const rint<K>& a, const rint<K>& b) {
        if (a.isNegative()) {
	    if (b.isNegative()) {
		div_q(q.Value, (-a).Value, (-b).Value);
	    }
	    else {
		div_q(q.Value, (-a).Value, b.Value);
		q = -q;
	    }
        }
        else {
	    if (b.isNegative()) {
		div_q(q.Value, a.Value, (-b).Value);
		q = -q;
	    }
	    else {
		div_q(q.Value, a.Value, b.Value);
	    }
        }
	
    	return q;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, rint<K>&) div_q(rint<K>& q, const rint<K>& a, const T& b) {
        // TODO
	std::cerr << "DIV with rint not implemented yet." << std::endl;
	return q;
    }
    
    // r = a mod b
    template <size_t K>
    inline rint<K>& div_r(rint<K>& r, const rint<K>& a, const rint<K>& b) {
	assert(b > 1);
	
        if (a.isNegative()) {
            div_r(r.Value, (-a).Value, b.Value);
            r = -r;
        }
        else {
            div_r(r.Value, a.Value, b.Value);
        }

    	return r;
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, T&) div_r(T& r, const rint<K>& a, const T& b) {
	assert(b > 1);
	
        if (a.isNegative()) {
            rint<K> aa(-a);
            div_r(r, aa.Value, b);
    	    return -r;
        }
        else {
            div_r(r, a.Value, b);
    	    return r;
        }
    }
}

#endif

