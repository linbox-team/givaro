/* rmint/arith.h - Arithmetic functions for rmint

   Copyright Université Grenoble Alpes
Contributors :
Alexis BREUST (alexis.breust@gmail.com 2014)
Christophe CHABOT (christophechabotcc@gmail.com 2011)


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


#ifndef RMINT_COMMON_ARITH_MUL_H
#define RMINT_COMMON_ARITH_MUL_H

/** NOTE : For this common file, either basic/reduc.h or mg/reduc.h
  has to be pre-included. The same goes for (basic|mg)/arith/mul.h **/

// --------------------------------------------------------------
// ----------------------- DEFINTIONS ---------------------------

namespace RecInt
{
    template <size_t K, size_t MG> rmint<K, MG>& operator*=(rmint<K, MG>&, const rmint<K, MG>&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, rmint<K, MG>&) operator*=(rmint<K, MG>&, const T&);

    template <size_t K, size_t MG> rmint<K, MG> operator*(const rmint<K, MG>&, const rmint<K, MG>&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, rmint<K, MG>) operator*(const rmint<K, MG>&, const T&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, rmint<K, MG>) operator*(const T&, const rmint<K, MG>&);

    template <size_t K, size_t MG> void mul(rmint<K, MG>&, const rmint<K, MG>&, const rmint<K, MG>&);
    template <size_t K, size_t MG> void mul(rmint<K, MG>&, const rmint<K, MG>&);

    //----------- square -----------//

    // a = b * b mod a.p
    template <size_t K, size_t MG> void square(rmint<K, MG>& a, const rmint<K, MG>& b);

    // a *= a mod a.p
    template <size_t K, size_t MG> void square(rmint<K, MG>& a);
}


// --------------------------------------------------------------
// ------------------------ Operators ---------------------------

namespace RecInt
{
    // Operator *=
    template <size_t K, size_t MG>
    inline rmint<K, MG>& operator*=(rmint<K, MG>& a, const rmint<K, MG>& b) {
        mul(a, a, b);
        return a;
    }

    template <size_t K, size_t MG, typename T>
    inline __RECINT_IS_ARITH(T, rmint<K, MG>&) operator*=(rmint<K, MG>& a, const T& b) {
        mul(a, a, b);
        return a;
    }

    // Operator *
    template <size_t K, size_t MG>
    inline rmint<K, MG> operator*(const rmint<K, MG>& b, const rmint<K, MG>& c) {
        rmint<K, MG> a;
        mul(a, b, c);
        return a;
    }

    template <size_t K, size_t MG, typename T>
    inline __RECINT_IS_ARITH(T, rmint<K, MG>) operator*(const rmint<K, MG>& b, const T& c) {
        rmint<K, MG> a;
        mul(a, b, c);
        return a;
    }

    template <size_t K, size_t MG, typename T>
    inline __RECINT_IS_ARITH(T, rmint<K, MG>) operator*(const T& c, const rmint<K, MG>& b) {
        rmint<K, MG> a;
        mul(a, b, c);
        return a;
    }
}


// --------------------------------------------------------------
// ---------------------- Multiplication ------------------------

namespace RecInt
{
    // a = b * c mod a.p
    template <size_t K, size_t MG>
    inline void mul(rmint<K, MG>& a, const rmint<K, MG>& b, const rmint<K, MG>& c) {
        ruint<K+1> resmul;
        lmul(resmul, b.Value, c.Value);
        reduction(a, resmul);
    }

    // a *= b mod a.p
    template <size_t K, size_t MG>
    inline void mul(rmint<K, MG>& a, const rmint<K, MG>& b) {
        ruint<K+1> resmul;
        lmul(resmul, a.Value, b.Value);
        reduction(a, resmul);
    }
}


// --------------------------------------------------------------
// -------------------------- Square ----------------------------

namespace RecInt
{
    // a = b * b mod a.p
    template <size_t K, size_t MG>
    inline void square(rmint<K, MG>& a, const rmint<K, MG>& b) {
        ruint<K+1> res;
        lsquare(res, b.Value);
        reduction(a, res);
    }

    // a *= a mod a.p
    template <size_t K, size_t MG>
    inline void square(rmint<K, MG>& a) {
        square(a, a);
    }
}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
