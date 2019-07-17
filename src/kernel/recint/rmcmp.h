/* rmint/cmp.h - Comparison functions for rmint

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


#ifndef RMINT_CMP_H
#define RMINT_CMP_H

/** NOTE : For this common file, either basic/reduc.h or mg/reduc.h
  has to be pre-included. **/

// --------------------------------------------------------------
// ----------------------- DEFINTIONS ---------------------------

namespace RecInt
{
    template <size_t K, size_t MG> bool operator==(const rmint<K, MG>&, const rmint<K, MG>&);
    template <size_t K, size_t MG> bool operator==(const rmint<K, MG>&, const ruint<K>&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, bool) operator==(const rmint<K, MG>&, const T&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, bool) operator==(const T&, const rmint<K, MG>&);

    template <size_t K, size_t MG> bool operator!=(const rmint<K, MG>&, const rmint<K, MG>&);
    template <size_t K, size_t MG> bool operator!=(const rmint<K, MG>&, const ruint<K>&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, bool) operator!=(const rmint<K, MG>&, const T&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, bool) operator!=(const T&, const rmint<K, MG>&);

    template <size_t K, size_t MG> bool operator>(const rmint<K, MG>&, const rmint<K, MG>&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, bool) operator>(const rmint<K, MG>&, const T&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, bool) operator>(const T&, const rmint<K, MG>&);

    template <size_t K, size_t MG> bool operator>=(const rmint<K, MG>&, const rmint<K, MG>&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, bool) operator>=(const rmint<K, MG>&, const T&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, bool) operator>=(const T&, const rmint<K, MG>&);

    template <size_t K, size_t MG> bool operator<(const rmint<K, MG>&, const rmint<K, MG>&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, bool) operator<(const rmint<K, MG>&, const T&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, bool) operator<(const T&, const rmint<K, MG>&);

    template <size_t K, size_t MG> bool operator<=(const rmint<K, MG>&, const rmint<K, MG>&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, bool) operator<=(const rmint<K, MG>&, const T&);
    template <size_t K, size_t MG, typename T> __RECINT_IS_ARITH(T, bool) operator<=(const T&, const rmint<K, MG>&);
}


// --------------------------------------------------------------
// --------------------- Implementation -------------------------

namespace RecInt
{
    // Operator ==
    template <size_t K, size_t MG>
    inline bool operator==(const rmint<K, MG>& a, const rmint<K, MG>& b) {
        return operator==(a.Value, b.Value);
    }
    template <size_t K>
    inline bool operator==(const rmint<K, MGI>& a, const ruint<K>& b) {
        return operator==(a.Value, b);
    }
    template <size_t K>
    inline bool operator==(const rmint<K, MGA>& a, const ruint<K>& b) {
        rmint<K, MGA> br(b);
        return operator==(a.Value, br.Value);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, bool) operator==(const rmint<K, MGI>& a, const T& b) {
        return operator==(a.Value, b);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, bool) operator==(const T& b, const rmint<K, MGI>& a) {
        return operator==(a.Value, b);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, bool) operator==(const rmint<K, MGA>& a, const T& b) {
        rmint<K, MGA> br(b);
        return operator==(a.Value, br.Value);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, bool) operator==(const T& b, const rmint<K, MGA>& a) {
        rmint<K, MGA> br(b);
        return operator==(a.Value, br.Value);
    }

    // Operator !=
    template <size_t K, size_t MG>
    inline bool operator!=(const rmint<K, MG>& a, const rmint<K, MG>& b) {
        return operator!=(a.Value, b.Value);
    }
    template <size_t K>
    inline bool operator!=(const rmint<K, MGI>& a, const ruint<K>& b) {
        return operator!=(a.Value, b);
    }
    template <size_t K>
    inline bool operator!=(const rmint<K, MGA>& a, const ruint<K>& b) {
        rmint<K, MGA> br(b);
        return operator!=(a.Value, br.Value);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, bool) operator!=(const rmint<K, MGI>& a, const T& b) {
        return operator!=(a.Value, b);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, bool) operator!=(const T& b, const rmint<K, MGI>& a) {
        return operator!=(a.Value, b);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, bool) operator!=(const rmint<K, MGA>& a, const T& b) {
        rmint<K, MGA> br(b);
        return operator!=(a.Value, br.Value);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, bool) operator!=(const T& b, const rmint<K, MGA>& a) {
        rmint<K, MGA> br(b);
        return operator!=(a.Value, br.Value);
    }

    // Operator >
    template <size_t K>
    inline bool operator>(const rmint<K, MGI>& a, const rmint<K, MGI>& b) {
        return operator>(a.Value, b.Value);
    }
    template <size_t K>
    inline bool operator>(const rmint<K, MGA>& a, const rmint<K, MGA>& b) {
        rmint<K, MGA> ar(a), br(b);
        return operator>(reduction(ar).Value, reduction(br).Value);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, bool) operator>(const rmint<K, MGI>& a, const T& b) {
        return operator>(a.Value, b);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, bool) operator>(const T& b, const rmint<K, MGI>& a) {
        return operator<(a.Value, b);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, bool) operator>(const rmint<K, MGA>& a, const T& b) {
        rmint<K, MGA> ar(a);
        return operator>(reduction(ar).Value, b);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, bool) operator>(const T& b, const rmint<K, MGA>& a) {
        rmint<K, MGA> ar(a);
        return operator<(reduction(ar).Value, b);
    }

    // Operator >=
    template <size_t K>
    inline bool operator>=(const rmint<K, MGI>& a, const rmint<K, MGI>& b) {
        return operator>=(a.Value, b.Value);
    }
    template <size_t K>
    inline bool operator>=(const rmint<K, MGA>& a, const rmint<K, MGA>& b) {
        rmint<K, MGA> ar(a), br(b);
        return operator>=(reduction(ar).Value, reduction(br).Value);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, bool) operator>=(const rmint<K, MGI>& a, const T& b) {
        return operator>=(a.Value, b);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, bool) operator>=(const T& b, const rmint<K, MGI>& a) {
        return operator<=(a.Value, b);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, bool) operator>=(const rmint<K, MGA>& a, const T& b) {
        rmint<K, MGA> ar(a);
        return operator>=(reduction(ar).Value, b);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, bool) operator>=(const T& b, const rmint<K, MGA>& a) {
        rmint<K, MGA> ar(a);
        return operator<=(reduction(ar).Value, b);
    }

    // Operator <
    template <size_t K>
    inline bool operator<(const rmint<K, MGI>& a, const rmint<K, MGI>& b) {
        return operator<(a.Value, b.Value);
    }
    template <size_t K>
    inline bool operator<(const rmint<K, MGA>& a, const rmint<K, MGA>& b) {
        rmint<K, MGA> ar(a), br(b);
        return operator<(reduction(ar).Value, reduction(br).Value);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, bool) operator<(const rmint<K, MGI>& a, const T& b) {
        return operator<(a.Value, b);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, bool) operator<(const T& b, const rmint<K, MGI>& a) {
        return operator>(a.Value, b);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, bool) operator<(const rmint<K, MGA>& a, const T& b) {
        rmint<K, MGA> ar(a);
        return operator<(reduction(ar).Value, b);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, bool) operator<(const T& b, const rmint<K, MGA>& a) {
        rmint<K, MGA> ar(a);
        return operator>(reduction(ar).Value, b);
    }

    // Operator <=
    template <size_t K>
    inline bool operator<=(const rmint<K, MGI>& a, const rmint<K, MGI>& b) {
        return operator<=(a.Value, b.Value);
    }
    template <size_t K>
    inline bool operator<=(const rmint<K, MGA>& a, const rmint<K, MGA>& b) {
        rmint<K, MGA> ar(a), br(b);
        return operator<=(reduction(ar).Value, reduction(br).Value);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, bool) operator<=(const rmint<K, MGI>& a, const T& b) {
        return operator<=(a.Value, b);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, bool) operator<=(const T& b, const rmint<K, MGI>& a) {
        return operator>=(a.Value, b);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, bool) operator<=(const rmint<K, MGA>& a, const T& b) {
        rmint<K, MGA> ar(a);
        return operator<=(reduction(ar).Value, b);
    }
    template <size_t K, typename T>
    inline __RECINT_IS_ARITH(T, bool) operator<=(const T& b, const rmint<K, MGA>& a) {
        rmint<K, MGA> ar(a);
        return operator>=(reduction(ar).Value, b);
    }
}

#endif

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
