// Copyright(c)'2020 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

#include <iostream>
#include <givaro/gfq.h>
#include <givaro/givmatrix.h>

typedef Givaro::GFqDom<int64_t> Field;
typedef Field::Element Element;


typedef Givaro::VectorDom<Field,Givaro::Dense> VDomain;
typedef VDomain::Element Vector;
typedef Givaro::MatrixDom<Field,Givaro::Dense> MDomain;
typedef MDomain::Element Matrix;



int main() {
    
    Field F4(2,2);
    VDomain V4(F4);
    MDomain M4(F4);
    
    const size_t n(3);

    Matrix A; M4.init(A, n, n);
    for(uint64_t i=0; i<n*n; ++i)
        F4.init(A(i/n,i%n), i); 
    M4.write(std::cout << "A:=", A) << ';' << std::endl;

    Vector x; V4.init(x, n);
    V4.write(std::cout << "z:=", x) << ';' << std::endl;
    for(uint64_t i=0; i<n; ++i)
        F4.init(x[i], i);
    V4.write(std::cout << "x:=", x) << ';' << std::endl;

    Vector b; V4.init(b, n);
 
    M4.mul(b, A, V4, x);
    
    V4.write(std::cout << "b:=", b) << ';' << std::endl;

    return 0;
    
}

