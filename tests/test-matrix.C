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
    
    Field F4(2,2); // modulo 2 and modulo 1+X+X^2
    VDomain V4(F4);
    MDomain M4(F4);
    
    const size_t n(3);

        // [[0,1,X],[1+X,0,1],[X,1+X,0]]
    Matrix A; M4.init(A, n, n);
    for(uint64_t i=0; i<n*n; ++i)
        F4.init(A(i/n,i%n), i);

        // [1,X,1+X]
    Vector x; V4.init(x, n);
    for(uint64_t i=0; i<n; ++i)
        F4.init(x[i], i+1);

    Vector b; V4.init(b, n);
    M4.mul(b, A, V4, x);
    
        // [1+X,0,1+X]
    Vector c; V4.init(c, n);
    F4.init(c[0],3);
    F4.init(c[1],0);
    F4.init(c[2],3);
    

    bool success = V4.areEqual(b,c);

    if (! success) {
        std::cerr << "Error: " << std::endl;
        M4.write(std::cerr << "A:=", A) << ';' << std::endl;
        V4.write(std::cerr << "x:=", x) << ';' << std::endl;
        V4.write(std::cerr << "c:=", c) << ';' << std::endl;
        V4.write(std::cerr << "But b:=", b) << ';' << std::endl;
    }

    return (! success);

}

