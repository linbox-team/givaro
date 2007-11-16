// ========================================================== //
// Time-stamp: <16 Nov 07 10:26:10 Jean-Guillaume.Dumas@imag.fr>
// Thanks to Dieter Schuster
// ========================================================== //
#include <givaro/givgfq.h>
#include <givaro/givtimer.h>

int main(int argc, char** argv)
{
    GFqDom<long> GF128(2, 7);
    GFqDom<long>::Element b, c, gen;
    GF128.init(b, 5UL);
    GF128.init(c, 3);
    GF128.write(std::cout, b) << std::endl;
    GF128.write(std::cout, c) << std::endl;

    GFqDom<long>::Element f,g,h;
    GFqDom<long> F2(2);
    Poly1Dom< GFqDom<long>, Dense> Pol2(F2);
    Poly1Dom< GFqDom<long>, Dense>::Element P, Q, R;
    Pol2.init(P,Degree(1));
    F2.init(P[0],1);
    F2.init(P[1],1);
    GF128.init(f, P);
    GF128.write(std::cout << "2-adic representation of 1+X is: ", f);
    std::cout << " ... while its internal representation is: " << f << std::endl;


    gen = 1;
    GF128.write(std::cout << "Indeed, we are in ") <<std::endl;
    GF128.write(std::cout << "In this field, the generator used is (in 2-adic): ", gen)
        << " whose internal representation is " << gen << std::endl;


    Poly1PadicDom< GFqDom<long>, Dense > Padic2(Pol2);
        // 

    std::cout << "Irreducible: " << GF128.irreducible() << std::endl;
    std::cout << "Now, the irreducible polynomial used to build this field is (in 2-adic)" << std::endl;


    GF128.init(g, Padic2.radix( Q, Integer(5) ));
    GF128.write(std::cout << "2-adic representation of 1+X^2 is: ", g) << std::endl;

    GF128.init(h);
    GF128.add(h, g, f);
    GF128.write(std::cout << "2-adic representation of X+X^2 is: ", h) << std::endl;   

    GF128.mul(h, g, f);
    GF128.write(std::cout << "2-adic representation of 1+X+X^2+X^3 is: ", h) << std::endl;    

    GF128.div(h, g, f);
    GF128.write(std::cout << "2-adic representation of 1+X is: ", h) << std::endl;
    
    return 0;
    
}
