#include <iostream>
#include <givaro/givgfq.h>
#include <givaro/givpoly1.h>

int main(int argc, char ** argv) {

 {
    GFqDom<int> Z13( 13, 1 );  // integers modulo 13

    // Polynomials over Z13, with X as indeterminate
    Poly1Dom< GFqDom<int>, Dense > DP13( Z13, "X" );
    Poly1Dom< GFqDom<int>, Dense>::Element P, Q, R, monomial;
    GFqDom<int>::Element tmp;

    DP13.assign( P, Z13.init(tmp,5) ); // P is degree 0 polynomial : 5 modulo 13
    DP13.init( monomial, Degree(1), -33) ; // -33 X
    DP13.addin( P, monomial ); // P += monomial
    DP13.init( monomial, Degree(2), 12) ; // 12 X^2
    DP13.addin( P, monomial ); // P is now 5-33*X+12*X^2

    // DP13.read( std::cin, P); // would read P as a succession of integers :
                                // deg leadcoeff (lead-1)coeff ... unitcoeff

    DP13.init( Q, Degree(0), 6 ); 
    DP13.init( monomial, Degree(4), 3);
    DP13.addin( Q, monomial) ;
    DP13.init( monomial, Degree(1), 75);
    DP13.addin( Q, monomial) ;
    DP13.init( monomial, Degree(3), 45);
    DP13.subin( Q, monomial) ;
    // Q is now 3*X^4+75*X-45*X^3+6

    DP13.mul ( R, P, Q); // R = P*Q;

    DP13.write( DP13.write( DP13.write(
         std::cout << "(" , P ) << ") * (", Q) << ") = ", R) << std::endl;

    DP13.gcd ( R, P, Q); // 

    DP13.write( DP13.write( DP13.write(
         std::cout << "gcd(", P ) << ",", Q) << ") = ", R) << std::endl;

    DP13.lcm ( R, P, Q); //
    DP13.write( DP13.write( DP13.write(
         std::cout << "lcm(", P ) << ",", Q) << ") = ", R) << std::endl;
    DP13.lcm ( R, Q, P); // 
    DP13.write( DP13.write( DP13.write(
         std::cout << "lcm(", Q ) << ",", P) << ") = ", R) << std::endl;

 }   

 return 0;
}
