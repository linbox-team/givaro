#include <iostream>
#include <givaro/givzpz.h>
#include <givaro/givpower.h>
#include <givaro/givtimer.h>

// Modular exponentiation
// argv[1] : a
// argv[2] : e
// argv[3] : p
// res ---> a^e % p
int main(int argc, char ** argv) {

 {
    ZpzDom<Integer>::Element a, b, p(argv[3]);
    ZpzDom<Integer> Zp( p );   
    unsigned long e = atoi(argv[2]) ;
    Zp.init(a, Integer(argv[1]));
    Zp.init(b);

    Timer tim;tim.clear();tim.start();
    dom_power(b, a, e, Zp);
    tim.stop();

    Zp.write( std::cout, a) << " ^ " << e << " % " << p << " = " << std::flush; 
    Zp.write( std::cerr, b) << std::endl;

    std::cerr << tim << std::endl;
    
 }
    
 


    return 0;
}

