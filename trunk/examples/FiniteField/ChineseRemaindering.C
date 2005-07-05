#include <iostream>
#include <givaro/givintprime.h>
#include <givaro/givmontg32.h>
#include <givaro/givzpz.h>
#include <givaro/givgfq.h>
#include <givaro/givcra.h>    // Chinese Remainder of two elements
#include <givaro/givrns.h>    // Chinese Remainder of an array of elements
#include <givaro/givtimer.h>
#include <givaro/givrandom.h>


typedef GFqDom<long> 		Field1;
typedef ZpzDom<Std16>           Field2; 
typedef ZpzDom<Log16>           Field3; 
typedef ZpzDom<Std32>  		Field4;
typedef ZpzDom<Std64>  		Field5;
typedef ZpzDom<Unsigned32>	Field6; 
typedef Montgomery<Std32>       Field7; 
typedef ZpzDom<Integer>         Field8;

template <typename Field> 
Integer tmain(int argc, char ** argv) {
    typedef RNSsystem<Integer, Field >  CRTSystem;
    typedef typename CRTSystem::domains	Domains;
    typedef typename CRTSystem::array	Elements;
    typedef typename CRTSystem::ring	Ring;
    
    IntPrimeDom ID; 
    GivRandom generator( argc>3 ? atoi(argv[3]):1234 );
    Integer a( generator() >>(argc>2?atoi(argv[2]):17) ), M(1);

    Domains Primes( argc>1 ? atoi(argv[1]):15);   
    Elements Moduli( Primes.size() );

    typename Domains::iterator i = Primes.begin();
    typename Elements::iterator e = Moduli.begin();
    for(; i != Primes.end(); ++i, ++e) {
        *i = Field( ID.nextprimein( a ) );
//         i->random( generator, *e );
        i->init(*e,  generator() );
        M *= a;
    }
    

    CRTSystem CRT( Primes );
    
    Timer tim; tim.clear(); tim.start();
    CRT.RnsToRing( a, Moduli );
    tim.stop();
    Field().write( std::cerr << tim << " using ") << std::endl;
    
    if (Primes.size() < 50) {
        i = Primes.begin();
        e = Moduli.begin();
        for( ; i != Primes.end(); ++i, ++e)
            i->write(std::cout << a << " mod " << i->characteristic() << " = ", *e) << ";" << std::endl;   
    }
    

   Elements Verifs( Primes.size() );    
   CRT.RingToRns( Verifs, a );
   typename Elements::const_iterator v = Verifs.begin();
   i = Primes.begin();
   e = Moduli.begin();
   for( ; i != Primes.end(); ++i, ++e, ++v)
       if (! i->areEqual(*e, *v) ) {
           i->write( std::cerr << "incoherency within ") << std::endl;
           break;
       }        
  

   
   Integer p( generator() >>(argc>2?atoi(argv[2]):17) ), res;
   Field F( ID.nextprimein(p) );
   typename Field::Element el;
   F.init(el, generator() );

   ChineseRemainder<IntPrimeDom, Field> CRA(ID, M, F);
   CRA( res, a, el);

   std::cout << res << " mod " << M << " = " << a << ";"  << std::endl;
   std::cout << res << " mod " << F.characteristic() << " = " << F.convert(a, el) << ";"  << std::endl;
   

   return  res;
}



int main(int argc, char ** argv) {
        // argv[1] : number of primes
        // argv[2] : 2^{32-j} is size of primes
        // argv[3] : seed for generator


    Integer a1 = tmain<Field1>(argc, argv);
    Integer a2 = tmain<Field2>(argc, argv);
    Integer a3 = tmain<Field3>(argc, argv);
    Integer a4 = tmain<Field4>(argc, argv);
    Integer a5 = tmain<Field5>(argc, argv);
    Integer a6 = tmain<Field6>(argc, argv);
    Integer a7 = tmain<Field7>(argc, argv);
    Integer a8 = tmain<Field8>(argc, argv);

    if (a1 != a2) std::cerr << "ERROR a1 != a2" << std::endl;
    if (a3 != a4) std::cerr << "ERROR a3 != a4" << std::endl;
    if (a6 != a5) std::cerr << "ERROR a5 != a6" << std::endl;
    if (a7 != a8) std::cerr << "ERROR a7 != a8" << std::endl;
    if (a1 != a3) std::cerr << "ERROR a1 != a3" << std::endl;
    if (a5 != a7) std::cerr << "ERROR a5 != a7" << std::endl;
    if (a1 != a5) std::cerr << "ERROR a1 != a5" << std::endl;



    


    return 0;
}
