#include <iostream>
#include <stdlib.h>
#include <givaro/givinteger.h>
#include <givaro/givintprime.h>
#include <givaro/givtimer.h>


int main(int argc, char** argv)
{
    Integer m, p;
    if (argc > 1) m = Integer(argv[1]);
    IntPrimeDom IP;
    
    {
        Timer tim; tim.clear(); tim.start();
        int a = isperfectpower(m);
        tim.stop();
        std::cout << a << std::endl;
        std::cerr << tim << std::endl;
    }
    {
        Timer tim; tim.clear(); tim.start();
        int a = IP.isprimepower(p, m);
        tim.stop();
        if (a) std::cout << "is " << p << "^" << a << std::endl;
        else   std::cout << a << std::endl;

        std::cerr << tim << std::endl;
    }
    return 0;
}




// #include <givaro/givintfactor.h>
// #include "container_ostream.C"

// bool factor_simple(std::vector<Integer>& Lf, std::vector<unsigned long>& Lo, const Integer& n, const unsigned long threshold = 0) {

//         // n = * Lf[i] ^ Lo[i]
//         // But Lf[i] might not be prime (cf. factor probability)
//     IntFactorDom<> IP;
//     GivRandom generator;
//     bool factocomplete = true;
    
//     Integer nn,g,r,u,k;
//     if (n<0) IP.neg(nn,n); else nn=n;
//     unsigned long c,nb=0;
//     while(nn > 1) { 
//            IP.Pollard(generator, g, nn, threshold);
//            if (g == 1) {
//              factocomplete = false;
//              g = nn;
//            } else if (! IP.isprimen(g)) {
//                 Integer gg(g);
//                 IP.factor(g, gg);
//            }
//            Lf.push_back(g);
//            c=0; r=IP.zero;
//            Integer::divexact(u, nn,g);
//             while(r == 0) {
//                   nn.copy(u); 
//                   Integer::divmod( u, r, nn,g );
//                   c++;
//         }
//         Lo.push_back( c );
//         nb++;
//     }
  
   
//    return factocomplete;
   

// } 

// double nthroot(const unsigned int k, const Integer& m, double l2) {
//     double dm=(double)m, xn, xnpu, xnkmu, dk=(double)(k), dkmu = dk-1.0;
// //     xnpu = (double) (Integer(1)<<((m.size()<<5)/k));
// //     std::cerr << "ms: " << m.size() << std::endl;
// //     std::cerr << "ms: " << ((m.size()<<5)/k) << std::endl;
// //     std::cerr << "ms: " << (Integer(1)<<((m.size()<<5)/k)) << std::endl;
// //     std::cerr << "x0: " << xnpu << std::endl;
//     xnpu = pow(2.0,l2/(dk));
// //     std::cerr << "x0: " << xnpu << std::endl;
//     do {
//         xn = xnpu;
//         xnkmu = pow(xn, dkmu);
//         xnpu = xnkmu;
//         xnpu *= xn;
//         xnpu -= dm;
//         xnpu /= k;
//         xnpu /= xnkmu;
//         xnpu = xn - xnpu;
//     } while (xn - xnpu > 0.5);
//     return xnpu;
// }

    
// int floatingipp(const Integer& m) {
//     double l2 = log2(m); unsigned int il2 = (unsigned int)(l2);
// //     std::cerr << "fipp: " << l2 << std::endl;
//     for (unsigned int i = 2; i<= il2; ++i) {
//         double t = nthroot(i,m, l2);
//         Integer it(t), x = pow( it, i);
// //         std::cerr << "i: " << i << ", t: "  << t << ", x:" << x << std::endl;
//         if (x == m) return i;
//         x = pow( ++it, i);
//         if (x == m) return i;
//     }
//     return 0;
// }


// int puissparfaite(const Integer& m) {
//     std::vector<Integer> Lf;
//     std::vector<unsigned long> Lo;
//     factor_simple(Lf, Lo, m, 1000);
//     if (Lo.size() > 1) {
// //         std::cerr << Lf << " ^ " << Lo << std::endl;
//         unsigned int e = Lo.front();
//         for(std::vector<unsigned long>::const_iterator it=Lo.begin(); ++it != Lo.end(); )
//             if (e != *it) return 0;
//         return e;
//     } else {
//         if (Lo.front() == 1)
//             return floatingipp(m);
//         else 
//             return Lo.front();
//     }
// }

    

    
