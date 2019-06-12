#include <iostream>
#include <cstdlib>
#include <omp.h>
#include <givaro/givtimer.h>
#include <givaro/givinteger.h>
using namespace Givaro;

    inline Givaro::Integer
    InfNorm (const size_t M, const size_t N, const Givaro::Integer* A, const size_t lda){
        Givaro::Integer max = 0;
        size_t log=0;
        for (size_t i=0; i<M; ++i)
            for (size_t j=0; j<N; ++j){
                const Givaro::Integer & x(A[i*lda+j]);
                if ((x.bitsize() >= log) && (abs(x) > max)){
                    max = abs(x);
                    log = x.bitsize();
                }
            }
        return max;
    }

    inline Givaro::Integer
    InfNorm2 (const size_t M, const size_t N, const Givaro::Integer* A, const size_t lda){
        Givaro::Integer max = 0;
        auto end=&A[M*N]; 
        for(auto * iter=A; iter != end; ++iter) 
            if (absCompare(*iter,max)>0) {
                max = *iter;
            }
        return abs(max);
    }

    inline Givaro::Integer
    InfNorp (const size_t M, const size_t N, const Givaro::Integer* A, const size_t lda){
        std::vector<Givaro::Integer> max(M,Givaro::Integer(0));

#pragma omp parallel for
        for (size_t i=0; i<M; ++i) {
            for (size_t j=0; j<N; ++j){
                const Givaro::Integer & x(A[i*lda+j]);
                if (absCompare(x,max[i])>0) {
                    max[i] = x;
                }
            }
            max[i] = abs(max[i]);
        }
        return *std::max_element(max.begin(),max.end());
    }


int main(int argc, char ** argv)
{
	size_t m(argc > 1? atoi(argv[1]) : 100);
	size_t n(argc > 2? atoi(argv[2]) : 100);
	size_t b(argc > 3? atoi(argv[3]) : 100);

    Givaro::Integer::seeding();
    Givaro::GivRandom generator;
    IntegerDom IP;

    Givaro::Timer ini,tim,ori,par;
    ini.clear(); ini.start();
	Givaro::Integer * A = new Givaro::Integer[m*n];
    const size_t mn(m*n);

    ini.stop();
    std::cout
        << "Alloc : " << ini << ' ' << m << 'x' << n << ':' << b 
        << std::endl;
    
    // Main loop
    ini.clear(); ini.start();
    for (size_t l = 0; l < mn; ++l) {
        IP.random(generator, A[l], b);
        if (A[l] & 1U) IP.negin(A[l]);
//         A[l] = Integer::random<false>(b); // signed
    }
    ini.stop();
    std::cout
        << "Random: " << ini << ' ' << m << 'x' << n << ':' << b 
        << std::endl;
    
    // Main loop
    ori.clear(); ori.start();
    Integer max=InfNorm(m,n,A,n);
    ori.stop();

    std::cout 
        << " InfNorm: " << ori
        << " Mflops: " << 1./ori.usertime()/1000.0/1000.0 
        << ' ' << max 
        << std::endl ;

    tim.clear(); tim.start();
    Integer m2=InfNorm2(m,n,A,n);
    tim.stop();

    std::cout 
        << " InfNor2: " << tim
        << " Mflops: " << 1./tim.usertime()/1000.0/1000.0 
        << ' ' << m2 
        << std::endl ;

    par.clear(); par.start();
    Integer mpa=InfNorp(m,n,A,n);
    par.stop();

    if ( (mpa != m2) || (mpa != max) || (m2 != max) )
        std::cerr << "ERROR: inconsistency" << std::endl;

    std::cout 
        << " InfNorp: " << par
        << " Mflops: " << 1./par.usertime()/1000.0/1000.0 
        << ' ' << mpa
        << std::endl ;
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
