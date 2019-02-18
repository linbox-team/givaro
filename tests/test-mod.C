// Copyright(c)'2010 by The Givaro group
// This file is part of Givaro.
// written by BB
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.

#include <iostream>

#include <givaro/givinteger.h>

using namespace Givaro;

#define AREEQUALVALUES(a,b)\
if ( (a) != (b) ) { \
    std::cout << "*** ERROR line " << __LINE__ << std::endl; \
    std::cout << a << "!=" << b << std::endl; \
    return -1 ; \
}

#define NB_ITERS 40

#include <cmath>

template<class U> inline bool IsNeg(const U p) { return (p<0); }
template<> inline bool IsNeg<uint64_t>(const uint64_t p) { return false; }

template<class T, class U>
int64_t ref_modulo(const T  m, const U p)
{
    Integer M(m);
    Integer P(p);
    Integer R ;
    mpz_mod(R.get_mpz() ,M.get_mpz() ,P.get_mpz());
    // if (R<0) R+=(P<0)?(-P):(P);
#if 0
    Integer Q = M/P ;

    if (Q*P+R != M){
        bool pass = false ;
        std::cout << Q << " * "  << P << " + " << R << " != " << M << std::endl;
        std::cout << "ref is no good..." << std::endl ;
        if (Q<0){
            if ((Q-1)*P+R != M)
                pass = true ;
        } else {
            if((Q+1)*P+R != M) pass = true ;
        }
    }
    if (R < 0 || R >= GMP__ABS(P)){
        std::cout << R << " _ "  << P << std::endl;
        std::cout << "ref is no good..." << std::endl ;
    }
#endif
    return R ;

}

template<class T, class U>
int64_t ref_modulobis(const T  m, const U p)
{
    Integer M(m);
    Integer P(p);
    Integer R ;
    mpz_tdiv_r(R.get_mpz() ,M.get_mpz() ,P.get_mpz());
    return R ;

}


// m and p should not be too large
template< class T, class U>
int test1( const T m, const U p)
{
    // double pi (p);
    int64_t r =  ref_modulo(m, p);

    const Integer M(m);
    const Integer P(p);
    Integer R ;

    Integer::mod(R,M,P);
    AREEQUALVALUES(r,R);

    Integer::mod(R,M,p);
    AREEQUALVALUES(r,R);


    Integer R1 = M ;
    Integer::modin(R1,P);
    AREEQUALVALUES(R,R1);

    R1 = M ;
    Integer::modin(R1,p);
    AREEQUALVALUES(R,R1);


    return 0;
}

// m and p should not be too large
template< class T, class U>
int test1bis( const T m, const U p)
{
    double pi = (double) p;
    int64_t r =  ref_modulobis(m, p);
    // if (r<0)  r += (IsNeg(p) )?(-p):(p); // r est dans [[0,p-1]]
    const Integer M(m);
    const Integer P(p);
    Integer R =r ;

    Integer R1 = M%p ;
    AREEQUALVALUES(R,R1);

    Integer R2 = M%P ;
    AREEQUALVALUES(R,R2);

    Integer R3 = M%pi ;
    AREEQUALVALUES(R,R3);


    R2 = M ;
    R2 %= P ;
    AREEQUALVALUES(R,R2);

    return 0;
}

int test2(Integer & M, Integer & P)
{
    Integer RR ;
    Integer MM = M ;
    Integer PP = P ;
    // reference
    mpz_mod(RR.get_mpz() ,MM.get_mpz() ,PP.get_mpz());

    Integer R = 0 ;

    //!@todo existe pas !
    //R = Integer:: imod(M,P) ;
    Integer:: mod(R,M,P) ;
    AREEQUALVALUES(RR,R);

    R = M ;
    Integer::modin(R,P);
    AREEQUALVALUES(RR,R);

    return 0;
}

int test2bis(Integer & M, Integer & P)
{
    Integer RR ;
    Integer MM = M ;
    Integer PP = P ;
    // reference
    mpz_tdiv_r(RR.get_mpz() ,MM.get_mpz() ,PP.get_mpz());

    Integer R = 0 ;

    R = M ;
    R %= P;
    AREEQUALVALUES(RR,R);

    // R = M ;
    R = M%P;
    AREEQUALVALUES(RR,R);

    return 0;
}

template< class T, class U>
int test3( const T m, const U p)
{
    int pi = int(p);
    int64_t q = (int64_t)(m / (T)p);
    const Integer M(m);
    const Integer P(p);
    Integer Q ;

    Integer::div(Q,M,P);
    AREEQUALVALUES(q,Q);

    Integer Q1 = M/p ;
    AREEQUALVALUES(Q,Q1);

    Integer Q2 = M/P ;
    AREEQUALVALUES(Q,Q2);

    Integer Q3 = M/pi ;
    AREEQUALVALUES(Q,Q3);

    Q1 = M ;
    Integer::divin(Q1,P);
    AREEQUALVALUES(Q,Q1);


    Q2 = M ;
    Q2 /= P ;
    AREEQUALVALUES(Q,Q2);

    return 0;
}

#include <cassert>

int main()
{
#if (SIZEOF_LONG==8)
    int64_t m = 1253345363665346363;
#else
    int64_t m = 1665346363;
#endif

    int64_t p = 78678675;
    uint64_t M((uint64_t)m);
    uint64_t P((uint64_t)p);

    Integer mOne(-1);
    // CONDITION: mpz_tdiv_ui does NOT consider the sign of gmp_rep
    assert(mpz_tdiv_ui( (mpz_ptr)&mOne, 3) == 1);

    // Integer r ;
    // Integer a = -6 ;
    // Integer b = 7 ;
    // mpz_tdiv_r( (mpz_ptr)&r, (mpz_ptr)&a,(mpz_ptr)&b) ;
    // std::cout << r << std::endl;
#if 0
    for (long i = -11 ; i < 11 ; ++i){
        long j = 2 ;
        std::cout << i/j << '=';
        Integer I = i;
        Integer J = j ;
        Integer::divin(I,j);
        std::cout << I << '=' ;
        I = i ;
        std::cout << I/J << '=' ;
        I /= J ;
        std::cout << I << " #  " ;

        I = i;
        j = -j ;
        J = j ;
        std::cout << i/j << '=';
        I = i;
        Integer::divin(I,j);
        I = i ;
        std::cout << I/J << "=" ;
        I /= J ;
        std::cout << I << "| @" ;




        std::cout << i << std::endl ;


    }
#endif

    int rez = 0;

    rez =  test1(m,p);   if (rez) return 1 ;
    rez =  test1(-m,p);  if (rez) return 2 ;
    rez =  test1(m,-p);  if (rez) return 3 ;
    rez =  test1(-m,-p); if (rez) return 4 ;

    rez =  test1(m,P);   if (rez) return 5 ;
    rez =  test1(m,-P);  if (rez) return 6 ; //
    rez =  test1(-m,P);  if (rez) return 7 ;
    // rez =  test1(-m,-P); if (rez) return 8 ;

    rez =  test1(M,P);   if (rez) return 9 ;
    rez =  test1(M,-P);  if (rez) return 10 ; //
    rez =  test1(-M,P);  if (rez) return 11 ; //
    // rez =  test1(-M,-P); if (rez) return 12 ;

    rez =  test1(-M,p);  if (rez) return 13 ; //
    rez =  test1(M,-p);  if (rez) return 14 ;
    rez =  test1(M,p);   if (rez) return 15 ;
    rez =  test1(-M,-p); if (rez) return 16 ; //

    rez =  test1(m,m);   if (rez) return 16 ;
    rez =  test1(-m,m);  if (rez) return 18 ;
    rez =  test1(m,-m);  if (rez) return 19 ;
    rez =  test1(-m,-m); if (rez) return 20 ;

    rez =  test1(P,P);   if (rez) return 21 ;
    rez =  test1(P,-P);  if (rez) return 22 ; //
    rez =  test1(-P,P);  if (rez) return 23 ; //
    // rez =  test1(-P,-P); if (rez) return 24 ;


    ////////////////////////////////////////
    rez =  test1bis(m,p);   if (rez) return 1 ;
    rez =  test1bis(-m,p);  if (rez) return 2 ;
    rez =  test1bis(m,-p);  if (rez) return 3 ;
    rez =  test1bis(-m,-p); if (rez) return 4 ;

    rez =  test1bis(m,P);   if (rez) return 5 ;
    rez =  test1bis(m,-P);  if (rez) return 6 ; //
    rez =  test1bis(-m,P);  if (rez) return 7 ;
    // rez =  test1bis(-m,-P); if (rez) return 8 ;

    rez =  test1bis(M,P);   if (rez) return 9 ;
    rez =  test1bis(M,-P);  if (rez) return 10 ; //
    rez =  test1bis(-M,P);  if (rez) return 11 ; //
    // rez =  test1bis(-M,-P); if (rez) return 12 ;

    rez =  test1bis(-M,p);  if (rez) return 13 ; //
    rez =  test1bis(M,-p);  if (rez) return 14 ;
    rez =  test1bis(M,p);   if (rez) return 15 ;
    rez =  test1bis(-M,-p); if (rez) return 16 ; //

    rez =  test1bis(m,m);   if (rez) return 16 ;
    rez =  test1bis(-m,m);  if (rez) return 18 ;
    rez =  test1bis(m,-m);  if (rez) return 19 ;
    rez =  test1bis(-m,-m); if (rez) return 20 ;

    rez =  test1bis(P,P);   if (rez) return 21 ;
    rez =  test1bis(P,-P);  if (rez) return 22 ; //
    rez =  test1bis(-P,P);  if (rez) return 23 ; //
    // rez =  test1bis(-P,-P); if (rez) return 24 ;


    ////////////////////////////////////////


    rez =  test3(p,m);   if (rez) return 11 ;
    rez =  test3(-p,m);  if (rez) return 12 ;
    rez =  test3(p,-m);  if (rez) return 13 ;
    rez =  test3(-p,-m); if (rez) return 14 ;

    rez =  test3(m,p);   if (rez) return 11 ;
    rez =  test3(-m,p);  if (rez) return 12 ;
    rez =  test3(m,-p);  if (rez) return 13 ;
    rez =  test3(-m,-p); if (rez) return 14 ;

    rez =  test3(m,P);   if (rez) return 15 ;
    // rez =  test3(m,-P);  if (rez) return 16 ;
    // rez =  test3(-m,P);  if (rez) return 17 ;
    // rez =  test3(-m,-P); if (rez) return 18 ;

    rez =  test3(M,P);   if (rez) return 19 ;
    // rez =  test3(M,-P);  if (rez) return 110 ;
    rez =  test3(-M,P);  if (rez) return 111 ;
    // rez =  test3(-M,-P); if (rez) return 112 ;

    rez =  test3(-M,p);  if (rez) return 113 ;
    // rez =  test3(M,-p);  if (rez) return 114 ;
    rez =  test3(M,p);   if (rez) return 115 ;
    // rez =  test3(-M,-p); if (rez) return 116 ;

    rez =  test3(m,m);   if (rez) return 116 ;
    rez =  test3(-m,m);  if (rez) return 118 ;
    rez =  test3(m,-m);  if (rez) return 119 ;
    rez =  test3(-m,-m); if (rez) return 120 ;

    rez =  test3(P,P);   if (rez) return 121 ;
    // rez =  test3(P,-P);  if (rez) return 122 ;
    rez =  test3(-P,P);  if (rez) return 123 ;
    // rez =  test3(-P,-P); if (rez) return 124 ;



    for (unsigned i = 0 ; i < NB_ITERS ; ++i)
    {
        Integer Mint = Integer::random_between(680,700);
        Integer Pint = Integer::random_between(134,198);
        rez = test2(Mint,Pint); if (rez) return 7 ;
        Integer::negin(Mint);
        rez = test2(Mint,Pint); if (rez) return 8 ;
        Integer::negin(Pint);
        rez = test2(Mint,Pint); if (rez) return 9 ;
        Integer::negin(Mint);
        rez = test2(Mint,Pint); if (rez) return 10 ;
    }
    for (unsigned i = 0 ; i < NB_ITERS ; ++i)
    {
        Integer Mint = Integer::random_between(680,700);
        Integer Pint = Integer::random_between(134,198);
        rez = test2bis(Mint,Pint); if (rez) return 7 ;
        Integer::negin(Mint);
        rez = test2bis(Mint,Pint); if (rez) return 8 ;
        Integer::negin(Pint);
        rez = test2bis(Mint,Pint); if (rez) return 9 ;
        Integer::negin(Mint);
        rez = test2bis(Mint,Pint); if (rez) return 10 ;
    }

    //std::cout << "ok6" << std::endl;

    rez =  test1(0,p);    if (rez) return 21 ;
    rez =  test1(0,-p);   if (rez) return 22 ;
    rez =  test1bis(0,p);    if (rez) return 21 ;
    rez =  test1bis(0,-p);   if (rez) return 22 ;


    rez =  test3(0,p);    if (rez) return 23 ;
    rez =  test3(0,-p);   if (rez) return 24 ;

    //std::cout << "ok7" << std::endl;



    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
