// =================================================================== //
// Givaro : Prime numbers
//              Modular powering,
//              Fermat numbers,
//              Primality tests, Factorization :
//                      (There are parameters to fix)
// Time-stamp: <08 Jun 04 17:32:35 Jean-Guillaume.Dumas@imag.fr> 
// =================================================================== //
#ifndef _GIVARO_INTEGERS_PRIME_H_
#define _GIVARO_INTEGERS_PRIME_H_

#ifndef _GIVARO_ISPRIMETESTS_
#define _GIVARO_ISPRIMETESTS_ 5
#endif

#include "givaro/givinteger.h"


// =================================================================== //
// Fermat numbers
// =================================================================== //
class FermatDom : public IntegerDom {
public:
    FermatDom() : IntegerDom() {}
    Rep& fermat (Rep&, const long)  const ;
    int pepin (const long) const ;
};


// =================================================================== //
// Primality tests and factorization algorithms
// =================================================================== //

// Those macros are parameters to fix

// primes known
// first array
#define LOGMAX 3512
#define TABMAX 32768
// second array
#define LOGMAX2 3031
#define TABMAX2 65536
// Bounds between big and small
#define BOUNDARY_isprime TABMAX    
#define BOUNDARY_2_isprime TABMAX2    

// =================================================================== //
// Primality tests 
// =================================================================== //
class IntPrimeDom : public IntegerDom {
public:
    IntPrimeDom() :  IntegerDom() {}

    int isprime(const Rep& n, int r=_GIVARO_ISPRIMETESTS_) const 
        {
/*
  return probab_prime(n);
*/
//             return ((n)<BOUNDARY_isprime ?  isprime_Tabule(n) : 
//                     (n)<BOUNDARY_2_isprime ? isprime_Tabule2(n) : 
//                     probab_prime(n));
            long l;
            return (islt(n,BOUNDARY_isprime) ?  isprime_Tabule(access(l,n)): 
                    islt(n,BOUNDARY_2_isprime) ? isprime_Tabule2(access(l,n)): 
                    local_prime(n,r));
        }

        // if p is a prime power, p = r^return
        // else return is 0 and r is undefined
    unsigned int isprimepower(Rep&, const Rep&) const ;

    template<class RandIter>
    unsigned int Miller(RandIter& g, const Rep& n=_GIVARO_ISPRIMETESTS_) const  ;

    template<class RandIter>
    Rep& test_Lehmann(RandIter& g, Rep&, const Rep& n=_GIVARO_ISPRIMETESTS_) const  ;

    template<class RandIter>
    int Lehmann(RandIter& g, const Rep& n=_GIVARO_ISPRIMETESTS_)  const ;

    int isprime_Tabule(const int n) const ;
    int isprime_Tabule2(const int n) const ;

    Rep& nextprime(Rep&, const Rep&, int r=_GIVARO_ISPRIMETESTS_) const ;
    Rep& prevprime(Rep&, const Rep&, int r=_GIVARO_ISPRIMETESTS_) const ;
    Rep& nextprimein(Rep&, int r=_GIVARO_ISPRIMETESTS_) const ;
    Rep& prevprimein(Rep&, int r=_GIVARO_ISPRIMETESTS_) const ;


// Using Integer
    int local_prime(const Rep& n, int r=_GIVARO_ISPRIMETESTS_) const { return probab_prime(n,r); }
    int& access(int& r, const Rep& a) const { return r=Integer2long(a); }
    long& access(long& r, const Rep& a) const { return r=Integer2long(a); }
    double& access(double& r, const Rep& a) const { return r=Integer2double(a); }

private:
    static int IP[LOGMAX+5];  // -- table for Tabule
    static const int * TP;    // -- shifted table 
    static int IP2[LOGMAX2+5]; // -- table for Tabule2
    static const int * TP2;    // -- shifted table 
/*
  static int Tabule2(const Integer& p) ;
  static int Tabule(const Integer& p) ;
  static int _memTab2[LOGMAX2+5];   // -- table for Tabule2
  static const int* _Tab2; // -- shifted _memTabule2
  static int _memTab[];    // -- table for Tabule
  static const int* _Tab;  // -- shifted _memTabule
*/
};

#include "givaro/givintprime.inl"
#endif
