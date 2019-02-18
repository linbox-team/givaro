// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/gmp++/gmp++_int_io.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama, T. Gautier
// $Id: gmp++_int_io.C,v 1.7 2009-09-17 14:28:22 jgdumas Exp $
// ==========================================================================
// Description:
/** @file gmp++/gmp++_int_io.C
 * ioing stuff.
 */

#ifndef __GIVARO_gmpxx_gmpxx_int_io_C
#define __GIVARO_gmpxx_gmpxx_int_io_C

#include <iostream>
#include <stdlib.h>
#ifndef __GIVARO_INLINE_ALL
#include "gmp++/gmp++.h"
#endif
#ifndef __GIVARO_GMP_NO_CXX
#include <sstream>
#endif


namespace Givaro {

    // Sortie nonsignee : 321321 meme si n = -321321, par exemple
    std::ostream& absOutput(std::ostream &o, const Integer&n)
    {
        int32_t base = 10;

        uint64_t strSize = mpz_sizeinbase((mpz_srcptr)&(n.gmp_rep), base) + 2;
        char *str = ::new char[strSize];
        mpz_get_str(str, base, (mpz_srcptr)&(n.gmp_rep));
        if (sign(n) < 0) {
            char *str1 = &(str[1]) ;
            o << str1;
        }
        else o << str;
        delete [] str ;
        return o;
    }

    // Sortie signee : +321321 ou -321321, par exemple
    std::ostream& Integer::print(std::ostream &o) const
    {
        // For some reason, ekopath does an undefined reference to the GMP-streams - A. Breust 2015-01-21
#if defined(__GIVARO_GMP_NO_CXX) || defined(__PATHCC__)
        int32_t base = 10;
        uint64_t strSize = mpz_sizeinbase((mpz_srcptr)&(gmp_rep), base) + 2;
        char *str = new char[strSize];
        mpz_get_str(str, base, (mpz_srcptr)&(gmp_rep));
        o << str;
        delete [] str ;
        return o;

#else
        return o << (mpz_srcptr)&gmp_rep;
#endif
    }

    // Entree au format de la sortie
    std::istream& operator>> (std::istream& inp, Integer& a)
    {
#if defined(__GIVARO_GMP_NO_CXX) || defined(__PATHCC__)
        static int64_t base[10] = {
            10,
            100,
            1000,
            10000,
            100000,
            1000000,
            10000000,
            100000000,
            1000000000
        } ;
        if (!inp) return inp ;
        // eat white
        inp >> std::ws  ;

        // Base : 10^9, we read by packet of length 9
        // the char.
        char Tmp[10] ;
        int32_t counter = 0 ;

        // Set the returned integer
        a = 0 ;
        char ch ;
        int32_t sign = 1 ;

        // find a sign:
        inp.get(ch) ;
        if ((ch != '+') && (ch != '-') && !((ch >= '0') && (ch <= '9')))
        {
            std::cerr << "Bad integer format: found: "<< ch ;
            std::cerr << ", in place of '+' '-' or a digit"<< std::endl ;
            return inp ;
        }
        switch (ch) {
        case '+' : break ;
        case '-' : sign = -1 ; break ;
        default  : inp.putback(ch) ; break ;
        }
        // eat white
        inp >> std::ws  ;

        int32_t noend = 1 ;
        while (noend)
        {
            counter = 0 ;

            // Read 9 digits or less
            while ((noend) && (counter < 9)) {
                inp.get(ch) ;
                if (inp.eof()) { noend = 0 ; }
                else if ((ch >= '0') && (ch <= '9')) Tmp[counter++] = ch ;
                else { noend = 0 ;  inp.putback(ch) ; }
            }
            if (counter >0) {
                int64_t l ;
                Tmp[counter] = '\0' ; // terminate the string
                l = atol(Tmp) ;
                a = a * base[counter-1] + l ;
            }
        }
        if (sign == -1) a = -a ;
        return inp ;

#else
        return inp >> (mpz_ptr)a.get_mpz();
#endif
    }

    //-------------------------------------------------inline >> & << operators
    std::ostream& operator<< (std::ostream& o, const Integer& a)
    {
        return a.print(o);
    }

}

#endif // __GIVARO_gmpxx_gmpxx_int_io_C

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
