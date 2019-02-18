// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/bstruct/givbits.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Author: T. Gautier
// $Id: givbits.C,v 1.2 2009-09-17 14:28:22 jgdumas Exp $
// ==========================================================================
// Description:
// - field of n bits, for any n.
// The bits field is and array of long: a_(k-1),...,a_0.
// each entry has sizeof(long) bits. Give an number i, to i-th
// bit in the representation is the r-th bit carried by entry a_q,
// where q = i quo sizeof(long) and r = i rem sizeof(long) (counting from 0).

#include <iostream>
#include <stddef.h>
#include "givaro/givbits.h"
#include "givaro/givmodule.h"

#define GIV_SIZE_LONG 32
#if (GIV_SIZE_LONG == 32)
#  define SIZE_IN_BYTE		4  					// size of the base
#  define SIZE_IN_BIT		32  					// size of the base
#  define QUO(x) 		( (x >>5)  ) 				// quotient by 32
#  define REM(x)		( (x & 0x1f) )				// remainder by 32
#  define MAX_WORD(x) 		( QUO(x) + (0!=REM(x) ? 1 : 0) )	// ceil( QUO(x) )
#else
#  define SIZE_IN_BYTE 		sizeof(Bits::base)
#  define SIZE_IN_BIT 		(sizeof(Bits::base)*8)
#  define QUO(x)		(x / SIZE_IN_BIT)
#  define REM(x)		(x % SIZE_IN_BIT)
#  define MAX_WORD(x) 	(QUO(x) + (0!=REM(x) ? 1 : 0))
#endif

namespace Givaro {

    // -- Table of 2^i = Table2pow[i], for i =0..sizeof(base)
    static Bits::base* Table2pow;


    Bits::Bits () :
        rep(0)
    { }

    Bits::Bits (const size_t n)
    {
        size_t len = MAX_WORD(n) ;
        rep.allocate( len );
    }

    Bits::Bits( const Bits& B ) :
        rep(B.rep, givWithCopy())
    {}

    Bits::Bits( const Bits::Rep& B ) :
        rep(B, givWithCopy())
    {}

    Bits::~Bits()
    { rep.destroy(); }

    // And
    const Bits Bits::operator& (const Bits& A) const
    {
        GIVARO_ASSERT( rep.size() == A.rep.size(), "[Bits]: invalide size in 'and'");
        size_t len = rep.size();
        Rep res( len );
        for (int i=0; i< int(len); ++i)
            res[i] = rep[i] & A.rep[i];
        return Bits(res);
    }

    Bits& Bits::andin( const Bits& A, const Bits B)
    {
        GIVARO_ASSERT( B.rep.size() == A.rep.size(), "[Bits]: invalide size in andin");
        GIVARO_ASSERT( rep.size() == A.rep.size(), "[Bits]: invalide size in andin");
        size_t len = rep.size();
        for (int i=0; i<int(len); ++i)
            rep[i] = A.rep[i] & B.rep[i];
        return *this;
    }

    Bits& Bits::operator&= (const Bits& A)
    {
        GIVARO_ASSERT( rep.size() == A.rep.size(), "[Bits]: invalide size in &=");
        size_t len = rep.size();
        for (int i=0; i< int(len); ++i)
            rep[i] &= A.rep[i];
        return *this;
    }

    // Or
    const Bits Bits::operator| (const Bits& A) const
    {
        GIVARO_ASSERT( rep.size() == A.rep.size(), "[Bits]: invalide size in or");
        size_t len = rep.size();
        Rep res( len );
        for (int i=0; i<int(len); ++i)
            res[i] = rep[i] | A.rep[i];
        return Bits(res);
    }

    Bits& Bits::orin( const Bits& A, const Bits B)
    {
        GIVARO_ASSERT( B.rep.size() == A.rep.size(), "[Bits]: invalide size in orin");
        GIVARO_ASSERT( rep.size() == A.rep.size(), "[Bits]: invalide size in orin");
        size_t len = rep.size();
        for (int i=0; i<int(len); ++i)
            rep[i] = A.rep[i] | B.rep[i];
        return *this;
    }

    Bits& Bits::operator|=( const Bits& A )
    {
        GIVARO_ASSERT( rep.size() == A.rep.size(), "[Bits]: invalide size in |=");
        size_t len = rep.size();
        for (int i=0; i<int(len); ++i)
            rep[i] |= A.rep[i];
        return *this;
    }

    // XOr
    const Bits Bits::operator^ (const Bits& A) const
    {
        GIVARO_ASSERT( rep.size() == A.rep.size(), "[Bits]: invalide size in xor");
        size_t len = rep.size();
        Rep res( len );
        for (int i=0; i<int(len); ++i)
            res[i] = rep[i] ^ A.rep[i];
        return Bits(res);
    }

    Bits& Bits::xorin( const Bits& A, const Bits B)
    {
        GIVARO_ASSERT( B.rep.size() == A.rep.size(), "[Bits]: invalide size in xorin");
        size_t len = rep.size();
        for (int i=0; i<int(len); ++i)
            rep[i] = A.rep[i] ^ B.rep[i];
        return *this;
    }

    Bits& Bits::operator^=( const Bits& A )
    {
        GIVARO_ASSERT( rep.size() == A.rep.size(), "[Bits]: invalide size in ^=");
        size_t len = rep.size();
        for (int i=0; i<int(len); ++i)
            rep[i] ^= A.rep[i];
        return *this;
    }

    // Not
    const Bits Bits::operator~() const
    {
        size_t len = rep.size();
        Rep res( len );
        for (int i=0; i<int(len); ++i) {
            res[i] = ~rep[i];
        }
        return res;
    }

    Bits& Bits::notin( const Bits& A )
    {
        GIVARO_ASSERT( rep.size() == A.rep.size(), "[Bits]: invalide size in notin");
        size_t len = rep.size();
        for (int i=0; i<int(len); ++i)
            rep[i] = ~(A.rep[i]);
        return *this;
    }

    // -- Returns the number of non zero bits :
    long Bits::numone() const
    {
        size_t l = rep.size();
        long num = 0;
        for (int i=0; i<int(l); ++i)
        {
            int quo = QUO(i);
            int rem = REM(i);
            if (((rep[quo] & Table2pow[rem]) >> rem) !=0) num++;
        }
        return num;
    }

    // -- Returns the index of non zero bits
    void Bits::indexofone( Array0<base>& index) const
    {
        size_t l = rep.size();
        index.allocate( (size_t) numone() );
        long num =0;
        for (int i=0; i<int(l); ++i)
        {
            int quo = QUO(i);
            int rem = REM(i);
            if (((rep[quo] & Table2pow[rem]) >> rem) !=0)
                index[int(num++)] = (base)i;
        }
    }


    // -- Returns the length (in byte) of this :
    size_t Bits::length() const
    {
        return rep.size();
    }


    // return the i-th bit of *this
    int Bits::operator[] (const int i) const
    {
        GIVARO_ASSERT( i>=0, "invalide index in Bits::operator[]");
        GIVARO_ASSERT( i<(int)length()*SIZE_IN_BIT, "invalide index in Bits::operator[]");
        int quo = QUO(i);
        int rem = REM(i);
        return int((rep[quo] & Table2pow[rem]) >> rem);
    }

    int Bits::operator[] (const size_t i) const
    {
        // GIVARO_ASSERT( i>=0, "invalide index in Bits::operator[]");
        GIVARO_ASSERT( i<length()*SIZE_IN_BIT, "invalide index in Bits::operator[]");
        int quo = QUO(int(i));
        int rem = REM(int(i));
        return (int) ((rep[quo] & Table2pow[rem]) >> rem);
    }


    int Bits::get (const int i) const
    {
        GIVARO_ASSERT( i>=0, "invalide index in Bits::get");
        GIVARO_ASSERT( i<(int)length()*SIZE_IN_BIT, "invalide index in Bits::get");
        int quo = QUO(i);
        int rem = REM(i);
        return int((rep[quo] & Table2pow[rem]) >> rem);
    }

    // Set all bits
    void Bits::set()
    {
        int len = (int) rep.size();
        for (int i=0; i< len; ++i)
            rep[i] = (base) ~0L;
    }

    // Set the i-th bit of *this
    void Bits::set(const int i)
    {
        GIVARO_ASSERT( i>=0, "invalide index in Bits::set");
        GIVARO_ASSERT( i<(int)length()*SIZE_IN_BIT, "invalide index in Bits::set");
        int quo = QUO(i);
        int rem = REM(i);
        rep[quo] |= Table2pow[rem];
    }

    void Bits::clear()
    {
        size_t len = rep.size();
        for (int i=0; i<int(len); ++i)
            rep[i] = 0;
    }

    // Clear the i-th bit of *this
    void Bits::clear(const int i)
    {
        GIVARO_ASSERT( i>=0, "invalide index in Bits::clear");
        GIVARO_ASSERT( i<(int)length()*SIZE_IN_BIT, "invalide index in Bits::clear");
        int quo = QUO(i);
        int rem = REM(i);
        rep[quo] ^= !(Table2pow[rem]); // is the true not ?
    }

    std::ostream& Bits::print( std::ostream& o ) const
    {
        //-  o << "";
        size_t len = rep.size();
        for (int i= int(len-1); i>=0; i--)
        {
            for (int j=SIZE_IN_BIT-1; j>=0; j--)
                if ((rep[i] & Table2pow[j]) !=0) o << '1';
                else o << '0';
        }
        return o;
    }


    // -- Module definition
    void Bits::Init(int*, char***)
    {
        Table2pow = new Bits::base[SIZE_IN_BIT];
        int i;
        Table2pow[0] = 1;
        for (i=1; i<SIZE_IN_BIT; ++i) {
            Table2pow[i] = Table2pow[i-1] << 1;
        }
    }
    void Bits::End()
    {
    }

    GivModule Bits::Module (Bits::Init, Bits::End, GivModule::DfltPriority, "[Bits]");

} // namespace Givaro
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
