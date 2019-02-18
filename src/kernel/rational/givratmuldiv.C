// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/rational/givratmuldiv.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama
// $Id: givratmuldiv.C,v 1.7 2010-04-14 16:20:30 jgdumas Exp $
// ==========================================================================
// Description:

#include "givaro/givrational.h"
#include "givaro/giverror.h"

namespace Givaro {

    // --------------------------------------- operator *
    Rational Rational::operator * (const Rational& r) const
    {

        if (isZero(r)) return Rational(0) ;
        if (isZero(*this)) return Rational(0) ;
        if (isOne(r)) return *this ;
        if (isOne(*this)) return r ;
        if (isInteger(*this) && isInteger(r))
            return Rational(num*r.num) ;

        if (absCompare(den, r.den) == 0)
            return Rational(num * r.num, den * r.den, 0) ;

        if (Rational::flags == Rational::NoReduce)
        {
            return Rational( num*r.num, den*r.den, 0) ;
        }
        Integer d1 = gcd(num, r.den);
        Integer d2 = gcd(den, r.num);

        return Rational( (num / d1) * (r.num / d2), (den / d2) * (r.den / d1), 0 );
    }


    // --------------------------------------- operator *=
    Rational& Rational::operator *= (const Rational& r)
    {

        if (isZero(r)) return *this=Rational(0) ;
        if (isZero(*this)) return *this ;
        if (isOne(r)) return *this ;
        if (isOne(*this)) return *this=r ;
        if (isInteger(*this) && isInteger(r)) {
            num *= r.num;
            return *this;
        }

        if ( (absCompare(den, r.den) == 0) || (Rational::flags == Rational::NoReduce) ) {
            num *= r.num;
            den *= r.den;
            return *this;
        }

        Integer d1 = gcd(num, r.den);
        Integer d2 = gcd(den, r.num);

        num /= d1;
        num *= (r.num / d2);
        den /= d2;
        den *= (r.den / d1);

        return *this;
    }


    // --------------------------------------- operator /
    Rational Rational::operator / (const Rational& r) const
    {
        if ( isZero(r) ) {
            throw GivMathDivZero("*** division by zero, in operator / (const Rational&)") ;
        }
        if (isZero(*this)) return Rational(0) ;
        if (isOne(r)) return *this ;
        if (isOne(*this))  {
            if (sign(r) < 0)
                return Rational(r.den, r.num, 0) ;
            else
                return Rational(-r.den, -r.num, 0) ;
        }

        if (absCompare(den, r.den) == 0)
            return Rational(num, r.num) ;

        if (Rational::flags == Rational::NoReduce)
            return Rational( num*r.den, den*r.num, 0) ;

        Integer d1 = gcd(num, r.num);
        Integer d2 = gcd(den, r.den);
        Integer resnum = (num / d1) * (r.den / d2);
        if (sign(r.num) < 0)
            resnum = -resnum ;
        Integer resden = (den / d2) * (r.num / d1);
        //rden can't be nul
        if (sign(resden) <0) resden = abs(resden) ;
        return Rational(resnum, resden, 0);
    }


    // --------------------------------------- operator /=
    Rational& Rational::operator /= (const Rational& r)
    {
        // std::cerr << "   BEGIN divin: " << *this << " /= " << r << std::endl;

        if ( isZero(r) ) {
            throw GivMathDivZero("*** division by zero, in operator / (const Rational&)") ;
        }
        if (isZero(*this)) return *this ;
        if (isOne(r)) return *this ;
        if (isOne(*this))  {
            if (sign(r.num)<0) {
                Integer::neg(this->num,r.den);
                Integer::neg(this->den,r.num);
            } else {
                this->num = r.den;
                this->den = r.num;
            }
            // std::cerr << "   DIVIN isone result: " << *this << std::endl;
            return *this;
        }

        if (compare(this->den, r.den) == 0) {
            if (sign(r.num)<0) {
                Integer::neg(this->den,r.num);
                Integer::negin( this->num );
            } else {
                this->den = r.num;
            }
            // std::cerr << "   DIVIN id result: " << this->num << " / " << this->den << std::endl;
            // std::cerr << "   DIVIN id den result: " << *this << std::endl;
            return this->reduce();
        }

        if (Rational::flags == Rational::NoReduce) {
            if (sign(r.num)<0) {
                this->num *= r.den;
                this->den *= r.num;
                Integer::negin( this->num );
                Integer::negin( this->den );
            } else {
                this->num *= r.den;
                this->den *= r.num;
            }
            return *this;
        }

        Integer d1 = gcd(this->num, r.num);
        Integer d2 = gcd(this->den, r.den);
        // std::cerr << "   d1: " << d1 << ", d2: " << d2 << std::endl;

        this->num /= d1;
        this->num *= (r.den / d2);
        //   if (sign(r.num) < 0)
        //     num = -num ;

        this->den /= d2;
        this->den *= (r.num / d1);

        //rden can't be neg
        if (sign(den) <0) {
            Integer::negin( this->num );
            Integer::negin( this->den );
        }

        // std::cerr << "   num: " << num << ", den: " << den << std::endl;
        // std::cerr << "   DIVIN result: " << *this << std::endl;

        return *this;
    }


    // --------------------------------------- operator /=
    Integer Rational::operator % (const Integer& r) const
    {
        if ( isZero(r) ) {
            throw GivMathDivZero("*** division by zero, in operator / (const Rational&)") ;
        }
        if (isZero(this->num)) return this->num ;

        Integer res(this->den);
        invin(res, r);
        return res *= this->num;
    }

} // namespace Givaro
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
