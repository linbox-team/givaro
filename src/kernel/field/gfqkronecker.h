// ==========================================================================
// Copyright(c)'1994-2010 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// file: gfqkronecker.h
// Time-stamp: <12 Apr 10 16:45:02 Jean-Guillaume.Dumas@imag.fr>
// date: 2007
// version:
// author: Jean-Guillaume.Dumas

/*! @file gfqkronecker.h
 * @ingroup zpz
 * @brief  Arithmetic on GF(p^k), with dynamic Kronecker substitution.
 * @pre  k(p-1)^2 < word size
 */

#ifndef __GIVARO_gfq_kronecker_H
#define __GIVARO_gfq_kronecker_H

#include "givaro/givzpz.h"
#include "givaro/givzpzInt.h"
#include "givaro/gfq.h"
#include "givaro/givpower.h"
#include <limits>
#include <vector>
#include <deque>

namespace Givaro {

	//! GFqKronecker
template<class TT, class Ints> struct GFqKronecker : public GFqDom<TT> {
protected:
    typedef typename Signed_Trait<TT>::unsigned_type UTT;
    typedef TT Rep;
    typedef GFqDom<TT> Father_t;

public:
    typedef GFqKronecker<TT,Ints> Self_t;

    typedef Rep Element;
    typedef Element* Element_ptr ;
  typedef const Element* ConstElement_ptr;


    typedef UTT Residu_t;

    typedef Rep* Array;
    typedef const Rep* constArray;

    typedef ModularRandIter< Father_t , Rep> RandIter;

    GFqKronecker(): Father_t() {}

        // Extension MUST be a parameter of the constructor
    GFqKronecker( const UTT P, const UTT e) : Father_t(P,e), _degree(e-1),_epmunsq(e*(P-1)*(P-1))  {
        buildsmalltables();
        setShift(std::numeric_limits<UTT>::digits/( (e<<1)-1) );
    }

    Ints getMaxn() const {
	return _sMAXN;
    }
    UTT getShift() const {
	return _SHIFTS;
    }
    UTT getBase() const {
	return _sBASE;
    }

        // Set shifts, returns maxn
    Ints setShift(const Ints& i) {
        _SHIFTS=i;
        _sBASE = Ints(1); _sBASE <<=_SHIFTS;
        _sMASK = _sBASE - 1;
        return _sMAXN = _sBASE / _epmunsq;
    }

        // Set maxn, returns shifts
    UTT setMaxn(const Ints& n) {
        _sMAXN = n;
        Ints m = _sMAXN * _epmunsq;
        _SHIFTS = 0;
        for(_sBASE = 1; _sBASE < m; ++_SHIFTS, _sBASE<<=1);
        _sMASK = _sBASE - 1;
        return _SHIFTS;
    }


    virtual ~GFqKronecker() {};

    using Father_t::init;
    using Father_t::convert;

    virtual Ints& convert(Ints& r, const Rep a) const {
            // First Step
            // 	from element to polynomial coefficient binary shifted
        UTT binpol=_log2bin[a];
            // Second step
            // 	from a0|a1|...|ad
            // 	to 0-0ad | ... | 0-0a1 | 0-0a0
        r = binpol & _pMASK; // a0
        for(size_t i=1; i<this->_exponent; ++i) {
            r <<= _SHIFTS;
            binpol >>= _pceil;
            r |= ( binpol & _pMASK);
        }
        return r;
    }



    virtual Rep& init(Rep& a, const Ints r) const {
            // WARNING: This could be speeded up with a REDQ transform

            // First Step lower part
            // 	from rd | ... | r1 | r0
            // 	to   a0|a1|...|ad
            // 	where ai = ri mod p
        Ints rs=r;
        UTT binpolLOW = (UTT)( (rs & _sMASK) % this->_characteristic);
        for(size_t i=1; i<this->_exponent; ++i) {
            binpolLOW <<= _pceil;
            rs >>= _SHIFTS;
            binpolLOW |= (UTT)( (rs & _sMASK) % this->_characteristic);
        }
            // First Step upper part
        rs >>= _SHIFTS;
        UTT binpolHIGH = (UTT)( (rs & _sMASK) % this->_characteristic);
        for(size_t i=1; i<this->_degree; ++i) {
            binpolHIGH <<= _pceil;
            rs >>= _SHIFTS;
            binpolHIGH |= (UTT)( (rs & _sMASK) % this->_characteristic);
        }
        binpolHIGH <<= _pceil;

            // Second step
            //  transform to element H*X^k+L
        return this->axpy(a, _bin2log[binpolHIGH], _Xk, _bin2log[binpolLOW]);
    }




protected:
    std::ostream& polywrite(std::ostream& out, const Element& a, const Indeter In= "B") const {
        static Modular<Integer> Zp(this->_characteristic);
        static Poly1PadicDom<Modular<Integer> > PAD(Zp ,In);
        static Poly1PadicDom<Modular<Integer> >::Element pol;
        Integer r;
        PAD.radixdirect(pol, this->convert(r, a), this->_exponent);
        return PAD.write(out, pol);
    }


    void buildsmalltables() {
        _log2bin.resize(this->_log2pol.size());
        _pceil = 1;
        for(unsigned long ppow = 2; ppow < this->_characteristic; ppow <<= 1,++_pceil) {}
        _pMASK = (1<<_pceil) - 1;

        Father_t Zp(this->_characteristic,1);
        typedef Poly1FactorDom< Father_t, Dense > PolDom;
        PolDom Pdom( Zp );
        typedef Poly1PadicDom< Father_t, Dense > PadicDom;
        PadicDom PAD(Pdom);

        typename std::vector<UTT>::iterator binit = _log2bin.begin();
        typename std::vector<UTT>::const_iterator polit = this->_log2pol.begin();
        for( ; polit != this->_log2pol.end(); ++polit, ++binit) {

            std::vector<unsigned long> vect;
            PAD.radixdirect( vect, (unsigned long)(*polit), this->_exponent);

            *binit = vect[0];
            for(size_t i =1; i<this->_exponent; ++i) {
                *binit <<= _pceil;
                *binit += vect[i];
            }
        }

        _bin2log.resize( 1<<(_pceil*this->_exponent) );
        for(size_t i=0; i<_log2bin.size(); ++i)
            _bin2log[ _log2bin[ i ] ] = i;


            // This is X
        _Xk = this->_pol2log[this->_characteristic];
//         polywrite(std::cerr << "Xk: " << _Xk << " rep ", _Xk) << std::endl;
            // This is X^{e}
        dom_power(_Xk,_Xk,this->_exponent,*this);
//         polywrite(std::cerr << "Xk: ", _Xk) << std::endl;
    }



    UTT _SHIFTS;
    Ints _sBASE,_sMASK,_sMAXN;

    UTT _pceil,_pMASK,_degree;
    Ints _epmunsq;

    std::vector<UTT> _log2bin;
    std::vector<UTT> _bin2log;

    Element _Xk;
};

} // namespace Givaro

#endif // __GIVARO_gfq_kronecker_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
