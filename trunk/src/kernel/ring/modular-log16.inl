// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/zpz/givzpz16table1.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: J.G. Dumas$
// Modified by Pascal Giorgi 2002/04/24
// $Id: givzpz16table1.inl,v 1.13 2011-02-04 14:11:46 jgdumas Exp $
// ==========================================================================
// Description:

// ---------
// -- normalized operations
// ---------

#ifndef __GIVARO_modular_log16_INL
#define __GIVARO_modular_log16_INL

#define __GIVARO_ZPZ16_LOG_MUL(r,p,a,b) ( (r)=  _tab_mul[(a) + (b)] )
#define __GIVARO_ZPZ16_LOG_DIV(r,p,a,b) ( (r)=  _tab_div[(a) - (b)] )
#define __GIVARO_ZPZ16_LOG_INV(r,p,b)   ( (r)=  _tab_div[ - (b)] )
#define __GIVARO_ZPZ16_LOG_SUB(r,p,a,b) ( (r)=  _tab_mul[(a) + _tab_subone[(b) - (a)] ] )
#define __GIVARO_ZPZ16_LOG_ADD(r,p,a,b) ( (r)=  _tab_mul[(a) + _tab_addone[(b) - (a)] ] )
#define __GIVARO_ZPZ16_LOG_NEG(r,p,a)   ( (r)=  _tab_neg[(a)] )

#define __GIVARO_ZPZ16_LOG_MUL_RES(r,p,a,b) ( (r)= (Residu_t) _tab_mul[(a) + (b)] )
#define __GIVARO_ZPZ16_LOG_DIV_RES(r,p,a,b) ( (r)= (Residu_t) _tab_div[(a) - (b)] )
#define __GIVARO_ZPZ16_LOG_INV_RES(r,p,b)   ( (r)= (Residu_t) _tab_div[ - (b)] )
#define __GIVARO_ZPZ16_LOG_SUB_RES(r,p,a,b) ( (r)= (Residu_t) _tab_mul[(a) + _tab_subone[(b) - (a)] ] )
#define __GIVARO_ZPZ16_LOG_ADD_RES(r,p,a,b) ( (r)= (Residu_t) _tab_mul[(a) + _tab_addone[(b) - (a)] ] )
#define __GIVARO_ZPZ16_LOG_NEG_RES(r,p,a)   ( (r)= (Residu_t) _tab_neg[(a)] )


/* Pascal Giorgi
   Changing the order of parameters.
   */
#define __GIVARO_ZPZ16_LOG_MULADD(r,p,a,b,c) \
{ __GIVARO_ZPZ16_LOG_MUL(r, p, a, b); __GIVARO_ZPZ16_LOG_ADD(r, p, r, c); }

// a*b-c
#define __GIVARO_ZPZ16_LOG_MULSUB(r,p,a,b,c) \
{ __GIVARO_ZPZ16_LOG_MUL(r, p, a, b); __GIVARO_ZPZ16_LOG_SUB(r, p, r, c); }

#define __GIVARO_ZPZ16_LOG_MULADD_RES(r,p,a,b,c) \
{ __GIVARO_ZPZ16_LOG_MUL_RES(r, p, a, b); __GIVARO_ZPZ16_LOG_ADD_RES(r, p, r, c); }
#define __GIVARO_ZPZ16_LOG_MULSUB_RES(r,p,a,b,c) \
{ __GIVARO_ZPZ16_LOG_MUL_RES(r, p, a, b); __GIVARO_ZPZ16_LOG_SUB_RES(r, p, r, c); }

namespace Givaro
{

    inline Modular<Log16, Log16>::Modular( Residu_t p ) :
        _p(p),_pmone(Residu_t(p-1)),zero(Rep(_pmone << 1)), one(0),mOne(Rep(_pmone>>1))
    {
        int32_t i,j;

        // tab value: Domain -> Rep, or something very similar
        _tab_value2rep = new Power_t[_p];
        // tab power: Rep -> Domain
        _tab_rep2value = new Residu_t[_p];

        _tab_rep2value[0] = 1;
        _tab_value2rep[1] = 0;

        int32_t fourp = ((int32_t)p) << 2, fourpmone= ((int32_t)_pmone)<<2;

        int not_found = 1;
        Residu_t accu = 1;
        Residu_t seed =2;

        // -- Find a generator of the multiplicative group
        while (_p > 2 && not_found == 1)
        {
            for(i=1; i<_p; i++)
            {
                accu = Residu_t((accu * seed) % _p);
                _tab_rep2value[i] = accu;
                if (accu == 1)
                    break;
                _tab_value2rep[accu] = Rep(i);
            }
            if (accu != 1){
                std::cerr << "attempted to build Log16 field with non-prime base "<<_p<<", halting\n";
                return;
            }
            if (i ==_p-1) not_found = 0;
            else {
                do {
                    seed = Residu_t(rand() % _p);
                } while ((seed ==0) && (seed !=1));
            }
        }
        // -- Set the zero at position 2 * _p - 2 in table
        _tab_value2rep[0] = zero;

        // -- Table for multiplication
        _tab_mul = new Power_t[(size_t)fourp];
        for(j=0; j<_pmone; j++)
            _tab_mul[j] = (Rep)j;
        for(j=_pmone; j< (int32_t)zero; j++)
                   _tab_mul[j] = Rep(j-_pmone);
        for(j=zero; j<= fourpmone; j++)
            _tab_mul[j] = zero;

        // -- Table for division and neg:
        _tab_div = &_tab_mul[_pmone];
        _tab_neg = &_tab_mul[_pmone/2];

        // -- Table for 1+value
        _tab_pone = new Power_t[(size_t)fourp];
        _tab_addone = &_tab_pone[(int32_t)(zero)];

        /* Pascal Giorgi 24/04/02
           Error between _tab_rep2value and _tab_value2rep
           corrected by inversing the array
           */
           
        for(j=0; j<_pmone; j++){
            if (_tab_rep2value[j] < _pmone)
                _tab_addone[j] = _tab_value2rep[ 1 + _tab_rep2value[j] ];
            else
                _tab_addone[j] = _tab_value2rep[0];
        }
        for(j=1-_pmone; j<0; j++){
            if (_tab_rep2value[j+_pmone] < _pmone)
                _tab_addone[j] = _tab_value2rep[ 1 + _tab_rep2value[j + _pmone] ];
            else
                _tab_addone[j] = _tab_value2rep[0];

        }
        for(j=_pmone; j<=(int32_t)zero; j++)
            _tab_addone[j] = 0;
        for(j=(int32_t)-zero; j<(int32_t)(1-_pmone); j++)
            _tab_addone[j] = (Rep)j;

        _tab_addone[_pmone / 2] = zero;
        _tab_addone[-_pmone / 2] = zero;


        // -- Table for 1-value
        _tab_mone = new Power_t[(size_t)fourp];
        _tab_subone = &_tab_mone[(int32_t)zero];

        for(j=_pmone; j<=(int32_t)zero; j++)
            _tab_subone[j] = 0;
        for(j=-zero; j<(int32_t)(1-3*_pmone/2); j++)
            _tab_subone[j] = Rep(j+_pmone/2);
        for(j=-3*_pmone/2; j<(1-_pmone); j++)
            _tab_subone[j] = Rep(j-_pmone/2);
        for(j=1-_pmone; j<(1-_pmone/2); j++)
            _tab_subone[j] = _tab_addone[j + _pmone/2 + _pmone];
        for(j=_pmone/2; j<_pmone; j++)
            _tab_subone[j] = _tab_addone[j - _pmone/2];

        for(j=-_pmone/2; j<_pmone/2; j++)
            _tab_subone[j] = _tab_addone[j+_pmone/2];

        numRefs = new int;
        (*numRefs) = 1;
    }

    inline Modular<Log16, Log16>::Modular(const Modular<Log16, Log16>& F) :
        _p ( F._p),
        _pmone ( F._pmone),
        _tab_value2rep ( F._tab_value2rep),
        _tab_rep2value ( F._tab_rep2value),
        _tab_mul ( F._tab_mul),
        _tab_div ( F._tab_div),
        _tab_neg ( F._tab_neg),
        _tab_addone ( F._tab_addone),
        _tab_subone ( F._tab_subone),
        _tab_mone ( F._tab_mone),
        _tab_pone ( F._tab_pone),
        numRefs ( F.numRefs),

        zero(F.zero), one(F.one),mOne(F.mOne)
    {
      (*numRefs)++;
    }

    Modular<Log16, Log16>& Modular<Log16, Log16>::operator=( const Modular<Log16, Log16>& F)
    {

        F.assign(const_cast<Element&>(one),F.one);
        F.assign(const_cast<Element&>(zero),F.zero);
        F.assign(const_cast<Element&>(mOne),F.mOne);


      if (this->numRefs) {
        (*(this->numRefs))--;
        if ((*(this->numRefs))==0) {
          delete [] _tab_value2rep;
          delete [] _tab_rep2value;
          delete [] _tab_mul;
          delete [] _tab_mone;
          delete [] _tab_pone;
          delete numRefs;
        }
      }

      this->_p = F.residu();
      this->_pmone = F._pmone;
      this->_tab_value2rep = F._tab_value2rep;
      this->_tab_rep2value = F._tab_rep2value;
      this->_tab_mul = F._tab_mul;
      this->_tab_div = F._tab_div;
      this->_tab_neg = F._tab_neg;
      this->_tab_mone = F._tab_mone;
      this->_tab_pone = F._tab_pone;
      this->_tab_addone = F._tab_addone;
      this->_tab_subone = F._tab_subone;
      this->numRefs = F.numRefs;
      (*(this->numRefs))++;

      return *this;
    }

    inline Modular<Log16, Log16>::~Modular()
    {
      (*numRefs)--;
      if (*numRefs == 0) {
        delete [] _tab_value2rep;
        delete [] _tab_rep2value;
        delete [] _tab_mul;
        delete [] _tab_mone;
        delete [] _tab_pone;
        delete numRefs;
      }
    }

    inline Modular<Log16, Log16>::Residu_t Modular<Log16, Log16>::residu( ) const
    {
        return _p;
    }

    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::mul (Rep& r, const Rep a, const Rep b) const
    {
        int32_t tmp;
        __GIVARO_ZPZ16_LOG_MUL(tmp,(int32_t)_p,(int32_t)a,(int32_t)b);
        return r= (Modular<Log16, Log16>::Rep)tmp;
    }

    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::div (Rep& r, const Rep a, const Rep b) const
    {
        __GIVARO_ZPZ16_LOG_DIV(r,_p,a,b);
        return r;
    }

    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::sub (Rep& r, const Rep a, const Rep b) const
    {
        __GIVARO_ZPZ16_LOG_SUB(r,_p,a,b);
        return r;
    }

    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::add (Rep& r, const Rep a, const Rep b) const
    {
        __GIVARO_ZPZ16_LOG_ADD(r,_p,a,b);
        return r;
    }

    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::neg (Rep& r, const Rep a) const
    {
        __GIVARO_ZPZ16_LOG_NEG(r,_p,a);
        return r;
    }

    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::inv (Rep& r, const Rep a) const
    {
        __GIVARO_ZPZ16_LOG_INV(r,_p,a);
        return r;
    }

    // -- inline array operations between Modular<Log16, Log16>::Rep
    inline void Modular<Log16, Log16>::mul (const size_t sz, Array r, constArray a, constArray b) const
    {
        for ( size_t i=sz ; --i ; ) {
            __GIVARO_ZPZ16_LOG_MUL_RES(r[i], _p,a[i], b[i]);
        }
    }

    inline void Modular<Log16, Log16>::mul (const size_t sz, Array r, constArray a, Rep b) const
    {
        for ( size_t i=sz ; --i ; ) {
            __GIVARO_ZPZ16_LOG_MUL_RES(r[i], _p, a[i], b);
        }
    }

    inline void Modular<Log16, Log16>::div (const size_t sz, Array r, constArray a, constArray b) const
    {
        for ( size_t i=sz ; --i ; ) {
            __GIVARO_ZPZ16_LOG_DIV_RES( r[i], _p, a[i], b[i]);
        }
    }

    inline void Modular<Log16, Log16>::div (const size_t sz, Array r, constArray a, Rep b) const
    {
        for ( size_t i=sz ; --i ; ) {
            __GIVARO_ZPZ16_LOG_DIV_RES( r[i], _p, a[i], b);
        }
    }

    inline void Modular<Log16, Log16>::add (const size_t sz, Array r, constArray a, constArray b) const
    {
        for ( size_t i=sz ; --i ; ) {
            __GIVARO_ZPZ16_LOG_ADD_RES(r[i], _p, a[i], b[i]);
        }
    }

    inline void Modular<Log16, Log16>::add (const size_t sz, Array r, constArray a, Rep b) const
    {
        for ( size_t i=sz ; --i ; ) {
            __GIVARO_ZPZ16_LOG_ADD_RES(r[i], _p, a[i], b);
        }
    }

    inline void Modular<Log16, Log16>::sub (const size_t sz, Array r, constArray a, constArray b) const
    {
        for ( size_t i=sz ; --i ; ) {
            __GIVARO_ZPZ16_LOG_SUB_RES(r[i], _p, a[i], b[i]);
        }
    }

    inline void Modular<Log16, Log16>::sub (const size_t sz, Array r, constArray a, Rep b) const
    {
        for ( size_t i=sz ; --i ; ) {
            __GIVARO_ZPZ16_LOG_SUB_RES(r[i], _p, a[i], b);
        }
    }

    inline void Modular<Log16, Log16>::neg (const size_t sz, Array r, constArray a) const
    {
        for ( size_t i=sz ; --i ; ) {
            __GIVARO_ZPZ16_LOG_NEG_RES(r[i], _p, a[i]);
        }
    }


    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::mulin (Rep& r, const Rep a) const
    {
        __GIVARO_ZPZ16_LOG_MUL(r,_p, r,a);
        return r;
    }

    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::divin (Rep& r, const Rep a) const
    {
        Modular<Log16, Log16>::Rep ia;
        inv(ia, a);
        mulin(r, ia);
        return r;
    }

    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::addin (Rep& r, const Rep a) const
    {
        __GIVARO_ZPZ16_LOG_ADD(r, _p, r,a);
        return r;
    }

    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::subin (Rep& r, const Rep a) const
    {
        __GIVARO_ZPZ16_LOG_SUB(r,_p, r,a);
        return r;
    }

    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::negin (Rep& r) const
    {
        __GIVARO_ZPZ16_LOG_NEG(r,_p,r);
        return r;
    }

    inline Modular<Log16, Log16>::Rep&  Modular<Log16, Log16>::invin (Rep& r) const
    {
        __GIVARO_ZPZ16_LOG_INV(r,_p,r);
        return r;
    }

    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::axpy
    (Rep& r, const Rep a, const Rep b, const Rep c) const
    {
        __GIVARO_ZPZ16_LOG_MULADD(r, _p, a, b, c);
        return r;
    }

    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::axpyin
    (Rep& r, const Rep a, const Rep b) const
    {
        return axpy(r,a,b,r);
    }


    inline void Modular<Log16, Log16>::axpy
    (const size_t sz, Array r, constArray a, constArray x, constArray y) const
    {
        for ( size_t i=sz ; --i ; ) {
            __GIVARO_ZPZ16_LOG_MULADD_RES(r[i], _p, a[i], x[i], y[i]);
        }
    }

    inline void Modular<Log16, Log16>::axpyin
    (const size_t sz, Array r, constArray a, constArray x) const
    {
        for ( size_t i=sz ; --i ; ) {
            __GIVARO_ZPZ16_LOG_MULADD_RES(r[i], _p, a[i], x[i], r[i]);
        }
    }

    // r <- a*b-c
    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::axmy
    (Rep& r, const Rep a, const Rep b, const Rep c) const
    {
        __GIVARO_ZPZ16_LOG_MULSUB(r,_p,a,b,c);
        return r;
    }

    // r <- r-a*b
    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::maxpyin
    (Rep& r, const Rep a, const Rep b) const
    {
        Rep t; __GIVARO_ZPZ16_LOG_MUL(t,_p,a,b);
        return this->subin(r,t);
    }

    // r <- c-a*b
    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::maxpy
    (Rep& r, const Rep a, const Rep b, const Rep c) const
    {
        Rep t; __GIVARO_ZPZ16_LOG_MUL(t,_p,a,b);
        return this->sub(r,c,t);
    }


    // r <- a*b-r
    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::axmyin (Rep& r,
                              const Rep a, const Rep b) const
    {
        Rep t; __GIVARO_ZPZ16_LOG_MUL(t,_p,a,b);
        return sub(r,t,r);
    }

    inline void Modular<Log16, Log16>::axmy
    (const size_t sz, Array r, constArray a, constArray x, constArray y) const
    {
        for ( size_t i=sz ; --i ; ) {
            __GIVARO_ZPZ16_LOG_MULSUB_RES(r[i], _p, a[i], x[i], y[i]);
        }
    }

    inline void Modular<Log16, Log16>::maxpyin (const size_t sz, Array r,
                        constArray a, constArray x) const
    {
        for ( size_t i=sz ; --i ; ) {
            __GIVARO_ZPZ16_LOG_MULSUB_RES(r[i], _p, a[i], x[i], r[i]);
            __GIVARO_ZPZ16_LOG_NEG_RES(r[i], _p, r[i]);
        }
    }

    // ------------------------- Miscellaneous functions

    inline int Modular<Log16, Log16>::iszero(const Rep a) const
    {
        return a >= _p;
    }

    inline int Modular<Log16, Log16>::isone(const Rep a) const
    {
        return a == Modular<Log16, Log16>::one;
    }

    inline int Modular<Log16, Log16>::ismone(const Rep a) const
    {
        return a == Modular<Log16, Log16>::mOne;
    }


    inline size_t Modular<Log16, Log16>::length(const Rep ) const
    {
        return Modular<Log16, Log16>::size_rep;
    }

    inline int Modular<Log16, Log16>::isZero( const Rep a ) const
    {
        return iszero(a);
    }
    inline int Modular<Log16, Log16>::isOne ( const Rep a ) const
    {
        return isone(a);
    }
    inline int Modular<Log16, Log16>::isMOne ( const Rep a ) const
    {
        return ismone(a);
    }


    // ---------
    // -- misc operations
    // ---------
#if 0
    inline void Modular<Log16, Log16>::assign
    ( const size_t sz, Array r, constArray a ) const
    {
        for ( size_t i=sz ; --i ; ) {
            if (a[i] <Modular<Log16, Log16>::zero) {
                r[i] = a[i] + _p;
                if (r[i] <Modular<Log16, Log16>::zero) r[i] = r[i] % _p;
            }
            else if (a[i] >_p) {
                r[i] = a[i] - _p;
                if (r[i] >_p) r[i] = r[i] % _p;
            }
            else r[i] = a[i];
        }
    }
#endif

    inline void Modular<Log16, Log16>::assign ( const size_t sz, Array r, constArray a ) const
    {
        for ( size_t i=sz ; --i ; )
            r[i] = a[i];
    }



    // initialized by a degree of the generator.
    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::init ( Rep& r ) const
    {
        return r = zero;
    }


    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::assign ( Rep& r, const Rep a ) const
    {
        return r = a;
    }

    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::init ( Rep& r, const long a ) const
    {
        int sign; unsigned long ua;
        if (a <0) {
            sign =-1;
            ua = (unsigned long)-a;
        }
        else {
            ua = (unsigned long)a;
            sign =1;
        }
        r = Rep( (ua >=_p) ? ua % _p : ua );
        if (sign ==-1)
            r = Rep(_p - r);
        assert(r < _p);
        return r = _tab_value2rep[r];
    }

    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::init ( Rep& r, const int a ) const
    {
        return Modular<Log16, Log16>::init( r, (long)a);
    }

    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::init ( Rep& r, const unsigned long a ) const
    {
        r = Rep((a >=_p) ? a % _p : a);
        assert(r < _p);
        return r= _tab_value2rep[r];
    }

    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::init ( Rep& r, const unsigned int a ) const
    {
        r = Rep((a >=_p) ? a % _p : a);
        assert(r < _p);
        return r= _tab_value2rep[r];
    }

    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::init ( Rep& r, const uint16_t a ) const
    {
        r = Rep((a >=_p) ? a % _p : a);
        assert(r < _p);
        return r= _tab_value2rep[r];
    }

    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::init ( Rep& r, const int16_t a ) const
    {
        return Modular<Log16, Log16>::init( r, (long)a);
    }

    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::init( Rep& a, const double i) const
    {
        return init(a,(long)i);
    }
    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::init( Rep& a, const float i) const
    {
        return init(a,(double)i);
    }

    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::init ( Rep& r, const Integer& Residu ) const
    {
        int16_t tr;
        if (Residu <0) {
            // -a = b [p]
            // a = p-b [p]
            if ( Residu <= (Integer)(-_p) ) tr = int16_t( (-Residu) % _p) ;
            else tr = int16_t(-Residu);
            if (tr){
                assert(_p -(uint16_t)tr < _p);
                return r = _tab_value2rep[ _p - (uint16_t)tr ];
            }
            else
                return r = (Rep) zero;
        } else {
            if (Residu >= (Integer)_p ) tr =   int16_t(Residu % _p) ;
            else tr = int16_t(Residu);
            assert(tr < _p);
            return r = _tab_value2rep[tr];
        }
    }



    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::dotprod
    ( Rep& r, const int bound, const size_t sz, constArray a, constArray b ) const
    {
        unsigned int stride = 1;
        if ((unsigned long)bound < GIVARO_MAXUINT16)
            stride = (unsigned int) ( GIVARO_MAXUINT32/((unsigned long)bound * (unsigned long)bound) );
        uint32_t dot = (uint32_t)zero;
        if ((sz <10) && (sz <stride)) {
            for(  size_t i= sz; i--; )
                dot += _tab_rep2value[a[i]] * _tab_rep2value[b[i]];
            if (dot > _p){
                assert( (Rep)(dot %_p) < _p);
                return r = _tab_value2rep[(Rep)(dot % _p)];
            }
            else {
                assert(dot < _p);
                return r = _tab_value2rep[dot];
            }

        }
        unsigned int i_begin=0;
        stride &= (unsigned int)~0x1;
        if (stride ==0) {
            for(  size_t i= sz-1; i>0; --i) {
                dot += _tab_rep2value[a[i]] * _tab_rep2value[b[i]];
                if (dot>_p) dot %= _p;
            }
            assert(dot < _p);
            return r = _tab_value2rep[dot];

        }
        do {
            size_t min_sz = ((sz-i_begin) < stride ? (sz-i_begin) : stride);
            if ( (min_sz & 0x1) !=0) {
                min_sz--; i_begin++;
                dot += _tab_rep2value[a++[min_sz]] * _tab_rep2value[b++[min_sz]];
            }
            if (min_sz > 1)
                for(  size_t i= min_sz; i>0; --i, --i, ++a, ++a, ++b, ++b )
                {
                    dot += _tab_rep2value[a[0]] * _tab_rep2value[b[0]];
                    dot += _tab_rep2value[a[1]] * _tab_rep2value[b[1]];
                }
            if (dot>_p) dot %= _p;
            i_begin += (unsigned int) min_sz;
        } while (i_begin <sz);
        assert(dot < _p);
        return r = _tab_value2rep[dot];
    }

    inline Modular<Log16, Log16>::Rep& Modular<Log16, Log16>::dotprod
    ( Rep& r, const size_t sz, constArray a, constArray b ) const
    {
        return Modular<Log16, Log16>::dotprod(r, _p, sz, a, b);
    }


    //  a -> r: int16_t to double
    inline void
    Modular<Log16, Log16>::i2d ( const size_t sz, double* r, constArray a ) const
    {
        for (size_t i=0; i<sz; ++i) r[i] = _tab_rep2value[a[i]];
    }

    //  a -> r: double to int16_t
    inline void
    Modular<Log16, Log16>::d2i ( const size_t sz, Array r, const double* a ) const
    {
        union d_2_l {
            double d;
            int32_t r[2];
        };
        static const double offset = 4503599627370496.0; // 2^52
        size_t i=sz-1;
        //warning todo while
        //do
label1:
        {
            d_2_l tmp;
            // - normalization: put fractional part at the end of the representation
            tmp.d = a[i] + offset;
            {
                assert((tmp.r[1] >=_p ? tmp.r[1] : tmp.r[1] % _p) < _p);
            }
            r[i--] = (Residu_t)_tab_value2rep[(tmp.r[1] >=_p ? tmp.r[1] : tmp.r[1] % _p)];
        }
        // while (i!=0)
        if (i >0) goto label1;
        //for (size_t i=sz-1; i>=0; --i)
    }


    // -- Input: (z, <_p>)
    inline std::istream& Modular<Log16, Log16>::read (std::istream& s)
    {
        char ch;
        s >> std::ws >> ch;
        //   if (ch != '(')
        //     GivError::throw_error( GivBadFormat("Modular<Log16, Log16>::read: syntax error: no '('"));
        if (ch != '(')
            std::cerr << "Modular<Log16, Log16>::read: syntax error: no '('" << std::endl;

        s >> std::ws >> ch;
        //   if (ch != 'z')
        //     GivError::throw_error( GivBadFormat("Modular<Log16, Log16>::read: bad domain object"));
        if (ch != 'z')
            std::cerr << "Modular<Log16, Log16>::read: bad domain object" << std::endl ;

        s >> std::ws >> ch;
        //   if (ch != ',')
        //     GivError::throw_error( GivBadFormat("Modular<Log16, Log16>::read: syntax error: no ','"));
        if (ch != ',')
            std::cerr << "Modular<Log16, Log16>::read: syntax error: no ','" << std::endl;


        s >> std::ws >> _p;

        s >> std::ws >> ch;
        //   if (ch != ')')
        //     GivError::throw_error( GivBadFormat("Modular<Log16, Log16>::read: syntax error: no ')'"));
        if (ch != ')')
            std::cerr << "Modular<Log16, Log16>::read: syntax error: no ')'" << std::endl;

        return s;
    }

    inline std::ostream& Modular<Log16, Log16>::write (std::ostream& s ) const
    {
        return s << "Log16 Givaro Z/pZ modulo " << residu();
    }

    inline std::istream& Modular<Log16, Log16>::read (std::istream& s, Rep& a) const
    {
        Integer tmp;
        s >> tmp;
        tmp %= _p;
        if (tmp < 0) tmp += _p;
        assert ( (uint)tmp < _p) ;
        a = _tab_value2rep[ (uint)tmp ];
        return s;
    }

    inline std::ostream& Modular<Log16, Log16>::write (std::ostream& s, const Rep a) const
    {
        if (a >= _p) return s << '0';
        return s << _tab_rep2value[a]; //dpritcha
    }

} // namespace Givaro

#endif // __GIVARO_modular_log16_INL

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
