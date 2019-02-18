// ========================================================================== //
// $Source: /var/lib/cvs/Givaro/src/kernel/bstruct/givarray0.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Author: T. Gautier
// $Id: givarray0.inl,v 1.7 2011-02-02 16:23:55 briceboyer Exp $
// ======================================================================= //
/** @internal
 * @file bstruct/givarray0.inl
 * implementation of operators of Array0<T>
 */

#ifndef __GIVARO__array0_INL
#define __GIVARO__array0_INL

namespace Givaro {

    // Default cstor : ctsor of s size array
    template<class T>
    inline void Array0<T>::build( size_t s, const T& t)
    {
        // JGD 24.02.2011 : size_t is unsigned
        //  GIVARO_ASSERT( s>=0, "[Array<T>::cstor(size_t)] must takes a >=0 parameter");
        _psz = _size = s;
        if (s !=0) {
            _d = GivaroMM<T>::allocate(s);
            GivaroMM<T>::initialize(_d, s, t);
            _cnt = GivaroMM<int>::allocate(1);
            *_cnt = 1;
        } else { _d =0; _cnt =0; }
    }

    template<class T>
    inline Array0<T>::Array0 ( size_t s )
    {
        build(s, T());
    }

    template<class T>
    inline Array0<T>::Array0 ( size_t s, const T& t)
    {
        build(s, t);
    }


    // Recopy cstor : logical copy
    template<class T>
    inline Array0<T>::Array0 (const Self_t& p, givNoCopy)
    {
        _psz = p._psz; _size = p._size;
        if (_size !=0)
        { // increment ref. counting
            _d = p._d;
            _cnt = p._cnt; (*_cnt) ++;
        }
        else { _d = 0; _cnt = 0; }
    }


    // Recopy cstor : physical copy
    template<class T>
    inline Array0<T>::Array0 (const Self_t& p, givWithCopy)
    {
        _psz = _size = p._size;
        if (_size !=0) {
            _d = GivaroMM<T>::allocate(_size);
            _cnt = GivaroMM<int>::allocate(1);
            *_cnt = 1;
            for (size_t i=0; i<_size; i++)
                GivaroMM<T>::initone(&_d[i], p._d[i]);
        } else { _d =0; _cnt =0; }
    }

    // Destroy of the array
    template<class T>
    inline void Array0<T>::destroy( )
    {
        if (_psz !=0) {
            if (--(*_cnt) ==0)
            {
                GivaroMM<T>::destroy(_d, _psz);
                GivaroMM<T>::desallocate(_d);
                GivaroMM<int>::desallocate(_cnt);
            }
        }
        _size = _psz = 0; _cnt = 0; _d = 0;
    }

    // Allocation of an array of s Elements
    template<class T>
    inline void Array0<T>::allocate( size_t s )
    {
        // JGD 24.02.2011 : size_t is unsigned
        //  GIVARO_ASSERT( s>=0, "[Array<T>::allocate]: must takes a >=0 parameter");
        if (_cnt !=0) {
            if (((*_cnt) ==1) && (_psz >= s)) { _size = s; return; }
            this->destroy();
        }
        if (s >0) {
            _d = GivaroMM<T>::allocate(s);
            GivaroMM<T>::initialize(_d, s);
            _cnt = GivaroMM<int>::allocate(1);
            *_cnt = 1;
        }
        else _cnt =0;
        _psz = _size = s;
    }

    // Reallocation of an array of s Elements
    // and recopy the min(_size,s) first Elements
    template<class T>
    inline void Array0<T>::reallocate( size_t s )
    {
        // JGD 24.02.2011 : size_t is unsigned
        //  GIVARO_ASSERT( s>=0, "[Array<T>::reallocate]: must takes a >=0 parameter");
        if (_cnt !=0) {
            if (*_cnt ==1) {
                if (_psz >=s) { _size = s; return; }
            }
            else (*_cnt) --;
        }
        if (s >0) {
            T* tmp = GivaroMM<T>::allocate(s);
            GivaroMM<T>::initialize(tmp+_size, s-_size);
            if (_cnt !=0) {
                for (size_t i=0; i<_size; i++)
                    GivaroMM<T>::initone(&(tmp[i]), _d[i]);
                this->destroy();
            }
            _cnt = GivaroMM<int>::allocate(1);
            *_cnt = 1;
            _d = tmp;
        } else _cnt =0;
        _psz = _size = s;
    }

    template<class T>
    inline void Array0<T>::push_back( const T& a )
    {
        this->reallocate(_size+1);
        this->back() = a;
    }

    // Logical destructor: identical to free
    template<class T>
    inline Array0<T>::~Array0 ()
    {
        this->destroy();
    }


    // Physical copy : recopy and assignement on each Element
    template <class T>
    inline Array0<T>& Array0<T>::copy (const Array0<T>& src)
    {
        if (src._d == _d) return *this;
        reallocate(src._size); // - try...
        // -- here we have a large enough array with refcount==1
        const T* baseP = src._d;
        T* baseThis = _d;
        for (size_t i=0; i<_size; i++)
            baseThis[i] = baseP[i];
        return *this;
    }

    // Logical copy
    template <class T>
    inline Array0<T>& Array0<T>::logcopy (const Array0<T>& src)
    {
        destroy();
        _psz = src._psz; _size = src._size;
        if (_psz !=0)
        {
            _d = src._d;
            _cnt = src._cnt; (*_cnt) ++;
        }
        else { _d = 0; _cnt = 0; }
        return *this;
    }


    // Physical copy
    template<class T>
    Array0<T>& Array0<T>::operator= (const Array0<T>& p)
    {
        //throw GivError("[Array0<T>::operator=] cannot be used" " File:" ##__FILE__ ", line:" ##__LINE__ );
        //   throw GivError("[Array0<T>::operator=] cannot be used");
        return this->copy(p);
    }

    template<class T>
    inline size_t Array0<T>::phsize() const { return _psz; }

    template<class T>
    inline T* Array0<T>::baseptr() { return _d; }

    template<class T>
    inline const T* Array0<T>::baseptr() const { return _d; }


    // This following functions directly access to protected
    template <class T>
    inline const T& Array0<T>::front () const
    {
        GIVARO_ASSERT(_size >0, "[Array<T>::[]]: try to access to an Element of null size Array0.");
        return _d[0];
    }

    // Access operator : Write access
    template <class T>
    inline T& Array0<T>::front ()
    {
        GIVARO_ASSERT(_size >0, "[Array<T>::[]]: try to access to an Element of null size Array0.");
        return _d[0];
    }


    // This following functions directly access to protected
    template <class T>
    inline const T& Array0<T>::back () const
    {
        GIVARO_ASSERT(_size >0, "[Array<T>::[]]: try to access to an Element of null size Array0.");
        return _d[_size-1];
    }

    // Access operator : Write access
    template <class T>
    inline T& Array0<T>::back ()
    {
        GIVARO_ASSERT(_size >0, "[Array<T>::[]]: try to access to an Element of null size Array0.");
        return _d[_size-1];
    }



    template <class T>
    inline const T& Array0<T>::operator[] (Indice_t i) const
    {
        GIVARO_ASSERT(_size >0, "[Array<T>::[]]: try to access to an Element of null size Array0.");
        GIVARO_ASSERT((i >=0)&&(i<(Indice_t)_size), "[Array<T>::[]]: index out of bounds.");
        return _d[i];
    }

    // Access operator : Write access
    template <class T>
    inline T& Array0<T>::operator[] (Indice_t i)
    {
        GIVARO_ASSERT(_size >0, "[Array<T>::[]]: try to access to an Element of null size Array0.");
        GIVARO_ASSERT((i >=0)&&(i<(Indice_t)_size), "[Array<T>]: index out of bounds.");
        return _d[i];
    }


    template <class T>
    inline void Array0<T>::write( Indice_t i, const T& val)
    {
        GIVARO_ASSERT(_size >0, "[Array<T>::write]: try to access to an Element of null size Array0.");
        GIVARO_ASSERT((i >=0)&&(i<(Indice_t)_size), "[Array<T>::write]: index out of bounds.");
        _d[i] = val;
    }

    template <class T>
    inline void Array0<T>::read ( Indice_t i, T& val ) const
    {
        GIVARO_ASSERT(_size >0, "[Array<T>::read]: try to access to an Element of null size Array0");
        GIVARO_ASSERT((i >=0)&&(i<(Indice_t)_size), "[Array<T>::read]: index out of bounds.");
        val = _d[i];
    }

    template <class T>
    inline typename Array0<T>::Iterator_t
    Array0<T>::begin()
    {
        return _d;
    }

    template <class T>
    inline typename Array0<T>::Iterator_t
    Array0<T>::end()
    {
        return _d + _size;
    }

    template <class T>
    inline typename Array0<T>::constIterator_t
    Array0<T>::begin() const
    {
        return _d;
    }

    template <class T>
    inline typename Array0<T>::constIterator_t
    Array0<T>::end() const
    {
        return _d + _size;
    }

} // namespace Givaro

#endif // __GIVARO__array0_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
