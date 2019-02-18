// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/bstruct/givlist0.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Author: T. Gautier
// $Id: givlist0.inl,v 1.3 2011-02-02 16:23:55 briceboyer Exp $
// ==========================================================================
// inline implementation of the List0 methods

#ifndef __GIVARO_list0_INL
#define __GIVARO_list0_INL

namespace Givaro {

    template <class T>
    inline List0<T>::List0()
    : _head(0), _queue(0), _size(0)
    {}

    template <class T>
    inline List0<T>::List0( const List0<T>& l)
    : _head(l._head), _queue(l._queue), _size(l._size)
    {
        node* curr = _head;
        while (curr !=0) {
            curr->cnt++ ;
            curr = curr->next ;
        }
    }

    template <class T>
    inline void List0<T>::destroy()
    {
        node* curr = _head;
        while (curr !=0) {
            node* tmp = curr->next ;
            if (--*(curr->cnt) ==0) {
                GivaroMM<node>::destroy(curr,1);
                GivaroMM<node>::desallocate(curr,1);
            }
            curr = tmp ;
        }
        _head = _queue = 0 ;
        _size = 0 ;
    }

    template <class T>
    inline List0<T>::~List0()
    { destroy() ; }

    template <class T>
    inline List0<T>& List0<T>::copy(const List0<T>& l)
    {
        destroy();
        node* curr = l._head ;
        while (curr !=0) {
            this->insertlast( *(curr->item) );
            curr = curr->next ;
        }
        return *this;
    }

    template <class T>
    inline List0<T>& List0<T>::logcopy(const List0<T>& l)
    {
        destroy();
        _head = l._head;
        _queue = l._queue;
        _size = l._size ;
        node* curr = _head;
        while (curr !=0) {
            curr->cnt++ ;
            curr = curr->next ;
        }
        return *this;
    }

    template <class T>
    inline List0<T>& List0<T>::operator=(const List0<T>& l)
    { return logcpy(l); }

    template <class T>
    inline  size_t List0<T>::size() const
    { return _size ; }

    template <class T>
    inline  int List0<T>::is_empty() const
    { return _size ==0; }

    template<class T>
    void List0<T>::insertlast( const T& val)
    {
        node* curr = GivaroMM<node>::allocate(1);
        curr->cnt =  GivaroMM<int>::allocate(1); *(curr->cnt) = 1 ;
        curr->item = GivaroMM<T>::allocate(1);
        GivaroMM<T>::initone(curr->item, val);
        curr->prev = _queue ;
        curr->next = 0;
        if (_queue !=0) _queue->next = curr ;
        _queue = curr ;
        if (_head ==0) _head = curr ;
        _size++ ;
    }

    template<class T>
    void List0<T>::insertfirst( const T& val)
    {
        node* curr = GivaroMM<node>::allocate(1);
        curr->cnt =  GivaroMM<int>::allocate(1); *(curr->cnt) = 1 ;
        curr->item = GivaroMM<T>::allocate(1);
        GivaroMM<T>::initone(curr->item, val);
        curr->prev = 0;
        curr->next = _head;
        if (_head !=0) _head->prev = curr ;
        _head = curr ;
        if (_queue ==0) _queue = curr ;
        _size++ ;
    }

    template<class T>
    int List0<T>::getlast(T& val) const
    {
        if (_queue ==0) return 0;
        val = *(_queue->item);
        return 1;
    }

    template<class T>
    int List0<T>::getfirst(T& val) const
    {
        if (_head ==0) return 0;
        val = *(_head->item);
        return 1;
    }

    template<class T>
    int List0<T>::getrmvlast(T& val)
    {
        if (_queue ==0) return 0;
        node* curr = _queue;
        val = *(curr->item);
        _queue = curr->prev;
        if (_queue !=0) _queue->next = 0;
        else _head = 0 ; // one item in the list !
        if (--*(curr->cnt) ==0) {
            GivaroMM<node>::destroy(curr,1);
            GivaroMM<node>::desallocate(curr,1);
        }
        --_size;
        return 1;
    }

    template<class T>
    int List0<T>::getrmvfirst(T& val)
    {
        if (_head ==0) return 0;
        node* curr = _head;
        val = *(curr->item);
        _head = curr->next;
        if (_head !=0) _head->prev = 0;
        else _queue = 0 ; // one item in the list !
        if (--*(curr->cnt) ==0) {
            GivaroMM<node>::destroy(curr,1);
            GivaroMM<node>::desallocate(curr,1);
        }
        --_size;
        return 1;
    }

} // namespace Givaro

#endif // __GIVARO_list0_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
