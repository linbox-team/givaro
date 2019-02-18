// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/bstruct/givhashtable.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Author: T. Gautier
// $Id: givhashtable.inl,v 1.3 2011-02-02 16:23:55 briceboyer Exp $
// ==========================================================================
// Description:
// - hash table implementation

#ifndef __GIVARO_hashtable_INL
#define __GIVARO_hashtable_INL
namespace Givaro {
    // Size of the table :
    template<class T, class Key>
    HashTable<T,Key>::HashTable ( int n )
    {
        allocate(n) ;
    }

    template<class T, class Key>
    void HashTable<T,Key>::allocate( int n )
    {
        num = n ;
        if (n ==0) {
            tabH = tabE = 0 ;
            return  ;
        }
        tabH = new _E*[n] ;
        tabE = new _E*[n] ;
        for (int i=0 ; i<num ; i++) tabH[i] = tabE[i] = 0 ;
    }

    template<class T, class Key>
    HashTable<T,Key>::~HashTable()
    {
        freeing() ;
    }

    template<class T, class Key>
    void HashTable<T,Key>::freeing()
    {
        if (tabH ==0) return ;
        for (int i=0 ; i<num ; i++)
        {
            _E* curr = tabH[i] ;
            while (curr !=0)
            {
                _E* tmp = curr->next ;
                delete curr ;
                curr = tmp ;
            }
        }
        delete [] tabH ;
        delete [] tabE ;
        tabH = 0 ;
        tabE = 0 ;
    }

    template<class T, class Key>
    void HashTable<T,Key>::removeall ( )
    {
        if (tabH ==0) return ;
        int i ;
        for (i=0 ; i<num ; i++)
        {
            _E* curr = tabH[i] ;
            while (curr !=0)
            {
                _E* tmp = curr->next ;
                delete curr ;
                curr = tmp ;
            }
        }
        for (i=0 ; i<num ; i++) tabH[i] = tabE[i] = 0 ;
    }

    template<class T, class Key>
    void HashTable<T,Key>::insert( const T& item, const Key& k )
    {
        if ((tabE ==0) || (tabH ==0)) return ;
        insertLast( item, k ) ;
    }

    template<class T, class Key>
    void HashTable<T,Key>::insertLast( const T& item, const Key& k )
    {
        if ((tabE ==0) || (tabH ==0)) return ;
        unsigned int e = ((unsigned long)k.godel()) % ((unsigned long)num) ;
        _E* cur = new _E ( item, k) ;
        if (tabE[e] == 0)
        {
            tabH[e] = tabE[e] = cur ;
            cur->next = cur->pred = 0 ;
            return ;
        }
        tabE[e]->next = cur ;
        cur->pred = tabE[e] ;
        cur->next = 0 ;
        tabE[e] = cur ;
    }


    template<class T, class Key>
    void HashTable<T,Key>::insertFront( const T& item, const Key& k )
    {
        if ((tabE ==0) || (tabH ==0)) return ;
        unsigned int e = ((unsigned long)k.godel()) % ((unsigned long)num) ;
        _E* cur = new _E ( item, k) ;
        if (tabH[e] == 0)
        {
            tabH[e] = tabE[e] = cur ;
            cur->next = cur->pred = 0 ;
            return ;
        }
        tabH[e]->pred = cur ;
        cur->next = tabE[e] ;
        cur->pred = 0 ;
        tabH[e] = cur ;
    }

    // Remove entry with key k
    template<class T, class Key>
    void HashTable<T,Key>::remove ( const Key& k)
    {
        if ((tabE ==0) || (tabH ==0)) return ;
        unsigned int e = ((unsigned long)k.godel()) % ((unsigned long)num) ;
        _E* cur = tabH[e] ;
        while (cur !=0)
        {
            if (k == cur->key)
            {
                if (cur->pred ==0) tabH[e] = cur->next ;
                else cur->pred->next = cur->next ;
                if (cur->next ==0) tabE[e] = cur->pred ;
                else cur->next->pred = cur->pred ;
                delete cur ;
            }
            cur = cur->next ;
        }
    }

    template<class T, class Key>
    int HashTable<T,Key>::get    ( T& item, const Key& k ) const
    {
        if ((tabE ==0) || (tabH ==0)) return 0 ;
        unsigned int e = ((unsigned long)k.godel()) % ((unsigned long)num) ;
        _E* cur = tabH[e] ;
        while (cur !=0)
        {
            if (k == cur->key)
            {
                item = cur->oneitem ;
                return 1 ;
            }
            cur = cur->next ;
        }
        return 0 ;
    }



    template<class T, class Key>
    int HashTable<T,Key>::getrmv ( T& item, const Key& k )
    {
        if ((tabE ==0) || (tabH ==0)) return 0 ;
        unsigned int e = ((unsigned long)k.godel()) % ((unsigned long)num) ;
        _E* cur = tabH[e] ;
        while (cur !=0)
        {
            if (k == cur->key)
            {
                item = cur->oneitem ;
                if (cur->pred ==0) tabH[e] = cur->next ;
                else cur->pred->next = cur->next ;
                if (cur->next ==0) tabE[e] = cur->pred ;
                else cur->next->pred = cur->pred ;
                delete cur ;
                return 1 ;
            }
            cur = cur->next ;
        }
        return 0 ;
    }

    // Returns the number of items of the same key :
    template<class T, class Key>
    int HashTable<T,Key>::num_item( const Key& k ) const
    {
        if ((tabE ==0) || (tabH ==0)) return 0 ;
        int count = 0 ;
        unsigned int e = ((unsigned long)k.godel()) % ((unsigned long)num) ;
        _E* cur = tabH[e] ;
        while (cur !=0)
        {
            if (k == cur->key) count++ ;
            cur = cur->next ;
        }
        return count ;
    }

} // namespace Givaro

#endif // __GIVARO_hashtable_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
