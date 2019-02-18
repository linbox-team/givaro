// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/bstruct/givhashtable.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Author: T. Gautier
// $Id: givhashtable.h,v 1.3 2011-02-02 16:23:55 briceboyer Exp $
// ==========================================================================
/*! @file givhashtable.h
 * @ingroup bstruct
 * @brief hash table
 */

#ifndef __GIVARO_hashtable_H
#define __GIVARO_hashtable_H

namespace Givaro {


    /*! @brief The class Key.
     * must have :
     * - default cstor
     * - operator== : (const Key&, const Key&) -> int
     * - godel      : void -> int
     * .
     */
    // Generic Key for class T which can be cast to and int
    template<class T>
    class Key {
    public:
        inline Key() {} ;
        inline Key(const T& aa) : a(aa) {} ;
        inline operator T& () { return a ; }
        inline operator const T& () const { return a ; }
        inline int godel() const { return (int)a ; }
        inline int operator== (const Key<T>& K) const { return a == K.a ; }
    private:
        T a ;
    } ;

    //! Hash table
    template<class T, class Key>
    class HashTable {
    public:
        // Size of the table :
        HashTable ( int n = 0 ) ;
        void allocate( int n ) ;
        ~HashTable ( ) ;
        void freeing() ;

        // Insert and remove functions
        void insert     ( const T& item, const Key& k) ;
        void insertLast ( const T& item, const Key& k) ;
        void insertFront( const T& item, const Key& k) ;
        // remove all entry with key k
        void remove     ( const Key& k) ;
        void removeall  ( ) ;

        // Get and "Get and Remove" functions, return 0 if no found
        // select the first item that they found
        int get    ( T& item, const Key& k ) const ;
        int getrmv ( T& item, const Key& k ) ;

        // Returns the number of items of the same key :
        int num_item( const Key& k ) const ;


    public: // ----- internal data structure
        // must be private but is access by Iterator class
        int num ;
        struct _E {
            _E( const T& item, const Key& k)
            : oneitem(item), key(k) {} ;
        public:
            T   oneitem ;
            Key key ;
            _E* next ;
            _E* pred ;
        } ;
        _E** tabH ;  // Head of list
        _E** tabE ;  // End of list


    public:
        // Iterator of all item of a  HashTable :
        class Iterator {
        public:
            Iterator ( HashTable<T,Key>& r)
            : ref(r)
            {
                // search the first non-void pointer in ref.tabH :
                for (e_curr =0; e_curr < ref.num ; e_curr++)
                    if ((curr =ref.tabH[e_curr]) !=0) break ;
            } ;

            operator int() const { return !isend() ; } ;
            int isend() const
            {
                if (e_curr >= ref.num) return 1 ;
                if (e_curr == ref.num-1)
                    if (curr ==0) return 1 ;
                return 0 ;
            }
            void next()
            {
                _E* tmp = curr->next ;
                if (tmp ==0) {
                    int i ;
                    for (i=e_curr+1 ; i<ref.num ; i++)
                    {
                        tmp = ref.tabH[i] ;
                        if (tmp !=0) break ;
                    }
                    if (i >= ref.num) { curr =0 ; e_curr = i ; }
                    else { curr = tmp ; e_curr = i ; }
                }
                else curr = tmp ;
            } ;
            T curr_item() { return curr->oneitem ; } ;
            Key curr_key() { return curr->key ; } ;
        private:
            HashTable<T,Key>& ref ;
            int e_curr ;
            _E* curr ;
        } ;

        // Iterator of all item for a same Key :
        class IteratorKey  {
        public:
            IteratorKey ( HashTable<T,Key>& r, const Key& k )
            : ref(r), key(k)
            {
                e_curr = key.godel() % r.num ;
                curr = ref.tabH[e_curr] ;
            } ;
            operator int() const { return !isend() ; } ;
            int isend() const { return curr == 0 ; }
            void next()
            {
                curr = curr->next ;
                while ( curr != 0 )
                {
                    if (curr->key == key) break ;
                    curr = curr->next ;
                }
            } ;
            T curr_item() { return curr->oneitem ; } ;
            Key curr_key() { return curr->key ; } ;
        private:
            HashTable<T,Key>& ref ;
            Key key ;
            int e_curr ;
            _E* curr ;
        } ;



    } ;

} // namespace Givaro

#include "givaro/givhashtable.inl"

#endif // __GIVARO_hashtable_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
