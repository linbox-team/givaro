// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/bstruct/givlist0.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Author: T. Gautier
// $Id: givlist0.h,v 1.3 2011-02-02 16:23:55 briceboyer Exp $
// ==========================================================================

/*! @file givlist0.h
 * @ingroup bstruct
 * @brief List of type T with double link and various insert/get/rmv method.
 * Used reference counting on each node of the list.
 */

#ifndef __GIVARO_list0_H
#define __GIVARO_list0_H

#include "givaro/givbasictype.h"
#include "givaro/giverror.h"


namespace Givaro {

    //! ListO
    template <class T>
    class List0 {
    public :
        // -- Default cstor : ctsor of s size array
        List0 ();

        // -- Recopy cstor : logical copy
        List0 ( const List0<T>& p ) ;

        // -- Destructor
        ~List0 ();

        // -- Destroy of the list
        void destroy( );

        // -- Physical copy operator
        List0<T>& copy(const List0<T>& p);

        // -- Logical recopy operator
        List0<T>& logcopy(const List0<T>& p);

        // -- Logical recopy
        List0<T>& operator= (const List0<T>& p);

        // -- Return 1 if the list is empty or 0
        int is_empty() const;

        // -- Return the occuped size of the list
        size_t size() const;

        // -- insertlast: add at the end of the list
        void insertlast( const T& );
        // -- insertfirst: add at the end of the list
        void insertfirst( const T& );
        // -- The four next methods return 1 is val is set (list not empty)
        // -- getlast: get the last item of the list, don't remove the item
        int getlast(T& item) const;
        // -- getfirst: get the first item of the list, don't remove the item
        int getfirst(T& item) const;
        // -- getrmvlast: get and remove the last item of the list
        int getrmvlast(T& item);
        // -- getrmvfirst: get and remove the first item of the list
        int getrmvfirst(T& item);

    public:
        struct node {
            ~node() {
                GivaroMM<T>::destroy(item,1);
                GivaroMM<T>::desallocate(item,1);
                GivaroMM<int>::desallocate(cnt,1);
            }
            node* next;
            node* prev;
            T* item;
            int* cnt;     // reference counter on item
        };

        friend class Iterator;
        class Iterator {
        public:
            Iterator( List0<T>& lst )
            : currnode(lst._head)
            { }
            int is_empty() const { return currnode =0; };
            void prev() { currnode = currnode->prev;};
            void next() { currnode = currnode->next;};
            T& get_item() { return *currnode->item;};
            const T& get_item() const { return *currnode->item; };
        private:
            List0<T>::node* currnode;
        };


    protected :  // --------------------- Public Internal representation
        size_t _size;  // actual size of the list
        node*  _head;   // head of the list
        node*  _queue;  // queue of the list
    };

} // namespace Givaro

#include "givaro/givlist0.inl"

#endif // __GIVARO_list0_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
