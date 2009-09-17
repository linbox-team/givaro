#ifndef _LIST0_H_
#define _LIST0_H_
// ========================================================================== 
// $Source: /var/lib/cvs/Givaro/src/kernel/bstruct/givlist0.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Author: T. Gautier
// $Id: givlist0.h,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ========================================================================== 
// Description:
// List of type T with double link and various insert/get/rmv method.
// Used reference counting on each node of the list.

#include "givaro/givbasictype.h"
#include "givaro/giverror.h"

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

#include "givaro/givlist0.inl"

#endif