// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/memory/givaromm.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Author: T. Gautier
// $Id: givaromm.h,v 1.7 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================

/** @file givaromm.h
 * @ingroup memory
 * @brief Memory management in Givaro
 * two memory managers:
 * - the first one handle a set on free-list of blocs ;
 * - the second one implement a reference mecanism on the bloc.
 * .
 * The latter used method of the former.
 */
#ifndef __GIVARO_mm_H
#define __GIVARO_mm_H

#include <cstring>            // gcc 4.3
#include <stddef.h>           // size_t definition
#include <iostream>         // ostream definition
#include <new>
#include "givaro/givmodule.h"
#include "givaro/givbasictype.h"
#ifdef __GIVARO_DEBUG
#include "givaro/giverror.h"
#endif

#include <givaro/givconfig.h>
namespace Givaro {

    // ==================================================================== //

    //!  Static informations of memory allocation
    class GivMMInfo {
    public:
        GivMMInfo();
        ~GivMMInfo();
        size_t physalloc; // size in bytes of physical allocated bloc
        size_t logalloc;  // size in bytes of "logical" allocated bloc
        size_t  sizetab;    // length of next arrays
        size_t* tabbloc;  // size of all the blocs
        size_t* tablog;     // number of each logical allocated bloc
        size_t* tabphy;     // number of each physical allocated bloc
        std::ostream& print( std::ostream& so ) const;
    };

    //! IO
    inline std::ostream& operator<<( std::ostream& o, const GivMMInfo& T)
    { return T.print(o); }

    // --------------------------------------------------------
    //! Data structure of a bloc.
    //! Each bloc in TabFree[id] has a data field of size TabSize[id]
    class BlocFreeList {
        union header {
            int index ;             // - index in free list
            BlocFreeList* nextfree; // - pointer to the next free bloc (of the same size)
            double dummy; 			// - here to force alignment on data on 8bytes boundary
        } u;
        int64_t data[1];

        // -- Array of list of free bloc
        static BlocFreeList* TabFree[];
        static const size_t TabSize[] ;
        static const int lenTables;
        static int search_binary( size_t sz ) ;

        friend class GivMMInfo;
        friend class GivMMFreeList;
        friend class GivMMRefCount;
    };


    // --------------------------------------------------------
    //! Implementation of a memory manager with free-lists.
    //! All members are static methods.
    class GivMMFreeList {
    public:

        // -- Free all internal data structures
        static void Destroy();

        // -- Allocation of a new bloc.
        // These two next functions returns a pointer or 0.
        // The frequently required bloc is treated in a special version
        // in order to speed-up the allocation.
        static BlocFreeList* _allocate (const size_t sz);
        inline static void* allocate (const size_t sz)
        {
#ifdef __GIVARO_DEBUG
            if (sz ==0) return 0 ;
#endif
            size_t index;
            BlocFreeList* tmp;
            if ((sz <= 32) && ((tmp=BlocFreeList::TabFree[index =sz-1]) !=0)) {
                BlocFreeList::TabFree[index] = tmp->u.nextfree;
                tmp->u.index = (int)index;
#ifdef GIVARO_STATMEM
                tablog[index] ++; logalloc += BlocFreeList::TabSize[index];
#endif
                return (void*) tmp->data;
            }
            tmp = _allocate(sz);
            return tmp->data ;
        }

        // -- Reallocation of a bloc.
        // The function returns a pointer to the possibly new
        // bloc. If the bloc must be moved then only oldsize are
        // recopied.
        static void* resize (void* p, const size_t oldsize, const size_t newsize);

        // -- Free the bloc pointed by p and allocated by the manager GivMMFreeList.
        inline static void desallocate(void* p, const size_t = 0)
        {
            if (p==0) return ;
            BlocFreeList* tmp = reinterpret_cast<BlocFreeList*>(((char*)p) -
                                                                (sizeof(BlocFreeList)-sizeof(int64_t)));
            int index = tmp->u.index;
#ifdef __GIVARO_DEBUG
            if ((index <0) || (index >= BlocFreeList::lenTables))
                throw GivError("[GivMMFreeList::desallocate]: bad pointer.");
#endif
            tmp->u.nextfree = BlocFreeList::TabFree[index];
            BlocFreeList::TabFree[index] = tmp;
        }

        // -- Recopy size bytes pointed by src into bytes pointer by dest
        static void memcpy( void* dest, const void* src, const size_t size );

        // -- Returns some memory usage informations
        static const GivMMInfo& Usage();

        //private:
        // -- bug (on CW10): must comment the private part to initialize next field
        static GivMMInfo info;

#ifdef GIVARO_STATMEM
        friend class GivMMInfo;
        friend class GivMMRefCount;
        static size_t& physalloc;  // total amount of physical allocated bloc
        static size_t& logalloc;   // total amoun of "logical" allocated bloc
        static size_t*& tablog;      // number of each logical allocated bloc
        static size_t*& tabphy;      // number of each physical allocated bloc
#endif

        // -- Initialization module
    public:
        static GivModule Module;
    private:
        static void Init(int* argc, char***argv);
        static void End();
        friend class GivModule;
    };


    //! Memory management with reference counter on allocated data.
    //! The memory manager uses the BlocFreeList data structure
    //! and stores the refcounter in the field data[0]
    class GivMMRefCount {
    public:
        // -- Allocation of a new bloc of size at least s
        inline static void* allocate (const size_t s)
        {

#ifdef __GIVARO_DEBUG
            if (s ==0) return 0 ;
#endif
            int index;
            BlocFreeList* tmp;
            size_t sz = s + sizeof(int64_t);
            if ((sz <= 32) && ((tmp=BlocFreeList::TabFree[index =int(sz-1)]) !=0)) {
                BlocFreeList::TabFree[index] = tmp->u.nextfree;
                tmp->u.index = index;
#ifdef GIVARO_STATMEM
                GivMMFreeList::tablog[index] ++; GivMMFreeList::logalloc += BlocFreeList::TabSize[index];
#endif
                tmp->data[0] = 1 ;
                return (void*) &(tmp->data[1]);
            }
            tmp = GivMMFreeList::_allocate(sz);
            tmp->data[0] = 1 ;


            return (void*) &(tmp->data[1]);
        }

        // -- Reallocation of a bloc. See description in GivMMFreeList's class.
        // Here, if ref count on p is >1 then a new bloc is allocated.
        static void* resize (void* p, const size_t oldsize, const size_t newsize);

        // -- Free the bloc allocated by the manager GivMMRefCount.
        inline static void desallocate(void* p, const size_t = 0)
        {
            if (p==0) return ;
            BlocFreeList* tmp = reinterpret_cast<BlocFreeList*>((char *) p - sizeof(BlocFreeList) ) ;
            if ( --(tmp->data[0]) ==0) {
                int index = tmp->u.index;
#ifdef __GIVARO_DEBUG
                if ((index <0) || (index >= BlocFreeList::lenTables))
                    GivError::throw_error(GivError("[GivMMRefCount::desallocate]: bad pointer."));
#endif
#ifdef GIVARO_JGD
                if ((index <0) || (index >= BlocFreeList::lenTables))
                    cerr << "[GivMMRefCount::desallocate]: bad pointer with index " << index << ", doing nothing ..." << endl;
                else {
#endif
                    tmp->u.nextfree = BlocFreeList::TabFree[index];
                    BlocFreeList::TabFree[index] = tmp;
#ifdef GIVARO_JGD
                }
#endif
            }
        }

        // -- Assignement of pointer:
        inline static void* assign (void** dest, void*src)
        {
            if (src == *dest) return *dest ;
            if (*dest !=0) GivMMRefCount::desallocate( *dest );
            if (src ==0) return *dest=src;
            BlocFreeList* s = reinterpret_cast<BlocFreeList*>(((char*)src) - sizeof(BlocFreeList));
#ifdef __GIVARO_DEBUG
            if ((s->u.index <0) || (s->u.index >= BlocFreeList::lenTables))
                GivError::throw_error(GivError("[GivMMRefCount::assign]: bad pointer 'src'."));
#endif
            ++(s->data[0]);
            return *dest = src ;
        }

        // -- Increment, decrement and getvalue of the refcount.
        // Returns the value after increment/decrment.
        inline static int incrc(void* p)
        {
            if (p ==0) return 0;
            BlocFreeList* bp = reinterpret_cast<BlocFreeList*>(((char*)p) - sizeof(BlocFreeList));
#ifdef __GIVARO_DEBUG
            if ((bp->u.index <0) || (bp->u.index >= BlocFreeList::lenTables))
                throw GivError("[GivMMRefCount::incrc]: bad pointer.");
#endif
            return int(++(bp->data[0]));
        }
        inline static int decrc(void* p)
        {
            if (p ==0) return 0;
            BlocFreeList* bp = reinterpret_cast<BlocFreeList*>(((char*)p) - sizeof(BlocFreeList));
#ifdef __GIVARO_DEBUG
            if ((bp->u.index <0) || (bp->u.index >= BlocFreeList::lenTables))
                throw GivError("[GivMMRefCount::incrc]: bad pointer.");
#endif
            return int(--(bp->data[0]));
        }
        inline static int getrc(void* p)
        {
            if (p ==0) return 0;
            BlocFreeList* bp = reinterpret_cast<BlocFreeList*>(((char*)p) - sizeof(BlocFreeList));
#ifdef __GIVARO_DEBUG
            if ((bp->u.index <0) || (bp->u.index >= BlocFreeList::lenTables))
                throw GivError("[GivMMRefCount::incrc]: bad pointer.");
#endif
            return int(bp->data[0]);
        }
    };


    //! Memory manager that allocates array of object of type T for
    template<class T>
    class GivaroMM {
    public:
        typedef T* ptT;
        // -- Allocation of a new bloc contains s unintialized elements
        static T* allocate (const size_t s);
        // -- Desallocation of bloc
        static void desallocate ( T* bloc, const size_t sz =0 );
        // -- Initialize each element pointed by bloc with value V
        static void initone( T*p, const T& V = T());
        static void initialize(T* bloc, const size_t s, const T& V =T());
        // -- Call destructor on each elements pointed by bloc
        static void destroy(T* bloc, const size_t s);
    };


    //! @bug implem does not belong here
    //@{
    template<class T>
    inline T* GivaroMM<T>::allocate (const size_t s)
    { return (T*)GivMMFreeList::allocate(s*sizeof(T)); }
    template<class T>
    inline void GivaroMM<T>::desallocate(T* bloc, const size_t sz)
    { GivMMFreeList::desallocate((void*)bloc,sz); }
    template<class T>
    inline void GivaroMM<T>::initone( T* p, const T& V)
    { new (p) T(V); }
    template<class T>
    inline void GivaroMM<T>::initialize(T* bloc, const size_t s, const T& V)
    { for (size_t i=0; i<s; i++) GivaroMM<T>::initone(&(bloc[i]), V); }
    template<class T>
    inline void GivaroMM<T>::destroy(T* bloc, const size_t s)
    { for (size_t i=0; i<s; i++) bloc[i].~T(); }
    //@}

    // -- specialized version: for basic C++ type and pointer on this type
#define GIVARO_MM_SPECIALIZED(TYPE) \
    template<> class GivaroMM<TYPE> {\
    public: \
            typedef TYPE* ptType;\
        typedef TYPE Type;\
        static ptType allocate (const size_t s);\
        static void desallocate ( ptType bloc, const size_t sz =0 );\
        static void initone( ptType p, const Type V = 0);\
        static void initialize(ptType bloc, const size_t s, const Type V =0);\
        static void destroy(GivaroMM<TYPE>::ptType bloc, const size_t s);\
    };\
    inline GivaroMM<TYPE>::ptType GivaroMM<TYPE>::allocate (const size_t s)\
    { return (GivaroMM<TYPE>::ptType)GivMMFreeList::allocate(s*sizeof(GivaroMM<TYPE>::Type)); }\
    inline void GivaroMM<TYPE>::desallocate(GivaroMM<TYPE>::ptType bloc, const size_t sz)\
    { GivMMFreeList::desallocate((void*)bloc,sz); }\
    inline void GivaroMM<TYPE>::initone( GivaroMM<TYPE>::ptType p, const GivaroMM<TYPE>::Type v) \
    { *p = v; }\
    inline void GivaroMM<TYPE>::initialize(GivaroMM<TYPE>::ptType bloc, const size_t s, const GivaroMM<TYPE>::Type V)\
    { for (size_t i=0; i<s; i++) bloc[i] = V; }\
    inline void GivaroMM<TYPE>::destroy(GivaroMM<TYPE>::ptType bloc, const size_t s){}

#if 0
    GIVARO_MM_SPECIALIZED(char)
    GIVARO_MM_SPECIALIZED(short)
    GIVARO_MM_SPECIALIZED(int)
    GIVARO_MM_SPECIALIZED(long)
    GIVARO_MM_SPECIALIZED(float)
    GIVARO_MM_SPECIALIZED(double)

    GIVARO_MM_SPECIALIZED(unsigned char)
    GIVARO_MM_SPECIALIZED(unsigned short)
    GIVARO_MM_SPECIALIZED(unsigned int)
    GIVARO_MM_SPECIALIZED(unsigned long)

    GIVARO_MM_SPECIALIZED(char*)
    GIVARO_MM_SPECIALIZED(short*)
    GIVARO_MM_SPECIALIZED(int*)
    GIVARO_MM_SPECIALIZED(long*)
    GIVARO_MM_SPECIALIZED(float*)
    GIVARO_MM_SPECIALIZED(double*)

    GIVARO_MM_SPECIALIZED(unsigned char*)
    GIVARO_MM_SPECIALIZED(unsigned short*)
    GIVARO_MM_SPECIALIZED(unsigned int*)
    GIVARO_MM_SPECIALIZED(unsigned long*)
#endif

} // namespace Givaro

#endif // __GIVARO_mm_H

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
