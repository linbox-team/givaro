// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/system/givmodule.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givmodule.C,v 1.3 2009-09-17 14:28:23 jgdumas Exp $
// ==========================================================================
// Description:
// Definition of initialization module. Each predefined
// module (a static variable in a compilation unit) is
// initialized after the main ([TN]: by the Givaro::Init).
// Each module can defined one precedence relation with
// and other module: let M1 and M2 are two modules, such that
// in the construction of the associated GivModule static
// variables the priority of M2 is defined as InitAfter(M1),
// then Init method of M2 is call after Init method of M1
// (and in the reverse order for End methods). This leads
// to a block oriented initialization mecanism, where each
// module can be see as a block: the Init method is call
// when enter into the block, and End when go out. Precedent
// relationships are equivalent to block inclusion.

#include "givaro/givmodule.h"
#include "givaro/giverror.h"
#include <iostream>

namespace Givaro {

    InitAfter::InitAfter( const GivModule& MI )
    : M(&MI), priority(GivModule::UndefPriority)
    {}

    InitAfter::InitAfter(int p )
    : M(0), priority(p)
    {}

    int InitAfter::operator < ( const InitAfter& Mi ) const
    {
        int pThis = ( this->M ==0 ? priority : M->priority ) ;
        int pMi   = (  Mi.M  ==0  ? Mi.priority : Mi.M->priority ) ;
        return (pThis < pMi) ;
    }


    InitAfter InitAfter::Default (GivModule::DfltPriority ) ;
    InitAfter InitAfter::First (GivModule::MaxPriority -1 ) ;
    InitAfter InitAfter::Last (GivModule::MinPriority ) ;

    // -- GivModule static object:
    static GivModule* All[1024] ;
    static int SortedAll[1024] ;
    static int counter = 0 ;

    // -- ObjectInit static object:
    static ObjectInit* head = 0;

    ObjectInit::ObjectInit( )
    {
#ifdef __GIVARO_DEBUG
        std::cout << "[GivModule] ObjectInit::cstor " << (void*)head << std::endl ;
        std::cout << "[GivModule] ObjectInit::cstor curr:" << (void*)this << std::endl ;
#endif
        // - link the new object:
        _next = head; head = this ;
#ifdef __GIVARO_DEBUG
        std::cout << "[GivModule] ObjectInit::cstor " << (void*)head << std::endl ;
#endif
    }

    ObjectInit::~ObjectInit(){}

    GivModule::GivModule ( ptFuncInit init, ptFuncEnd end, const int p, const char* n)
    : priority(UndefPriority), which(p), f_init(init), f_end(end), name(n)
    {
        All[counter++] = this ;
    }

    GivModule::GivModule ( ptFuncInit init, ptFuncEnd end, const InitAfter& M, const char* n)
    : priority(UndefPriority), which(M), f_init(init), f_end(end), name(n)
    {
        All[counter++] = this ;
    }

    void GivModule::SortGivModule()
    {
#ifdef __GIVARO_DEBUG
        {
            for (int j=0; j<counter; j++)
                std::cout << j << ':' << All[j]->which.priority << ' ' << All[j]->name << std::endl ;
            std::cout <<std::endl ;
        }
#endif

        int curr ;

        // -- Set priority to each module

        // - Set undefined priority of constraints
        for (curr=0; curr<counter; curr++) {
            if (All[curr]->which.priority == GivModule::UndefPriority)
            {
                All[curr]->which.priority = GivModule::MinPriority ;
                All[curr]->which.M = 0 ;
                SortedAll[curr] = curr ;
            }
        }

        // - Set priority of modules
        int isundef ; // - ==1 if exists undef priority
        int cpt = 0 ; // - count the number of tentative to resolve constraint
        do {
            isundef = 0 ;
            for (curr=0; curr<counter; curr++)
                if (All[curr]->priority == GivModule::UndefPriority) {
                    if (All[curr]->which.M ==0)
                        // -- Set the priority field with those of which
                        All[curr]->priority = All[curr]->which.priority+1 ;
                    else {
                        // -- Set the priority field with the value pointed by which
                        All[curr]->priority = All[curr]->which.M->priority+1 ;
                        if (All[curr]->priority == GivModule::UndefPriority) isundef = 1 ;
                    }
                }
        } while ((isundef) && (cpt <= counter*counter));

        if (cpt > counter*counter)
            throw GivError("*** Can't resolve constraint of initialization of modules");

        // -- Sort Initialization Module by priority
        SortedAll[0] = 0 ; // - the first module is ordered
        SortedAll[1] = 1 ; // - the first module is ordered
        for (curr=1; curr<counter; curr++)
        {
            GivModule* currM = All[curr] ;
#ifdef __GIVARO_DEBUG
            std::cout << " Search for : " << currM->name << std::endl ;
#endif

            // -- Search the first GivModule with strict less priority
            int k;
            for (k=0; k<curr; k++)
                if ( currM->priority < All[SortedAll[k]]->priority) break ;

#ifdef __GIVARO_DEBUG
            std::cout << currM->name << " inserted at : " << k << std::endl ;
#endif

            // -- either k == curr, and this is inserted at position curr
            // -- either k != curr and this is inserted at position k
            if (k == curr) SortedAll[curr] = curr ;
            else {
                // -- SortedAll is shift by one from curr downto k
                for (int i=curr-1; i>=k ; i--) SortedAll[i+1] = SortedAll[i] ;
                // -- this is inserted :
                SortedAll[k] = curr ;
            }
#ifdef __GIVARO_DEBUG
            {
                for (int j=0; j<=curr; j++) std::cout << All[SortedAll[j]]->name << ',' ;
                std::cout <<std::endl << std::endl ;
            }
#endif
        }

#ifdef __GIVARO_DEBUG
        {
            for (int j=0; j<counter; j++) std::cout << j << ':' << All[SortedAll[j]]->name << std::endl ;
            std::cout <<std::endl ;
        }
#endif
    }

    GivModule::~GivModule ()
    {
    }

    void GivModule::InitApp(int* argc, char***argv)
    {
        // -- Computation of a topological sort from All[]
        SortGivModule() ;

        // -- Initialization of the modules
        for (int i=0; i<counter; i++)
        {
            GivModule* curr = All[SortedAll[i]] ;
#ifdef __GIVARO_DEBUG
            if (curr->name !=0) std::cout << curr->name << " initializing..."  ;
#endif
            if (curr->f_init !=0) curr->f_init(argc, argv) ;
#ifdef __GIVARO_DEBUG
            if (curr->name !=0) std::cout << "done !" << std::endl ;
#endif
        }

        // -- Init of object:
        ObjectInit* curr = head;
#ifdef __GIVARO_DEBUG
        std::cout << "[GivModule] ObjectInit::Head: " << (void*)head << std::endl;
#endif
        while (curr !=0) {
            curr->objinit();
            curr = curr->_next;
#ifdef __GIVARO_DEBUG
            std::cout << "[GivModule] ObjectInit::next one: " << (void*)curr << std::endl;
#endif
            if (curr ==head) break ; // -- on MacOS I can make circular list !!! (shared lib)
        }
    }

    void GivModule::EndApp()
    {
        for (int i=counter-1; i>=0; i--)
        {
            GivModule* curr = All[SortedAll[i]] ;
#ifdef __GIVARO_DEBUG
            if (curr->name !=0) std::cout << curr->name << " freeing..."  ;
#endif
            if (curr->f_end !=0) curr->f_end() ;
#ifdef __GIVARO_DEBUG
            if (curr->name !=0) std::cout << "done !" << std::endl ;
#endif
        }
    }

} // namespace Givaro
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
