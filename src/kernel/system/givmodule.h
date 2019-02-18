// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/system/givmodule.h,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givmodule.h,v 1.4 2011-02-02 16:23:56 briceboyer Exp $
// ==========================================================================
/** @file givmodule.h
 * @ingroup system
 * @brief NO DOC
 */
#ifndef __GIVARO_module_init_H
#define __GIVARO_module_init_H

namespace Givaro {

    // -------------------------------------------------------- Fwd declaration

    class GivModule;


    //! GivaroNoInit.
    /**
     * Purpose: used to delay construction of object in the init function
     of a module definition.
     */
    class GivaroNoInit {};


    //! InitAfter.
    //! Purpose: define a precedence relation between two modules.

    class InitAfter {
    public:
        InitAfter( const GivModule& MI );
        static InitAfter Default;
        static InitAfter First;
        static InitAfter Last;
    private:
        InitAfter( int p );
        const GivModule* M;
        int priority;
        friend class GivModule;
        int operator < ( const InitAfter& M ) const;
    };


    //! GivModule.
    //! Purpose: definition of module with precedence relation use to initialize
    //! them between different units compilation.

    class GivModule {
    public:
        typedef void (*ptFuncInit)(int* argc, char** *argv);
        typedef void (*ptFuncEnd)();
        enum {
            MaxPriority  = -100000 ,          // - maximum priority
            MinPriority  = -MaxPriority,      // - minimum priority
            DfltPriority = 0,                 // - default priority
            UndefPriority = MaxPriority-1     // - use to build depedences
        };
        // - Cstor of a module with a priority
        GivModule ( ptFuncInit init, ptFuncEnd end,
                    const int p, const char* n=0 );

        // - Cstor of a module with precedence relation between an other module
        GivModule ( ptFuncInit init, ptFuncEnd end,
                    const InitAfter& M, const char* n=0 );
        ~GivModule ();

    private:
        // - Call by the Givaro::Init and Givaro::End functions
        static void InitApp(int* argc, char***argv);
        static void EndApp();
        friend class GivaroMain;

    private:
        // - Internal data of a module
        int priority;
        InitAfter which;
        ptFuncInit f_init;
        ptFuncEnd f_end;
        const char* name;

        friend class InitAfter;
        static void SortGivModule();
    };


    //! GivModule.
    //! Purpose: definition of object to be initialized after all modules
    //! initialization

    class ObjectInit {
    public:
        // -- when call: link in a global list, then ...
        virtual ~ObjectInit();
        ObjectInit();
        // -- ... call init during the initialization phase
        virtual void objinit() {};
    private:
        ObjectInit* _next;
        friend class GivModule;
    };

} // namespace Givaro
#endif // __GIVARO_module_init_H
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
