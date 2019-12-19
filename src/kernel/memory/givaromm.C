// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/kernel/memory/givaromm.C,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: M. Samama, T. Gautier
// $Id: givaromm.C,v 1.5 2010-10-15 13:43:48 briceboyer Exp $
// ==========================================================================
// Description:

#include <iostream>
#include <fstream>
#include <string>
#include "givaro/givmodule.h"
#include "givaro/givaromm.h"
#include "givaro/giverror.h"

namespace Givaro {

#ifdef GIVARO_MAPMEM
    ofstream memlog;
#endif

    // -- Current information of the memory manager:
    GivMMInfo GivMMFreeList::info;

#ifdef GIVARO_STATMEM
    size_t& GivMMFreeList::physalloc = GivMMFreeList::info.physalloc;
    size_t& GivMMFreeList::logalloc = GivMMFreeList::info.logalloc;
    size_t*& GivMMFreeList::tablog = GivMMFreeList::info.tablog;
    size_t*& GivMMFreeList::tabphy = GivMMFreeList::info.tabphy;
#endif


    // -- for debug #define

    // ======================================================================= //
    // Short description of the Givaro memory mechanism:
    // each request of sz bytes returns the smallest bloc that contains
    // at least sz bytes. Each bloc begins by a private data used to store
    // size information or a pointer. All free blocs are linked in list.
    // The allocation algorithm search in the corresponding list a free bloc,
    // before call the system "malloc" routine.
    // 512 different sizes of bloc are predefined.
    // The maximum bloc size is 8054880 bytes.
    // For higher bloc, the manager directly call malloc/free of the underlaying system.

    const int  BlocFreeList::lenTables = 512 ;
    BlocFreeList* BlocFreeList::TabFree[512];

    const size_t BlocFreeList::TabSize[] = { // -- array of different sizes of bloc (bytes)
        // --- 8 values per line (I hope !)
        1, 2, 3, 4, 5, 6, 7, 8,
        9, 10, 11, 12, 13, 14, 15, 16,
        17, 18, 19, 20, 21, 22, 23, 24,
        25, 26, 27, 28, 29, 30, 31, 32,
        64, 96, 128, 160, 192, 224, 256, 288,
        320, 352, 384, 416, 448, 480, 512, 544,
        576, 608, 640, 672, 704, 736, 768, 800,
        832, 864, 896, 928, 960, 992, 1024, 1056,
        1088, 1120, 1152, 1184, 1216, 1248, 1280, 1312,
        1344, 1376, 1408, 1440, 1472, 1504, 1536, 1568,
        1600, 1632, 1664, 1696, 1728, 1760, 1792, 1824,
        1856, 1920, 1952, 1984, 2016, 2048, 2112, 2144,
        2176, 2240, 2272, 2336, 2368, 2400, 2464, 2496,
        2560, 2624, 2656, 2720, 2784, 2816, 2880, 2944,
        3008, 3072, 3104, 3168, 3232, 3296, 3360, 3456,
        3520, 3584, 3648, 3712, 3808, 3872, 3936, 4032,
        4096, 4192, 4288, 4352, 4448, 4544, 4640, 4704,
        4800, 4896, 4992, 5120, 5216, 5312, 5408, 5536,
        5632, 5760, 5856, 5984, 6112, 6208, 6336, 6464,
        6592, 6720, 6880, 7008, 7136, 7296, 7424, 7584,
        7744, 7872, 8032, 8192, 8384, 8544, 8704, 8864,
        9056, 9248, 9408, 9600, 9792, 9984, 10208, 10400,
        10624, 10816, 11040, 11264, 11488, 11712, 11936, 12192,
        12416, 12672, 12928, 13184, 13440, 13728, 13984, 14272,
        14560, 14848, 15136, 15456, 15744, 16064, 16384, 16704,
        17056, 17408, 17728, 18112, 18464, 18816, 19200, 19584,
        19968, 20384, 20800, 21184, 21632, 22048, 22496, 22944,
        23392, 23872, 24352, 24832, 25344, 25824, 26368, 26880,
        27424, 27968, 28512, 29088, 29664, 30272, 30880, 31488,
        32128, 32768, 33408, 34080, 34784, 35456, 36160, 36896,
        37632, 38400, 39168, 39936, 40736, 41536, 42368, 43232,
        44096, 44960, 45888, 46784, 47712, 48672, 49664, 50656,
        51648, 52704, 53760, 54816, 55904, 57024, 58176, 59328,
        60512, 61728, 62976, 64224, 65504, 66816, 68160, 69536,
        70912, 72320, 73760, 75264, 76768, 78304, 79840, 81440,
        83072, 84736, 86432, 88160, 89920, 91712, 93568, 95424,
        97344, 99296, 101280, 103296, 105376, 107456, 109632, 111808,
        114048, 116320, 118656, 121024, 123456, 125920, 128448, 131008,
        133632, 136288, 139008, 141792, 144640, 147520, 150464, 153472,
        156544, 159680, 162880, 166144, 169472, 172864, 176320, 179840,
        183424, 187104, 190848, 194656, 198560, 202528, 206560, 210688,
        214912, 219200, 223584, 228064, 232640, 237280, 242016, 246848,
        251808, 256832, 261984, 267200, 272544, 278016, 283552, 289248,
        295008, 300928, 306944, 313088, 319328, 325728, 332224, 338880,
        345664, 352576, 359616, 366816, 374144, 381632, 389280, 397056,
        404992, 413088, 421344, 429792, 438368, 447136, 456096, 465216,
        474528, 484000, 493696, 503552, 513632, 523904, 534368, 545056,
        555968, 567072, 578432, 589984, 601792, 613824, 626112, 638624,
        651392, 664416, 677728, 691264, 705088, 719200, 733568, 748256,
        763200, 778464, 794048, 809920, 826112, 842656, 859488, 876704,
        894240, 912096, 930336, 948960, 967936, 987296, 1007040, 1027168,
        1047712, 1068672, 1090048, 1111840, 1134080, 1156768, 1179904, 1203488,
        1227584, 1252128, 1277152, 1302720, 1328768, 1355328, 1382432, 1410080,
        1438304, 1467072, 1496416, 1526336, 1556864, 1588000, 1619744, 1652160,
        1685184, 1718880, 1753280, 1788320, 1824096, 1860576, 1897792, 1935744,
        1974464, 2013952, 2054240, 2095328, 2137216, 2179968, 2223552, 2268032,
        2313408, 2359648, 2406848, 2455008, 2504096, 2554176, 2605248, 2657376,
        2710496, 2764704, 2820000, 2876416, 2933952, 2992608, 3052480, 3113536,
        3175776, 3239296, 3304096, 3370176, 3437568, 3506336, 3576448, 3647968,
        3720928, 3795360, 3871264, 3948704, 4027680, 4108224, 4190400, 4274176,
        4359680, 4446880, 4535808, 4626528, 4719040, 4813440, 4909696, 5007904,
        5108064, 5210208, 5314400, 5420704, 5529120, 5639712, 5752480, 5867552,
        5984896, 6104576, 6226688, 6351200, 6478240, 6607808, 6739968, 6874752,
        7012256, 7152512, 7295552, 7441472, 7590272, 7742080, 7896928, 8054880
    };



    // ======================================================================= //
    // Class GivMMFreeList : allocation/destroy methods

    inline int BlocFreeList::search_binary( size_t sz )
    {
        if (sz <= 32)
            return int(sz-1);
        int max = BlocFreeList::lenTables-1; // -- last element in TabSize
        if (sz > BlocFreeList::TabSize[max])
            throw GivError("[GivaroMM]: unable to allocate this size of memory");
        int min = 0;
        int med;
        unsigned int curr;
        med = 8; // may be value < TabSize[8]
        do {
            curr = (unsigned int) BlocFreeList::TabSize[med];
            if (curr == sz) return med;
            if (curr < sz) { min = med; }
            else { max = med; }
            med = (max+min)>>1; // /2
        } while  ( min != med );
        return max;
    }


    // // assume 4 bytes in a long:
    // #define ALIGN(s) (s>>2 + 1)  // may one more ?

    BlocFreeList* GivMMFreeList::_allocate (const size_t s)
    {
        int index = BlocFreeList::search_binary( s );
#ifdef __GIVARO_DEBUG
        if ((index <0) || (index >= BlocFreeList::lenTables))
            throw GivError("[GivMMFreeList::_allocate]: index error, invalid bloc size.");
#endif
        BlocFreeList* tmp;
        if (BlocFreeList::TabFree[index] !=0) {
            tmp = BlocFreeList::TabFree[index];
            BlocFreeList::TabFree[index] = tmp->u.nextfree;
        }
        else {
            tmp = (BlocFreeList*) malloc( (BlocFreeList::TabSize[index]+sizeof(BlocFreeList)-sizeof(int64_t)) );

#ifdef GIVARO_STATMEM
            tabphy[index] ++; physalloc += BlocFreeList::TabSize[index];
#endif
        }
#ifdef GIVARO_STATMEM
        tablog[index] ++; logalloc += BlocFreeList::TabSize[index];
#endif
        tmp->u.index = index;
#ifdef GIVARO_MAPMEM
        memlog << "alloc: in:" << (void*) tmp << ", user:" << (void*)tmp->data << std::endl;
#endif
        return tmp;
    }


    void* GivMMFreeList::resize (void* src, const size_t oldsize, const size_t newsize)
    {
        if (src ==0) return _allocate(newsize) ;
        if (newsize <= oldsize) return src;
        BlocFreeList* tmp = reinterpret_cast<BlocFreeList*>(((char*)src)-sizeof(BlocFreeList)+sizeof(int64_t));
#ifdef __GIVARO_DEBUG
        if ((tmp->u.index <0) || (tmp->u.index >= BlocFreeList::lenTables))
            throw GivError("[GivMMFreeList::resize]: bad pointer 'src'");
#endif
#ifdef GIVARO_MAPMEM
        memlog << "reall: in:" << (void*) tmp << ", user:" << (void*)src << std::endl;
#endif
        int index = tmp->u.index;
        if (BlocFreeList::TabSize[index] >= newsize) return src;
        tmp = GivMMFreeList::_allocate( newsize );
        if (oldsize !=0) ::memcpy( tmp->data, src, oldsize );
        return tmp->data;
    }

    void GivMMFreeList::memcpy( void* dest, const void* src, const size_t size )
    {
        BlocFreeList* tmp1 = reinterpret_cast<BlocFreeList*>(((char*)dest) - sizeof(BlocFreeList)+sizeof(int64_t));
        void * toto = const_cast<void*>( src );
        BlocFreeList* tmp2 = reinterpret_cast<BlocFreeList*>(((char*)toto) - sizeof(BlocFreeList)+sizeof(int64_t));
#ifdef __GIVARO_DEBUG
        if ((tmp1->u.index <0) || (tmp1->u.index >= BlocFreeList::lenTables))
            throw GivError("[GivMMFreeList::memcpy]: bad pointer 'dest'");
        if ((tmp2->u.index <0) || (tmp2->u.index >= BlocFreeList::lenTables))
            throw GivError("[GivMMFreeList::memcpy]: bad pointer 'src'");
#endif
        ::memcpy ( tmp1->data, tmp2->data, size );
    }

    void GivMMFreeList::Destroy()
    {
        for (int i=0; i<BlocFreeList::lenTables; i++)
        {
            BlocFreeList* curr = BlocFreeList::TabFree[i];
            while (curr !=0) {
                BlocFreeList* tmp = curr->u.nextfree;
                free((char*)curr);
                curr = tmp;
            }
        }
    }

    void* GivMMRefCount::resize (void* p, const size_t oldsize, const size_t newsize )
    {

        if (p ==0)
            return &(GivMMFreeList::_allocate(newsize+sizeof(int64_t))->data[1]) ;


        BlocFreeList* tmp = reinterpret_cast<BlocFreeList*>(((char*)p)-sizeof(BlocFreeList));
#ifdef __GIVARO_DEBUG
        if ((tmp->u.index <0) || (tmp->u.index >= BlocFreeList::lenTables))
            throw GivError("[GivMMRefCount::resize]: bad pointer");
#endif
        if (tmp->data[0] ==1) { // -- one pointer on the bloc, use standard optimization
            if (newsize <= oldsize) return p;
#ifdef GIVARO_MAPMEM
            memlog << "reall: in:" << (void*) tmp << ", user:" << (void*)p << std::endl;
#endif
            int index = tmp->u.index;
            if (BlocFreeList::TabSize[index] >= sizeof(int64_t)+newsize) return p;
            GivMMRefCount::desallocate(p);
        }
        else --(tmp->data[0]);  // -- two pointer on the bloc:
        tmp = GivMMFreeList::_allocate( newsize+sizeof(int64_t) );
        tmp->data[0] = 1 ;
        if (oldsize !=0) {
            if (newsize <= oldsize) ::memcpy( &(tmp->data[1]), p, newsize );
            else ::memcpy( &(tmp->data[1]), p, oldsize );
        }
        return &(tmp->data[1]);
    }


    // -------- Static informations of memory allocation
    GivMMInfo::GivMMInfo()
    {
        tabbloc = ::new size_t[BlocFreeList::lenTables];
        tablog = ::new size_t[BlocFreeList::lenTables];
        tabphy = ::new size_t[BlocFreeList::lenTables];
        sizetab = BlocFreeList::lenTables;
        for (size_t i=0; i<sizetab; i++) {
            tabbloc[i] = BlocFreeList::TabSize[i];
            tabphy[i] = 0;
            tablog[i] = 0;
        }
    }

    GivMMInfo::~GivMMInfo()
    {
        ::delete[] tabbloc;
        ::delete[] tablog;
        ::delete[] tabphy;
    }

    std::ostream& GivMMInfo::print( std::ostream& so ) const
    {
        so << "--- Memory usage" << std::endl;
        so << "- physical allocated memory (in bytes):" << physalloc << std::endl;
        so << "- logical  allocated memory (in bytes):" << logalloc << std::endl;
        so << "- details for each bloc size:\n";
        so.width(7); so << "index" << ": "; so.width(9); so << "bytes" << " | ";
        so.width(9); so << "#physical" << " | "; so.width(9); so << "#logical" << std::endl;
        for (size_t i=0; i<sizetab; i++) {
            if (tabphy[i] !=0) {
                so.width(7); so << i << ": ";
                so.width(9); so << tabbloc[i] << " | " ;
                so.width(9); so << tabphy[i] << " | ";
                so.width(9); so << tablog[i] << std::endl;
            }
        }
        return so;
    }

    // -- Returns some memory usage informations
    const GivMMInfo& GivMMFreeList::Usage()
    {
        return GivMMFreeList::info;
    }



    // ------ Initialization module:
    GivModule GivMMFreeList::Module (GivMMFreeList::Init,
                                     GivMMFreeList::End,
                                     GivModule::MaxPriority,
                                     "Givaro Memory Manager") ;

    void GivMMFreeList::Init(int* , char***)
    {
#ifdef GIVARO_MAPMEM
        memlog.open("mem.log", ios::out);
#endif
        for (int i=0; i<BlocFreeList::lenTables; i++)
            BlocFreeList::TabFree[i] = 0;
    }

    void GivMMFreeList::End()
    {
#ifdef GIVARO_MAPMEM
        memlog.close();
#endif
        Destroy();

    }

    /*
       \documentstyle[12pt,a4]{article}
       \pagestyle{empty}
       \title{Mousse au chocolat}
       \date{Pour 5 personnes}
       \textwidth 12cm
       \begin{document}
       \maketitle
       \thispagestyle{empty}
       \vspace*{3cm}

       {\it \noindent {\large Ingr\'edients} \\
       250 g chocolat noir \\
       5 oeufs \\
       80 g beurre (ou moins)\\
       50 g sucre\\
       }

       \vspace*{1.5cm}
       \noindent {\large Pr\'eparation} \\
       S\'eparer les jaunes des blancs.
       M\'elanger les jaunes avec le sucre jusqu'\`a ce que le m\'elange mousse et blanchisse. \\
       Faire fondre le chocolat doucement avec le beurre. L'incorporer au m\'elange
       jaunes/sucre. \\
       Monter les blancs en neige tr\`es ferme. Les incorporer au m\'elange jaunes/chocolat. \\

       \noindent Un truc : incorporer le premier tiers des blancs en remuant \'ener\-giquement pour que le m\'elange soit moins ferme, puis incorporer d\'elicatement le reste des blancs.

       \end{document}
       */

} // namespace Givaro

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
