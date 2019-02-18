// ==========================================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <22 Mar 00 20:06:18 Jean-Guillaume.Dumas@imag.fr>
// ==========================================================================
// Tabulation of factors of cyclotomic polynomials
// of degree expo modulo mod
// By vectors. P = v[0] + v[1] X + ... + v[n] X^n
#ifndef __GIVARO_poly1_cyclo_table_INL
#define __GIVARO_poly1_cyclo_table_INL

#include "givaro/givpoly1tabcycl.h"


namespace Givaro {

    template<class Domain, class Tag> inline void CyclotomicTable<Domain,Tag>::table_0 (const typename Domain::Residu_t mod, const long expo) {
        switch (mod) {
        case 2:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],1); _domain.read(v[1],1); _domain.read(v[2],1); _Irreductible =  v; }; break; // 2^2 : 1 + X + X^2
            case 3: { Element v(4); _domain.read(v[0],1); _domain.read(v[1],0); _domain.read(v[2],1); _domain.read(v[3],1); _Irreductible =  v; }; break; // 2^3 : 1 + X^2 + X^3
            case 4: { Element v(5); _domain.read(v[0],1); _domain.read(v[1],1); _domain.read(v[2],0); _domain.read(v[3],0); _domain.read(v[4],1); _Irreductible =  v; }; break; // 2^4 : 1 + X + X^4
            case 5: { Element v(6); _domain.read(v[0],1); _domain.read(v[1],0); _domain.read(v[2],1); _domain.read(v[3],1); _domain.read(v[4],1); _domain.read(v[5],1); _Irreductible =  v; }; break; // 2^5 : 1 + X^2 + X^3 + X^4 + X^5
            case 6: { Element v(7); _domain.read(v[0],1); _domain.read(v[1],1); _domain.read(v[2],1); _domain.read(v[3],0); _domain.read(v[4],0); _domain.read(v[5],1); _domain.read(v[6],1); _Irreductible =  v; }; break; // 2^6 : 1 + X + X^2 + X^5 + X^6
            case 7: { Element v(8); _domain.read(v[0],1); _domain.read(v[1],1); _domain.read(v[2],1); _domain.read(v[3],1); _domain.read(v[4],1); _domain.read(v[5],1); _domain.read(v[6],0); _domain.read(v[7],1); _Irreductible =  v; }; break; // 2^7 : 1 + X + X^2 + X^3 + X^4 + X^5 + X^7
            case 8: { Element v(9); _domain.read(v[0],1); _domain.read(v[1],0); _domain.read(v[2],0); _domain.read(v[3],0); _domain.read(v[4],1); _domain.read(v[5],1); _domain.read(v[6],1); _domain.read(v[7],0); _domain.read(v[8],1); _Irreductible =  v; }; break; // 2^8 : 1 + X^4 + X^5 + X^6 + X^8
            case 9: { Element v(10); _domain.read(v[0],1); _domain.read(v[1],0); _domain.read(v[2],1); _domain.read(v[3],0); _domain.read(v[4],1); _domain.read(v[5],1); _domain.read(v[6],1); _domain.read(v[7],1); _domain.read(v[8],0); _domain.read(v[9],1); _Irreductible =  v; }; break; // 2^9 : 1 + X^2 + X^4 + X^5 + X^6 + X^7 + X^9
            case 10: { Element v(11); _domain.read(v[0],1); _domain.read(v[1],1); _domain.read(v[2],0); _domain.read(v[3],1); _domain.read(v[4],0); _domain.read(v[5],1); _domain.read(v[6],1); _domain.read(v[7],0); _domain.read(v[8],1); _domain.read(v[9],0); _domain.read(v[10],1); _Irreductible =  v; }; break; // 2^10 : 1 + X + X^3 + X^5 + X^6 + X^8 + X^10
            case 11: { Element v(12); _domain.read(v[0],1); _domain.read(v[1],0); _domain.read(v[2],0); _domain.read(v[3],1); _domain.read(v[4],1); _domain.read(v[5],1); _domain.read(v[6],0); _domain.read(v[7],0); _domain.read(v[8],0); _domain.read(v[9],1); _domain.read(v[10],1); _domain.read(v[11],1); _Irreductible =  v; }; break; // 2^11 : 1 + X^3 + X^4 + X^5 + X^9 + X^10 + X^11
            case 12: { Element v(13); _domain.read(v[0],1); _domain.read(v[1],1); _domain.read(v[2],1); _domain.read(v[3],1); _domain.read(v[4],0); _domain.read(v[5],0); _domain.read(v[6],0); _domain.read(v[7],0); _domain.read(v[8],1); _domain.read(v[9],0); _domain.read(v[10],1); _domain.read(v[11],0); _domain.read(v[12],1); _Irreductible =  v; }; break; // 2^12 : 1 + X + X^2 + X^3 + X^8 + X^10 + X^12
            case 13: { Element v(14); _domain.read(v[0],1); _domain.read(v[1],1); _domain.read(v[2],1); _domain.read(v[3],0); _domain.read(v[4],0); _domain.read(v[5],0); _domain.read(v[6],0); _domain.read(v[7],0); _domain.read(v[8],0); _domain.read(v[9],1); _domain.read(v[10],0); _domain.read(v[11],1); _domain.read(v[12],1); _domain.read(v[13],1); _Irreductible =  v; }; break; // 2^13 : 1 + X + X^2 + X^9 + X^11 + X^12 + X^13
            case 14: { Element v(15); _domain.read(v[0],1); _domain.read(v[1],0); _domain.read(v[2],1); _domain.read(v[3],1); _domain.read(v[4],0); _domain.read(v[5],1); _domain.read(v[6],1); _domain.read(v[7],1); _domain.read(v[8],1); _domain.read(v[9],1); _domain.read(v[10],0); _domain.read(v[11],0); _domain.read(v[12],0); _domain.read(v[13],0); _domain.read(v[14],1); _Irreductible =  v; }; break; // 2^14 : 1 + X^2 + X^3 + X^5 + X^6 + X^7 + X^8 + X^9 + X^14
            case 15: { Element v(16); _domain.read(v[0],1); _domain.read(v[1],1); _domain.read(v[2],0); _domain.read(v[3],0); _domain.read(v[4],0); _domain.read(v[5],1); _domain.read(v[6],1); _domain.read(v[7],1); _domain.read(v[8],1); _domain.read(v[9],0); _domain.read(v[10],0); _domain.read(v[11],1); _domain.read(v[12],1); _domain.read(v[13],1); _domain.read(v[14],1); _domain.read(v[15],1); _Irreductible =  v; }; break; // 2^15 : 1 + X + X^5 + X^6 + X^7 + X^8 + X^11 + X^12 + X^13 + X^14 + X^15
            case 16: { Element v(17); _domain.read(v[0],1); _domain.read(v[1],0); _domain.read(v[2],1); _domain.read(v[3],1); _domain.read(v[4],0); _domain.read(v[5],1); _domain.read(v[6],0); _domain.read(v[7],0); _domain.read(v[8],1); _domain.read(v[9],0); _domain.read(v[10],1); _domain.read(v[11],0); _domain.read(v[12],1); _domain.read(v[13],0); _domain.read(v[14],1); _domain.read(v[15],0); _domain.read(v[16],1); _Irreductible =  v; }; break; // 2^16 : 1 + X^2 + X^3 + X^5 + X^8 + X^10 + X^12 + X^14 + X^16
            }; break;
        case 3:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],1); _domain.read(v[1],1); _domain.read(v[2],2); _Irreductible =  v; }; break; // 3^2 : 1 + X + (2)*X^2
            case 3: { Element v(4); _domain.read(v[0],1); _domain.read(v[1],2); _domain.read(v[2],1); _domain.read(v[3],1); _Irreductible =  v; }; break; // 3^3 : 1 + (2)*X + X^2 + X^3
            case 4: { Element v(5); _domain.read(v[0],1); _domain.read(v[1],2); _domain.read(v[2],0); _domain.read(v[3],0); _domain.read(v[4],2); _Irreductible =  v; }; break; // 3^4 : 1 + (2)*X + (2)*X^4
            case 5: { Element v(6); _domain.read(v[0],1); _domain.read(v[1],2); _domain.read(v[2],0); _domain.read(v[3],0); _domain.read(v[4],1); _domain.read(v[5],1); _Irreductible =  v; }; break; // 3^5 : 1 + (2)*X + X^4 + X^5
            case 6: { Element v(7); _domain.read(v[0],1); _domain.read(v[1],1); _domain.read(v[2],2); _domain.read(v[3],1); _domain.read(v[4],0); _domain.read(v[5],1); _domain.read(v[6],2); _Irreductible =  v; }; break; // 3^6 : 1 + X + (2)*X^2 + X^3 + X^5 + (2)*X^6
            case 7: { Element v(8); _domain.read(v[0],1); _domain.read(v[1],0); _domain.read(v[2],1); _domain.read(v[3],2); _domain.read(v[4],2); _domain.read(v[5],0); _domain.read(v[6],1); _domain.read(v[7],1); _Irreductible =  v; }; break; // 3^7 : 1 + X^2 + (2)*X^3 + (2)*X^4 + X^6 + X^7
            case 8: { Element v(9); _domain.read(v[0],2); _domain.read(v[1],0); _domain.read(v[2],2); _domain.read(v[3],1); _domain.read(v[4],2); _domain.read(v[5],2); _domain.read(v[6],2); _domain.read(v[7],2); _domain.read(v[8],1); _Irreductible =  v; }; break; // 3^8 : (2) + (2)*X^2 + X^3 + (2)*X^4 + (2)*X^5 + (2)*X^6 + (2)*X^7 + X^8
            case 9: { Element v(10); _domain.read(v[0],1); _domain.read(v[1],2); _domain.read(v[2],2); _domain.read(v[3],1); _domain.read(v[4],2); _domain.read(v[5],0); _domain.read(v[6],1); _domain.read(v[7],0); _domain.read(v[8],0); _domain.read(v[9],1); _Irreductible =  v; }; break; // 3^9 : 1 + (2)*X + (2)*X^2 + X^3 + (2)*X^4 + X^6 + X^9
            case 10: { Element v(11); _domain.read(v[0],2); _domain.read(v[1],2); _domain.read(v[2],0); _domain.read(v[3],1); _domain.read(v[4],1); _domain.read(v[5],1); _domain.read(v[6],0); _domain.read(v[7],1); _domain.read(v[8],0); _domain.read(v[9],1); _domain.read(v[10],1); _Irreductible =  v; }; break; // 3^10 : (2) + (2)*X + X^3 + X^4 + X^5 + X^7 + X^9 + X^10
            }; break;
        case 5:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],1); _domain.read(v[1],1); _domain.read(v[2],2); _Irreductible =  v; }; break; // 5^2 : 1 + X + (2)*X^2
            case 3: { Element v(4); _domain.read(v[0],4); _domain.read(v[1],4); _domain.read(v[2],0); _domain.read(v[3],3); _Irreductible =  v; }; break; // 5^3 : (4) + (4)*X + (3)*X^3
            case 4: { Element v(5); _domain.read(v[0],3); _domain.read(v[1],1); _domain.read(v[2],0); _domain.read(v[3],3); _domain.read(v[4],4); _Irreductible =  v; }; break; // 5^4 : (3) + X + (3)*X^3 + (4)*X^4
            case 5: { Element v(6); _domain.read(v[0],4); _domain.read(v[1],1); _domain.read(v[2],3); _domain.read(v[3],0); _domain.read(v[4],4); _domain.read(v[5],2); _Irreductible =  v; }; break; // 5^5 : (4) + X + (3)*X^2 + (4)*X^4 + (2)*X^5
            case 6: { Element v(7); _domain.read(v[0],2); _domain.read(v[1],3); _domain.read(v[2],2); _domain.read(v[3],4); _domain.read(v[4],0); _domain.read(v[5],1); _domain.read(v[6],1); _Irreductible =  v; }; break; // 5^6 : (2) + (3)*X + (2)*X^2 + (4)*X^3 + X^5 + X^6
            }; break;
        case 7:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],2); _domain.read(v[1],4); _domain.read(v[2],6); _Irreductible =  v; }; break; // 7^2 : (2) + (4)*X + (6)*X^2
            case 3: { Element v(4); _domain.read(v[0],1); _domain.read(v[1],6); _domain.read(v[2],5); _domain.read(v[3],4); _Irreductible =  v; }; break; // 7^3 : 1 + (6)*X + (5)*X^2 + (4)*X^3
            case 4: { Element v(5); _domain.read(v[0],4); _domain.read(v[1],1); _domain.read(v[2],1); _domain.read(v[3],6); _domain.read(v[4],5); _Irreductible =  v; }; break; // 7^4 : (4) + X + X^2 + (6)*X^3 + (5)*X^4
            case 5: { Element v(6); _domain.read(v[0],2); _domain.read(v[1],6); _domain.read(v[2],2); _domain.read(v[3],6); _domain.read(v[4],2); _domain.read(v[5],4); _Irreductible =  v; }; break; // 7^5 : (2) + (6)*X + (2)*X^2 + (6)*X^3 + (2)*X^4 + (4)*X^5
            }; break;
        case 11:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],6); _domain.read(v[1],2); _domain.read(v[2],9); _Irreductible =  v; }; break; // 11^2 : (6) + (2)*X + (9)*X^2
            case 3: { Element v(4); _domain.read(v[0],3); _domain.read(v[1],6); _domain.read(v[2],3); _domain.read(v[3],5); _Irreductible =  v; }; break; // 11^3 : (3) + (6)*X + (3)*X^2 + (5)*X^3
            case 4: { Element v(5); _domain.read(v[0],7); _domain.read(v[1],10); _domain.read(v[2],4); _domain.read(v[3],8); _domain.read(v[4],3); _Irreductible =  v; }; break; // 11^4 : (7) + (10)*X + (4)*X^2 + (8)*X^3 + (3)*X^4
            }; break;
        case 13:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],8); _domain.read(v[1],3); _domain.read(v[2],9); _Irreductible =  v; }; break; // 13^2 : (8) + (3)*X + (9)*X^2
            case 3: { Element v(4); _domain.read(v[0],10); _domain.read(v[1],6); _domain.read(v[2],7); _domain.read(v[3],7); _Irreductible =  v; }; break; // 13^3 : (10) + (6)*X + (7)*X^2 + (7)*X^3
            case 4: { Element v(5); _domain.read(v[0],1); _domain.read(v[1],12); _domain.read(v[2],10); _domain.read(v[3],4); _domain.read(v[4],6); _Irreductible =  v; }; break; // 13^4 : 1 + (12)*X + (10)*X^2 + (4)*X^3 + (6)*X^4
            }; break;
        case 17:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],9); _domain.read(v[1],2); _domain.read(v[2],12); _Irreductible =  v; }; break; // 17^2 : (9) + (2)*X + (12)*X^2
            case 3: { Element v(4); _domain.read(v[0],15); _domain.read(v[1],13); _domain.read(v[2],7); _domain.read(v[3],5); _Irreductible =  v; }; break; // 17^3 : (15) + (13)*X + (7)*X^2 + (5)*X^3
            }; break;
        case 19:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],4); _domain.read(v[1],18); _domain.read(v[2],3); _Irreductible =  v; }; break; // 19^2 : (4) + (18)*X + (3)*X^2
            case 3: { Element v(4); _domain.read(v[0],2); _domain.read(v[1],11); _domain.read(v[2],17); _domain.read(v[3],15); _Irreductible =  v; }; break; // 19^3 : (2) + (11)*X + (17)*X^2 + (15)*X^3
            }; break;
        case 23:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],15); _domain.read(v[1],1); _domain.read(v[2],3); _Irreductible =  v; }; break; // 23^2 : (15) + X + (3)*X^2
            case 3: { Element v(4); _domain.read(v[0],12); _domain.read(v[1],6); _domain.read(v[2],11); _domain.read(v[3],18); _Irreductible =  v; }; break; // 23^3 : (12) + (6)*X + (11)*X^2 + (18)*X^3
            }; break;
        case 29:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],13); _domain.read(v[1],17); _domain.read(v[2],19); _Irreductible =  v; }; break; // 29^2 : (13) + (17)*X + (19)*X^2
            case 3: { Element v(4); _domain.read(v[0],5); _domain.read(v[1],10); _domain.read(v[2],10); _domain.read(v[3],21); _Irreductible =  v; }; break; // 29^3 : (5) + (10)*X + (10)*X^2 + (21)*X^3
            }; break;
        case 31:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],24); _domain.read(v[1],1); _domain.read(v[2],2); _Irreductible =  v; }; break; // 31^2 : (24) + X + (2)*X^2
            case 3: { Element v(4); _domain.read(v[0],20); _domain.read(v[1],29); _domain.read(v[2],5); _domain.read(v[3],1); _Irreductible =  v; }; break; // 31^3 : (20) + (29)*X + (5)*X^2 + X^3
            }; break;
        case 37:
            switch (expo) {
            }; break;
        };
        if (isZero(_Irreductible)) table_50(mod,expo);
    }
    template<class Domain, class Tag> inline void CyclotomicTable<Domain,Tag>::table_50 (const typename Domain::Residu_t mod, const long expo) {
        switch (mod) {
        case 37:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],10); _domain.read(v[1],32); _domain.read(v[2],20); _Irreductible =  v; }; break; // 37^2 : (10) + (32)*X + (20)*X^2
            case 3: { Element v(4); _domain.read(v[0],34); _domain.read(v[1],16); _domain.read(v[2],35); _domain.read(v[3],35); _Irreductible =  v; }; break; // 37^3 : (34) + (16)*X + (35)*X^2 + (35)*X^3
            }; break;
        case 41:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],27); _domain.read(v[1],37); _domain.read(v[2],21); _Irreductible =  v; }; break; // 41^2 : (27) + (37)*X + (21)*X^2
            }; break;
        case 43:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],24); _domain.read(v[1],18); _domain.read(v[2],22); _Irreductible =  v; }; break; // 43^2 : (24) + (18)*X + (22)*X^2
            }; break;
        case 47:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],43); _domain.read(v[1],44); _domain.read(v[2],3); _Irreductible =  v; }; break; // 47^2 : (43) + (44)*X + (3)*X^2
            }; break;
        case 53:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],43); _domain.read(v[1],52); _domain.read(v[2],22); _Irreductible =  v; }; break; // 53^2 : (43) + (52)*X + (22)*X^2
            }; break;
        case 59:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],11); _domain.read(v[1],31); _domain.read(v[2],15); _Irreductible =  v; }; break; // 59^2 : (11) + (31)*X + (15)*X^2
            }; break;
        case 61:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],15); _domain.read(v[1],60); _domain.read(v[2],17); _Irreductible =  v; }; break; // 61^2 : (15) + (60)*X + (17)*X^2
            }; break;
        case 67:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],9); _domain.read(v[1],54); _domain.read(v[2],63); _Irreductible =  v; }; break; // 67^2 : (9) + (54)*X + (63)*X^2
            }; break;
        case 71:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],39); _domain.read(v[1],63); _domain.read(v[2],15); _Irreductible =  v; }; break; // 71^2 : (39) + (63)*X + (15)*X^2
            }; break;
        case 73:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],23); _domain.read(v[1],1); _domain.read(v[2],29); _Irreductible =  v; }; break; // 73^2 : (23) + X + (29)*X^2
            }; break;
        case 79:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],39); _domain.read(v[1],32); _domain.read(v[2],36); _Irreductible =  v; }; break; // 79^2 : (39) + (32)*X + (36)*X^2
            }; break;
        case 83:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],79); _domain.read(v[1],8); _domain.read(v[2],30); _Irreductible =  v; }; break; // 83^2 : (79) + (8)*X + (30)*X^2
            }; break;
        case 89:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],88); _domain.read(v[1],13); _domain.read(v[2],76); _Irreductible =  v; }; break; // 89^2 : (88) + (13)*X + (76)*X^2
            }; break;
        case 97:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],16); _domain.read(v[1],39); _domain.read(v[2],67); _Irreductible =  v; }; break; // 97^2 : (16) + (39)*X + (67)*X^2
            }; break;
        case 101:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],93); _domain.read(v[1],61); _domain.read(v[2],31); _Irreductible =  v; }; break; // 101^2 : (93) + (61)*X + (31)*X^2
            }; break;
        case 103:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],16); _domain.read(v[1],34); _domain.read(v[2],22); _Irreductible =  v; }; break; // 103^2 : (16) + (34)*X + (22)*X^2
            }; break;
        case 107:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],90); _domain.read(v[1],80); _domain.read(v[2],6); _Irreductible =  v; }; break; // 107^2 : (90) + (80)*X + (6)*X^2
            }; break;
        case 109:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],41); _domain.read(v[1],101); _domain.read(v[2],20); _Irreductible =  v; }; break; // 109^2 : (41) + (101)*X + (20)*X^2
            }; break;
        case 113:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],99); _domain.read(v[1],111); _domain.read(v[2],76); _Irreductible =  v; }; break; // 113^2 : (99) + (111)*X + (76)*X^2
            }; break;
        case 127:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],60); _domain.read(v[1],113); _domain.read(v[2],54); _Irreductible =  v; }; break; // 127^2 : (60) + (113)*X + (54)*X^2
            }; break;
        case 131:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],128); _domain.read(v[1],90); _domain.read(v[2],84); _Irreductible =  v; }; break; // 131^2 : (128) + (90)*X + (84)*X^2
            }; break;
        case 137:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],120); _domain.read(v[1],85); _domain.read(v[2],127); _Irreductible =  v; }; break; // 137^2 : (120) + (85)*X + (127)*X^2
            }; break;
        case 139:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],131); _domain.read(v[1],16); _domain.read(v[2],88); _Irreductible =  v; }; break; // 139^2 : (131) + (16)*X + (88)*X^2
            }; break;
        case 149:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],36); _domain.read(v[1],84); _domain.read(v[2],97); _Irreductible =  v; }; break; // 149^2 : (36) + (84)*X + (97)*X^2
            }; break;
        case 151:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],54); _domain.read(v[1],149); _domain.read(v[2],29); _Irreductible =  v; }; break; // 151^2 : (54) + (149)*X + (29)*X^2
            }; break;
        case 157:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],93); _domain.read(v[1],8); _domain.read(v[2],62); _Irreductible =  v; }; break; // 157^2 : (93) + (8)*X + (62)*X^2
            }; break;
        case 163:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],38); _domain.read(v[1],112); _domain.read(v[2],122); _Irreductible =  v; }; break; // 163^2 : (38) + (112)*X + (122)*X^2
            }; break;
        case 167:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],35); _domain.read(v[1],70); _domain.read(v[2],108); _Irreductible =  v; }; break; // 167^2 : (35) + (70)*X + (108)*X^2
            }; break;
        case 173:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],154); _domain.read(v[1],119); _domain.read(v[2],81); _Irreductible =  v; }; break; // 173^2 : (154) + (119)*X + (81)*X^2
            }; break;
        case 179:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],131); _domain.read(v[1],103); _domain.read(v[2],74); _Irreductible =  v; }; break; // 179^2 : (131) + (103)*X + (74)*X^2
            }; break;
        case 181:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],156); _domain.read(v[1],135); _domain.read(v[2],164); _Irreductible =  v; }; break; // 181^2 : (156) + (135)*X + (164)*X^2
            }; break;
        case 191:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],91); _domain.read(v[1],15); _domain.read(v[2],49); _Irreductible =  v; }; break; // 191^2 : (91) + (15)*X + (49)*X^2
            }; break;
        case 193:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],65); _domain.read(v[1],7); _domain.read(v[2],82); _Irreductible =  v; }; break; // 193^2 : (65) + (7)*X + (82)*X^2
            }; break;
        case 197:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],179); _domain.read(v[1],72); _domain.read(v[2],132); _Irreductible =  v; }; break; // 197^2 : (179) + (72)*X + (132)*X^2
            }; break;
        case 199:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],35); _domain.read(v[1],157); _domain.read(v[2],179); _Irreductible =  v; }; break; // 199^2 : (35) + (157)*X + (179)*X^2
            }; break;
        case 211:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],68); _domain.read(v[1],122); _domain.read(v[2],43); _Irreductible =  v; }; break; // 211^2 : (68) + (122)*X + (43)*X^2
            }; break;
        case 223:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],63); _domain.read(v[1],61); _domain.read(v[2],123); _Irreductible =  v; }; break; // 223^2 : (63) + (61)*X + (123)*X^2
            }; break;
        case 227:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],73); _domain.read(v[1],60); _domain.read(v[2],15); _Irreductible =  v; }; break; // 227^2 : (73) + (60)*X + (15)*X^2
            }; break;
        case 229:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],171); _domain.read(v[1],75); _domain.read(v[2],79); _Irreductible =  v; }; break; // 229^2 : (171) + (75)*X + (79)*X^2
            }; break;
        case 233:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],224); _domain.read(v[1],55); _domain.read(v[2],176); _Irreductible =  v; }; break; // 233^2 : (224) + (55)*X + (176)*X^2
            }; break;
        case 239:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],58); _domain.read(v[1],70); _domain.read(v[2],221); _Irreductible =  v; }; break; // 239^2 : (58) + (70)*X + (221)*X^2
            }; break;
        case 241:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],48); _domain.read(v[1],194); _domain.read(v[2],219); _Irreductible =  v; }; break; // 241^2 : (48) + (194)*X + (219)*X^2
            }; break;
        case 251:
            switch (expo) {
            case 2: { Element v(3); _domain.read(v[0],108); _domain.read(v[1],61); _domain.read(v[2],163); _Irreductible =  v; }; break; // 251^2 : (108) + (61)*X + (163)*X^2
            }; break;
        };
        if (isZero(_Irreductible)) set_random_irreducible(_domain,expo);
    }

} // Givaro

#endif // __GIVARO_poly1_cyclo_table_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
