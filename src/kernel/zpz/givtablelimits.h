#ifndef __GIVARO_TABLESIZE_MAX__
#define __GIVARO_TABLESIZE_MAX__

// 2^23 ---> 2^23*4*3 = 100K
// #define FF_TABLE_MAX 8388608UL
// 2^20 ---> 2s on 735MHz
//#define FF_TABLE_MAX 1048576UL
// Now 2^21+1 seems OK
#define FF_TABLE_MAX 2097153UL

#endif
