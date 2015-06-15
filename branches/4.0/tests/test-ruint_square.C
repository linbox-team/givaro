/* ruint_square.cpp - Square functions of RecInt library test file

Return value.
    0    No error
    != 0 Bad result for an operation

The following constants have to be defined.
    LOOPS           number of loops of randomized tests
*/

#include <recint/ruint.h>

#if not defined(LOOPS)
#define LOOPS 10000
#endif

using namespace RecInt;

int main(void)
{
    ruint<STD_RECINT_SIZE> a, ma, sa;
    ruint<STD_RECINT_SIZE+1> m, s;
    
    RecInt::srand(limb(time(NULL)));
    
    // Loop
    for (UDItype i = 1; i < LOOPS; i++) {
        rand(a);
        
        //------ lsquare -------
        lmul(m, a, a);
        lsquare(s, a);
        if (s != m) return 1;
    }

    return 0; 
}

