#include <iostream>
using namespace std;
#define GIVARO_LENSTRA
#include "givaro/givintrsa.h"
#include "givaro/givtimer.h"



int main(int argc, char** argv)
{
    Timer tim;
    tim.clear();
    
    IntRSADom<>::Element m,k,u;
    if (argc > 1)
        m = IntRSADom<>::Element( argv[1] );
    else 
        cin >> m;
    if (argc > 2)
        k = IntRSADom<>::Element( argv[2] );
    else 
        cin >> k;
    
    IntRSADom<> IR(m,k);
    tim.start();
    IR.point_break(u);
    tim.stop();
    
    /* For a factored output :

    IR.write( cerr << "m=pq: ", IR.getm()) ;
    IR.write( cerr << ", cipher k: " , IR.getk()) << "   ----> decipering key: " ;
    IR.write( cout, u ) << endl;

    */

    // Unfactored output
    cerr << "n=pq: " << IR.getn() << ", cipher key: " << IR.gete() << "   ----> decipering key: ";
    cout << u << endl;



    cerr << tim << endl;

    return 0;
}

