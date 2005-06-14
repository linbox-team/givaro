// ========================================================== //
// Time-stamp: <04 Sep 02 18:10:22 Jean-Guillaume.Dumas@imag.fr>
// ========================================================== //
#include <iostream>
using namespace std;
#include <stdlib.h>
#include <givaro/givintnumtheo.h>
#include <givaro/givtimer.h>


int main(int argc, char ** argv)
{
    IntNumTheoDom<> IP;
    IntNumTheoDom<>::Element a,q,o;
    if (argc > 1) a = Integer(argv[1]); else cin >> a;
    if (argc > 2) q = Integer(argv[2]); else cin >> q;

    Timer tim; tim.clear(); tim.start();
	// Ordre de a dans GF(q)
    IP.order(o, a, q);
    tim.stop();

    cout << o << endl;
    cerr << tim << endl;

    return 0;
}


