#include <iostream>
using namespace std;
#include <givaro/givinit.h>
#include <givaro/givintfactor.h>
#include <givaro/givtimer.h>


int main(int argc, char** argv)
{
    IntFactorDom<> IP;
    Integer m;
    if (argc > 1)
       m = Integer(argv[1]);
    else
        cin >> m;
    if (IP.islt(m,0) ) {
        cerr << "-";
        IP.negin(m);
   } 
    if (IP.islt(m,4))
        IP.write(cerr,m) << endl;
    else {
        Timer tim; tim.clear(); tim.start();
        IP.write(cerr,m) << endl;
        tim.stop();
        cerr << tim << endl;
    }
    return 0;
}

