#include <iostream>
using namespace std;
#include <stdlib.h>
#include <givaro/givintprime.h>
#include <givaro/givtimer.h>
#include <givaro/givinit.h>         // Givaro initialization



int main(int argc, char** argv)
{
//  Givaro::Init(&argc, &argv);


  IntPrimeDom IP;
  IntPrimeDom::element m, ff;
  if (argc > 1) m = Integer(argv[1]);
  
        Timer tim; tim.clear(); tim.start();
        bool a = IP.isprime(m);
        tim.stop();
        cout << (a?"true":"false") << endl;
        cerr << tim << endl;

//  Givaro::End();

  return 0;
}

