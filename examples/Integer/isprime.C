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
  IntPrimeDom::element m, ff(5);
  if (argc > 1) m = Integer(argv[1]);
  if (argc > 2) ff = Integer(argv[2]);
  
        Timer tim; tim.clear(); tim.start();
        bool a = IP.isprime(m,ff);
        tim.stop();
        cout << (a?"true":"false") << endl;
        cerr << tim << endl;

//  Givaro::End();

  return 0;
}

