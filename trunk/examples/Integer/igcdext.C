#include <iostream>
using namespace std;
#include <stdlib.h>
#include <givaro/givintprime.h>
#include <givaro/givtimer.h>
#include <givaro/givinit.h>         // Givaro initialization



int main(int argc, char** argv)
{
//  Givaro::Init(&argc, &argv);


  IntegerDom IP;
  IntegerDom::Element g,a,b,u,v;
  if (argc > 1) a = Integer(argv[1]); else cin >> a;
  if (argc > 2) b = Integer(argv[2]); else cin >> b;
  
        Timer tim; tim.clear(); tim.start();
        IP.gcd(g,u,v,a,b);
        tim.stop();
        cout << g << " = " << u << " * " << a << " + " << v << " * " << b << endl;
        cerr << tim << endl;

//  Givaro::End();

  return 0;
}

