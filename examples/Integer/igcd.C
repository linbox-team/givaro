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
  IntegerDom::Element GG, g,a,b;
  int offset = 0;
  if (argc > ++offset) a = Integer(argv[offset]); else cin >> a;
  if (argc > ++offset) b = Integer(argv[offset]); else cin >> b;
  
        Timer tim; tim.clear(); tim.start();
        IP.gcd(GG,a,b);
	for ( ; argc > ++offset; ) {
		a = Integer(argv[offset]);
		IP.gcd(g, GG, a);
		GG = g;
	}	
        tim.stop();
        cout << GG << endl;
        cerr << tim << endl;

//  Givaro::End();

  return 0;
}

