#include <iostream>
#include <cstdlib>

#include <recint/ruint.h>

#if not defined(STD_RECINT_SIZE)
#define STD_RECINT_SIZE 9
#endif

using namespace RecInt;

int main(void)
{
    ruint<STD_RECINT_SIZE> n;
    ruint<STD_RECINT_SIZE>::cr_iterator itr;
    
    RecInt::srand(time(NULL));
    rand(n);

    std::cout << "Our " << NBBITS<STD_RECINT_SIZE>::value << "-bit-wide number is: " << std::endl << std::hex << n << std::endl;

    const ruint<STD_RECINT_SIZE> n2(n);
    std::cout << std::endl << "Reverse:" << std::endl;
    for (itr = n2.rbegin(); itr != n2.rend(); itr++)
        std::cout << std::hex << *itr << std::endl;

    return 0;
}

