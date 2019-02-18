/* Copyright Université Grenoble Alpes
Contributors :
Jean-Guillaume DUMAS (Jean-Guillaume.Dumas@imag.fr)

Time-stamp: <27 Jul 16 12:06:56 Jean-Guillaume.Dumas@imag.fr>

This software is a computer program whose purpose is to provide a
fixed precision arithmetic library.

This software is governed by the CeCILL-B license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL-B
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL-B license and that you accept its terms.
*/
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

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
