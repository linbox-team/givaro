/* Copyright Université Grenoble Alpes
Contributors :
Jean-Guillaume DUMAS (Jean-Guillaume.Dumas@imag.fr)

Time-stamp: <27 Jul 16 15:23:32 Jean-Guillaume.Dumas@imag.fr>

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
#include <recint/ruint.h>
#include <recint/rint.h>

int main(void)
{
    RecInt::ruint64  a(1234567890);
    RecInt::ruint128 b("12345678901234567890123456789012345678901234567890123456789012345678901234567890");
    RecInt::ruint256 c("12345678901234567890123456789012345678901234567890123456789012345678901234567890");
    RecInt::ruint512 d("12345678901234567890123456789012345678901234567890123456789012345678901234567890");
    RecInt::rint64  r(1234567890);
    RecInt::rint128 s("12345678901234567890123456789012345678901234567890123456789012345678901234567890");
    RecInt::rint256 t("12345678901234567890123456789012345678901234567890123456789012345678901234567890");
    RecInt::rint512 u("12345678901234567890123456789012345678901234567890123456789012345678901234567890");

    std::cerr << " 64 bits; " << sizeof(a) << " octets: " <<  a << std::endl;
    std::cerr << "128 bits;" << sizeof(b) << " octets: " <<  b << std::endl;
    std::cerr << "256 bits;" << sizeof(c) << " octets: " <<  c << std::endl;
    std::cerr << "512 bits;" << sizeof(d) << " octets: " <<  d << std::endl;
    std::cerr << " 64 bits; " << sizeof(r) << " octets: " <<  r << std::endl;
    std::cerr << "128 bits;" << sizeof(s) << " octets: " <<  s << std::endl;
    std::cerr << "256 bits;" << sizeof(t) << " octets: " <<  t << std::endl;
    std::cerr << "512 bits;" << sizeof(u) << " octets: " <<  u << std::endl;

    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
