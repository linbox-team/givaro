// ==================================================================== //
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <05 Apr 11 17:28:53 Jean-Guillaume.Dumas@imag.fr>        //
// ==================================================================== //
// Givaro replacement for ssh-keygen: generated keys use strong primes  //
// Random generator is seeded by					 //
//	- true randomness in file, e.g. 				 //
//	  dd if=/dev/urandom of=.rand bs=1k count=16			 //
//	- otherwise gettimeofday is used, see givrandom.h		 //
// File formats are then managed by openssl and openssh, thus this file //
// requires to be compiled with -lssl -lssh -lopenbsd-compat            //
// Warning: Sometimes putting -lssh twice is required                   //
// You also need the openssl and openssh development headers below:     //
// openssh/buffer.h openssh/key.h openssh/uuencode.h openssh/xmalloc.h  //
// The latter libraries/heqders are avaible from openssl and openssh    //
// openssl > 0.9.8  and  openssh > 5.2p1  are expected                  //
// ==================================================================== //
/*! @file examples/Integer/givaro-ssh-keygen.C
 * @ingroup examples
 * @ingroup integers
 * @example examples/Integer/givaro-ssh-keygen.C
 * @brief NO DOC
 */
#include <iostream>
#include <fstream>
#include <sys/stat.h>

#include <openssl/rsa.h>
#include <openssl/pem.h>

extern "C" {
#include "openssh/key.h"
}

#include <givaro/givintrsa.h>
#include <givaro/givtimer.h>

using namespace Givaro;

/* Strong pseudo-prime generation using Gordon's algorithm */
template<class RandIter=GivRandom>
struct Givaro_keygen {
    void operator()(Integer& n, Integer& e, Integer& d,
                   Integer& p, Integer& q,
                   Integer& dmp1, Integer& dmq1,
                   Integer& iqmp, long size, unsigned long seed) {

        typename IntRSADom<RandIter>::random_generator gen(seed);
        Integer::seeding(gen.seed());
        IntRSADom<RandIter> IRD(false,gen);

        IRD.keys_gen(gen, (size>>1)+1, (size>>1)-1, n, e, d, p, q);

        Integer phim,p1,q1; IRD.mul(phim, IRD.sub(p1,p,IRD.one), IRD.sub(q1,q,IRD.one));

        Integer v, g;

        IRD.mod(e, 65537, phim);
        IRD.gcd(g,d,v,e,phim);

        IRD.modin(d,phim);
        if ( IRD.islt(d,IRD.zero) ) IRD.addin(d,phim);


        IRD.mod(dmp1,d,p1);
        IRD.mod(dmq1,d,q1);

        IRD.gcd(g,iqmp,v,q,p);
        if ( IRD.islt(iqmp,IRD.zero) ) IRD.addin(iqmp,p);
    }
};

/* Converts a givaro integer to an openssl BIGNUM */
BIGNUM* Integer2BN(BIGNUM * n, const Integer& a) {
    std::string str(a);
    BN_dec2bn(&n,str.c_str());
    return n;
}

int mymain(FILE* fileout, FILE* filepub, long s, unsigned long seed) {

    RSA *rsa= new RSA();

    Integer in, ie, id, ip, iq, idmp1, idmq1, iiqmp;
    Timer tim;tim.clear();tim.start();
    Givaro_keygen<>()(in,ie,id,ip,iq,idmp1,idmq1,iiqmp, s, seed);
    tim.stop();

/*
    std::cerr << "n: " << in << std::endl;
    std::cerr << "e: " << ie << std::endl;
    std::cerr << "d: " << id << std::endl;
    std::cerr << "p: " << ip << std::endl;
    std::cerr << "q: " << iq << std::endl;
    std::cerr << "dmp1: " << idmp1 << std::endl;
    std::cerr << "dmq1: " << idmq1 << std::endl;
    std::cerr << "iqmp: " << iiqmp << std::endl;
*/

    std::cerr << tim << std::endl;

    rsa->n = BN_new(); Integer2BN(rsa->n, in);
    rsa->e = BN_new(); Integer2BN(rsa->e, ie);
    rsa->d = BN_new(); Integer2BN(rsa->d, id);
    rsa->p = BN_new(); Integer2BN(rsa->p, ip);
    rsa->q = BN_new(); Integer2BN(rsa->q, iq);
    rsa->dmp1 = BN_new(); Integer2BN(rsa->dmp1, idmp1);
    rsa->dmq1 = BN_new(); Integer2BN(rsa->dmq1, idmq1);
    rsa->iqmp = BN_new(); Integer2BN(rsa->iqmp, iiqmp);


    Key rsakey;
    rsakey.type=KEY_RSA;
    rsakey.rsa = rsa;

/*
    std::cerr << "print: " << std::endl;
    RSA_print_fp(stdout, rsa, 0);

    std::cerr << "PEM Write: " << std::endl;
    PEM_write_RSAPublicKey(stdout,rsa);
*/

    std::cerr << "key's randomart: \n" << key_fingerprint(&rsakey, SSH_FP_MD5,SSH_FP_RANDOMART) << std::endl;

        // Write Private Key in ssl PEM format
    PEM_write_RSAPrivateKey(fileout,rsa,NULL,NULL,0,NULL,NULL);

        // Write Public key in ssh b64 format
    key_write(&rsakey, filepub);
    fprintf(filepub,"\n");

    return 0;

}

unsigned long seedfromfile(char * filename) {
     std::ifstream filrand(filename);
     unsigned long seed=0;
     for(int i=0; i<sizeof(unsigned long); ++i) {
         unsigned char t; filrand >> t;
         seed <<= 8;
         seed |= t;
     }
std::cerr << "Generated seed: " << seed << ", using " << filename << std::endl;
     return seed;
}




int main(int argc, char** argv)
{
    if (argc > 10) {
        std::cerr << "Usage: givaro-ssh-keygen [-b bits] [-f private-key-file] [-p public-key-file] [-r randomness-file]" << std::endl;
        return 0;
    }

    long s = 4096;
    unsigned long seed = 0;
    long files = 0;
    std::string filprivname, filpubname;


    for (long i = 1; i < argc; i++) {
        if (argv[i][0] == '-') {
            switch(argv[i][1]) {
                case 'b':; case 'B': {
                    s = atoi(argv[++i]);
                    break;
                }
                case 'f':; case 'F': {
                    filprivname = std::string(argv[++i]);
                    ++files;
                    break;
                }
                case 'p':; case 'P': {
                    filpubname = std::string(argv[++i]);
                    ++files;
                    break;
                }
                case 'r':; case 'R': {
                    seed = seedfromfile(argv[++i]);
                    break;
                }
            }
        }
    }


    FILE * filpriv, *filpub;
    if (files > 1) {
        filpriv = fopen(filprivname.c_str(),"w");
        if (argc>3) {
            filpub = fopen(filpubname.c_str(),"w");
            mymain(filpriv,filpub,s,seed);
            fclose(filpub);
        } else
            mymain(filpriv,stdout,s,seed);
        fclose(filpriv);
        chmod(filprivname.c_str(),(S_IRUSR|S_IWUSR));
    } else
        mymain(stdout,stdout,s,seed);
    return 0;
}
