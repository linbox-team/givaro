// ==================================================================== //
// Givaro replacement for ssh-keygen: generated keys use strong primes  //
// File formats are then managed by openssl and openssh, thus this file //
// requires to be compiled with -lssl -lssh -lopenbsd-compat            //
// The latter libraries are avaible from openssl and openssh            //
// Time-stamp: <12 Jun 09 14:45:41 Jean-Guillaume.Dumas@imag.fr>        //
// ==================================================================== //
#include <iostream>
#include <sys/stat.h>

#include <openssl/rsa.h>
#include <openssl/pem.h>

#include <givaro/givintrsa.h>
#include <givaro/givtimer.h>

extern "C" {
#include "openssh/key.h"
}

/* Strong pseudo-prime generation using Gordon's algorithm */
void Givaro_keygen(Integer& n, Integer& e, Integer& d, 
                   Integer& p, Integer& q,
                   Integer& dmp1, Integer& dmq1,
                   Integer& iqmp, long size) {

    IntRSADom<> IRD; IntRSADom<>::random_generator gen;
    
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

};

/* Converts a givaro integer to an openssl BIGNUM */
BIGNUM* Integer2BN(BIGNUM * n, const Integer& a) {
    std::string str(a);
    BN_dec2bn(&n,str.c_str());
    return n;
}

int mymain(FILE* fileout, FILE* filepub, long s) {

    RSA *rsa= new RSA();

    Integer in, ie, id, ip, iq, idmp1, idmq1, iiqmp;
    Timer tim;tim.clear();tim.start();
    Givaro_keygen(in,ie,id,ip,iq,idmp1,idmq1,iiqmp, s);
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

int main(int argc, char** argv) 
{
    long s =  argc>1? atoi(argv[1]) : 4096;


    FILE * filpriv, *filpub;
    if (argc>2) {
        filpriv = fopen(argv[2],"w");
        if (argc>3) {
            filpub = fopen(argv[3],"w");
            mymain(filpriv,filpub,s);
            fclose(filpub);
        } else
            mymain(filpriv,stdout,s);
        fclose(filpriv);
    } else
        mymain(stdout,stdout,s);
    return 0;
}
