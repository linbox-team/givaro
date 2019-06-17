// ===============================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Author: J-G. Dumas
// Time-stamp: <17 Jun 19 14:48:01 Jean-Guillaume.Dumas@imag.fr>
// Description: generic rational fraction reconstruction
// ===============================================================
#ifndef __GIVARO_poly1_ratrecon_INL
#define __GIVARO_poly1_ratrecon_INL


namespace Givaro {

    template <class Domain>
    bool Poly1Dom<Domain,Dense>::ratrecon(typename Poly1Dom<Domain,Dense>::Rep& N, typename Poly1Dom<Domain,Dense>::Rep& D, const typename Poly1Dom<Domain,Dense>::Rep& P, const typename Poly1Dom<Domain,Dense>::Rep& M, const Degree& dk) const {

        Degree degU, degV;
        this->degree(degU,P); this->degree(degV,M);
        if ((degU < dk) || (degV == 0)) { this->assign(N,P); this->assign(D,one); return true; }
        if ((degV < 0) || (degU == 0)) { this->assign(N,one); this->assign(D,one); return false; }

        typename Poly1Dom<Domain,Dense>::Rep U;
        this->assign(N, M);
        this->assign(U, P);

        Degree degN;
        typename Poly1Dom<Domain,Dense>::Rep Q, D0;
        this->assign(D0,this->zero);
        this->assign(D,this->one);
        //   do {
        //       this->divmod(Q,N,V,U);

        //       this->assign(V,U);
        //       this->assign(U,N);
        //       this->maxpy(TMP2,Q,D,D0);
        //       this->assign(D0,D);
        //       this->assign(D,TMP2);
        //       this->degree(degN, N);
        //       if (degN <= dk) break;
        //   } while (degN>=0);

        do {
            this->divmodin(Q,N,U);
            this->maxpyin(D0,Q,D);
            this->degree(degN, N);
            if ((degN <= dk) || (degN <0)) {
                this->assign(D,D0);
                break;
            }

            this->divmodin(Q,U,N);
            this->maxpyin(D,Q,D0);

            this->degree(degN, U);
            if (degN <= dk) {
                this->assign(N,U);
                break;
            }

        } while (degN>=0);


        //   do {
        //       this->divmod(Q,N,V,U);
        //       this->maxpyin(D0,Q,D);
        //       this->degree(degN, N);
        //       if ((degN <= dk) || (degN <0)) {
        //           this->assign(D,D0);
        //           break;
        //       }
        //       this->assign(V,N);

        //       this->divmod(Q,N,U,V);
        //       this->maxpyin(D,Q,D0);

        //       this->degree(degN, N);
        //       if (degN <= dk) break;

        //       this->assign(U,N);
        //   } while (degN>=0);

        return (degN <= dk);

    }


    template <class Domain>
    bool Poly1Dom<Domain,Dense>::ratreconcheck(typename Poly1Dom<Domain,Dense>::Rep& N, typename Poly1Dom<Domain,Dense>::Rep& D, const typename Poly1Dom<Domain,Dense>::Rep& P, const typename Poly1Dom<Domain,Dense>::Rep& M, const Degree& dk) const {
        bool pass=ratrecon(N,D,P,M,dk);


        typename Poly1Dom<Domain,Dense>::Rep G;
        Degree degG;
        if ( degree(degG, gcd(G, N, D)) > 0)
            return false;

        Type_t r; leadcoef(r, D);
        if (! _domain.isOne(r)) {
            this->divin(D, r);
            this->divin(N, r);
        }

        return pass;
    }

    template <class Domain>
    bool Poly1Dom<Domain,Dense>::ratrecon(typename Poly1Dom<Domain,Dense>::Rep& N, typename Poly1Dom<Domain,Dense>::Rep& D, const typename Poly1Dom<Domain,Dense>::Rep& P, const typename Poly1Dom<Domain,Dense>::Rep& M, const Degree& dk, const bool forcereduce) const {
        if (forcereduce) 
            return ratreconcheck(N,D,P,M,dk);
        else
            return ratrecon(N,D,P,M,dk);
    }
    
} // Givaro

#endif // __GIVARO_poly1_ratrecon_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
