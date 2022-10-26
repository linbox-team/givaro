 // ==========================================================================
// Copyright(c)'1994-2022 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: B. Grenet
//
// Generalized middle product:
// For polynomials P of size m+n-1 and Q of size n, define their
// middle product MP(P,Q) as the central m coefficients of P*Q
// that is MP(P,Q) = (PQ div X^(n-1)) mod X^m
//
// Balanced case: m = n
// ==========================================================================
#ifndef __GIVARO_poly1_midmul_INL
#define __GIVARO_poly1_midmul_INL

namespace Givaro {

#ifndef KARA_THRESHOLD
#define KARA_THRESHOLD 50
#endif

    // Generic middle product with dynamic recursive choices between stdmidmul and karamidmul
    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::midmul(
        Rep& R, const Rep& P, const Rep& Q ) const
    {
        size_t sR = R.size();
        size_t sP = P.size();
        size_t sQ = Q.size();
        if ((sQ ==0) || (sP ==0)) { R.resize(0); return R; }
        if (sR != sP-sQ+1) R.resize(sR = sP-sQ+1);

        // Generic middle product handler
        // Can use e.g. Karatsuba multiplication
        midmul(R, R.begin(), R.end(),
               P, P.begin(), P.end(),
               Q, Q.begin(), Q.end());

        return setdegree(R);

    }

    // forces standard middle product
    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::stdmidmul(
        Rep& R, const Rep& P, const Rep& Q ) const
    {
        const size_t sP = P.size();
        const size_t sQ = Q.size();
        const size_t sR = sP-sQ+1;
        if (sR != R.size()) R.resize(sR);

        stdmidmul(R, R.begin(), R.end(),
               P, P.begin(), P.end(),
               Q, Q.begin(), Q.end());

        return setdegree(R);
    }

    // forces FIRST recursive level with Karatsuba algorithm
    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::karamidmul(
        Rep& R, const Rep& P, const Rep& Q ) const
    {
        // Assumes sP = 2*sQ-1; otherwise undefined behavior
        const size_t sP = P.size();
        const size_t sQ = Q.size();
        const size_t sR = sP-sQ+1;
        if (sR != R.size()) R.resize(sR);

        karamidmul(R, R.begin(), R.end(),
               P, P.begin(), P.end(),
               Q, Q.begin(), Q.end());

        return setdegree(R);
    }

    // Generic midmul that dispatches to balanced midmul
    // Multiplies between the iterator bounds.
    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::midmul(
        Rep& R, const RepIterator Rbeg, const RepIterator Rend,
        const Rep& P, const RepConstIterator Pbeg, const RepConstIterator Pend,
        const Rep& Q, const RepConstIterator Qbeg, const RepConstIterator Qend)
        const {

            const ssize_t n = (ssize_t) (Qend-Qbeg);
            const ssize_t m = (ssize_t) (Pend-Pbeg) - n + 1;

        //std::cerr << "Entering GenMidMul with "
        //          << "P[" << (Pbeg-P.begin()) << ":" << (Pend-P.begin()) << "/" << (P.end()-P.begin()) << "] × "
        //          << "Q[" << (Qbeg-Q.begin()) << ":" << (Qend-Q.begin()) << "/" << (Q.end()-Q.begin()) << "] → "
        //          << "R[" << (Rbeg-R.begin()) << ":" << (Rend-R.begin()) << "/" << (R.end()-R.begin()) << "] : "
        //          << "m=" << m << ", n=" << n << std::endl;

            if ( std::min(m,n) <= KARA_THRESHOLD)
                return stdmidmul(R, Rbeg, Rend,
                                 P, Pbeg, Pend,
                                 Q, Qbeg, Qend);

            if (m == n)
                return karamidmul(R, Rbeg, Rend,
                                  P, Pbeg, Pend,
                                  Q, Qbeg, Qend);

            if (m > n) {
                ssize_t i = 0;
                for (; i <= m-n; i+=n)
                    karamidmul(R, Rbeg+i, Rbeg+i+n,
                               P, Pbeg+i, Pbeg+i+2*n-1,
                               Q, Qbeg, Qend);

                if (i < m)
                    midmul(R, Rbeg + i, Rend,
                           P, Pbeg + i, Pend,
                           Q, Qbeg, Qend);

                return R;
            }

            // m < n

            RepConstIterator Pe = Pend, Qb = Qbeg;

            karamidmul(R, Rbeg, Rend,
                       P, Pe-(2*m-1), Pe,
                       Q, Qb, Qb+m);

            Pe -= m; Qb += m;

            Rep Tmp; Tmp.resize(m);

            for(;Qb <= Qend-m; Pe-=m, Qb+=m) {
                karamidmul(Tmp, Tmp.begin(), Tmp.end(),
                           P, Pe-(2*m-1), Pe,
                           Q, Qb, Qb+m);
                RepIterator Ri = Rbeg;
                for( RepConstIterator Ti=Tmp.begin(); Ti != Tmp.end(); ++Ti,++Ri)
                    _domain.addin(*Ri, *Ti);
            }

            if (Qb < Qend) {
                midmul(Tmp, Tmp.begin(), Tmp.end(),
                       P, Pbeg, Pe,
                       Q, Qb, Qend);
                RepIterator Ri = Rbeg;
                for( RepConstIterator Ti=Tmp.begin(); Ti != Tmp.end(); ++Ti,++Ri)
                    _domain.addin(*Ri, *Ti);
            }

            return R;

    }

    // Standard middle product between iterator bounds
    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::stdmidmul(
        Rep& R, const RepIterator Rbeg, const RepIterator Rend,
        const Rep& P, const RepConstIterator Pbeg, const RepConstIterator Pend,
        const Rep& Q, const RepConstIterator Qbeg, const RepConstIterator Qend )
        const {

        //std::cerr << "Entering StdMidMul with "
        //          << "P[" << (Pbeg-P.begin()) << ":" << (Pend-P.begin()) << "/" << (P.end()-P.begin()) << "] × "
        //          << "Q[" << (Qbeg-Q.begin()) << ":" << (Qend-Q.begin()) << "/" << (Q.end()-Q.begin()) << "] → "
        //          << "R[" << (Rbeg-R.begin()) << ":" << (Rend-R.begin()) << "/" << (R.end()-R.begin()) << "]"
        //          << std::endl;

        if (Rbeg == Rend) return R;

        RepConstIterator ai=Pbeg,aig=Pbeg,bi=Qend; --bi;
        RepIterator ri=Rbeg;
        if (_domain.isZero(*bi))
            for(; (ai != Pend) && (ri != Rend); ++ai,++ri)
                _domain.assign(*ri,_domain.zero);
        else
            for(; (ai!=Pend) && (ri!=Rend);++ai,++ri)
                if (_domain.isZero(*ai))
                    _domain.assign(*ri, _domain.zero);
                else
                    _domain.mul(*ri,*ai,*bi);

        for(;ri!=Rend;++ri)
            _domain.assign(*ri,_domain.zero);

        for(--bi,++aig; (aig!=Pend) && (bi>=Qbeg);++aig,--bi)
            if (! _domain.isZero(*bi))
                for(ri=Rbeg,ai=aig; (ai!=Pend) && (ri!=Rend);++ai,++ri)
                    _domain.axpyin(*ri,*ai,*bi);

        return R;
    }

    // Karastuba balanced middle product between iterator bounds
    template <class Domain>
    inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::karamidmul(
        Rep& R, const RepIterator Rbeg, const RepIterator Rend,
        const Rep& P, const RepConstIterator Pbeg, const RepConstIterator Pend,
        const Rep& Q, const RepConstIterator Qbeg, const RepConstIterator Qend )
        const {

        // Assumes Pend-Pbeg == 2*(Qend-Qbeg)-1
        // Otherwise: undefined behavior!
        if (Rbeg == Rend) return R;

        if (Qend-Qbeg == 1) {
            _domain.mul(*Rbeg, *Pbeg, *Qbeg);
            return R;
        }

        // Initialize R to zero
        for(RepIterator ri=Rbeg; ri!= Rend; ++ri) _domain.assign(*ri,_domain.zero);

        const ssize_t n = (ssize_t) (Qend-Qbeg);
        const ssize_t n0 = n>>1;
        const ssize_t n1 = n0 + (n&1);
        const RepConstIterator P0end=Pbeg+2*n1-1, // end of P0
                               P1beg=Pbeg+n1, // beg of P1
                               P1plus=Pbeg+3*n1-1, // end of P1+
                               P1minus=Pbeg+n1+2*n0-1, // end of P1-
                               P2beg=Pbeg+2*n1, // beg of P2
                               Qmid=Qbeg+n0; // mid of Q = end of Q0 = beg of Q1
        const RepIterator      Rmid=Rbeg+n1; // mid of R = end of R0 = beg of R1

        //std::cerr << "Entering KaraMidMul " << (Pend-Pbeg) << "×" << (Qend-Qbeg) << "→" << (Rend-Rbeg) << std::endl;
        //std::cerr << "P0 = slice(P," << (Pbeg-P.begin()) << "," << (P0end-P.begin()) << ")" << std::endl;
        //std::cerr << "P1m= slice(P," << (P1beg-P.begin()) << "," << (P1minus-P.begin()) << ")" << std::endl;
        //std::cerr << "P1p= slice(P," << (P1beg-P.begin()) << "," << (P1plus-P.begin()) << ")" << std::endl;
        //std::cerr << "P2 = slice(P," << (P2beg-P.begin()) << "," << (Pend-P.begin()) << ")" << std::endl;
        //std::cerr << "Q0 = slice(Q," << (Qbeg-Q.begin()) << "," << (Qmid-Q.begin()) << ")" << std::endl;
        //std::cerr << "Q1 = slice(Q," << (Qmid-Q.begin()) << "," << (Qend-Q.begin()) << ")" << std::endl;
        //std::cerr << "R0 = slice(R," << (Rbeg-R.begin()) << "," << (Rmid-R.begin()) << ")" << std::endl;
        //std::cerr << "R1 = slice(R," << (Rmid-R.begin()) << "," << (Rend-R.begin()) << ")" << std::endl;


        Rep TMP, Rtmp;

        // R0 ← S0 = MP(P1+ + P0, Q1)
        TMP.resize(2*n1-1);
        RepIterator TMPi=TMP.begin();
        for(RepConstIterator P0i=Pbeg, P1i = P1beg;(P0i!=P0end) && (P1i!=P1plus); ++P0i, ++P1i, ++TMPi)
            _domain.add(*TMPi, *P0i, *P1i);
        midmul(R, Rbeg, Rmid,
               TMP, TMP.begin(), TMP.end(),
               Q, Qmid, Qend);

        // R1 ← S1 = MP(P1- + P2, Q0)
        TMP.resize(2*n0-1);
        TMPi = TMP.begin();
        for (RepConstIterator P1i = P1beg, P2i = P2beg; (P1i!=P1minus) && (P2i!=Pend); ++P1i, ++P2i, ++TMPi)
            _domain.add(*TMPi,*P1i,*P2i);
        midmul(R, Rmid, Rend,
               TMP, TMP.begin(), TMP.end(),
               Q, Qbeg, Qmid);

        // Rtmp ← S2 = MP(P1+, Q1- X^(n%2) Q0)
        TMP.resize(n1);
        RepConstIterator Q0i = Qbeg, Q1i = Qmid;
        TMPi = TMP.begin();
        if (n0 == n1) _domain.sub(*TMPi++,*Q1i++,*Q0i++);
        else          _domain.assign(*TMPi++,*Q1i++);
        for (;(Q0i!=Qmid) && (Q1i!=Qend); ++Q0i, ++Q1i, ++TMPi)
            _domain.sub(*TMPi,*Q1i,*Q0i);
        Rtmp.resize(n1);
        midmul(Rtmp, Rtmp.begin(), Rtmp.end(),
               P, P1beg, P1plus,
               TMP, TMP.begin(), TMP.end());

        // R0 ← R0 - S2
        RepConstIterator S2i = Rtmp.begin();
        for (RepIterator R0i = Rbeg; (R0i!=Rmid)&&(S2i!=Rtmp.end()); ++R0i, ++S2i)
            _domain.subin(*R0i, *S2i);
        // R1 ← R1 + S2*
        S2i = Rtmp.begin();
        for (RepIterator R1i = Rmid; (R1i!=Rend)&&(S2i!=Rtmp.end()); ++R1i, ++S2i)
            _domain.addin(*R1i, *S2i);

        return R;
    }


}
#endif // __GIVARO_poly1_mid_INL
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
