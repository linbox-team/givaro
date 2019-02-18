// ===============================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <05 May 10 09:02:07 Jean-Guillaume.Dumas@imag.fr>
// Author: J-G. Dumas
// Description: Pieces of polynomials as defined in
// [Arne Storjohann, High-Order Lifting
//  ISSAC'2002, pp. 246-254, ACM Press, July 2002].
// ===============================================================

#ifndef __GIVARO_trunc_domain_INL
#define __GIVARO_trunc_domain_INL

#include <givaro/givtruncdomain.h>

namespace Givaro {

    template <class Domain>
    inline typename TruncDom<Domain>::Rep& TruncDom<Domain>::truncin(Rep& p, const Degree& v, const Degree& d) const {
        Degree dP;degree(dP,p);
        Degree vP;val(vP,p);
        if (vP<v) {
            if (dP<v) {
                p.first.resize(0);
                p.second = 0;
                dP=Degree::deginfty;
            } else {
                size_t diVV=(size_t)value(v-vP);
                p.first.erase(p.first.begin(),p.first.begin()+(ssize_t)diVV);
                p.second = v;
            }
            vP=v;
        }
        if (d<dP) {
            if (d<vP) {
                p.first.resize(0);
                p.second = 0;
            } else {
                p.first.resize((size_t)value(d-vP)+1);
            }
        }
        return p;
    }


    template <class Domain>
    inline typename TruncDom<Domain>::Rep& TruncDom<Domain>::addin ( Rep& R, const Rep& P) const {
        Degree vP; val(vP, P);
        Degree vR; val(vR, R);
        size_t sP = P.first.size();
        if (sP == 0) return R;
        size_t sR = R.first.size();
        if (sR == 0) { return assign(R,P); }


        if (vR > vP) {
            expand(R, vP);
            sR = R.first.size();
            size_t diRP = (size_t)value(vR-vP);
            size_t i=0;
            for( ; (i<diRP) && (i<sP); ++i)
                this->_domain.assign(R.first[i],P.first[i]);
            if(sR < sP) {
                for(i=diRP; i<sR; ++i)
                    this->_domain.addin(R.first[i], P.first[i]);
                R.first.resize(sP);
                for(; i<sP; ++i)
                    this->_domain.assign(R.first[i],P.first[i]);
            } else {
                for(i=diRP; i<sP; ++i)
                    this->_domain.addin(R.first[i], P.first[i]);
            }
        } else {
            size_t diRP = (size_t)value(vP-vR);
            size_t newS = sP+diRP;
            if (newS>sR) R.first.resize(newS);
            size_t i = diRP, j=0;
            for( ; (i<sR) && (j<sP); ++i, ++j)
                this->_domain.addin(R.first[i], P.first[j]);
            for( ; (j<sP); ++i,++j)
                this->_domain.assign(R.first[i],P.first[j]);
        }
        return R;
    }

    template <class Domain>
    inline typename TruncDom<Domain>::Rep& TruncDom<Domain>::addin ( Rep& R, const Rep& P, const Degree& v, const Degree& d) const {
        Degree vP; val(vP, P);
        Degree vR; val(vR, R);
        size_t sP = P.first.size();
        if (sP == 0) return this->truncin(R,v,d);
        size_t sR = R.first.size();
        if (sR == 0) { return this->trunc(R,P,v,d); }
        // Degree dP; degree(dP, P); Degree dR; degree(dR, R); std::cout << "vR: " << vR << ", dR: "<< dR << ",sR: " << sR << ", vP: "<< vP << ", dP: "<< dP << ", sP: " << sP << ", v: "<< v << ", d: "<< d << std::endl;
        if (vP<=d) {
            if (vR<=d) {
                if (v<=vP) {
                    if (vP < vR) {
                        // v <= vP <= d  &&   v <= vP < vR
                        // v < vP <= vR
                        expand(R, vP);
                        // v < vP = vR
                        sR = R.first.size();
                        size_t diDV=(size_t)value(d-vP)+1;
                        sR = (sR>diDV?diDV:sR);
                        R.first.resize(sR);
                        size_t diDVP=(size_t)value(d-vP)+1;
                        sP=(sP>diDVP?diDVP:sP);
                        size_t i=0;
                        size_t diRP = (size_t)value(vR-vP);
                        if (sR < sP) R.first.resize(sP);
                        for( ; (i<diRP) && (i<sP); ++i)
                            this->_domain.assign(R.first[i],P.first[i]);
                        for( ; (i<sR) && (i<sP); ++i)
                            this->_domain.addin(R.first[i], P.first[i]);
                        for(; i<sP; ++i)
                            this->_domain.assign(R.first[i],P.first[i]);
                    } else {
                        // v <= vP <= d  &&   vR <= vP
                        if (vR < v) {
                            // vR < v <= vP <= d
                            if (static_cast<size_t>(value(v-vR))<sR) {
                                R.first.erase(R.first.begin(),R.first.begin()+value(v-vR));
                                sR = R.first.size();
                            } else {
                                R.first.resize(0);
                                sR = 0;
                            }
                            R.second = v;
                            vR=v;
                        }
                        // v = vR <= vP <= d
                        long diRP = value(vP-vR);
                        size_t diDV=(size_t)value(d-vP)+1;
                        sP=(sP>diDV?diDV:sP);
                        long newS = (long)(sP)+diRP;
                        long newsR = (long) R.first.size();
                        Degree dR(vR+newsR-1);
                        newsR = (dR>d?value(d-vR)+1:newsR);
                        newS = (newsR>newS?newsR:newS);
                        R.first.resize((size_t)newS);
                        size_t i = (size_t)diRP, j=0;
                        for( ; (i<sR) && (j<sP); ++i, ++j)
                            this->_domain.addin(R.first[i], P.first[j]);
                        for( ; (j<sP); ++i,++j)
                            this->_domain.assign(R.first[i],P.first[j]);
                    }
                } else {
                    // vP < v <= d
                    if (vP < vR) {
                        // vP < vR   &&   vP < v
                        if (vR < v) {
                            // vP < vR < v
                            Degree lastR(vR+(long)sR-1);
                            if (lastR < v) {
                                // vP <= vR < lastR < v
                                Degree lastP(vP+(long)sP-1);
                                if (lastP<v) {
                                    R.first.resize(0);
                                    R.second = 0;
                                } else {
                                    // vP < vR < lastR < v <= lastP
                                    lastP=(lastP>d?d:lastP);
                                    R.first.resize((size_t)value(lastP-v)+1);
                                    R.second = v;
                                    size_t i=0;
                                    size_t j=(size_t)value(v-vP);
                                    for( ;(j<sP) && (i<R.first.size());++i,++j)
                                        this->_domain.assign(R.first[i],P.first[j]);
                                }
                            }
                            else {
                                // vP <= vR < v <= lastR
                                R.first.erase(R.first.begin(),R.first.begin()+value(v-vR));
                                R.second = v;
                                sR=R.first.size();
                                size_t i = 0;
                                size_t j = (size_t)value(v-vP);
                                Degree lastP(vP+(long)sP-1);
                                lastP = (lastP>d?d:lastP);
                                Degree lastRR(v+(long)sR-1);
                                lastRR = (lastRR>d?d:lastRR);
                                size_t endi=(size_t)value(lastRR-v)+1;
                                size_t endj=(size_t)value(lastP-vP)+1;
                                lastRR = (lastP>lastRR?lastP:lastRR);
                                R.first.resize((size_t)value(lastRR-v)+1);
                                for( ; (j<endj) && (i<endi); ++i,++j)
                                    this->_domain.addin(R.first[i], P.first[j]);
                                for( ; j<endj; ++i,++j)
                                    this->_domain.assign(R.first[i],P.first[j]);
                            }
                        } else {
                            // vP < v <= vR <= d
                            if ((vP+(long)sP-1)<v) {
                                size_t diRD = (size_t)value(d-vR)+1;
                                if (diRD<sR) R.first.resize(diRD);
                            } else {
                                expand(R, v);
                                // vP < v = vR
                                sR = R.first.size();
                                size_t diDV=(size_t)value(d-v)+1;
                                sR = (sR>diDV?diDV:sR);
                                R.first.resize(sR);
                                size_t diDVP=(size_t)value(d-vP)+1;
                                sP=(sP>diDVP?diDVP:sP);
                                size_t diRP = (size_t)value(vR-v);
                                size_t i=0;
                                size_t j=(size_t)value(v-vP);
                                if (sR < (sP-j)) R.first.resize(sP-j);
                                for( ; (i<diRP) && (j<sP); ++i,++j)
                                    this->_domain.assign(R.first[i],P.first[j]);
                                for( ; (i<sR) && (j<sP); ++i,++j)
                                    this->_domain.addin(R.first[i], P.first[j]);
                                for(; j<sP; ++i,++j)
                                    this->_domain.assign(R.first[i],P.first[j]);
                            }
                        }
                    } else {
                        // vR <= vP < v <= d
                        Degree lastR(vR+(long)sR-1);
                        if (lastR < v) {
                            // vR <= vP < lastR < v
                            Degree lastP(vP+(long)sP-1);
                            if (lastP<v) {
                                R.first.resize(0);
                                R.second = 0;
                            } else {
                                // vR <= vP < lastR < v <= lastP
                                lastP=(lastP>d?d:lastP);
                                R.first.resize((size_t)value(lastP-v)+1);
                                R.second = v;
                                size_t i=0;
                                size_t j=(size_t)value(v-vP);
                                for( ;(j<sP) && (i<R.first.size());++i,++j)
                                    this->_domain.assign(R.first[i],P.first[j]);
                            }
                        } else {
                            // vR <= vP < v <= lastR
                            R.first.erase(R.first.begin(),R.first.begin()+value(v-vR));
                            R.second = v;
                            sR=R.first.size();
                            size_t i = 0;
                            size_t j = (size_t)value(v-vP);
                            Degree lastP(vP+(long)sP-1);
                            lastP = (lastP>d?d:lastP);
                            Degree lastRR(v+(long)sR-1);
                            lastRR = (lastRR>d?d:lastRR);
                            size_t endi=(size_t)value(lastRR-v)+1;
                            size_t endj=(size_t)value(lastP-vP)+1;
                            lastRR=(lastP>lastRR?lastP:lastRR);
                            R.first.resize((size_t)value(lastRR-v)+1);
                            for( ; (j<endj) && (i<endi); ++i,++j)
                                this->_domain.addin(R.first[i], P.first[j]);
                            for( ; j<endj; ++i,++j)
                                this->_domain.assign(R.first[i],P.first[j]);
                        }
                    }
                }
            } else {
                // vP <= d < vR
                size_t j=0;
                if (vP<v) {
                    size_t diPV = (size_t)value(v-vP)+1;
                    if (sP<diPV) {
                        R.first.resize(0);
                        R.second = 0;
                        return R;
                    }
                    j+=(size_t)value(v-vP);
                    vP=v;
                }
                size_t i=0;
                sR=(size_t)value(d-vP)+1;
                size_t inP = sP-j;
                sR=(inP>sR?sR:inP);
                R.first.resize(sR);
                for( ; i<sR; ++i,++j)
                    this->_domain.assign(R.first[i],P.first[j]);
                R.second = vP;
            }
        } else {
            if (d<vR) {
                // d < vR,vP
                R.first.resize(0);
                R.second = 0;
            } else {
                // vR <= d < vP
                if (vR<v) {
                    if (static_cast<size_t>(value(v-vR))<sR) {
                        R.first.erase(R.first.begin(),R.first.begin()+value(v-vR));
                        R.second = v;
                        R.first.resize((size_t)value(d-v)+1);
                    } else {
                        R.first.resize(0);
                        R.second = 0;
                    }
                } else {
                    R.first.resize((size_t)value(d-vR)+1);
                }
            }
        }
        return R;
    }

    template <class Domain>
    inline typename TruncDom<Domain>::Rep& TruncDom<Domain>::subin ( Rep& R, const Rep& P) const {
        size_t sP = P.first.size();
        if (sP == 0) return R;
        size_t sR = R.first.size();
        if (sR == 0) { return assign(R,P); }
        Degree vP; val(vP, P);
        Degree vR; val(vR, R);
        if (vR > vP) {
            expand(R, vP);
            sR = R.first.size();
            size_t diRP = (size_t)value(vR-vP);
            size_t i=0;
            for( ; (i<diRP) && (i<sP); ++i)
                this->_domain.neg(R.first[i],P.first[i]);
            if(sR < sP) {
                for(i=diRP; i<sR; ++i)
                    this->_domain.subin(R.first[i], P.first[i]);
                R.first.resize(sP);
                for(; i<sP; ++i)
                    this->_domain.neg(R.first[i],P.first[i]);
            } else {
                for(i=diRP; i<sP; ++i)
                    this->_domain.subin(R.first[i], P.first[i]);
            }
        } else {
            size_t diRP = (size_t)value(vP-vR);
            size_t newS = sP+diRP;
            if (newS>sR) R.first.resize(newS);
            size_t i = diRP, j=0;
            for( ; (i<sR) && (j<sP); ++i, ++j)
                this->_domain.subin(R.first[i], P.first[j]);
            for( ; (j<sP); ++i,++j)
                this->_domain.neg(R.first[i],P.first[j]);
        }
        return R;
    }

    template <class Domain>
    inline typename TruncDom<Domain>::Rep& TruncDom<Domain>::mul( Rep& r, const Rep& u, const Rep& v, const Degree& Val, const Degree& deg) const {
        Degree vU; this->val(vU, u);
        Degree vV; this->val(vV, v);
        r.second = u.second+v.second;
        if (Val >= r.second) {
            Father_t::mul(r.first,u.first,v.first,Val-r.second,deg-r.second);
            r.second=Val.value();
        } else if (deg >= r.second) {
            Father_t::mul(r.first,u.first,v.first,Degree(0),deg-r.second);
        } else {
            //         r.first=this->zero;
            //         r.second=0;
            r=this->zero;
        }
        return r;
    }

} // Givaro

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
