// ===============================================================
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <15 Apr 11 16:37:12 Jean-Guillaume.Dumas@imag.fr>
// Author: J-G. Dumas
// Description: Storjohann's high-order lifting
// Reference:   A. Storjohann. High-order lifting. ISSAC 2002.
// ===============================================================
#ifndef __GIVARO_highorder_H
#define __GIVARO_highorder_H

#ifndef _GIVARO_HIGHORDER_THRESHOLD
#define _GIVARO_HIGHORDER_THRESHOLD 30
#endif


#include <givaro/givpoly1.h>
#include <givaro/givfractiondomain.h>
#include <givaro/givtruncdomain.h>

#include <givaro/givpower.h>
#include <givaro/givquotientdomain.h>


namespace Givaro {

    template<class Domain>
    struct HighOrder {

        // -- Self_t
        typedef HighOrder<Domain>			Self_t;

        // -- Exported types
        typedef FracDom< Poly1Dom<Domain,Dense> >	Father_t;
        typedef FracDom< Poly1Dom<Domain,Dense> >	Frac_t;

        typedef Poly1Dom<Domain,Dense>		Ring_t;
        typedef Poly1Dom<Domain,Dense>		Poly_t;
        typedef TruncDom<Domain> 			Trunc_t;
        typedef typename Trunc_t::Element	Truncated;
        typedef typename Ring_t::Element	Ring_E;
        typedef typename Ring_t::Element	Polynomial;
        typedef Frac<Ring_E, Ring_E>		Element;
        typedef Frac<Ring_E, Ring_E>		Rep;

        typedef Domain				Domain_t;
        typedef typename Domain::Element	Type_t;

        const Poly_t	_poldom;
        const Domain& _dom;
        const Trunc_t _truncdom;

        const Poly_t& getpoldom() const { return _poldom; }
        const Domain& getdomain() const { return _dom; }
        const Trunc_t& gettruncdom() const { return _truncdom; }



        HighOrder(const Domain& d, const Indeter& Z = Indeter() ) : _poldom(d, Z), _dom(_poldom.getdomain()), _truncdom(_poldom) {}

        Polynomial& taylor(Polynomial& Tay, const Rep& Fra, Degree order) const {
            Degree d; _poldom.degree(d,Fra._den);
            Tay.resize((size_t)order.value()+1);
            size_t i = 0;
            for( ; (i<Fra._num.size()) && (order>=(long)i); ++i) {
                _dom.assign(Tay[i],Fra._num[i]);
                for(size_t j = 1; (j<Fra._den.size()) && (j<=i); ++j) {
                    _dom.maxpyin(Tay[i],Fra._den[j],Tay[i-j]);
                }
                _dom.divin(Tay[i], Fra._den.front());
            }
            for( ; (order>=(long)i); ++i) {
                _dom.assign(Tay[i], _dom.zero);
                for(size_t j = 1; (j<Fra._den.size()) && (j<=i); ++j) {
                    _dom.maxpyin(Tay[i],Fra._den[j],Tay[i-j]);
                }
                _dom.divin(Tay[i], Fra._den.front());
            }
            return _poldom.setdegree(Tay);
        }


        Truncated& Fiduccia(Truncated& F, const Rep& Fra, Degree a, Degree b) const {
            Polynomial Tay;
            Degree dA; _poldom.degree(dA, Fra._den);
            this->taylor(Tay, Fra, dA);
            return Fiduccia(F, Tay, Fra._den, a, b);
        }

        Truncated& Fiduccia(Truncated& F, const Polynomial& Tay, const Polynomial& FraDen, Degree a, Degree b) const {
            Degree dA; _poldom.degree(dA, FraDen);
            Polynomial Rev; _poldom.init(Rev,dA);
            for(int i=0;i<dA.value();++i)
                _dom.div(Rev[(size_t)i],(FraDen)[(size_t)dA.value()-(size_t)i],FraDen.front());
            _dom.assign(Rev[(size_t)dA.value()],_dom.one);
            return FiducciaReversed(F, Tay, Rev, a, b);
        }

        Truncated& Fiduccia(Truncated& F, const Polynomial& Tay, const Polynomial& FraDen, Degree b) const {
            return Fiduccia(F, Tay, FraDen, b, b);
        }

        Truncated& Fiduccia(Truncated& F, const Rep& Fra, Degree b) const {
            return Fiduccia(F, Fra, b, b);
        }


        Type_t& shifteddotproduct(Type_t& dp, const Polynomial& P, const Polynomial& Q, const long shift) const {
            Degree dl; _poldom.degree(dl, P);

            GIVARO_STATE( Degree dq; _poldom.degree(dq, Q);  );
            GIVARO_ASSERT( (dl.value()+shift) <= (dq.value()+1), " in dotproduct HighOrder dP: " << dl << ", dQ: " << dq << ", shift: " << shift );

            _dom.assign(dp,_dom.zero);
            for(long i=0;i<=dl.value();++i)
                _dom.axpyin(dp, P[(size_t)i], Q[(size_t)(i+shift)]);

            return dp;
        }



        Truncated& FiducciaReversed(Truncated& F, const Polynomial& Tay, const Polynomial& FraDen, Degree a, Degree b) const {

            Degree dT; _poldom.degree(dT, Tay);
            //         std::cout << "a: " << a << ", b: " << b << ", dT: " << dT << std::endl;
            if (b > dT) {
                Degree dA; _poldom.degree(dA, FraDen);

                Polynomial Xl; _poldom.init(Xl);
                Polynomial Xone; _poldom.init(Xone,Degree(1));
                QuotientDom<Poly_t> Qdom(_poldom, FraDen);
                Qdom.assign(Xone);

                const long bonus = dT.value()-dA.value()+1;

                Degree dR(b-a);
                Polynomial Res; _poldom.init(Res, dR);
                long iterT=a.value();
                long iterR=0;

                for( ; iterT < dT.value(); ++iterR,++iterT)
                    _dom.assign(Res[(size_t)iterR],Tay[(size_t)iterT]);

                dom_power(Xl, Xone, iterT-bonus, Qdom);

                Degree dl; _poldom.degree(dl, Xl);
                shifteddotproduct( Res[(size_t)iterR], Xl, Tay, bonus);


                for( ++iterR; iterR<=dR.value(); ++iterR) {
                    Qdom.mulin(Xl, Xone);
                    shifteddotproduct( Res[(size_t)iterR], Xl, Tay, bonus);
                }
                _truncdom.assign(F, Res);
                return _truncdom.mulin(F, a);
            } else {
                return _truncdom.assign(F, Tay, Degree(a), Degree(b));
            }
        }

        Truncated& FiducciaReversed(Truncated& F, const Polynomial& Tay, const Polynomial& FraDen, Degree b) const {
            return FiducciaReversed(F, Tay, FraDen, b, b);
        }


        //     Truncated& FiducciaReversed(Truncated& F, const Polynomial& Tay, const Polynomial& FraDen, Degree b) const {
        // 	Degree dT; _poldom.degree(dT, Tay);
        // 	if (b > dT) {
        // 	  Degree dA; _poldom.degree(dA, FraDen);

        // 	  Polynomial Xl; _poldom.init(Xl);
        // 	  Polynomial Xone; _poldom.init(Xone,Degree(1));
        // 	  QuotientDom<Poly_t> Qdom(_poldom, FraDen);
        // 	  Qdom.assign(Xone);

        // 	  const long bonus = dT.value()-dA.value()+1;

        // 	  dom_power(Xl, Xone, b.value()-bonus, Qdom);

        //  	  Degree dl; _poldom.degree(dl, Xl);

        //           Type_t Tl; _dom.init(Tl);
        //           shifteddotproduct( Tl, Xl, Tay, bonus);

        // 	  return _truncdom.assign(F, Degree(b), Tl);
        // 	} else {
        // 	  return _truncdom.assign(F, Degree(b), Tay[b.value()]);
        // 	}
        //     }


        Truncated& FracDevel(Truncated& F, const Rep& Fra, Degree a, Degree b) const {

            // Precondition Degree(Fra._num)<Degree(Fra._den)
            if (b < _GIVARO_HIGHORDER_THRESHOLD) {
                Polynomial Tay;
                this->taylor(Tay, Fra, b);
                return _truncdom.assign(F, Tay, a, b);
            } else {
                std::vector<Truncated> Gam, T;
                std::vector<Degree> Deg;
                Polynomial Tay; Degree dT;
                Degree dA; _poldom.degree(dA, Fra._den);

#ifdef GIVARO_HIGHORDER_TIMER
                Timer GivHOTimer; GivHOTimer.clear(); GivHOTimer.start();
#endif
                this->highorder(Gam, T, Deg, Tay, dT, a, b, Fra._den, dA);
#ifdef GIVARO_HIGHORDER_TIMER
                GivHOTimer.stop();
                std::cerr << "HighOrder_" << a << '^' << b << " : " << GivHOTimer << std::endl;
#endif

                Truncated A, B;
                _truncdom.assign(A, Fra._den);
                _truncdom.assign(B, Fra._num);
                return this->FracDevel(F, B, A, dA, a, b, Tay, dT, Gam, T, Deg);
            }
        }

        Truncated& GammaId(Truncated& Gam, const Polynomial& Tay, const Degree dT, const Degree k0, const Polynomial& A, const Degree dA) const {
            Truncated U0,TA,One;
            _truncdom.assign(TA, A);
            _truncdom.assign(One,_poldom.one );
            _truncdom.assign(U0, Tay, k0-dA,k0-1);
            _truncdom.maxpy(Gam, TA, U0, One, k0, k0+dA-1);
            _truncdom.divin(Gam, k0);
            return Gam;
        }

        Truncated& doubleorder(Truncated& Gam, Truncated& T, const Truncated& Gamp, const Truncated& Tp, const Truncated& S, const Truncated& TA, const Degree dA, const Degree ke) const {
            //         std::cerr << "Double order: " << ke << std::endl;
            // ke = 2^e
            Truncated AL;
            _truncdom.mul(AL, Tp, Gamp, ke-dA, ke-1);
            _truncdom.mulin(AL, ke-dA);
            //write(std::cout << "AL" << ke.value() << ":=", AL) << ';' << std::endl;

            _truncdom.mul(Gam, TA, AL, ke*2-dA, ke*2-1);
            _truncdom.negin(Gam);
            _truncdom.divin(Gam, ke*2-dA);
            //         _truncdom.setval(Gam);
            // write(std::cout << "Gam" << (ke*2-dA) << ":=", Gam) << ';' << std::endl;
            // std::cerr << "[G]_" << (ke*2-dA) << std::endl;

            Truncated AH;
            _truncdom.mul(AH, Gam, S, 0, dA-1);
            _truncdom.mulin(AH,ke*2-dA);
            //write(std::cout << "AH" << ke.value() << ":=", AH) << ';' << std::endl;

            // Forget the last term
            _truncdom.truncin(AL, ke*2-dA*2+1, ke*2-dA-1);
            _truncdom.add(T, AL, AH);
            // std::cerr << "[I]_" << (ke*2-dA*2+1) << '^' << (ke*2-1) << std::endl;

            // write(std::cout << "T" << (ke*2-dA) << ":=", T) << ';' << std::endl;

            return Gam;
        }

        std::vector<Truncated>& highorder(std::vector<Truncated>& Gam, std::vector<Truncated>& T, std::vector<Degree>& Deg, Polynomial& Tay, Degree& dT, Degree a, Degree order, const Polynomial& A, const Degree dA) const {
            Gam.resize(0); T.resize(0); Deg.resize(0);
            size_t e;
            for(e=0; (1<<e)<dA.value(); ++e) {}
            ++e; // 2^{e-2} < d <= 2^{e-1}
        int dt = int (1<<e);

        Degree k0 = Degree(1<<e);
        Deg.push_back(k0-dA);

        Degree dif = order-a;
        dt = int(dif.value()>dt? dif.value() : dt);
        // std::cout << "BEG HighOrder" << std::endl;
        // std::cout << "d: " << dA << std::endl;
        // std::cout << "e: " << e << std::endl;
        // std::cout << "k0: " << k0 << std::endl;
        // std::cout << "k0-d: " << Deg.back() << std::endl;
        // std::cout << "a: " << a << std::endl;
        // std::cout << "b: " << order << std::endl;
        // std::cout << "dif: " << dif << std::endl;
        // std::cout << "dt: " << dt << std::endl;


        Rep Fra; Fra._num=_poldom.one; Fra._den=A;
        this->taylor(Tay, Fra, dt);
        // std::cerr << "[I]_" << 0 << '^' << dt << std::endl;
        // write(std::cout << "Tay" << dt << ":=", Tay) << ';' << std::endl;
        dT = dt;

        Truncated G0; this->GammaId(G0, Tay, dT, Deg.back(), A, dA);
        Gam.push_back(G0);
        // std::cerr << "[G]_" << Deg.back() << std::endl;

        // write(std::cout << "Gam" << Deg.back() << ":=", G0) << ';' << std::endl;

        Truncated S; _truncdom.assign(S, Tay, Degree(0), dA-1);
        // write(std::cout << "S:=", S) << ';' << std::endl;


        Truncated T0; _truncdom.assign(T0,Tay, k0-dA*2+1,k0-1);
        T.push_back(T0);
        // write(std::cout << "T" << e << ":=", T0) << ';' << std::endl;

        Truncated TA; _truncdom.assign(TA, A);
        size_t ordero2 = (size_t)order.value()/2;


        for( ; k0<(long)ordero2; ) {
#ifdef GIVARO_HIGHORDER_TIMER
            Timer GivHOTimer; GivHOTimer.clear(); GivHOTimer.start();
#endif
            doubleorder(G0, T0, Gam.back(), T.back(), S, TA, dA, k0);
            Gam.push_back(G0);
            T.push_back(T0);
            k0*=2; Deg.push_back(k0-dA);
#ifdef GIVARO_HIGHORDER_TIMER
            GivHOTimer.stop();
            std::cerr << "DoubleOrder[" << dA << "]^" << k0 << " out of " << ordero2 << " : " << GivHOTimer << std::endl;
#endif
        }
        return Gam;
        }

        Truncated& Betta(Truncated& B, const Truncated& TB, const Truncated& TA, const Degree dA, const Degree a, const Polynomial& Tay, const Degree dT, const std::vector<Truncated>& Gam, const std::vector<Truncated>& T, const std::vector<Degree>& Deg) const {
            // this->write(std::cerr << "BEG Betta" << a << "TB: [", TB) << "]_" << a << std::endl;
#ifdef GIVARO_HIGHORDER_TIMER
            Timer GivHOTimer; GivHOTimer.clear(); GivHOTimer.start();
#endif

            Degree a0=dA*2-1;
            a0 = (a0>a? 0 : a-a0);
            Truncated S;
            this->Inverse(S, a0, a-1, TA, dA, Tay, dT, Gam, T, Deg);
            Truncated U;
            _truncdom.mul(U, S, TB, a-dA, a-1);
            //   _truncdom.setval(U);

            _truncdom.maxpy(B, TA, U, TB, a, a+dA-1);
            // this->write(std::cerr << "END Betta" << a << " B: [", B) << "]/" << a << std::endl;
#ifdef GIVARO_HIGHORDER_TIMER
            GivHOTimer.stop();
            std::cerr << "Betta[" << dA << "]_" << a << " : " << GivHOTimer << std::endl;
#endif
            return _truncdom.divin(B, a);
        }


        Truncated& FracDevel(Truncated& F, const Truncated& TB, const Truncated& TA, const Degree dA, const Degree a, const Degree b, const Polynomial& Tay, const Degree dT, const std::vector<Truncated>& Gam, const std::vector<Truncated>& T, const std::vector<Degree>& Deg) const {
            // this->write(this->write(std::cerr << "BEG FracDevel" << a << '-' << b << ": [", TB) << " x ", TA)<< "]_" << a << '^' << b << std::endl;
            if (b <= dT) {
                Truncated tTay; _truncdom.assign(tTay, Tay, a-dA-1, b);
                //             this->write(std::cerr << "END FracDevel" << a << '-' << b << " b<" << dT << ": ", _truncdom.mul(F, TB, tTay, a, b)) << std::endl;
                return _truncdom.mul(F, TB, tTay, a, b);
            } else if (a == 0) {
                std::cerr << "ERROR: This should not happen" << std::endl;
                std::cerr <<"[FracDevel" << a << '-' << b << "]_0^" << b.value() << " with degree(taylor)=" << dT << std::endl;
                return _truncdom.assign(F,_poldom.zero );
            } else {
                Truncated Beta;
                this->Betta(Beta, TB, TA, dA, a, Tay, dT, Gam, T, Deg);
                //             _truncdom.setval(Beta);

                //             this->FracDevel(F, Beta, TA, dA, Degree(0), b-a, Tay, dT, Gam, T, Deg);
                Truncated tTay; _truncdom.assign(tTay, Tay, 0, b-a);
                _truncdom.mul(F, Beta, tTay, 0, b-a);
                //             _truncdom.setval(F);
                //             this->write(std::cerr << "END FracDevel" << a << '-' << b << " : [", F) << "]/" << a << std::endl;
                return _truncdom.mulin(F,a);
            }
        }


        Truncated& Inverse(Truncated& I, const Degree a, const Degree b, const Truncated& TA, const Degree dA, const Polynomial& Tay, const Degree dT, const std::vector<Truncated>& Gam, const std::vector<Truncated>& T, const std::vector<Degree>& Deg) const {
            // std::cerr << "[I]_" << a << '^' << b << std::endl;

            // this->write(std::cerr << "BEG Inverse" << a << '-' << b << ": [", TA) << "]_" << a << '^' << b << std::endl;
            if (b <= dT) {
                // this->write(std::cerr << "END Inverse" << a << '-' << b << " b<" << dT << ": ", _truncdom.assign(I, Tay, a, b) ) << std::endl;
                return _truncdom.assign(I, Tay, a, b);
            } else {
                Truncated Gama;
                this->Gamma(Gama, a, TA, dA, Tay, dT, Gam, T, Deg);
                //             _truncdom.setval(Gama);
                //             this->FracDevel(I, Gama, TA, dA, Degree(0), b-a, Tay, dT, Gam, T, Deg);
                Truncated tTay; _truncdom.assign(tTay, Tay, 0, b-a);
                _truncdom.mul(I, Gama, tTay, 0, b-a);

                //             _truncdom.setval(I);
                // this->write(std::cerr << "END Inverse" << a << '-' << b << ": [", I ) << "]*" << a << std::endl;
                return _truncdom.mulin(I, a);
            }
        }


        Truncated& Gamma(Truncated& G, const Degree a, const Truncated& TA, const Degree dA, const Polynomial& Tay, const Degree dT, const std::vector<Truncated>& Gam, const std::vector<Truncated>& T, const std::vector<Degree>& Deg) const {
            // std::cerr << "[G]_" << a << std::endl;
            // this->write(std::cerr << "BEG Gamma" << a << ": [", TA) << "]_" << a << std::endl;
            //         if ((a+dA)>dT) {
            if (a>dT) {
                size_t i=0;
                for( ; i<Deg.size(); ++i) {
                    if (Deg[i] > a) { --i; break; }
                }
                if (i>=Deg.size()) --i;

                this->Betta(G, Gam[i], TA, dA, a-Deg[i], Tay, dT, Gam, T, Deg);
                // this->write(std::cerr << "END Gamma" << a << " B: [", G) << "]" << std::endl;
                return _truncdom.setval(G);
            } else if (a == 0) {
                // this->write(std::cerr << "END Gamma" << a << " O: ", _truncdom.one) << std::endl;
                return _truncdom.assign(G, _truncdom.one);
            } else {
                Truncated One; _truncdom.assign(One,_poldom.one );
                Truncated tTay; _truncdom.assign(tTay, Tay, 0, a-1);
                _truncdom.maxpy(G, TA, tTay, One,a,a+dA-1);
                // this->write(std::cerr << "END Gamma" << a << " D: [", G) << "]/" << a << std::endl;
                return _truncdom.divin(G, a);
            }
        }



        std::ostream& write( std::ostream& o) const {
            return _truncdom.write(o << "HighOrder<") << '>';
        }
        std::ostream& write( std::ostream& o, const Truncated& n) const {
            return _truncdom.write(o, n);
        }
        std::ostream& write( std::ostream& o, const Polynomial& n) const {
            return _poldom.write(o, n);
        }
        std::ostream& write( std::ostream& o, const Rep& n) const {
            return _poldom.write(_poldom.write(o << '(',n._num) << ")/(", n._den) << ')';
        }



        };

    } // Givaro

#endif
    /* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
    // vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
