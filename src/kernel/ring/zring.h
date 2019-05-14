// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: W. J. Turner <wjturner@acm.org>
//          Bradford Hovinen <hovinen@cis.udel.edu>
//          Clement Pernet <clement.pernet@gmail.com> (inserted into FFLAS-FFPACK)
//          A. Breust (taken from FFLAS-FFPACK)
// ==========================================================================


/*! @file field/zring.h
 * @ingroup field
 * @brief  representation of a field of characteristic 0.
 */

#ifndef __GIVARO_ring_zring_H
#define __GIVARO_ring_zring_H

#include <algorithm>
#include <typeinfo>
#include <math.h>

#include "givaro/unparametric-operations.h"
#include "givaro/givranditer.h"
#include "givaro/givinteger.h"

namespace Givaro
{
    template<class _Element> class ZRing;

    template<typename Domain> struct DomainRandIter {
        typedef GeneralRingRandIter<Domain> RandIter;
    };

    template<> struct DomainRandIter<ZRing<Integer>> {
        typedef IntegerDom::RandIter RandIter;
    };


    /** Class ZRing.
     *  Ring of integers, using the _Element base type.
     */
    template<class _Element>
    class ZRing : public UnparametricOperations<_Element>
    {
    public:

        // ----- Exported Types and constantes
        using Element = _Element;
        using Rep = _Element;
        using Self_t = ZRing<Element>;
        using Residu_t = int64_t;
        using Element_ptr = Element*;
        using ConstElement_ptr = const Element*;
        enum { size_rep = sizeof(Element) };

        const Element one  = 1;
        const Element zero = 0;
        const Element mOne = -1;

        //----- Constructors
        ZRing() {}
        ZRing(const ZRing& F) {}
        // Needed in FFLAS, when ZRing is used as delayed field.
        template<class T> ZRing(const T&) {}

        //----- Access
        Residu_t residu() const { return static_cast<Residu_t>(0); }
        Residu_t size() const { return static_cast<Residu_t>(0); }
        Residu_t cardinality() const { return static_cast<Residu_t>(0); }
        Residu_t characteristic() const { return static_cast<Residu_t>(0); }
        template<typename T> T& cardinality(T& c) const { return c = static_cast<T>(0); }
        template<typename T> T& characteristic(T& c) const { return c = static_cast<T>(0); }

        static inline Residu_t maxCardinality() { return -1; }
        static inline Residu_t minCardinality() { return 2; }

        //----- Ring-wise operations
        inline bool operator==(const Self_t& F) const { return true; }
        inline bool operator!=(const Self_t& F) const { return false; }
        inline ZRing<Element>& operator=(const ZRing<Element>&) { return *this; }
        // Ring tests
        bool isZero(const Element& a) const { return a == zero; }
        bool isOne (const Element& a) const { return a == one; }
        bool isMOne(const Element& a) const { return a == mOne; }
        bool isUnit(const Element& a) const { return isOne(a) || isMOne(a); }

        Element& abs(Element& x, const Element& a) const {return x= (a>0)? a: -a;}
        Element abs(const Element& a) const {return (a>0)? a: -a;}
        long compare (const Element& a, const Element& b) const {return (a>b)? 1: ((a<b)? -1 : 0);}
        Element& gcd (Element& g, const Element& a, const Element& b) const {return Givaro::gcd(g,a,b);}
        Element& gcdin (Element& g, const Element& b) const{return gcd(g, g, b);}
        Element& gcd (Element& g, Element& s, Element& t, const Element& a, const Element& b) const{return Givaro::gcd(g,s,t,a,b);}
        Element &dxgcd(Element &g, Element &s, Element &t, Element &u, Element &v, const Element &a, const Element &b) const
        {
            gcd(g,s,t,a,b);
            div(u,a,g);
            div(v,b,g);
            return g;
        }
        Element& lcm (Element& c, const Element& a, const Element& b) const
        {
            if ((a==Element(0)) || (b==Element(0))) return c = Element(0);
            else {
                Element g;
                gcd (g, a, b);
                c= a*b;
                c /= g;
                c=abs (c);
                return c;
            }
        }

        Element& lcmin (Element& l, const Element& b) const
        {
            if ((l==Element(0)) || (b==Element(0))) return l = Element(0);
            else {
                Element g;
                gcd (g, l, b);
                l*= b;
                l/= g;
                l=abs (l);
                return l;
            }
        }

        void reconstructRational (Element& a, Element& b, const Element& x, const Element& m) const
        {this->RationalReconstruction(a,b, x, m, Givaro::sqrt(m), true, true);}
        void reconstructRational (Element& a, Element& b, const Element& x, const Element& m, const Element& bound) const
        {this->RationalReconstruction(a,b, x, m, bound, true, true);}
        bool reconstructRational (Element& a, Element& b, const Element& x, const Element& m, const Element& a_bound, const Element& b_bound) const
        {
            Element bound = x/b_bound;
            this->RationalReconstruction(a,b,x,m, (bound>a_bound?bound:a_bound), true, false);
            return b <= b_bound;
        }
        Element& quo (Element& q, const Element& a, const Element& b) const {return Integer::floor(q, a, b);}
        Element& rem (Element& r, const Element& a, const Element& b) const {return Integer::mod(r,a,b);}
        Element& quoin (Element& a, const Element& b) const{return quo(a,a,b);}
        Element& remin (Element& a, const Element& b) const {return rem(a,a,b);}
        void quoRem (Element& q, Element& r, const Element& a, const Element& b) const{Integer::divmod(q,r,a,b);}
        bool isDivisor (const Element& a, const Element& b) const
        {
            Element r;
            return rem(r,a,b)==Element(0);
        }

        inline Element& sqrt(Element& x, const Element& y) const{return Givaro::sqrt(x,y);}
        inline  Element powtwo(Element& z, const Element& x) const
        {
            z = 1;
            Element max; init(max, (int64_t)(1<<30));
            if (x < 0) return z;
            //if (x < (Element)max-1) {
            if (x < max) {
                z<<=(int32_t)x;
                return z;
            }
            else {
                Element n,m;
                quoRem(n,m,x,max);
                //quoRem(n,m,x,(Element)(LONG_MAX-1));
                for (int i=0; i < n; ++i) {
                    z <<=(int32_t)max;
                    //z <<=(long int)(LONG_MAX-1);
                }
                z <<= (int32_t)m;
                return z;
            }
            //for (Element i=0; i < x; ++i) {
            //      z <<= 1;
            //}
            //return z; // BB peut pas !
        }

        inline  Element logtwo(Element& z, const Element& x) const {return z = x.bitsize() - 1;}


        //----- Initialisation
        Element& init(Element& x) const { return x; }
        template <typename T> Element& init(Element& x, const T& s) const
        { return Caster(x,s); }

        Element& assign(Element& x, const Element& y) const { return x = y; }

        //----- Convert
        template <typename T> T& convert(T& x, const Element& y) const
        { return Caster(x,y); }

        Element& reduce (Element& x, const Element& y) const { return x = y; }
        Element& reduce (Element& x) const { return x; }

        // To ensure interface consistency
        Element minElement() const { return 0; }
        Element maxElement() const { return 0; }


        // ----- Random generators
        typedef typename DomainRandIter<Self_t>::RandIter RandIter;
        typedef GeneralRingNonZeroRandIter<Self_t> NonZeroRandIter;
        template< class Random > Element& random(const Random& g, Element& r) const
        { return init(r, g()); }
        template< class Random > Element& nonzerorandom(const Random& g, Element& a) const
        { while (isZero(init(a, g())))
            ;
            return a; }

        //----- IO
        std::ostream& write(std::ostream &os) const
        {
            return os << "ZRing<" << typeid(Element).name() << ">";
        }
        std::ostream& write(std::ostream &os, const Element& a) const
        {
            return os << a;
        }
        std::istream& read(std::istream &is, Element& a) const
        {
            return is >> a;
        }
        protected:
        /*! Rational number reconstruction.
         * \f$\frac{n}{d} \equiv f \mod m\f$, with \f$\vert n
         \vert <k\f$ and \f$0 < \vert d \vert \leq \frac{f}{k}\f$.
         * @bib
         * - von zur Gathen & Gerhard, <i>Modern Computer Algebra</i>,
         *      5.10, Cambridge Univ. Press 1999
         */
        inline void RationalReconstruction( Element& a, Element& b,
                                            const Element& f, const Element& m,
                                            const Element& k,
                                            bool forcereduce, bool recursive ) const
        {
            Element x(f);
            if (x<0) {
                if ((-x)>m)
                    x %= m;
                if (x<0)
                    x += m;
            }
            else {
                if (x>m)
                    x %= m;
            }

            if (x == 0) {
                a = 0;
                b = 1;
            }
            else {
                bool res = ratrecon(a,b,x,m,k, forcereduce, recursive);
                if (recursive)
                    for( Element newk = k + 1; (!res) && (newk<f) ; ++newk)
                        res = ratrecon(a,b,x,m,newk,forcereduce, true);
            }
        }

        // Precondition f is suppposed strictly positive and strictly less than m
        inline  bool ratrecon( Element& num, Element& den,
                               const Element& f, const Element& m,
                               const Element& k,
                               bool forcereduce, bool recursive ) const
        {

            //std::cerr << "RatRecon : " << f << " " << m << " " << k << std::endl;
            Element  r0, t0, q, u;
            r0=m;
            t0=0;
            num=f;
            den=1;
            while(num>=k)
            {
                q = r0;
                q /= num;   // r0/num
                u = num;
                num = r0;  	// num <-- r0
                r0 = u;	// r0 <-- num
                //maxpyin(num,u,q);
                Integer::maxpyin(num,u,q);
                if (num == 0) return false;

                u = den;
                den = t0;  	// num <-- r0
                t0 = u;	// r0 <-- num
                //maxpyin(den,u,q);
                Integer::maxpyin(den,u,q);
            }

            if (forcereduce) {
                // [GG, MCA, 1999] Theorem 5.26

                // (ii)
                Element gg;
                if (gcd(gg,num,den) != 1) {

                    Element ganum, gar2;
                    for( q = 1, ganum = r0-num, gar2 = r0 ; (ganum < k) && (gar2>=k); ++q ) {
                        ganum -= num;
                        gar2 -= num;
                    }

                    //maxpyin(r0,q,num);
                    Integer::maxpyin(r0,q,num);
                    //maxpyin(t0,q,den);
                    Integer::maxpyin(t0,q,den);

                    if (t0 < 0) {
                        num = -r0;
                        den = -t0;
                    }
                    else {
                        num = r0;
                        den = t0;
                    }

                    // if (t0 > m/k)
                    if (den > m/k) {
                        if (!recursive)
                            std::cerr
                            << "*** Error *** No rational reconstruction of "
                            << f
                            << " modulo "
                            << m
                            << " with denominator <= "
                            << (m/k)
                            << std::endl;
                    }
                    if (gcd(gg,num,den) != 1) {
                        if (!recursive)
                            std::cerr
                            << "*** Error *** There exists no rational reconstruction of "
                            << f
                            << " modulo "
                            << m
                            << " with |numerator| < "
                            << k
                            << std::endl
                            << "*** Error *** But "
                            << num
                            << " = "
                            << den
                            << " * "
                            << f
                            << " modulo "
                            << m
                            << std::endl;
                        return false;
                    }
                }
            }
            // (i)
            if (den < 0) {
                Integer::negin(num);
                Integer::negin(den);
            }

            // std::cerr << "RatRecon End " << num << "/" << den << std::endl;
            return true;
        }
        };

        typedef ZRing<float> FloatDomain;
        typedef ZRing<double> DoubleDomain;
        typedef ZRing<Integer> IntegerDomain;
    }

#endif
    /* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
    // vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
