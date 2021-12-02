// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: Clement Pernet <clement.pernet@gmail.com>
//          A. Breust (taken from FFLAS-FFPACK)
// ==========================================================================

#ifndef __GIVARO_modular_balanced_double_INL
#define __GIVARO_modular_balanced_double_INL

#include <cmath> // fmod

#define NORMALISE(x)				\
{						\
    if (x < _mhalfp) x += _p;		\
    else if (x > _halfp) x -= _p;		\
}

#define NORMALISE_HI(x)				\
{						\
    if (x > _halfp) x -= _p;		\
}

namespace Givaro
{

    //----- Classic arithmetic

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::mul(Element& r, const Element& a, const Element& b) const
    {
        return reduce(r = a * b);
    }

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::div(Element& r, const Element& a, const Element& b) const
    {
        Element tmp;
        return mul (r, a, inv(tmp, b));
    }

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::add(Element& r, const Element& a, const Element& b) const
    {
        r = a + b;
        NORMALISE(r);
        return r;
    }

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::sub(Element& r, const Element& a, const Element& b) const
    {
        r = a - b;
        NORMALISE(r);
        return r;
    }

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::neg(Element& r, const Element& a) const
    {
        return r = -a;
    }

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::inv(Element& r, const Element& a) const
    {
        r = invext(a, _p);
        NORMALISE(r);
        return r;
    }

    inline bool ModularBalanced<double>::isUnit(const Element& a) const
    {
        Element u,d;
        extended_euclid(u,d,a,_p);
        return isOne(d) || isMOne(d);
    }


    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::mulin(Element& r, const Element& a) const
    {
        return mul(r, r, a);
    }

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::divin(Element& r, const Element& a) const
    {
        return div(r, r, a);
    }

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::addin(Element& r, const Element& a) const
    {
        return add(r, r, a);
    }

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::subin(Element& r, const Element& a) const
    {
        return sub(r, r, a);
    }

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::negin(Element& r) const
    {
        return neg(r, r);
    }

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::invin(Element& r) const
    {
        return inv(r, r);
    }

    //----- Special arithmetic

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::axpy(Element& r, const Element& a, const Element& x, const Element& y) const
    {
        return reduce(r = a * x + y);
    }

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::axpyin(Element& r, const Element& a, const Element& x) const
    {
        return reduce(r += a * x);
    }

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::axmy(Element& r, const Element& a, const Element& x, const Element& y) const
    {
        return reduce(r = a * x - y);
    }

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::axmyin(Element& r, const Element& a, const Element& x) const
    {
        return reduce(r = a * x - r);
    }

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::maxpy(Element& r, const Element& a, const Element& x, const Element& y) const
    {
        return reduce(r = y - a * x);
    }

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>:: maxpyin(Element& r, const Element& a, const Element& x) const
    {
        return reduce(r -= a * x);
    }

    //----- Init

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::init(Element& x) const
    {
        return x = 0;
    }

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::init(Element& x, const float y) const
    {
        x = static_cast<Element>(fmod(y, double(_p)));
        NORMALISE(x);
        return x;
    }

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::init(Element& x, const double y) const
    {
        x = static_cast<Element>(fmod(y, _p));
        NORMALISE(x);
        return x;
    }

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::init(Element& x, const int64_t y) const
    {
        x = static_cast<Element>(y % static_cast<int64_t>(_up));
        NORMALISE(x);
        return x;
    }

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::init(Element& x, const uint64_t y) const
    {
        x = static_cast<Element>(y % static_cast<uint64_t>(_up));
        NORMALISE_HI(x);
        return x;
    }

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::init(Element& x, const Integer& y) const
    {
        x = static_cast<Element>(y % _p);
        NORMALISE_HI(x);
        return x;
    }

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::assign(Element& x, const Element& y) const
    {
        return x = y;
    }

    //----- Reduce

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::reduce(Element& x, const Element& y) const
    {
        x = fmod(y, _p);
        NORMALISE(x);
        return x;
    }

    inline ModularBalanced<double>::Element&
    ModularBalanced<double>::reduce(Element& x) const
    {
        x = fmod(x, _p);
        NORMALISE(x);
        return x;
    }

    //----- IO

    inline std::ostream&
    ModularBalanced<double>::write(std::ostream& os) const
    {
        return os << "ModularBalanced<double> modulo " << static_cast<uint64_t>(_p); // _p in ModularBalanced<double> is always < 2**64
    }

    inline std::ostream&
    ModularBalanced<double>::write(std::ostream& os, const Element& x) const
    {
        return os << x;
    }

    inline std::istream&
    ModularBalanced<double>::read(std::istream& is, Element& x) const
    {
        Element tmp;
        is >> tmp;
        init(x, tmp);
        return is;
    }

}

#undef NORMALISE
#undef NORMALISE_HI

#endif
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
