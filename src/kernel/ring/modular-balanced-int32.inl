// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: Clement Pernet <clement.pernet@gmail.com>
//          A. Breust (taken from FFLAS-FFPACK)
// ==========================================================================

#ifndef __GIVARO_modular_balanced_int32_INL
#define __GIVARO_modular_balanced_int32_INL

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

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::mul(Element& r, const Element& a, const Element& b) const
    {
        Element q = static_cast<Element>(double(a) * double(b) * _dinvp);
        r = static_cast<Element>(a * b - q * _p);
        NORMALISE(r);
        return r;
    }

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::div(Element& r, const Element& a, const Element& b) const
    {
        Element tmp;
        return mul (r, a, inv(tmp, b));
    }

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::add(Element& r, const Element& a, const Element& b) const
    {
        r = a + b;
        NORMALISE(r);
        return r;
    }

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::sub(Element& r, const Element& a, const Element& b) const
    {
        r = a - b;
        NORMALISE(r);
        return r;
    }

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::neg(Element& r, const Element& a) const
    {
        return r = -a;
    }

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::inv(Element& r, const Element& a) const
    {
        r = invext(r, (a < 0)? a + _p : a, _p);
        NORMALISE(r);
        return r;
    }

    inline bool ModularBalanced<int32_t>::isUnit(const Element& a) const
    {
        Element u,d;
        extended_euclid(u,d,a,_p);
        return isOne(d) || isMOne(d);
    }

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::mulin(Element& r, const Element& a) const
    {
        return mul(r, r, a);
    }

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::divin(Element& r, const Element& a) const
    {
        return div(r, r, a);
    }

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::addin(Element& r, const Element& a) const
    {
        return add(r, r, a);
    }

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::subin(Element& r, const Element& a) const
    {
        return sub(r, r, a);
    }

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::negin(Element& r) const
    {
        return neg(r, r);
    }

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::invin(Element& r) const
    {
        return inv(r, r);
    }

    //----- Special arithmetic

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::axpy(Element& r, const Element& a, const Element& x, const Element& y) const
    {
        // q could be off by (+/-) 1
        Element q = static_cast<Element>(((((double) a) * ((double) x)) + (double) y) * _dinvp);
        r = static_cast<Element>(a * x + y - q * _p);
        NORMALISE(r);
        return r;
    }

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::axpyin(Element& r, const Element& a, const Element& x) const
    {
        // q could be off by (+/-) 1
        Element q = static_cast<Element>(((((double) a) * ((double) x)) + (double) r) * _dinvp);
        r = static_cast<Element>(a * x + r - q * _p);
        NORMALISE(r);
        return r;
    }

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::axmy(Element& r, const Element& a, const Element& x, const Element& y) const
    {
        // q could be off by (+/-) 1
        Element q = static_cast<Element>(((((double) a) * ((double) x)) - (double) y) * _dinvp);
        r = static_cast<Element>(a * x - y - q * _p);
        NORMALISE(r);
        return r;
    }

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::axmyin(Element& r, const Element& a, const Element& x) const
    {
        // q could be off by (+/-) 1
        Element q = static_cast<Element>(((((double) a) * ((double) x)) - (double) r) * _dinvp);
        r = static_cast<Element>(a * x - r - q * _p);
        NORMALISE(r);
        return r;
    }

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::maxpy(Element& r, const Element& a, const Element& x, const Element& y) const
    {
        return negin(axmy(r, a, x, y));
    }

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>:: maxpyin(Element& r, const Element& a, const Element& x) const
    {
        return negin(axmyin(r, a, x));
    }

    //----- Init

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::init(Element& x) const
    {
        return x = 0;
    }

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::init(Element& x, const float y) const
    {
        x = static_cast<Element>(fmod(y, double(_p)));
        NORMALISE(x);
        return x;
    }

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::init(Element& x, const double y) const
    {
        x = static_cast<Element>(fmod(y, double(_p)));
        NORMALISE(x);
        return x;
    }

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::init(Element& x, const int64_t y) const
    {
        x = static_cast<Element>(y % static_cast<int64_t>(_p));
        NORMALISE(x);
        return x;
    }

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::init(Element& x, const uint64_t y) const
    {
        x = static_cast<Element>(y % static_cast<uint64_t>(_p));
        NORMALISE_HI(x);
        return x;
    }

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::init(Element& x, const Integer& y) const
    {
        x = static_cast<Element>(y % _p);
        NORMALISE_HI(x);
        return x;
    }

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::assign(Element& x, const Element& y) const
    {
        return x = y;
    }

    //----- Reduce

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::reduce(Element& x, const Element& y) const
    {
        x = y % _p;
        NORMALISE(x);
        return x;
    }

    inline ModularBalanced<int32_t>::Element&
    ModularBalanced<int32_t>::reduce(Element& x) const
    {
        x %= _p;
        NORMALISE(x);
        return x;
    }

    //----- IO

    inline std::ostream&
    ModularBalanced<int32_t>::write(std::ostream& os) const
    {
        return os << "ModularBalanced<int32_t> modulo " << _p;
    }

    inline std::ostream&
    ModularBalanced<int32_t>::write(std::ostream& os, const Element& x) const
    {
        return os << x;
    }

    inline std::istream&
    ModularBalanced<int32_t>::read(std::istream& is, Element& x) const
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
