// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: Clement Pernet <clement.pernet@gmail.com>
//          A. Breust (taken from FFLAS-FFPACK)
// ==========================================================================

#ifndef __GIVARO_modular_balanced_float_INL
#define __GIVARO_modular_balanced_float_INL

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

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::mul(Element& r, const Element& a, const Element& b) const
    {
        return reduce(r = a * b);
    }

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::div(Element& r, const Element& a, const Element& b) const
    {
        Element tmp;
        return mul (r, a, inv(tmp, b));
    }

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::add(Element& r, const Element& a, const Element& b) const
    {
        r = a + b;
        NORMALISE(r);
        return r;
    }

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::sub(Element& r, const Element& a, const Element& b) const
    {
        r = a - b;
        NORMALISE(r);
        return r;
    }

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::neg(Element& r, const Element& a) const
    {
        return r = -a;
    }

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::inv(Element& r, const Element& a) const
    {
        r = invext(a, _p);
        NORMALISE(r);
        return r;
    }

    inline bool ModularBalanced<float>::isUnit(const Element& a) const
    {
        Element u,d;
        extended_euclid(u,d,a,_p);
        return isOne(d) || isMOne(d);
    }

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::mulin(Element& r, const Element& a) const
    {
        return mul(r, r, a);
    }

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::divin(Element& r, const Element& a) const
    {
        return div(r, r, a);
    }

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::addin(Element& r, const Element& a) const
    {
        return add(r, r, a);
    }

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::subin(Element& r, const Element& a) const
    {
        return sub(r, r, a);
    }

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::negin(Element& r) const
    {
        return neg(r, r);
    }

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::invin(Element& r) const
    {
        return inv(r, r);
    }

    //----- Special arithmetic

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::axpy(Element& r, const Element& a, const Element& x, const Element& y) const
    {
        return reduce(r = a * x + y);
    }

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::axpyin(Element& r, const Element& a, const Element& x) const
    {
        return reduce(r += a * x);
    }

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::axmy(Element& r, const Element& a, const Element& x, const Element& y) const
    {
        return reduce(r = a * x - y);
    }

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::axmyin(Element& r, const Element& a, const Element& x) const
    {
        return reduce(r = a * x - r);
    }

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::maxpy(Element& r, const Element& a, const Element& x, const Element& y) const
    {
        return reduce(r = y - a * x);
    }

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>:: maxpyin(Element& r, const Element& a, const Element& x) const
    {
        return reduce(r -= a * x);
    }

    //----- Init

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::init(Element& x) const
    {
        return x = 0;
    }

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::init(Element& x, const float y) const
    {
        x = static_cast<Element>(fmodf(y, _p));
        NORMALISE(x);
        return x;
    }

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::init(Element& x, const double y) const
    {
        x = static_cast<Element>(fmod(y, double(_p)));
        NORMALISE(x);
        return x;
    }

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::init(Element& x, const int32_t y) const
    {
        x = static_cast<Element>(y % static_cast<int32_t>(_up));
        NORMALISE(x);
        return x;
    }

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::init(Element& x, const uint32_t y) const
    {
        x = static_cast<Element>(y %  static_cast<uint32_t>(_up));
        NORMALISE_HI(x);
        return x;
    }

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::init(Element& x, const int64_t y) const
    {
        x = static_cast<Element>(y % static_cast<int64_t>(_up));
        NORMALISE(x);
        return x;
    }

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::init(Element& x, const uint64_t y) const
    {
        x = static_cast<Element>(y % static_cast<uint64_t>(_up));
        NORMALISE_HI(x);
        return x;
    }

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::init(Element& x, const Integer& y) const
    {
        x = static_cast<Element>(y % _p);
        NORMALISE_HI(x);
        return x;
    }

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::assign(Element& x, const Element& y) const
    {
        return x = y;
    }

    //----- Reduce

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::reduce(Element& x, const Element& y) const
    {
        x = fmodf(y, _p);
        NORMALISE(x);
        return x;
    }

    inline ModularBalanced<float>::Element&
    ModularBalanced<float>::reduce(Element& x) const
    {
        x = fmodf(x, _p);
        NORMALISE(x);
        return x;
    }

    //----- IO

    inline std::ostream&
    ModularBalanced<float>::write(std::ostream& os) const
    {
        return os << "ModularBalanced<float> modulo " << static_cast<uint32_t>(_p); // _p in ModularBalanced<double> is always < 2**32 
    }

    inline std::ostream&
    ModularBalanced<float>::write(std::ostream& os, const Element& x) const
    {
        return os << x;
    }

    inline std::istream&
    ModularBalanced<float>::read(std::istream& is, Element& x) const
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
