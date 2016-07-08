#define SUITABLE_INTEGRAL_TYPES(ELEMENT, COMPUTE_T, s) \
        (2 * sizeof(ELEMENT) >= sizeof(COMPUTE_T)) && (sizeof(COMPUTE_T) == s) \
        && (std::is_integral<ELEMENT>::value) && (std::is_integral<COMPUTE_T>::value)

#define SUITABLE_INTEGRAL_TYPES_AND_HALF_SIZE(ELEMENT, COMPUTE_T, s) \
        (2 * sizeof(ELEMENT) == sizeof(COMPUTE_T)) && (sizeof(COMPUTE_T) == s) \
        && (std::is_integral<ELEMENT>::value) && (std::is_integral<COMPUTE_T>::value)

#define SUITABLE_INTEGRAL_TYPES_AND_FULL_SIZE(ELEMENT, COMPUTE_T, s) \
        (sizeof(ELEMENT) == sizeof(COMPUTE_T)) && (sizeof(COMPUTE_T) == s) \
        && (std::is_integral<ELEMENT>::value) && (std::is_integral<COMPUTE_T>::value)

template<typename ELEMENT = Element, typename COMPUTE_T = Compute_t, int s = sizeof(COMPUTE_T)>
inline
typename std::enable_if<SUITABLE_INTEGRAL_TYPES (ELEMENT, COMPUTE_T, s)>::type
precomp_p
(Compute_t& invp) const
{
        assert( _bitsizep <= (4*s-2));
        invp = (static_cast<Compute_t>(1) << (4*s + _bitsizep - 1)) / static_cast<Compute_t>(_p);
}

template<typename ELEMENT = Element, typename COMPUTE_T = Compute_t, int s = sizeof(COMPUTE_T)>
inline
typename std::enable_if<SUITABLE_INTEGRAL_TYPES (ELEMENT, COMPUTE_T, s), Element&>::type
mul_precomp_p
(Element& r, const Element& a, const Element& b, const Compute_t& invp) const
{
        Compute_t prod = static_cast<Compute_t>(static_cast<Residu_t>(a))*static_cast<Compute_t>(static_cast<Residu_t>(b));
        Compute_t prodhi = (prod >> (_bitsizep - 2)); // Could fit into an Element but no use
        Residu_t c = static_cast<Residu_t>((prodhi * invp) >> (4*s+1));
        Residu_t rr = static_cast<Residu_t>(prod) - c * _p;
        rr -= (rr >= _p)?_p:0;
        return r = static_cast<Element>(rr);
}

template<typename ELEMENT = Element, typename COMPUTE_T = Compute_t, int s = sizeof(COMPUTE_T)>
inline
typename std::enable_if<SUITABLE_INTEGRAL_TYPES (ELEMENT, COMPUTE_T, s)>::type
precomp_b
(Compute_t& invb, const Element& b) const{
        assert( _bitsizep <= (4*s-1));
        invb = (static_cast<Compute_t>(1) << (4*s)) * static_cast<Compute_t>(static_cast<Residu_t>(b)) / static_cast<Compute_t>(_p);
}

template<typename ELEMENT = Element, typename COMPUTE_T = Compute_t, int s = sizeof(COMPUTE_T)>
inline
typename std::enable_if<SUITABLE_INTEGRAL_TYPES_AND_HALF_SIZE (ELEMENT, COMPUTE_T, s)>::type
precomp_b
(Compute_t& invb, const Element& b, const Compute_t& invp) const{
        assert( _bitsizep <= (4*s-2));

        invb = (static_cast<Compute_t>(static_cast<Residu_t>(b)) * invp) >> (_bitsizep - 1);
        Residu_t r = - static_cast<Residu_t>(invb) * _p;

        // Quotient can only be off by two since (2^n * b) / 2^(r-1) is an exact division
        bool flag = (r >= (_p));
        invb += (flag)?1:0;
        r -= (flag)?_p:0;
        invb += (r >= (_p))?1:0;
}

template<typename ELEMENT = Element, typename COMPUTE_T = Compute_t, int s = sizeof(COMPUTE_T)>
inline
typename std::enable_if<SUITABLE_INTEGRAL_TYPES_AND_FULL_SIZE (ELEMENT, COMPUTE_T, s)>::type
precomp_b
(Compute_t& invb, const Element& b, const Element& invp) const{
        assert( _bitsizep <= (4*s-2));
        // Prod = 2^(4*s) * b
        invb = (static_cast<Compute_t>(b) * static_cast<Compute_t>(static_cast<Residu_t>(invp))) >> (_bitsizep - 1);

        // The problem is that Residu_t should be unsigned half size
        // Warning : left shift count >= width of type when Element is half size of Compute_t
        // and then << does nothing instead of mapping to zero
        Residu_t r = (static_cast<Residu_t>(b) << (4*s)) - static_cast<Residu_t>(invb) * _p;

        // Quotient can only be off by two since (2^n * b) / 2^(r-1) is an exact division
        bool flag = (r >= (_p));
        invb += (flag)?1:0;
        r -= (flag)?_p:0;
        invb += (r >= (_p))?1:0;
}

template<typename ELEMENT = Element, typename COMPUTE_T = Compute_t, int s = sizeof(COMPUTE_T)>
inline
typename std::enable_if<SUITABLE_INTEGRAL_TYPES (ELEMENT, COMPUTE_T, s), Element&>::type
mul_precomp_b
(Element& r, const Element& a, const Element& b, const Compute_t& invb) const {
        Residu_t c = static_cast<Residu_t> ((static_cast<Compute_t>(a) * invb) >> (4*s));
        Residu_t rr = static_cast<Residu_t>(a) * static_cast<Residu_t>(b) - c * _p;
        rr -= (rr >= _p)?_p:0;
        return r = static_cast<Residu_t>(rr);
}

#undef SUITABLE_INTEGRAL_TYPES
