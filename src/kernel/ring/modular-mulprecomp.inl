#define SUITABLE_INTEGRAL_TYPES(ELEMENT, COMPUTE_T, s) \
	(2 * sizeof(ELEMENT) >= sizeof(COMPUTE_T)) && (sizeof(COMPUTE_T) == s) \
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
	Compute_t prod = static_cast<Compute_t>(a)*static_cast<Compute_t>(b);
	Compute_t prodhi = (prod >> (_bitsizep - 2)); // Could fit into an Element but no use
	Element c = static_cast<Element>((prodhi * invp) >> (4*s+1));
	r = static_cast<Element>(prod) - c * static_cast<Element>(_p);
	r -= (r >= static_cast<Element>(_p))?_p:0;
	return r;
}


template<typename ELEMENT = Element, typename COMPUTE_T = Compute_t, int s = sizeof(COMPUTE_T)>
inline
typename std::enable_if<SUITABLE_INTEGRAL_TYPES (ELEMENT, COMPUTE_T, s)>::type
precomp_b
(Compute_t& invb, const Element& b) const{
	assert( _bitsizep <= (4*s-1));
	invb = (static_cast<Compute_t>(1) << (4*s)) * static_cast<Compute_t>(b) / static_cast<Compute_t>(_p);
}

template<typename ELEMENT = Element, typename COMPUTE_T = Compute_t, int s = sizeof(COMPUTE_T)>
inline
typename std::enable_if<SUITABLE_INTEGRAL_TYPES (ELEMENT, COMPUTE_T, s), Element&>::type
mul_precomp_b
(Element& r, const Element& a, const Element& b, const Compute_t& invb) const {
	Element c = static_cast<Element> ((static_cast<Compute_t>(a) * invb) >> (4*s));
	r = a * b - c * static_cast<Element> (_p);
	r -= (r >= static_cast<Element>(_p))?_p:0;
	return r;
}

#undef SUITABLE_INTEGRAL_TYPES
