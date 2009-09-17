// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givpoly1cstor.inl,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givpoly1cstor.inl,v 1.9 2009-09-17 11:30:08 jgdumas Exp $
// ==========================================================================

template<class Domain>
inline Poly1Dom<Domain,Dense>::Poly1Dom(const Domain& d, const Indeter& X )
  : _domain(d), _x(X) ,zero(0), one(1)
{ 
	Type_t _one;
	_domain.init( _one, 1.0);
	_domain.assign( one[0], _one);
}

template<class Domain>
inline Poly1Dom<Domain,Dense>::Poly1Dom(const Self_t& P)
  : _domain(P._domain), _x(P._x) ,zero(P.zero), one(P.one)
{}

template<class Domain>
inline Poly1Dom<Domain,Dense>::~Poly1Dom()
{ 
}


template<class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::init( Rep& P) const
{ P.reallocate(0); return P; }

/*
template<class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::init( Rep& P, const Rep& Q ) const
{ 
  Degree degQ; 
  degree(degQ,Q);
  if (degQ <0) { 
    P.reallocate(0);
    return P;
  }
  P.reallocate(++degQ);
  for (int i=0; degQ>i; ++i)
    _domain.init(P[i], Q[i]);
  return P;
}
*/

template<class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::assign( Rep& P, const Rep& Q ) const
{
  Degree degQ; 
  degree(degQ,Q);
  if (degQ <0) { 
    P.reallocate(0);
    return P;
  }
  P.reallocate(++degQ);
  for (int i=0; degQ>i; ++i)
    _domain.assign(P[i], Q[i]);
  return P;
}


template<class Domain> template<class XXX>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::init( Rep& P, const XXX& val ) const
{ 
    P.reallocate(1); 
    _domain.init(P[0], val);
    return P;
}


template<class Domain> 
inline typename Poly1Dom<Domain,Dense>::Type_t& Poly1Dom<Domain,Dense>::convert(typename Poly1Dom<Domain,Dense>::Type_t & val, const typename Poly1Dom<Domain,Dense>::Rep& P ) const
{ 
    if (P.size())
        return _domain.assign(val, P[0]);
    else
        return _domain.init(val, 0UL);
}

template<class Domain> template<class XXX>
inline XXX& Poly1Dom<Domain,Dense>::convert( XXX& val, const typename Poly1Dom<Domain,Dense>::Rep& P ) const
{ 
    if (P.size())
        return _domain.convert(val, P[0]);
    else
        return _domain.convert(val, 0UL);
}

template<class Domain> template<class UU, template<class XX> class Vect>
inline Vect<UU>& Poly1Dom<Domain,Dense>::convert( Vect<UU>& val, const typename Poly1Dom<Domain,Dense>::Rep& P ) const
{ 
    val.resize( P.size() );
    typename Vect<UU>::iterator vit = val.begin();
    typename Rep::const_iterator        pit = P.begin();
    for ( ; pit != P.end(); ++pit, ++vit)
        _domain.convert(*vit, *pit);
    return val;
}



template<class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::init( Rep& P, const Degree deg ) const
{ 
  P.reallocate(value(deg+1)); 
	Type_t _one,_zero;
	_domain.init( _one, 1.0);
	_domain.init( _zero, 0.0);
	
	size_t sz = P.size();
	for (unsigned int i=0; i<sz-1; ++i)
		_domain.assign(P[i], _zero);
	_domain.assign(P[sz-1], _one);
	return P;
}

template<class Domain> template<class XXX>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::init
 ( Rep& P, const Degree d, const XXX& val ) const
{
    Type_t _zero;
    _domain.init( _zero, 0.0);
    long deg = value(d);
    P.reallocate(deg+1);
    for (int i=0; i<deg; ++i)
        _domain.assign(P[i], _zero);
    _domain.init(P[deg], val);
    if (_domain.isZero(P[deg])) { 
        P.reallocate(0);
    } 
    return P;
}


template<class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::assign
 ( Rep& P, const Degree d, const Type_t& lcoeff ) const
{
    Type_t _zero;
    _domain.init( _zero, 0.0);
    long deg = value(d);
    if (_domain.isZero(lcoeff)) { 
        P.reallocate(0);
    } else {
        P.reallocate(deg+1);
        for (int i=0; i<deg; ++i)
            _domain.assign(P[i], _zero);
        _domain.assign(P[deg], lcoeff);
    }
    return P;
}

