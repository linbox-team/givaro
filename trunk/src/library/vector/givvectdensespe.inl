// ==========================================================================
// $Source
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id
// ==========================================================================


template<>
inline void VectorDom<ZpzDom<Std16>,Dense>::dot
  ( Type_t& res, const Rep& op1, const Rep& op2) const
{
  size_t sz = dim(op1);
  const ZpzDom<Std16>& domain = subdomain();
  domain.dotprod( res, sz, op1.baseptr(), op2.baseptr() );
}

