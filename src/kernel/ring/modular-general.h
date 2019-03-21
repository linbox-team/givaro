// ========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: A. Breust
// Time-stamp: <28 May 15 09:40:00 Alexis.Breust@imag.fr>
// ========================================================================
// Description:
// Forward declarations for Givaro::Modular and associated functions

#pragma once

namespace Givaro
{
    /*! Forward declaration for Givaro::Modular.
     *  Elements will be stored in the storage type.
     *  While arithmetics will occur on the unsigned version of COMP type.
     *  Example: Modular<int32_t, uint64_t>
     */
    //template<typename Storage_t, typename COMP = Storage_t> class Modular;
    template<typename Storage_t, typename Compute_t = Storage_t, typename Enable = void> class Modular;

    //! Generalized extended GCD used by specialized Modular.
    template<typename Storage_t>
    inline Storage_t& gcdext(Storage_t& d,  Storage_t& u, Storage_t& v, const Storage_t a, const Storage_t b);

    //! Extended Euclidean algorithm computing only the  Bezout coefficient for a
    //  2 different versions, for floating or not floating point elements
    //  both require a different forward declaration with std::enable_if
    template<typename Storage_t>
    inline typename std::enable_if<std::is_floating_point<Storage_t>::value, Storage_t&>::type
    extended_euclid (Storage_t& x, Storage_t& d, const Storage_t a, const Storage_t b);

    template<typename Storage_t>
    inline typename std::enable_if<!std::is_floating_point<Storage_t>::value, Storage_t&>::type
    extended_euclid (Storage_t& x, Storage_t& d, const Storage_t a, const Storage_t b);

    //! Generalized inversion used by specialized Modular.
    //  uses extended_euclid, adding a zero division protection in DEBUG mode
    template<typename Storage_t>
      inline Storage_t& invext (Storage_t& x, Storage_t& d, const Storage_t a, const Storage_t b);

    template<typename Storage_t>
      inline Storage_t& invext(Storage_t& x, const Storage_t a, const Storage_t b);

    template<typename Storage_t>
      inline Storage_t invext(const Storage_t a, const Storage_t b);
}

#include "givaro/modular-general.inl"
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
