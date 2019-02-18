// ==========================================================================
// Copyright(c)'1994-2015 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: A. Breust <alexis.breust@gmail.com>
// ==========================================================================

#pragma once

#include <cstdint>

/*! User-defined literals.
 *  Allows cross-platform literals with compile-time cast.
 */

constexpr uint64_t operator"" _ui64(unsigned long long x)
{
    return static_cast<uint64_t>(x);
}

constexpr uint32_t operator"" _ui32(unsigned long long x)
{
    return static_cast<uint32_t>(x);
}

constexpr uint16_t operator"" _ui16(unsigned long long x)
{
    return static_cast<uint16_t>(x);
}

constexpr uint8_t operator"" _ui8(unsigned long long x)
{
    return static_cast<uint8_t>(x);
}

constexpr int64_t operator"" _i64(unsigned long long x)
{
    return static_cast<int64_t>(x);
}

constexpr int32_t operator"" _i32(unsigned long long x)
{
    return static_cast<int32_t>(x);
}

constexpr int16_t operator"" _i16(unsigned long long x)
{
    return static_cast<int16_t>(x);
}

constexpr int8_t operator"" _i8(unsigned long long x)
{
    return static_cast<int8_t>(x);
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
