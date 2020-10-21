// ==========================================================================
// Copyright(c)'1994-2020 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: C. Bouvier
// ==========================================================================

#ifndef __GIVARO_typestring_H
#define __GIVARO_typestring_H

#include <string>
#include <type_traits>

namespace Givaro {

    template <typename T>
    class HasTypeString
    {
    private:
        typedef char YesType[1];
        typedef char NoType[2];

        template <typename C> static YesType& test (decltype(&C::type_string));
        template <typename C> static NoType& test (...);

    public:
        static constexpr bool value = sizeof(test<T>(0)) == sizeof(YesType);
    };


    template <class T>
    struct TypeString
    {
        template <bool B, class R = void>
        using enable_if_t = typename std::enable_if<B, R>::type;

        /* base type, integral: char, int8_t, uint32_t, ... */
        template <class C = T, enable_if_t<std::is_integral<C>::value>* = nullptr>
        static std::string get () {
            std::string s = std::string(std::is_signed<C>::value ? "int" : "uint");
            return s + std::to_string(8*sizeof(C)) + "_t";
        }

        /* base type: bool */
        template <class C = T, enable_if_t<std::is_same<C, bool>::value>* = nullptr>
        static std::string get () {
            return "bool";
        }

        /* base type: float */
        template <class C = T, enable_if_t<std::is_same<C,float>::value>* = nullptr>
        static std::string get () {
            return "float";
        }

        /* base type: double */
        template <class C = T, enable_if_t<std::is_same<C,double>::value>* = nullptr>
        static std::string get () {
            return "double";
        }

        /* class with type_string static method */
        template <class U = T, enable_if_t<HasTypeString<U>::value>* = nullptr>
        static std::string get () {
            return T::type_string();
        }
    };

}

#endif /* __GIVARO_typestring_H */
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
