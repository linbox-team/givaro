// bb: c'est quoi ce fichier sans guardes, sans licence, sans auteur, sans modeline etc ?

#pragma once

// This file is necessary as default template parameters cannot
// be redeclared. - A.B. 2015-18-03

namespace Givaro
{
    // Elements will be stored in the storage type.
    // While arithmetics will occur on the unsigned version of COMP type.
    // Example: Modular<int32_t, uint64_t>
    template<typename Storage_t, typename COMP = Storage_t> class Modular;

    //---- GCD

    template<typename Storage_t>
    inline Storage_t& gcdext(Storage_t& d,  Storage_t& u, Storage_t& v, const Storage_t a, const Storage_t b);

    template<typename Storage_t>
    inline Storage_t& invext(Storage_t& u, const Storage_t a, const Storage_t b);

    template<typename Storage_t>
    inline Storage_t invext(const Storage_t a, const Storage_t b);
}

#include "givaro/modular-general.inl"

