// bb: c'est quoi ce fichier sans guardes, sans licence, sans auteur, sans modeline , etc ?

#pragma once

namespace Givaro
{
    //---- GCD

    template<typename Storage_t>
    inline Storage_t& gcdext(Storage_t& d, Storage_t& u, Storage_t& v, const Storage_t a, const Storage_t b)
    {
        using Compute_t = typename std::make_signed<Storage_t>::type;

        Compute_t u1 = 1, u2 = 0, u3 = static_cast<Compute_t>(a);
        Compute_t v1 = 0, v2 = 1, v3 = static_cast<Compute_t>(b);

        while (v3 != 0)
        {
            Compute_t q = u3 / v3;
            Compute_t t1, t2 , t3;

            t1 = u1 - q * v1;
            t2 = u2 - q * v2;
            t3 = u3 - q * v3;

            u1 = v1;
            u2 = v2;
            u3 = v3;

            v1 = t1;
            v2 = t2;
            v3 = t3;
        }

        u = static_cast<Storage_t>(u1);
        v = static_cast<Storage_t>(u2);
        return d = static_cast<Storage_t>(u3);
    }

    template<typename Storage_t>
    inline Storage_t& invext(Storage_t& x, const Storage_t a, const Storage_t b)
    {
        using Compute_t = typename std::make_signed<Storage_t>::type;

        Compute_t u1 = 1, u3 = static_cast<Compute_t>(a);
        Compute_t v1 = 0, v3 = static_cast<Compute_t>(b);

        while (v3 != 0)
        {
            Compute_t q = u3 / v3;
            Compute_t t;

            t = v1;
            v1 = u1 - q * v1;
            u1 = t;

            t = v3;
            v3 = u3 - q * v3;
            u3 = t;
        }

        return x = static_cast<Storage_t>(u1);
    }

    template<typename Storage_t>
    inline Storage_t invext(const Storage_t a, const Storage_t b)
    {
        Storage_t r;
        return invext(r, a, b);
    }
}

