// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/matrix/givmatdenseops.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givmatdenseops.inl,v 1.3 2011-01-19 18:29:09 briceboyer Exp $
// ==========================================================================

namespace Givaro {

    template<class Domain>
    int MatrixDom<Domain,Dense>::areEqual ( const Element& A, const Element& B ) const
    {
        if (ncol(A) != ncol(B)) return 0;
        if (nrow(A) != nrow(B)) return 0;
        return _supportdomain.areEqual(A,B);
    }

    template<class Domain>
    int MatrixDom<Domain,Dense>::areNEqual ( const Element& A, const Element& B ) const
    {
        return !areEqual(A,B);
    }

    template<class Domain>
    int MatrixDom<Domain,Dense>::isZero ( const Element& A ) const
    {
        return _supportdomain.isZero(A);
    }

    // res <- A + B; aliases of operands are allowed
    template<class Domain>
    void MatrixDom<Domain,Dense>::add( Element& R, const Element& A, const Element& B) const
    {
        GIVARO_ASSERT((ncol(R) == ncol(A)) && (ncol(A) == ncol(B)), "Bad size");
        GIVARO_ASSERT((nrow(R) == nrow(A)) && (nrow(A) == nrow(B)), "Bad size");
        _supportdomain.add(R,A,B);
    }

    // res <- res + A; aliases of operands are allowed
    template<class Domain>
    void MatrixDom<Domain,Dense>::addin( Element& R, const Element& A ) const
    {
        GIVARO_ASSERT((ncol(R) == ncol(A)), "Bad size");
        GIVARO_ASSERT((nrow(R) == nrow(A)), "Bad size");
        _supportdomain.addin(R,A);
    }

    // res <- A + val; aliases of operands are allowed
    template<class Domain>
    void MatrixDom<Domain,Dense>::add( Element& R, const Element& A, const Type_t& val ) const
    {
        GIVARO_ASSERT((ncol(R) == ncol(A)), "Bad size");
        GIVARO_ASSERT((nrow(R) == nrow(A)), "Bad size");
        long i;
        size_t min = (ncol(R) < nrow(R) ? ncol(R) : nrow(R));
        for (i=min; --i>=0; ) _domain.add(R(i,i), A(i,i), val);
    }

    // res <- val + A; aliases of operands are allowed
    template<class Domain>
    void MatrixDom<Domain,Dense>::add( Element& R, const Type_t& val, const Element& A) const
    {
        GIVARO_ASSERT((ncol(R) == ncol(A)), "Bad size");
        GIVARO_ASSERT((nrow(R) == nrow(A)), "Bad size");
        long i;
        size_t min = (ncol(R) < nrow(R) ? ncol(R) : nrow(R));
        for (i=min; --i>=0; ) _domain.add(R(i,i), A(i,i), val);
    }



    // res <- A - B; aliases of operands are allowed
    template<class Domain>
    void MatrixDom<Domain,Dense>::sub( Element& R, const Element& A, const Element& B) const
    {
        GIVARO_ASSERT((ncol(R) == ncol(A)) && (ncol(A) == ncol(B)), "Bad size");
        GIVARO_ASSERT((nrow(R) == nrow(A)) && (nrow(A) == nrow(B)), "Bad size");
        _supportdomain.sub(R,A,B);
    }

    // res <- res - A; aliases of operands are allowed
    template<class Domain>
    void MatrixDom<Domain,Dense>::subin( Element& R, const Element& A ) const
    {
        GIVARO_ASSERT((ncol(R) == ncol(A)), "Bad size");
        GIVARO_ASSERT((nrow(R) == nrow(A)), "Bad size");
        _supportdomain.subin(R,A);
    }

    // res <- res - A; aliases of operands are allowed
    template<class Domain>
    void MatrixDom<Domain,Dense>::sub( Element& R, const Element& A, const Type_t& val ) const
    {
        GIVARO_ASSERT((ncol(R) == ncol(A)), "Bad size");
        GIVARO_ASSERT((nrow(R) == nrow(A)), "Bad size");
        long i;
        size_t min = (ncol(R) < nrow(R) ? ncol(R) : nrow(R));
        for (i=min; --i>=0; ) _domain.sub(R(i,i), A(i,i), val);
    }

    // res <- res - A; aliases of operands are allowed
    template<class Domain>
    void MatrixDom<Domain,Dense>::sub( Element& R, const Type_t& val, const Element& A) const
    {
        GIVARO_ASSERT((ncol(R) == ncol(A)), "Bad size");
        GIVARO_ASSERT((nrow(R) == nrow(A)), "Bad size");
        long i;
        size_t min = (ncol(R) < nrow(R) ? ncol(R) : nrow(R));
        for (i=min; --i>=0; ) _domain.sub(R(i,i), val, A(i,i));
    };

    // res <- - A; aliases of operands are allowed
    template<class Domain>
    void MatrixDom<Domain,Dense>::neg( Element& R, const Element& A) const
    {
        GIVARO_ASSERT((ncol(R) == ncol(A)), "Bad size");
        GIVARO_ASSERT((nrow(R) == nrow(A)), "Bad size");
        _supportdomain.neg(R,A);
    }

    template<class Domain>
    void MatrixDom<Domain,Dense>::negin( Element& R ) const
    {
        _supportdomain.negin(R);
    }


    // R <- A * B, R cannot be an alias for A or B
    template<class Domain>
    void MatrixDom<Domain,Dense>::mul( Element& A, const Element& B, const Element& C) const
    {
        GIVARO_ASSERT((ncol(A) == ncol(B)) && (ncol(B) == ncol(C)), "Bad size");
        GIVARO_ASSERT((nrow(A) == nrow(B)) && (nrow(B) == nrow(C)), "Bad size");
        Indice_t k, i, j;
        int startBi, startCj, startAi;
        Indice_t nrA = nrow(A);
        Indice_t ncA = ncol(A);
        Indice_t ncB = ncol(B);

        // -- the algorithm used the fact that matrices are stored by row
        Type_t tmp;
        for (i=0; i < nrA; ++i)
        {
            startAi = i*ncA;
            for (j=0; j < ncA; ++j, ++startAi)
            {
                startBi = i*ncB; startCj = j;
                _domain.assign(tmp, _domain.zero);
                for (k=0; k<ncB; ++k, ++startBi, startCj+=ncA) // ncA = ncC
                    _domain.axpy(tmp, B[startBi], C[startCj], tmp);
                _domain.assign( A[startAi], tmp );
            }
        }
    }

    template<class Domain>
    void MatrixDom<Domain,Dense>::mulin
    ( Element& R, const Type_t& val) const
    {
        _supportdomain.mulin(R, val);
    }

    template<class Domain>
    void MatrixDom<Domain,Dense>::mul
    ( typename VectorDom<Domain,Dense>::Element& res,
      const Element& M,
      const VectorDom<Domain,Dense>& VD,
      const typename VectorDom<Domain,Dense>::Element& u ) const
    {
        GIVARO_ASSERT( (nrow(M) == VD.dim(res)), "Bad size");
        GIVARO_ASSERT( (ncol(M) == VD.dim(u)), "Bad size");
        Indice_t i, j;
        int startMi;
        Indice_t nrM = nrow(M);
        Indice_t ncM = ncol(M);

        // -- the algorithm used the fact that matrices are stored by row
        for (i=0; i < nrM; ++i)
        {
            _domain.assign(res[i], _domain.zero);
            startMi = i*ncM;
            for (j=0; j < ncM; ++j, ++startMi)
                _domain.axpyin(res[i], M[startMi], u[j]);
        }
    }


    // res <- tr(u) * tr(A)
    template<class Domain>
    void MatrixDom<Domain,Dense>::multrans
    ( typename VectorDom<Domain,Dense>::Element& res,
      const Element& M,
      const VectorDom<Domain,Dense>& VD,
      const typename VectorDom<Domain,Dense>::Element& u ) const
    {
        GIVARO_ASSERT( (nrow(M) == VD.dim(u)), "Bad size");
        GIVARO_ASSERT( (ncol(M) == VD.dim(res)), "Bad size");
        int startMi;
        Indice_t nrM = nrow(M);
        Indice_t ncM = ncol(M);

        // -- the algorithm used the fact that matrices are stored by row
        for (Indice_t j=0; j < ncM; ++j)
            _domain.assign(res[j], _domain.zero);
        for (Indice_t i=0; i < nrM; ++i)
        {
            startMi = i*ncM;
            for (Indice_t j=0; j < ncM; ++j, ++startMi)
                _domain.axpyin(res[j], u[i], M[startMi]);
        }
    }



    template<class Domain>
    void MatrixDom<Domain,Dense>::mul
    ( Element& R, const Type_t& val, const Element& A) const
    {
        GIVARO_ASSERT((ncol(R) == ncol(A)), "Bad size");
        GIVARO_ASSERT((nrow(R) == nrow(A)), "Bad size");
        _supportdomain.mul(R, val, A);
    }

    template<class Domain>
    void MatrixDom<Domain,Dense>::mul
    ( Element& R, const Element& A, const Type_t& val) const
    {
        GIVARO_ASSERT((ncol(R) == ncol(A)), "Bad size");
        GIVARO_ASSERT((nrow(R) == nrow(A)), "Bad size");
        _supportdomain.mul(R, A, val);
    }

    // r <- a*x+y
    template<class Domain>
    void MatrixDom<Domain,Dense>::axpy
    ( Element& R, const Type_t& a, const Element& X, const Element& Y ) const
    {
        GIVARO_ASSERT((ncol(R) == ncol(X)) && (ncol(X) == ncol(Y)), "Bad size");
        GIVARO_ASSERT((nrow(R) == nrow(X)) && (nrow(X) == nrow(Y)), "Bad size");
        _supportdomain.axpy(R, a, X, Y);
    }
    // r <- r+a*x
    template<class Domain>
    void MatrixDom<Domain,Dense>::axpyin
    ( Element& R, const Type_t& a, const Element& X ) const
    {
        GIVARO_ASSERT((ncol(R) == ncol(X)), "Bad size");
        GIVARO_ASSERT((nrow(R) == nrow(X)), "Bad size");
        _supportdomain.axpyin(R, a, X);
    }
    // r <- a*x-y
    template<class Domain>
    void MatrixDom<Domain,Dense>::axmy ( Element& R,
                                         const Type_t& a, const Element& X, const Element& Y ) const
    {
        GIVARO_ASSERT((ncol(R) == ncol(X)) && (ncol(X) == ncol(Y)), "Bad size");
        GIVARO_ASSERT((nrow(R) == nrow(X)) && (nrow(X) == nrow(Y)), "Bad size");
        _supportdomain.axmy(R, a, X, Y);
    }
    // r <- a*x-r
    template<class Domain>
    void MatrixDom<Domain,Dense>::axmyin ( Element& R,
                                           const Type_t& a, const Type_t& X ) const
    {
        GIVARO_ASSERT((ncol(R) == ncol(X)), "Bad size");
        GIVARO_ASSERT((nrow(R) == nrow(X)), "Bad size");
        _supportdomain.axmy(R, a, X);
    }


    // -- axpy operations:
    // A*X + Y
    template<class Domain>
    void MatrixDom<Domain,Dense>::axpy
    ( Element& R, const Element& A, const Element& X, const Element& Y ) const
    {
        int k, i, j;
        int startAi, startXj, startRi;
        int nrR = nrow(R);
        int ncR = ncol(R);
        int ncA = ncol(A);

        // -- the algorithm used the fact that matrices are stored by row
        Type_t tmp;
        for (i=0; i < nrR; ++i)
        {
            startRi = i*ncR;
            for (j=0; j < ncR; ++j, ++startRi)
            {
                startAi = i*ncA; startXj = j;
                _domain.assign(tmp, _domain.zero);
                for (k=0; k<ncA; ++k, ++startAi, startXj+=ncR) // ncR = ncX
                    _domain.axpy(tmp, A[startAi], X[startXj], tmp);
                _domain.add( R[startRi], tmp, Y[startRi] );
            }
        }
    }

    // a*A*X - b*Y
    template<class Domain>
    void MatrixDom<Domain,Dense>::axpy
    ( Element& R, const Type_t& a, const Element& A,
      const Element& X, const Type_t& b, const Element& Y ) const
    {
        Indice_t k, i, j;
        int startAi, startXj, startRi;
        Indice_t nrR = nrow(R);
        Indice_t ncR = ncol(R);
        Indice_t ncA = ncol(A);

        // -- the algorithm used the fact that matrices are stored by row
        Type_t tmp;
        for (i=0; i < nrR; ++i)
        {
            startRi = i*ncR;
            for (j=0; j < ncR; ++j, ++startRi)
            {
                startAi = i*ncA; startXj = j;
                _domain.assign(tmp, _domain.zero);
                for (k=0; k<ncA; ++k, ++startAi, startXj+=ncR) // ncR = ncX
                    _domain.axpy(tmp, A[startAi], X[startXj], tmp);
                _domain.mul( tmp, a, tmp );
                _domain.axpy( R[startRi], b, Y[startRi], tmp);
            }
        }
    }

    // Map: used for constant operation !, not that OP
    template<class Domain>
    template<class OP>
    void MatrixDom<Domain,Dense>::map
    ( Element& R, OP& op ) const
    {
        const Indice_t max = _supportdomain.dim(R);
        for (Indice_t i(0); i<max; ++i) op(R[i]);
    }

    template<class Domain>
    template<class OP>
    void MatrixDom<Domain,Dense>::map
    ( Element& R, OP& op, const Element& A ) const
    {
        const Indice_t max = _supportdomain.dim(R);
        for (Indice_t i(0); i<max; ++i) op(R[i], A[i]);
    }

    template<class Domain>
    template<class OP>
    void MatrixDom<Domain,Dense>::map
    ( Element& R, OP& op, const Element& A, const Element& B ) const
    {
        const Indice_t max = _supportdomain.dim(R);
        for (Indice_t i(0); i<max; ++i) op(R[i], A[i], B[i]);
    }


    template <class Domain>
    inline std::ostream& MatrixDom<Domain,Dense>::write( std::ostream& sout ) const
    {
        return _domain.write(sout << '(') << ",Dense)";
    }


    // -- read the domain
    template<class Domain>
    std::istream& MatrixDom<Domain,Dense>::read( std::istream& sin )
    {
        char ch;
        sin >> std::ws >> ch;
        if (ch != '(')
            GivError::throw_error(
                                  GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no '('"));

        // -- read subdomain:
        _domain.read(sin);

        // -- read ,
        sin >> std::ws >> ch;
        if (ch != ',')
            GivError::throw_error(
                                  GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no ','"));

        // -- read dense:
        char name[6];
        sin >> std::ws; sin.read(name, 5); name[5] = (char)0;
        if (strcmp(name, "Dense") !=0)
            GivError::throw_error(
                                  GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no 'Dense'"));

        // -- read )
        sin >> std::ws >> ch;
        if (ch != ')')
            GivError::throw_error(
                                  GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no ')'"));
        return sin;
    }


    template <class Domain>
    std::ostream& MatrixDom<Domain,Dense>::write( std::ostream& sout, const Element& A) const
    {
        Indice_t i,j;
        sout << '[' << nrow(A) << ',' << ncol(A) << ",[";
        for (i=0; i < nrow(A); i++)
        {
            sout << "[";
            for (j=0; j < ncol(A); j++)
            {
                _domain.write(sout,A(i,j));
                if (j < ncol(A) - 1) sout << ",";
            }
            if (i < nrow(A) - 1) sout << "]," << std::endl;
            else sout << "]";
        }
        return sout << "]]";
    }

    template <class Domain>
    std::istream& MatrixDom<Domain,Dense>::read (std::istream& sin, Element& R) const
    {
        char ch;
        Indice_t i,j;
        long nr,nc;

        //  -- read [
        sin >> std::ws >> ch;
        if (ch != '[') {
            std::cerr << "Read: " << ch << std::endl;
            GivError::throw_error(
                                  GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no '['"));
        }

        // -- read dimensions: nr , nc
        // -- read nrow
        sin >> std::ws >> nr;
        if (nr < 0)
            GivError::throw_error(
                                  GivBadFormat("MatrixDom<Domain,Dense>::read: bad row dimension"));

        //  -- read ,
        sin >> std::ws >> ch;
        if (ch != ',')
            GivError::throw_error(
                                  GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no ','"));

        //  -- read ncol
        sin >> std::ws >> nc;
        if (nc < 0)
            GivError::throw_error(
                                  GivBadFormat("MatrixDom<Domain,Dense>::read: bad col dimension"));

        //  -- read ,
        sin >> std::ws >> ch;
        if (ch != ',')
            GivError::throw_error(
                                  GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no ','"));

        //  -- read [
        sin >> std::ws >> ch;
        if (ch != '[')
            GivError::throw_error(
                                  GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no '[' before entries"));

        // -- Read the matrix:
        R.resize( nr, nc );
        for (i=0; i<nr; ++i) {
            if (i != 0) {
                sin >> std::ws >> ch;
                if (ch != ',')
                    GivError::throw_error(
                                          GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no ','"));
            }

//             cout << "Row: " << i << std::endl;
            //  -- read [
            sin >> std::ws >> ch;
            if (ch != '[')
                GivError::throw_error(
                                      GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no '['"));

            // read nc-1 times: val ,
            for (j=0; j< nc-1; ++j) {
//                 cout << "Read: i=" << i << ", j=" << j << std::endl;
                _domain.read(sin, R(i,j));
                //  -- read [
                sin >> std::ws >> ch;
                if (ch != ',')
                    GivError::throw_error(
                                          GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no ','"));
            }

            // read 1 : val ]
//             cout << "Read: i=" << i << ", j=" << j << std::endl;
            _domain.read(sin, R(i,j));
            sin >> std::ws >> ch;
            if (ch != ']')
                GivError::throw_error(
                                      GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no ']'"));
        }

        //  -- read ] ]
        sin >> std::ws >> ch;
        if (ch != ']')
            GivError::throw_error(
                                  GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no ']'"));
        sin >> std::ws >> ch;
        if (ch != ']')
            GivError::throw_error(
                                  GivBadFormat("MatrixDom<Domain,Dense>::read: syntax error no ']'"));
        return sin;
    }
} // Givaro

//#include "givaro/givmatdenseops.f.spe"
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
