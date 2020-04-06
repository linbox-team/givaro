// ==========================================================================
// $Source
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id
// ==========================================================================
namespace Givaro {

    // -- map of a unary operator, with operator()( Type_t& res )
    // res and u could be aliases if OP permits it
    template<class Domain>
    template<class UNOP>
    inline void VectorDom<Domain,Dense>::
    map ( Element& res, UNOP& op ) const
    {
        size_t sz = dim(res);
        for (size_t i=0; i<sz; ++i) op(res[i]);
    }


    template<class Domain>
    template<class BINOP>
    inline void VectorDom<Domain,Dense>::
    map( Element& res, const BINOP& OP, const Element& u, const Element& v ) const
    {
        GIVARO_ASSERT((dim(res) == dim(u)) && (dim(u) == dim(v)), "Bad size");
        size_t sz = dim(res);
        for (size_t i=0; i<sz; ++i)
            OP(res[i], u[i], v[i]);
    }

    template<class Domain>
    template<class UNOP>
    inline void VectorDom<Domain,Dense>::
    map( Element& res, UNOP& OP, const Element& u ) const
    {
        GIVARO_ASSERT((dim(res) == dim(u)), "Bad size");
        size_t sz = dim(u);
        for (size_t i=0; i<sz; ++i)
            OP(res[i], u[i]);
    }

    // --- mul
    template<class Domain>
    inline void VectorDom<Domain,Dense>::mulin( Element& res, const Type_t& val ) const
    {
        Curried2<MulOp<Domain> > opcode(_domain, (Type_t&)val);
        map( res, opcode);
    }

    template<class Domain>
    inline void VectorDom<Domain,Dense>::mul
    ( Element& res, const Element& u, const Type_t& val ) const
    {
        GIVARO_ASSERT((dim(res) == dim(u)), "Bad size");
        Curried2<MulOp<Domain> > opcode(_domain, (Type_t&)val);
        map( res, opcode, u);
    }

    template<class Domain>
    inline void VectorDom<Domain,Dense>::mul
    ( Element& res, const Type_t& val, const Element& v ) const
    {
        GIVARO_ASSERT((dim(res) == dim(v)), "Bad size");
        Curried1<MulOp<Domain> > opcode(_domain, val);
        map( res, opcode, v);
    }



    // --- add
    template<class Domain>
    inline void VectorDom<Domain,Dense>::addin( Element& res, const Element& u ) const
    {
        GIVARO_ASSERT((dim(res) == dim(u)), "Bad size");
        AddOp<Domain> opcode( _domain );
        map( res, opcode, u );
    }

    template<class Domain>
    inline void VectorDom<Domain,Dense>::add
    ( Element& res, const Element& u, const Element& v ) const
    {
        GIVARO_ASSERT((dim(res) == dim(u)) && (dim(u) == dim(v)), "Bad size");
        AddOp<Domain> opcode (_domain);
        map( res, opcode, u, v);
    }

    template<class Domain>
    inline void VectorDom<Domain,Dense>::add
    ( Element& res, const Element& u, const Type_t& val ) const
    {
        GIVARO_ASSERT((dim(res) == dim(u)), "Bad size");
        Curried2<AddOp<Domain> > opcode(_domain, (Type_t&)val);
        map( res, opcode, u);
    }

    template<class Domain>
    inline void VectorDom<Domain,Dense>::add
    ( Element& res, const Type_t& val, const Element& v ) const
    {
        GIVARO_ASSERT((dim(res) == dim(v)), "Bad size");
        Curried1<AddOp<Domain> > opcode(_domain, val);
        map( res, opcode, v);
    }

    template<class Domain>
    inline void VectorDom<Domain,Dense>::add
    ( Element& res, const SparseVector& u, const Element& v ) const
    {
        GIVARO_ASSERT((dim(res) == dim(v)) && (dim(v) == u.dim()), "Bad size");
        // -- here assume : subdomains of res, u and v are equal
        assign(res, v);
        addin(res, u);
    }

    template<class Domain>
    inline void VectorDom<Domain,Dense>::add
    ( Element& res, const Element& v, const SparseVector& u ) const
    {
        GIVARO_ASSERT((dim(res) == dim(v)) && (dim(v) == u.dim()), "Bad size");
        // -- here assume : subdomains of res, u and v are equal
        assign(res, v);
        addin(res,u);
    }

    template<class Domain>
    inline void VectorDom<Domain,Dense>::addin
    ( Element& res, const SparseVector& u ) const
    {
        GIVARO_ASSERT((dim(res) == u.dim()), "Bad size");
        // -- here assume : subdomains of res, u and v are equal
        typedef VectorDom<Domain,Sparse> VSparseDom;
        typedef typename VSparseDom::constIterator_t cIterator_t;
        typedef typename VSparseDom::IndiceIterator_t cIndiceIterator_t;
        cIterator_t u_curr = u.begin_data();
        cIterator_t u_end = u.end_data();
        cIndiceIterator_t u_indice = u.begin_indice();
        for ( ; u_curr != u_end; ++u_curr, ++u_indice )
            _domain.addin(res[*u_indice], *u_curr);
    }



    // --- sub
    template<class Domain>
    inline void VectorDom<Domain,Dense>::subin( Element& res, const Element& u ) const
    {
        GIVARO_ASSERT((dim(res) == dim(u)), "Bad size");
        SubOp<Domain> opcode ( _domain );
        map( res, opcode, u );
    }

    template<class Domain>
    inline void VectorDom<Domain,Dense>::sub( Element& res, const Element& u, const Element& v ) const
    {
        GIVARO_ASSERT((dim(res) == dim(v)) && (dim(v) == dim(u)), "Bad size");
        SubOp<Domain> opcode ( _domain );
        map( res, opcode, u, v);
    }

    template<class Domain>
    inline void VectorDom<Domain,Dense>::sub( Element& res, const Element& u, const Type_t& val ) const
    {
        GIVARO_ASSERT((dim(res) == dim(u)), "Bad size");
        Curried2<SubOp<Domain> > opcode( _domain, (Type_t&)val);
        map( res, opcode, u);
    }

    template<class Domain>
    inline void VectorDom<Domain,Dense>::sub( Element& res, const Type_t& val, const Element& v ) const
    {
        GIVARO_ASSERT((dim(res) == dim(v)), "Bad size");
        Curried1<SubOp<Domain> > opcode( _domain, val);
        map( res, opcode, v);
    }

    template<class Domain>
    inline void VectorDom<Domain,Dense>::subin
    ( Element& res, const SparseVector& u ) const
    {
        GIVARO_ASSERT((dim(res) == u.dim()), "Bad size");
        // -- here assume : subdomains of res, u and v are equal
        typedef VectorDom<Domain,Sparse> VSparseDom;
        typedef typename VSparseDom::constIterator_t cIterator_t;
        typedef typename VSparseDom::IndiceIterator_t cIndiceIterator_t;
        cIterator_t u_curr = u.begin_data();
        cIterator_t u_end = u.end_data();
        cIndiceIterator_t u_indice = u.begin_indice();
        for ( ; u_curr != u_end; ++u_curr, ++u_indice )
            _domain.subin(res[*u_indice], *u_curr);
    }

    template<class Domain>
    inline void VectorDom<Domain,Dense>::sub
    ( Element& res, const SparseVector& u, const Element& v ) const
    {
        GIVARO_ASSERT((dim(res) == dim(v)) && (dim(v) == u.dim()), "Bad size");
        // -- here assume : subdomains of res, u and v are equal
        neg(res, v);
        addin(res, u);
    }

    template<class Domain>
    inline void VectorDom<Domain,Dense>::sub
    ( Element& res, const Element& v, const SparseVector& u ) const
    {
        GIVARO_ASSERT((dim(res) == dim(v)) && (dim(v) == u.dim()), "Bad size");
        // -- here assume : subdomains of res, u and v are equal
        assign(res,v);
        subin(res,u);
    }


    // --- sub

    // --- neg
    template<class Domain>
    inline void VectorDom<Domain,Dense>::negin( Element& res ) const
    {
        NegOp<Domain> opcode ( _domain );
        map( res, opcode, res );
    }

    template<class Domain>
    inline void VectorDom<Domain,Dense>::neg( Element& res, const Element& u ) const
    {
        GIVARO_ASSERT((dim(res) == dim(u)), "Bad size");
        NegOp<Domain> opcode ( _domain );
        map( res, opcode, u );
    }

    template<class Domain>
    void VectorDom<Domain,Dense>::axpy
    ( Element& res, const Type_t& a, const Element& x, const Element& y ) const
    {
        GIVARO_ASSERT((dim(res) == dim(x)) && (dim(x) == dim(y)), "Bad size");
        size_t sz = dim(res);
        for ( size_t i=0; i<sz; ++i )
            _domain.axpy(res[i], a, x[i], y[i]);
    }

    template<class Domain>
    void VectorDom<Domain,Dense>::axpyin
    ( Element& res, const Type_t& a, const Element& x ) const
    {
        GIVARO_ASSERT((dim(res) == dim(x)), "Bad size");
        size_t sz = dim(res);
        for ( size_t i=0; i<sz; ++i )
            _domain.axpyin(res[i], a, x[i]);
    }

    template<class Domain>
    void VectorDom<Domain,Dense>::axmy
    ( Element& res, const Type_t& a, const Element& x, const Element& y ) const
    {
        GIVARO_ASSERT((dim(res) == dim(x)) && (dim(x) == dim(y)), "Bad size");
        size_t sz = dim(res);
        for ( size_t i=0; i<sz; ++i )
            _domain.axmy(res[i], a, x[i], y[i]);
    }

    template<class Domain>
    void VectorDom<Domain,Dense>::axmyin
    ( Element& res, const Type_t& a, const Element& x ) const
    {
        GIVARO_ASSERT((dim(res) == dim(x)), "Bad size");
        size_t sz = dim(res);
        for ( size_t i=0; i<sz; ++i )
            _domain.axmyin(res[i], a, x[i]);
    }


    template<class Domain>
    void VectorDom<Domain,Dense>::axpy
    ( Element& res, const Element& a, const Element& b, const Element& y ) const
    {
        GIVARO_ASSERT((dim(res) == dim(a)) && (dim(a) == dim(b)) && (dim(b) == dim(y)), "Bad size");
        size_t sz = dim(res);
        for ( size_t i=0; i<sz; ++i )
            _domain.axpy(res[i], a[i], b[i], y[i]);
    }

    template<class Domain>
    void VectorDom<Domain,Dense>::axpyin
    ( Element& res, const Element& a, const Element& b ) const
    {
        GIVARO_ASSERT((dim(res) == dim(a)) && (dim(a) == dim(b)), "Bad size");
        size_t sz = dim(res);
        for ( size_t i=0; i<sz; ++i )
            _domain.axpyin(res[i], a[i], b[i]);
    }

    template<class Domain>
    void VectorDom<Domain,Dense>::axmy
    ( Element& res, const Element& a, const Element& b, const Element& y ) const
    {
        GIVARO_ASSERT((dim(res) == dim(a)) && (dim(a) == dim(b)) && (dim(b) == dim(y)), "Bad size");
        size_t sz = dim(res);
        for ( size_t i=0; i<sz; ++i )
            _domain.axmy(res[i], a[i], b[i], y[i]);
    }

    template<class Domain>
    void VectorDom<Domain,Dense>::axmyin
    ( Element& res, const Element& a, const Element& b ) const
    {
        GIVARO_ASSERT((dim(res) == dim(a)) && (dim(a) == dim(b)), "Bad size");
        size_t sz = dim(res);
        for ( size_t i=0; i<sz; ++i )
            _domain.axmyin(res[i], a[i], b[i]);
    }



    template<class Domain>
    inline void VectorDom<Domain,Dense>::dot
    ( Type_t& res, const Element& op1, const Element& op2) const
    {
        GIVARO_ASSERT((dim(op1) == dim(op2)), "Bad size");
        size_t sz = dim(op1);
        _domain.assign(res, _domain.zero);
        for ( size_t i=0; i<sz; ++i ) _domain.axpyin(res, op1[i], op2[i] );
    }



    // ==========================================================================
    //
    // -- Write the domain
    //
    template<class Domain>
    std::ostream& VectorDom<Domain, Dense>::write( std::ostream& o ) const
    {
        return _domain.write(o << "(") << ",Dense)";
    }

    // -- read the domain
    template<class Domain>
    std::istream& VectorDom<Domain, Dense>::read( std::istream& sin )
    {
        char ch;
        sin >> std::ws >> ch;
        if (ch != '(')
            GivError::throw_error(
                                  GivBadFormat("VectorDom<Domain,Dense>::read: syntax error no '('"));

        // -- read subdomain:
        _domain.read(sin);

        // -- read ,
        sin >> std::ws >> ch;
        if (ch != ',')
            GivError::throw_error(
                                  GivBadFormat("VectorDom<Domain,Dense>::read: syntax error no ','"));

        // -- read dense:
        char name[6];
        sin >> std::ws; sin.read(name, 5); name[5] = (char)0;
        if (strcmp(name, "Dense") !=0)
            GivError::throw_error(
                                  GivBadFormat("VectorDom<Domain,Dense>::read: syntax error no 'Dense'"));

        // -- read )
        sin >> std::ws >> ch;
        if (ch != ')')
            GivError::throw_error(
                                  GivBadFormat("VectorDom<Domain,Dense>::read: syntax error no ')'"));
        return sin;
    }


    // ==========================================================================
    //
    // -- Don't write the domain !, only write a rep / domain
    template<class Domain>
    std::ostream& VectorDom<Domain, Dense>::write(std::ostream& o, const Element& V) const
    {
        if (dim(V) ==0) return o << "[ ]";
        o << '[';
        size_t i;
        for (i=0; i<dim(V)-1; i++) _domain.write(o, V[i]) << ',';
        return _domain.write(o, V[i]) << ']';
    }


    //
    // Read a vector given by a Maple-like format:
    // the grammar is :
    //   s  --->  '[' list_of_elt ']'
    //   list_of_elt --->   elt
    //                    | list_of_elt ',' elt
    // The contraints are :
    //   All lines are the same number of Elements.
    //   The separators a those of the C lexical-convention, i.e.
    //    ' ', '\n', '\t', '\f' .

    template<class Domain>
    std::istream&  VectorDom<Domain,Dense>::read (std::istream& fin, Element& A) const
    {
        char ch;

        // -- Skip the first "white":
        fin >> std::ws; fin.get(ch);
        if (ch != '[')
            GivError::throw_error(
                                  GivBadFormat("VectorDom<Domain,Dense>::read: syntax error no '['"));

        // -- Read the line and count the nb of elts
        int i = 0;
        Element rep;
        rep.allocate(1);
        Type_t Tmp;
        fin >> Tmp;
        fin >> std::ws; fin.get(ch);
        rep[0] = Tmp;
        while (ch != ']')
        {
            if (ch != ',')
                GivError::throw_error(
                                      GivBadFormat("VectorDom<T,Dense>::read: syntax error no ','"));
            i++;
            fin >> std::ws >> Tmp;
            fin >> std::ws; fin.get(ch);

            // resize the vector :
            rep.resize(i+1);
            rep[i] = Tmp;
        }
        A.logcopy( rep );
        return fin;
    }
} // Givaro
#include "givaro/givvectdensespe.inl"
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
