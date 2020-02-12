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

    template<class Domain>
    template<class UNOP>
    inline void VectorDom<Domain,Sparse>::
    map( Element& res, const UNOP& OP, const Element& u ) const
    {
        res.copy(u);
        size_t sz = res.size();
        for (size_t i=0; i<sz; ++i)
            OP(res._data[i], u._data[i]);
    }

    template<class Domain>
    template<class UNOP>
    inline void VectorDom<Domain,Sparse>::
    map( Element& res, UNOP& OP, const Element& u ) const
    {
        res.copy(u);
        size_t sz = res.size();
        for (size_t i=0; i<sz; ++i)
            OP(res._data[i], u._data[i]);
    }

    template<class Domain>
    int VectorDom<Domain,Sparse>::areEqual ( const Element& P, const Element& Q) const
    {
        size_t d;
        if ((d =dim(P)) != dim(Q)) return 0;
        if (P._index.size() != Q._index.size()) return 0;
        for (size_t i=0; i<d; ++i) {
            if (P._index[i] != Q._index[i]) return 0;
            if (!_domain.areEqual(P._data[i], Q._data[i])) return 0;
        }
        return 1;
    }

    template<class Domain>
    int VectorDom<Domain,Sparse>::areNEqual( const Element& P, const Element& Q) const
    {
        return !areEqual(P,Q);
    }

    template<class Domain>
    int VectorDom<Domain,Sparse>::isZero  ( const Element& P ) const
    {
        size_t d;
        if ((d =dim(P)) == 0) return 1;
        for (size_t i=0; i<d; ++i)
            if (!_domain.isZero(P._data[i])) return 0;
        return 1;
    }


    // --
    // -- Compression method to compact a dense vector
    // --
    template<class Domain>
    void VectorDom<Domain,Sparse>::compact (
                                            Element& u,
                                            const VectorDom<Domain, Dense>& VDom,
                                            const typename VectorDom<Domain, Dense>::Element& v ) const
    {
        size_t dim = VDom.dim(v);
        u.allocate(dim, 0);
        Indice_t pos_next =0;
        for (size_t i=0; i< dim; ++i)
        {
            if ( !_domain.isZero(v[i]) ) {
                u.resize( dim, pos_next + 1 );
                u._index[pos_next] = i;
                _domain.assign(u._data[pos_next], v[i]);
                ++pos_next;
            }
        }
    }

    // --
    // -- Compression method to compact a sparse vector
    // --
    template<class Domain>
    void VectorDom<Domain,Sparse>::compact (
                                            Element& u,
                                            const VectorDom<Domain, Sparse>& VDom,
                                            const typename VectorDom<Domain, Sparse>::Element& v ) const
    {
        u.copy(v);
    }


    template<class Domain>
    inline void VectorDom<Domain,Sparse>::dot
    ( Type_t& res, const Element& op1, const Element& op2) const
    {
        const Domain_t& domain = subdomain();
        domain.assign(res, domain.zero);
        size_t op1size =op1._index.size(), op2size =op2._index.size();
        for ( long i =op1size-1, j=op2size-1; (i >=0) && (j >=0); ) {
            long diff = op1._index[i] - op2._index[j];
            if (diff >0) --i;
            else if (diff <0) --j;
            else domain.axpy(res, op1._data[i--], op2._data[j--], res );
        }
    }

    template<class Domain>
    inline void VectorDom<Domain,Sparse>::add
    ( Element& res, const Element& op1, const Element& op2) const
    {
        long i,j;
        size_t curr =0;
        size_t op1size =op1._index.size(), op2size =op2._index.size();
        res.resize(dim(op1), op1size+op2size);
        for (i =0, j=0; (i <op1size) && (j <op2size); ) {
            long diff = op1._index[i] - op2._index[j];
            if (diff <0) {
                res._index[curr] = op1._index[i];
                _domain.assign(res._data[curr],op1._data[i]);
                ++i;
                if (!_domain.isZero(res._data[curr])) ++curr;
            }
            else if (diff >0) {
                res._index[curr] = op2._index[j];
                _domain.assign(res._data[curr],op2._data[j]);
                ++j;
                if (!_domain.isZero(res._data[curr])) ++curr;
            }
            else {
                res._index[curr] = op1._index[i];
                _domain.add(res._data[curr], op1._data[i], op2._data[j]);
                ++i; ++j;
                if (!_domain.isZero(res._data[curr])) ++curr;
            }
        }
        // -- Test if i != op1size or j != op2size, then complete the res
        if (i == op1size) {
            for ( ; j < op2size; ++j, ++curr) {
                res._index[curr] = op2._index[j];
                _domain.assign(res._data[curr],op2._data[j]);
            }
        }
        else if (j == op2size) {
            for ( ; i < op1size; ++i, ++curr) {
                res._index[curr] = op1._index[i];
                _domain.assign(res._data[curr],op1._data[i]);
            }
        }
        // -- Set the correct size:
        res.resize( dim(op1), curr );
    }


    template<class Domain>
    inline void VectorDom<Domain,Sparse>::sub
    ( Element& res, const Element& op1, const Element& op2) const
    {
        long i,j;
        size_t curr =0;
        size_t op1size =op1._index.size(), op2size =op2._index.size();
        res.resize(dim(op1), op1size+op2size);
        for (i =0, j=0; (i <op1size) && (j <op2size); ) {
            long diff = op1._index[i] - op2._index[j];
            if (diff <0) {
                res._index[curr] = op1._index[i];
                _domain.assign(res._data[curr],op1._data[i]);
                ++i;
                if (!_domain.isZero(res._data[curr])) ++curr;
            }
            else if (diff >0) {
                res._index[curr] = op2._index[j];
                _domain.neg(res._data[curr],op2._data[j]);
                ++j;
                if (!_domain.isZero(res._data[curr])) ++curr;
            }
            else {
                res._index[curr] = op1._index[i];
                _domain.sub(res._data[curr], op1._data[i], op2._data[j]);
                ++i; ++j;
                if (!_domain.isZero(res._data[curr])) ++curr;
            }
        }
        // -- Test if i != op1size or j != op2size, then complete the res
        if (i == op1size) {
            for ( ; j < op2size; ++j, ++curr) {
                res._index[curr] = op2._index[j];
                _domain.assign(res._data[curr],op2._data[j]);
            }
        }
        else if (j == op2size) {
            for ( ; i < op1size; ++i, ++curr) {
                res._index[curr] = op1._index[i];
                _domain.neg(res._data[curr],op1._data[i]);
            }
        }
        // -- Set the correct size:
        res.resize( dim(op1), curr );
    }


    template<class Domain>
    inline void VectorDom<Domain,Sparse>::addin
    ( Element& res, const Element& u ) const
    { Element tmp; init( tmp ); add(tmp, res, u); assign(res, tmp); }

    template<class Domain>
    inline void VectorDom<Domain,Sparse>::add
    ( Element& res, const Element& u, const Type_t& val ) const
    {
        Curried2<AddOp<Domain> > opcode(_domain, val);
        map( res, opcode, u);
    }

    template<class Domain>
    inline void VectorDom<Domain,Sparse>::add ( Element& res, const Type_t& val, const Element& u ) const
    {
        Curried1<AddOp<Domain> > opcode(_domain, val);
        map( res, opcode, u);
    }

    template<class Domain>
    inline void VectorDom<Domain,Sparse>::subin
    ( Element& res, const Element& u ) const
    { Element tmp; init( tmp ); sub(tmp, res, u); assign(res, tmp); }

    template<class Domain>
    inline void VectorDom<Domain,Sparse>::sub ( Element& res, const Element& u, const Type_t& val ) const
    {
        Curried2<SubOp<Domain> > opcode(_domain, val);
        map( res, opcode, u);
    }

    template<class Domain>
    inline void VectorDom<Domain,Sparse>::sub ( Element& res, const Type_t& val, const Element& u ) const
    {
        Curried2<SubOp<Domain> > opcode(_domain, val);
        map( res, opcode, u);
    }

    template<class Domain>
    inline void VectorDom<Domain,Sparse>::negin ( Element& res ) const
    {
        NegOp<Domain> opcode ( _domain );
        map( res, opcode, res );
    }

    template<class Domain>
    inline void VectorDom<Domain,Sparse>::neg ( Element& res, const Element& u ) const
    {
        NegOp<Domain> opcode ( _domain );
        map( res, opcode, u );
    }




    // ==========================================================================
    //
    // -- Write the domain
    template<class Domain>
    std::ostream& VectorDom<Domain, Sparse>::write( std::ostream& o ) const
    {
        return _domain.write(o << '(') << ",Sparse)";
    }

    // -- read the domain
    template<class Domain>
    std::istream& VectorDom<Domain, Sparse>::read( std::istream& sin )
    {
        char ch;
        sin >> std::ws >> ch;
        if (ch != '(')
            GivError::throw_error(
                                  GivBadFormat("VectorDom<Domain,Sparse>::read: syntax error no '('"));

        // -- read the subdomain
        _domain.read(sin);

        // -- read ,
        sin >> std::ws >> ch;
        if (ch != ',')
            GivError::throw_error(
                                  GivBadFormat("VectorDom<Domain,Sparse>::read: syntax error no ','"));

        // -- read dense:
        char name[7];
        sin >> std::ws; sin.read(name, 6); name[6] = (char)0;
        if (strcmp(name, "Sparse") !=0)
            GivError::throw_error(
                                  GivBadFormat("VectorDom<Domain,Sparse>::read: syntax error no 'Sparse'"));

        // -- read )
        sin >> std::ws >> ch;
        if (ch != ')')
            GivError::throw_error(
                                  GivBadFormat("VectorDom<Domain,Dense>::read: syntax error no ')'"));
        return sin;
    }


    // ==========================================================================
    //
    // Write an Element of the domain, not the domain itself.
    // grammar: list of pairs:
    //   s  --->  '[' size, '[' list_of_elt ']]'
    //   list_of_elt --->   (index, value)
    //                    | list_of_elt ',' (index, value)
    template<class Domain>
    std::ostream& VectorDom<Domain,Sparse>::write (std::ostream& o, const Element& V) const
    {
        if (dim(V) ==0) return o << "[0,[]]";
        size_t sz = V.size();
        o << '[' << sz << ",[";
        if (sz != 0) {
            Pair<Indice_t, Type_t> p (V._index[0], V._data[0]);
            o << p;
        }
        for (size_t i=1; i< sz; ++i) {
            Pair<Indice_t, Type_t> p (V._index[i], V._data[i]);
            o << ',' << p;
        }
        return o << "]]";
    }

    //
    // Read an Element of the domain, not the domain itself.
    // Read a sparse vector given by the grammar:
    //   s  --->  '[' size, '[' list_of_elt ']]'
    //   list_of_elt --->   (index, value)
    //                    | list_of_elt ',' (index, value)
    // The contraints are :
    //   All lines are the same number of Elements.
    //   The separators a those of the C lexical-convention, i.e.
    //    ' ', '\n', '\t', '\f' .
    template<class Domain>
    std::istream&  VectorDom<Domain,Sparse>::read (std::istream& fin, Element& V) const
    {
        char ch;

        // -- Skip the first blanks:
        fin >> std::ws; fin.get(ch);
        if (ch != '[')
            GivError::throw_error(
                                  GivBadFormat("VectorDom<Domain,Sparse>::read: syntax error no '['"));

        // -- Read the size of != 0 Element of the rep
        size_t size;
        fin >> std::ws >> size;

        // -- read ,
        fin >> std::ws; fin.get(ch);
        if (ch != ',')
            GivError::throw_error(
                                  GivBadFormat("VectorDom<Domain,Sparse>::read: syntax error no ','"));

        // -- read [
        fin >> std::ws; fin.get(ch);
        if (ch != '[')
            GivError::throw_error(
                                  GivBadFormat("VectorDom<Domain,Sparse>::read: syntax error no ']'"));

        // -- read the size pairs
        size_t i;
        V.allocate(size, size);
        for (i=0; i<size-1; ++i)
        {
            Pair<Indice_t, Type_t> p;
            fin >> p;
            V._index[i] = p.first();
            V._data[i] = p.second();

            // -- read ,
            fin >> std::ws; fin.get(ch);
            if (ch != ',')
                GivError::throw_error(
                                      GivBadFormat("VectorDom<Domain,Sparse>::read: syntax error no ','"));
        }

        {
            Pair<Indice_t, Type_t> p;
            fin >> p;
            V._index[size-1] = p.first();
            V._data[size-1] = p.second();
        }

        // -- read ]
        fin >> std::ws; fin.get(ch);
        if (ch != ']')
            GivError::throw_error(
                                  GivBadFormat("VectorDom<Domain,Sparse>::read: syntax error no ']'"));

        // -- read ]

        // - read the last characters:
        fin >> std::ws; fin.get(ch);
        if (ch != ']')
            GivError::throw_error(
                                  GivBadFormat("operator>><Vector<T,Sparse> >: syntax error no ']'"));
        return fin;
    }

} // Givaro
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
