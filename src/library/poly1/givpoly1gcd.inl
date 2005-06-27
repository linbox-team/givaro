// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givpoly1gcd.inl,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Author: J-L. Roch, T. Gautier, J-G. Dumas
// $Id: givpoly1gcd.inl,v 1.5 2005-06-27 09:07:21 jgdumas Exp $
// ==========================================================================
// friend void bezout (const Poly1<T> &P,
//                     const Poly1<T> &Q,
//                     Poly1<T> &d,
//                     Poly1<T> &u,
//                     Poly1<T> &v);
// computes d = unitary gcd(P,Q) and u,v such that :
//   u*P+v*Q = d
// u and v are the unique polynomisals such that :
//   deg(u) <= deg(Q)-deg(d) and deg(v) <= deg(P) - deg(d)
// friend Poly1<T> gcd (const Poly1<T> &P, const Poly1<T> &Q );
// computes d = unitary gcd(P,Q) 
// ===========================================================================
// A coder par PRS : cf Geddes page 283
template <class Domain> 
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::gcd ( Rep& G, const Rep& P, const Rep& Q) const
{
  Rep U,V; 
  Degree degU, degV;
  degree(degU,P); degree(degV,Q);
  if ((degU < 0) || (degV == 0)) return assign(G, Q);
  if ((degV < 0) || (degU == 0)) return assign(G, P);
  
  if (degU >= degV) {
    assign(U, P);
    assign(V, Q);
  }
  else {
    assign(U, Q);
    assign(V, P);
  }
  // -- PRS (U,V) using pmod:
  Type_t g;
  _domain.assign(g, _domain.one);
  Degree degR;
  Rep R;
  do {
// write(cout << "gcd: U:", U) << endl;
// write(cout << "gcd: V:", V) << endl;
    mod( R, U, V); 
    setdegree(R);
    degree(degR, R);
// write(cout << "mod: R:", R) << endl;
    if (degR < 0) break;
    assign(U,V); 
    assign(V,R); 
  } while (1);
  
  degree(degV, V);
  G.logcopy(V);
// JGD 15.12.1999
//   if (degV <= 1) assign(G,one);
  if (degV <= 0) assign(G,one);
// write(cout << "gcd: G:", G) << endl;
  return G;
}

template <class Domain> 
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::gcd ( Rep& F, Rep& S0, Rep& T0, const Rep& A, const Rep& B) const
{
  Type_t r0, r1, tt;
  Rep G; 
  Degree degF, degG;
  degree(degF,A); degree(degG,B);
  if ((degF < 0) || (degG == 0)) {
      assign(T0, Degree(0), _domain.inv( tt, leadcoef(r0,B)));
      init(S0, 0);
      assign(F, B);
      return mulin(F,tt);
  }
  if ((degG < 0) || (degF == 0)) {
      assign(S0, Degree(0), _domain.inv( tt, leadcoef(r0,A)));
      init(T0, 0);
      assign(F, A);
      return mulin(F,tt);
  }
  
  
//   if (degF >= degG) {
    assign(F, A);
    assign(G, B);
//   }
//   else {
//     assign(F, B);
//     assign(G, A);
//   }

  leadcoef(r0,  F);
  leadcoef(r1, G);

  divin(F,r0);
  divin(G,r1);

  Rep S1,R1,T1,Q,TMP, TMP2;

  assign(S0, 0, _domain.inv(tt,r0) );
  assign(S1,zero);
  assign(T0,zero);
  assign(T1, 0, _domain.inv(tt,r1) );
  
  while ( ! isZero(G) ) {
      divmod(Q,R1,F,G);
 leadcoef(r1, R1); if (_domain.isZero(r1)) _domain.assign(r1,_domain.one);
      assign(F,G);
      div(G,R1,r1);
      mul(TMP,Q,S1); sub(TMP2,S0,TMP); assign(S0,S1); div(S1,TMP2,r1);
      mul(TMP,Q,T1); sub(TMP2,T0,TMP); assign(T0,T1); div(T1,TMP2,r1);
  }
  
  return F;
}


template <class Domain> 
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::invmod ( Rep& S0, const Rep& A, const Rep& B) const
{
  Type_t r0, r1, tt;
  Rep F, G; 
  Degree degF, degG;
  degree(degF,A); degree(degG,B);
  if ((degF <= 0) || (degG <= 0) ) {
      return assign(S0, Degree(0), _domain.inv( tt, leadcoef(r0,A)));
  }
  
//   if (degF >= degG) {
    assign(F, A);
    assign(G, B);
//   }
//   else {
//     assign(F, B);
//     assign(G, A);
//   }

  leadcoef(r0, F);
  leadcoef(r1, G);

  divin(F,r0);
  divin(G,r1);

  Rep S1,R1,Q,TMP, TMP2;

  assign(S0, 0, _domain.inv(tt,r0) );
  assign(S1,zero);
  
  while ( ! isZero(G) ) {
      divmod(Q,R1,F,G);
      leadcoef(r1, R1); if (_domain.isZero(r1)) _domain.assign(r1,_domain.one);
      assign(F,G);
      div(G,R1,r1);
      mul(TMP,Q,S1); sub(TMP2,S0,TMP); assign(S0,S1); div(S1,TMP2,r1);
  }
  
  return S0;
}


template <class Domain> 
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::lcm ( Rep& F, const Rep& A, const Rep& B) const
{
  Rep G, S0, T0; 
  Degree degA, degB, degG;
  degree(degA,A); degree(degB,B);
  if ((degA < 0) || (degB == 0)) return assign(F, B);
  if ((degB < 0) || (degA == 0)) return assign(F, A);
  
  if (degA >= degB) {
    assign(F, A);
    assign(G, B);
  }
  else {
    assign(F, B);
    assign(G, A);
  }

  Type_t r0, r1, tt;
  leadcoef(r0, F);
  leadcoef(r1, G);

  divin(F,r0);
  divin(G,r1);

  Rep S1,R1,T1,Q,TMP, TMP2;

  assign(S0, 0, _domain.inv(tt,r0) );
  assign(S1,zero);
  assign(T0,zero);
  assign(T1, 0, _domain.inv(tt,r1) );
  
  while ( ! isZero(G) ) {
      divmod(Q,R1,F,G);
 leadcoef(r1, R1); if (_domain.isZero(r1)) _domain.assign(r1,_domain.one);
      assign(F,G);
      div(G,R1,r1);
      mul(TMP,Q,S1); sub(TMP2,S0,TMP); assign(S0,S1); div(S1,TMP2,r1);
      mul(TMP,Q,T1); sub(TMP2,T0,TMP); assign(T0,T1); div(T1,TMP2,r1);


//   std::cerr << "Degree Q : " << degree(degF,Q) << std::endl;
//   std::cerr << "Degree R1 : " << degree(degF,R1) << std::endl;
//   std::cerr << "Degree F : " << degree(degF,F) << std::endl;
//   std::cerr << "Degree G : " << degree(degF,G) << std::endl;
//   std::cerr << "Degree S0 : " << degree(degF,S0) << std::endl;
//   std::cerr << "Degree T0 : " << degree(degF,T0) << std::endl;
//   std::cerr << "Degree S1 : " << degree(degF,S1) << std::endl;
//   std::cerr << "Degree T1 : " << degree(degF,T1) << std::endl;

//   write( write( write( write( write( std::cerr, G ) << " = ((", T1) << ") * (", A) << ")) + ((", S1) << ") * (", B) << "))" << std::endl;

  }
  
  degree(degG, G);

    
  if ( degG <= 0) {
// If normalisation is needed
//       if (degA >= degB) {
//           return mul(G, S1, A);
//       } else {
//           return mul(G, S1, B);
//       }
//       leadcoef(tt, G);
//       return div(F,G,tt);
      if (degA >= degB) {
          return mul(F, S1, A);
      } else {
          return mul(F, S1, B);
      }
  } else {
      return mul(F, A, B);
  }
}





























#if 0
// A coder par PRS : cf Geddes page 283
template <class Domain> 
void Poly1Dom<Domain,Dense>::gcd 
 ( Poly1<T> &d, Poly1<T> &u, Poly1<T> &v,
   const Poly1<T> &P, const Poly1<T> &Q ) 
{
  // computes d = unitary gcd(P,Q) and u,v such that :
  //   u*P+v*Q = d
  // u and v are the unique polynomisals such that :
  //   deg(u) <= deg(Q)-deg(d) and deg(v) <= deg(P) - deg(d)
 
  // Auxilliary matrix : [oldu oldv] which is left multiplies by [0  1]
  //                     [newu newv]                             [1 -q]
  // at each step. 
  // At the end, u = oldu, v = old v

  Rep A;
  Rep B;
  Rep quot;
  Rep rem;
  int permuter;

  if (degree(P) < degree(Q)) { 
    assign(A, Q); 
    assign(B, P); 
    permuter = 1; 
  }
  else { 
    assign(A, P); 
    assign(B, Q); 
    permuter = 0; 
  }

  if (isZero(B) {
    Type_t inv_lcoeff_A;
    if (isZero(A)) { 
      _domain.assign(inv_lcoeff_A, _domain.one); 
    }
    else { 
      _domain.inv(inv_lcoeff_A, A[degree(A)]; 
    }

    d = A * inv_lcoeff_A;
    Poly1<T> cste(0, inv_lcoeff_A);
    if (permuter) { u = Poly1<T>::Zero; v = cste; }
    else { u = cste; v = Poly1<T>::Zero; }
  }
  else {
    // -- On rend B unitaire 
    T inv_lcoeff_B = csteT1 / B[B.degree()]; 
    B = B * inv_lcoeff_B;

    Poly1<T> oldu(0, csteT1);
    Poly1<T> oldv(0, csteT0 );
    Poly1<T> newu(0, csteT0 );
    Poly1<T> newv(0, inv_lcoeff_B);

    int cont;
    do {
      Poly1<T>::divide(quot, rem, A, B);
      A = B;
      cont = !isZero(rem);
      if (cont) { 
        inv_lcoeff_B = csteT1/rem[rem.degree()]; 
        B = rem*inv_lcoeff_B;

        Poly1<T> tmpu = (oldu - quot*newu) * inv_lcoeff_B;
        Poly1<T> tmpv = (oldv - quot*newv) * inv_lcoeff_B;
        oldu = newu;
        oldv = newv;
        newu = tmpu;
        newv = tmpv;
      }
    } while (cont); 

    
    d = A;

    if (permuter) { u = newv; v = newu; }
    else { u = newu; v = newv; }
  }

#ifdef GIVARO_ASSERT
  if (!isZero(u*P+v*Q-d)) {
    cout << "Erreur dans la verif de bezout. " << endl
         << "   P  =" << P << endl
         << "   Q  =" << Q << endl
         << "   d  =" << d  << endl
         << "   u  =" << u << endl 
         << "   v  =" << v << endl;
  }
#endif  
}

template <class T>
Poly1<T>& Poly1<T>::gcd1 (Poly1<T>& G, Poly1<T>& u, Poly1<T>& v, 
                          const Poly1<T>& P, const Poly1<T>& Q) 
{ bezout(P,Q,G,u,v); return G;}


template <class T>
Poly1<T>& Poly1<T>::gcd1 (Poly1<T>& G, const Poly1<T> &P, const Poly1<T> &Q) 
{
  Poly1<T> A;
  Poly1<T> B;

  if (P.degree() < Q.degree()) { A = Q; B = P; }
  else { A = P; B = Q; }

  if (isZero(B)) {
    if (isZero(A)) { return G.logcopy(A); }
    else { return G.logcopy(A/A[A.degree()]); }
  }

  B = B/B[B.degree()]; // B is made unitary
  Poly1<T> rem = A % B;
  while (!isZero(rem = A % B)) {
    A = B;
    B = rem / rem.leadcoef();
  }
  return G.logcopy(B);
}


template<class T>
inline const Poly1<T> prime_part(const Poly1<T>& P, const Poly1<T>& Q)
{ 
  Poly1<T> E,D;
  E = P;
  D = Q;
  while (D.degree() >0)
  {
    E = E / D;
    D = gcd(E,D);
  }  
  return E;
}

// Split the argument into two factors prime_Q and divisor_Q such that :
// divisor_Q = gcd(*this, Q)
// There exists a constant k such that P divides (prime_Q * divisor_Q)^k
template<class T>
inline void decomposition(const Poly1<T> &Q, Poly1<T> & prime_Q, Poly1<T> & divisor_Q)
{
    divisor_Q = gcd(*this, Q);
    prime_Q= *this / gcd( ::pow(divisor_Q,((long)(degree() - divisor_Q.degree()))),
                        *this );
}

#endif // 0
