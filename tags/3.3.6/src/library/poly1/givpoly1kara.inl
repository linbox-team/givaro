// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givpoly1kara.inl,v $
// Copyright(c)'1994-2009 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: T. Gautier
// $Id: givpoly1kara.inl,v 1.3 2011-02-02 16:23:56 bboyer Exp $
// ==========================================================================
// Description: To be rewrite using Array0 and subarray0 or generator!!!
#ifndef __GIVARO_poly1_kara_INL
#define __GIVARO_poly1_kara_INL

#define KARA_THRESHOLD 10

template<class T>
void Poly1<T>::stdmul(Array0<T>& R,
                         const Array0<T>& P, const Array0<T>& Q,
                         Degree tP, Degree tQ)
{
   int i,j ;
   for (i=0 ; i <= tP+tQ ; i++) R[i] = T(zero) ;
   for (i=0 ; i <= tP ; i++)
     for (j=0 ; j <= tQ ; j++)
           R[i+j] += P[i]*Q[j] ;
}

  // Tailles requises par les entrees.
  // P  : polynome de degre tP
  // Q  :   "       "       tQ
  // R  :   "       "       tP+tQ
  // Tmp1 : "       "       2*([1+max(tP,tQ)] div 2)
  // Tmp2 : "       "       ([1+max(tP,tQ)] div 2)
  // Au sein de l'algo :
  // P0 :   "       "       k <= tP
  // P1 :   "       "       tP - k <= k
  // Q0 :   "       "       k <= tQ
  // Q1 :   "       "       tQ - k <= k
template<class T>
void Poly1<T>::karatsuba(Array0<T>& P, const Array0<T>& Q,
                Array0<T>& R, Array0<T>& Tmp1,
                Degree tP, Degree tQ)
{
   if ((tP <=KARA_THRESHOLD) || (tQ <=KARA_THRESHOLD))
    {
      Poly1<T>::stdmul(P, Q, R, tP, tQ) ;
      return ;
    }

   Degree k ;
   if (tQ <= tP)
    {
       int i ;
       k = tP/2 ;
       if (tQ <= k)
        {
//cout << "heheh : tP:" << tP << " tQ: " << tQ << endl  ;
          Poly1<T>::karatsuba(P, Q, R, Tmp1, k, tQ) ;
          Poly1<T>::karatsuba(P+k+1, Q, Tmp1, P, tP-k-1, tQ) ;
          for (i=0 ; i <= tP+tQ-k-1 ; i++) R[i+k+1] += Tmp1[i] ;
          return ;
        }
       // Debut des tableaux :
       T* R00 = R ;
       T* R01 = R + k + 1 ;
       T* R0 = R ;
       T* R1 = R + 2*k + 2 ;
       T* P0 = P ;
       T* P1 = P + k+1 ;
       const T* Q0 = Q ;
       const T* Q1 = Q + k+1 ;

       // P1-P0 -> R00
       for (i=0 ; i < tP-k ; i++)     R00[i] = P1[i] - P0[i] ;
       for (i=tP-k ; i <= k ; i++)    R00[i] = -P0[i] ;

       // Q0-Q1 -> R01, degQ0 = k > degQ1 = tQ-k
       for (i=0 ; i < tQ-k ; i++)     R01[i] = Q0[i] - Q1[i] ;
       for (i=tQ-k ; i <= k ; i++)    R01[i] = Q0[i] ;

       // Recursive calls
       _karatsuba(R00, R01, Tmp1, R1, k, k) ;
       _karatsuba(P0,  Q0,  R0,   R1, k, k) ;
       _karatsuba(P1,  Q1,  R1,   P0, tP-k-1, tQ-k-1) ;

       // Merge des termes
       // Correction des termes du milieu:
       for (i=0 ; i <=2*k ; i++) Tmp1[i] += R0[i] ;
       for (i=0 ; i <=tP+tQ-2*k-2 ; i++) Tmp1[i] += R1[i] ;

       // Mis a jour du resultat
       R[2*k+1] = T(zero) ;
       for (i=0 ; i <=2*k ; i++) R[i+k+1] += Tmp1[i] ;
    }
   else
    {
       int i ;
       k = tQ/2 ;
       if (tP <= k)
        {
//cout << "heheh : tP:" << tP << " tQ: " << tQ << endl  ;
          _karatsuba(P, Q, R, Tmp1, tP, k) ;
          _karatsuba(P, Q+k+1, Tmp1, P, tP, tQ-k-1) ;
          for (i=0 ; i <= tP+tQ-k-1 ; i++) R[i+k+1] += Tmp1[i] ;
          return ;
        }
       // Debut des tableaux :
       T* R00 = R ;
       T* R01 = R + k + 1 ;
       T* R0 = R ;
       T* R1 = R + 2*k + 2 ;
       T* P0 = P ;
       T* P1 = P + k+1 ;
       const T* Q0 = Q ;
       const T* Q1 = Q + k+1 ;

       // P1-P0 -> R00
       for (i=0 ; i < tP-k ; i++)     R00[i] = P1[i] - P0[i] ;
       for (i=tP-k ; i <= k ; i++)    R00[i] = -P0[i] ;

       // Q0-Q1 -> R01, degQ0 = k > degQ1 = tQ-k
       for (i=0 ; i < tQ-k ; i++)     R01[i] = Q0[i] - Q1[i] ;
       for (i=tQ-k ; i <= k ; i++)    R01[i] = Q0[i] ;

       // Recursive calls
       _karatsuba(R00, R01, Tmp1, R1, k, k) ;
       _karatsuba(P0,  Q0,  R0,   R1, k, k) ;
       _karatsuba(P1,  Q1,  R1,   P0, tP-k-1, tQ-k-1) ;

       // Merge des termes
       // Correction des termes du milieu:
       for (i=0 ; i <=2*k ; i++) Tmp1[i] += R0[i] ;
       for (i=0 ; i <=tP+tQ-2*k-2 ; i++) Tmp1[i] += R1[i] ;

       // Mis a jour du resultat
       R[2*k+1] = T(zero) ;
       for (i=0 ; i <=2*k ; i++) R[i+k+1] += Tmp1[i] ;
    }
}

template<class T>
const Poly1<T> Karatsuba(const Poly1<T>& P, const Poly1<T>& Q)
{
   Degree degP = P.degree() ;
   Degree degQ = Q.degree() ;

   T* Res = new T[degP+degQ+1];
   T* PP  = new T[degP+1];
   T* Tmp1, *Tmp2 ;
   if (degQ <= degP)
    {
        Tmp1 = new T[degP+1];
//        Tmp2 = new T[degP+1];
    }
   else
    {
        Tmp1 = new T[degQ+1];
//        Tmp2 = new T[degQ+1];
    }
   int i ;
   for (i=0 ; i<= degQ+degP ; i++) Res[i] = T(zero) ;
   for (i=0 ; i<= degP ; i++) PP[i] = P[i] ;
   _karatsuba(PP, Q.baseptr(), Res, Tmp1, degP, degQ) ;
   Poly1<T> PRes(P.variable(), degP+degQ, Res) ;
   delete [] Tmp1 ; delete [] PP ; delete [] Res ;
   return PRes ;
}
#endif // __GIVARO_poly1_kara_INL
