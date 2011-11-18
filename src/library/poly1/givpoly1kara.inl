// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/poly1/givpoly1kara.inl,v $
// Copyright(c)'1994-2011 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: J-G Dumas
// $Id: givpoly1kara.inl,v 1.3 2011-11-08 10:38:00 jgdumas Exp $
// ==========================================================================
#ifndef __GIVARO_poly1_kara_INL
#define __GIVARO_poly1_kara_INL

namespace Givaro {

#ifndef KARA_THRESHOLD
#define KARA_THRESHOLD 50
#endif

#ifndef GIVMIN
#define GIVMIN(a,b) ((a)<(b)?(a):(b))
#endif

// forces standard multiplication
template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::stdmul( Rep& R, const Rep& P, const Rep& Q ) const
{
	size_t sR = R.size();
	size_t sP = P.size();
	size_t sQ = Q.size();
	if ((sQ ==0) || (sP ==0)) { R.reallocate(0); return R; }
	if (sR != sQ+sP) R.reallocate(sR = sP+sQ-1);

 	stdmul(R, R.begin(), R.end(),
            P, P.begin(), P.end(),
            Q, Q.begin(), Q.end());
        
        return setdegree(R);
}

// forces FIRST recursive level with Karatsuba multiplication
template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::karamul( Rep& R, const Rep& P, const Rep& Q ) const
{
	size_t sR = R.size();
	size_t sP = P.size();
	size_t sQ = Q.size();
	if ((sQ ==0) || (sP ==0)) { R.reallocate(0); return R; }
	if (sR != sQ+sP) R.reallocate(sR = sP+sQ-1);

 	karamul(R, R.begin(), R.end(),
            P, P.begin(), P.end(),
            Q, Q.begin(), Q.end());
        
        return setdegree(R);
}

// Generic mul with choices between standard and Karatsuba multiplication
// Multiplies between the iterator bounds.
template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::mul( 
    Rep& R, const RepIterator Rbeg, const RepIterator Rend, 
    const Rep& P, const RepConstIterator Pbeg, const RepConstIterator Pend, 
    const Rep& Q, const RepConstIterator Qbeg, const RepConstIterator Qend ) const {

    if ( ( (Pend-Pbeg)> KARA_THRESHOLD ) &&
         ( (Qend-Qbeg)> KARA_THRESHOLD) )
        return karamul(R, Rbeg, Rend,
                       P, Pbeg, Pend,
                       Q, Qbeg, Qend);
    else 
        return stdmul(R, Rbeg, Rend,
                      P, Pbeg, Pend,
                      Q, Qbeg, Qend);
}




// Standard multiplication between iterator bounds
template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::stdmul( 
    Rep& R, const RepIterator Rbeg, const RepIterator Rend, 
    const Rep& P, const RepConstIterator Pbeg, const RepConstIterator Pend, 
    const Rep& Q, const RepConstIterator Qbeg, const RepConstIterator Qend ) const {

	RepConstIterator ai=Pbeg,bi=Qbeg;
	RepIterator ri=Rbeg, rig=Rbeg;
	if (_domain.isZero(*ai))
		for(;bi!=Qend;++bi,++ri)
			*ri = _domain.zero;
	else
		for(;bi!=Qend;++bi,++ri)
			if (_domain.isZero(*bi))
				*ri = _domain.zero;
			else
				_domain.mul(*ri,*ai,*bi);
        
	for(;ri!=Rend;++ri)
		*ri = _domain.zero;
	for(++ai,++rig;ai!=Pend;++ai,++rig)
		if (! _domain.isZero(*ai))
			for(ri=rig,bi=Qbeg;bi!=Qend;++bi,++ri)
				_domain.axpyin(*ri,*ai,*bi);
	return R;
}

template <class Domain>
inline typename Poly1Dom<Domain,Dense>::Rep& Poly1Dom<Domain,Dense>::karamul( Rep& R, const RepIterator Rbeg, const RepIterator Rend, const Rep& P, const RepConstIterator Pbeg, const RepConstIterator Pend, const Rep& Q, const RepConstIterator Qbeg, const RepConstIterator Qend ) const
{
    // Initialize R to zero
    for(RepIterator ri=Rbeg; ri!= Rend; ++ri) _domain.assign(*ri,_domain.zero);
    

    size_t halfP = (Pend-Pbeg)>>1;
    size_t halfQ = (Qend-Qbeg)>>1;
    size_t half = GIVMIN(halfP, halfQ);
    size_t halfR = half<<1;

    RepConstIterator Pmid=Pbeg+half;		// cut P in halves
    RepConstIterator Qmid=Qbeg+half;		// cut Q in halves
    RepIterator Rmid=Rbeg+halfR;		// cut R in halves
    
    mul(R, Rbeg, Rmid, 				// Recursive dynamic choice
        P, Pbeg, Pmid,
        Q, Qbeg, Qmid);				// PlQl in first storage part of R
    
    mul(R, Rmid, Rend,				// Recursive dynamic choice
        P, Pmid, Pend,
        Q, Qmid, Qend);				// PhQh in second storage part of R
    
    Rep PHPL;
    for(RepConstIterator PHi=Pmid; PHi!=Pend; ++PHi)
        PHPL.push_back(*PHi);
    subin(PHPL, PHPL.begin(), P, Pbeg, Pmid);	// Ph - Pl
    setdegree(PHPL);
    
    Rep QHQL;
    for(RepConstIterator QHi=Qmid; QHi!=Qend; ++QHi)
        QHQL.push_back(*QHi);
    subin(QHQL, QHQL.begin(), Q, Qbeg, Qmid);	// Qh - Ql
    setdegree(QHQL);
    
    Rep M; 
    mul(M, 					// Recursive dynamic choice
        PHPL, 
        QHQL);					// (Ph-Pl)(Qh-Ql)
    setdegree(M);
    
    subin(M, M.begin(), M.end(), R, Rbeg, Rmid);// -= PlQl
    setdegree(M);

    subin(M, M.begin(), M.end(), R, Rmid, Rend);// -= PhQh
    setdegree(M);


    RepIterator ri=Rbeg+half;
    RepConstIterator mi=M.begin();		// update R with mid product
    for( ; mi != M.end(); ++ri, ++mi) _domain.subin(*ri, *mi);

    return R;
}



}
#endif // __GIVARO_poly1_kara_INL
