// Copyright(c)'1994-2011 by The Givaro group 
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Time-stamp: <24 Feb 11 17:20:55 Jean-Guillaume.Dumas@imag.fr>
// ==========================================================================
// Description:
//  Chinese Remainder Algorithm.
#ifndef __GIVARO_arithmodu_fixedprimes_INL
#define __GIVARO_arithmodu_fixedprimes_INL





template<class Ints>
Ints& RNSsystemFixed<Ints>::RnsToRingLeft( Ints& I, const typename RNSsystemFixed<Ints>::array& residues, const int level, const int col ) const {
    if (level) {
        int lowercol=col<<1, lowercolnext=lowercol+1;
        int lowerlevel=level-1;
        Ints u0;
        RnsToRingLeft(u0, residues, lowerlevel, lowercol);
        RnsToRingRight(I, residues, lowerlevel, lowercolnext);
        
        I -= u0; 				// u1-u0
        I *= _primes[lowerlevel][lowercolnext];	// (u1-u0) * (p0(p0^{-1} mod p1))
        I += u0;				// (u1-u0)*M01+u0
        Integer::modin(I, _primes[level][col]);// (u1-u0)*M01 +u0 mod p0p1, between 0 and p0p1-1
        return I;
//         return Integer::modin(I, _primes[level][col]);// (u1-u0)*M01 +u0 mod p0p1, between 0 and p0p1-1
    } else {
        return I=residues[col];
    }
}

template<class Ints>
Ints& RNSsystemFixed<Ints>::RnsToRingRight( Ints& I, const typename RNSsystemFixed<Ints>::array& residues, const int level, const int col ) const {
    if (level) {
        int lowercol=col<<1, lowercolnext=lowercol+1;
        int lowerlevel=level-1;
        Ints u0;
        RnsToRingLeft(u0, residues, lowerlevel, lowercol);
        RnsToRingRight(I, residues, lowerlevel, lowercolnext);
        
        I -= u0; 				// u1-u0
        I *= _primes[lowerlevel][lowercolnext];	// (u1-u0) * (p0(p0^{-1} mod p1))
        I += u0;				// (u1-u0)*M01+u0
        return I;
    } else {
        return I=residues[col];
    }
}



  // Convert to an Ints:
template<class Ints>
Ints& RNSsystemFixed<Ints>::RnsToRing( Ints& I, const RNSsystemFixed<Ints>::array& rns) const
{
    int ir = _RNS.Primes().size();
    typename RNS_t::array Reds( ir );
    for(int i = _primes.size(); (--i)>=0;) {
        int is = _primes[i].size();
        if (is & 1) {
            RnsToRingLeft(Reds[--ir], rns, i, --is);
        }
    }
   
    return _RNS.RnsToRing( I, Reds);
}


// -- free memory allocated in array !
template<class Ints>
RNSsystemFixed<Ints>::~RNSsystemFixed()
{}

template<class Ints>
RNSsystemFixed<Ints>::RNSsystemFixed ()
 : _primes(0)
{}

template<class Ints>
RNSsystemFixed<Ints>::RNSsystemFixed (const RNSsystemFixed<Ints>& R)
 : _primes(R._primes, givWithCopy()),
{}


  // -- Array of primes are given
template<class Ints>
RNSsystemFixed<Ints>::RNSsystemFixed( const RNSsystemFixed<Ints>::array& inprimes)
 : _primes(0)
{
   GIVARO_ASSERT( inprimes.size()>0, "[RNSsystemFixed<Ints>::RNSsystemFixed] bad size of array");
   _primes.resize(1);
   
   for(typename array::const_iterator pit = inprimes.begin(); pit != inprimes.end(); ++pit) {
       _primes.front().push_back( *pit );
       for(int i = 1; i< _primes.size(); ++i) {
           int s = _primes[i-1].size();
           if (s & 1) break;
           else {
               Ints& p0(_primes[i-1][s-2]), & p1(_primes[i-1][s-1]);
               Ints prod(p0*p1);
               inv(p1, p0, p1) *= p0;
               _primes[i].push_back( prod );
           }
       }

       int lastp = _primes.size()-1;

       if (! (_primes.back().size() & 1)) {

//            std::cout << '[';
//            for(int j=0; j<_primes.back().size(); ++j)
//                std::cout << ' ' << _primes.back()[j];
//            std::cout << ']'  << std::endl;
           


           array newlevel;
           int s = _primes.back().size();
           
           Ints& p0(_primes[lastp][s-2]), & p1(_primes[lastp][s-1]);
           newlevel.push_back( p0*p1  );
           
           inv(p1, p0, p1) *= p0;
           _primes.push_back( newlevel );
       }

//        std::cout  << std::endl << "----------------" << std::endl;
//        std::cout << '[';
//        for(int i=0; i<_primes.size(); ++i) {
//            std::cout << '[';
//            for(int j=0; j<_primes[i].size(); ++j)
//                std::cout << ' ' << _primes[i][j];
//            std::cout << ']'  << std::endl;
//        }
//        std::cout << ']' << std::endl;
//        std::cout << "----------------" << std::endl;
       
           

   }

   int numodd=0;
   for(int i = _primes.size(); (--i)>=0; )
       if (_primes[i].size() & 1) ++numodd;
   

   typename RNS_t::domains Mods(numodd);
   for(int i = _primes.size(); (--i)>=0; ) {
       if (_primes[i].size() & 1)
           Mods[--numodd] = _primes[i].back();
   }
   _RNS.setPrimes( Mods );
}


template<class Ints>
const typename RNSsystemFixed<Ints>::tree& RNSsystemFixed<Ints>::Primes() const
{
  return _primes;
}


template<class Ints>
const Ints RNSsystemFixed<Ints>::ith(const size_t i) const
{
  return _primes.front()[i];
}

#endif

