#ifndef _GIV_ELEM_H_
#define _GIV_ELEM_H_
// ==========================================================================
// $Source
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id
// ==========================================================================
// Description:
// definition of a reference to an object.
// 

template<class T>
struct ElemRef {
  typedef T Type_t;
  T& _ref;
  ElemRef( Type_t& ref ) : _ref(ref) {}
  operator Type_t& () { return _ref; }
  operator const Type_t& () const { return _ref; }
  ElemRef<T> operator= (const Type_t& v) { _ref = v; return *this; }
};

template<class T>
struct ElemConstRef {
  typedef T Type_t;
  const T& _ref;
  ElemConstRef( const Type_t& ref ) : _ref(ref) {}
  operator const Type_t& () const { return _ref; }
};

template<class T1, class T2>
struct Pair {
  T1 _val1;
  T2 _val2;
  Pair() {}
  Pair( const T1& v1, const T2& v2) : _val1(v1), _val2(v2) {}
  T1& first() { return _val1; }
  const T1& first() const { return _val1; }
  T2& second() { return _val2; }
  const T2& second() const { return _val2; }
};

template<class T1, class T2>
ostream& operator<< (ostream& o, const Pair<T1,T2>& p )
{ return o << '(' << p._val1 << ',' << p._val2 << ')'; }

template<class T1, class T2>
istream& operator>> (istream& fin, Pair<T1,T2>& p )
{ 
  char ch;
  // Skip the first blanks:
  fin >> std::ws; fin.get(ch);
  if (ch != '(') 
    GivError::throw_error(
      GivBadFormat("operator>><Pair<T1,T2> >: syntax error no '('"));
  
  // - read a T1 + ','
  fin >> std::ws >> p._val1;
  fin >> std::ws; fin.get(ch);
  if (ch != ',') 
    GivError::throw_error(
      GivBadFormat("operator>><Pair<T1,T2> >: syntax error no ','"));

  // - read a T1 + ')'
  fin >> p._val2;
  fin >> std::ws; fin.get(ch);
  if (ch != ')') 
    GivError::throw_error(
      GivBadFormat("operator>><Pair<T1,T2> >: syntax error no ')'"));
  return fin;
}


#endif
