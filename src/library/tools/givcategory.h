#ifndef _GIV_CATEGORY_H_
#define _GIV_CATEGORY_H_
// ==========================================================================
// $Source: /var/lib/cvs/Givaro/src/library/tools/givcategory.h,v $
// Copyright(c)'94-97 by Givaro Team
// see the copyright file.
// Authors: T. Gautier
// $Id: givcategory.h,v 1.1.1.1 2004-05-12 16:08:24 jgdumas Exp $
// ==========================================================================
// Description:
// - Definition of traits.
// It's a beta-beta version


// -- nothing
class Undefined{}; 


// ==========================================================================
// --
// -- Test is A is equal to B, A and B are two types.
// -- ISQUAL(A,B) or IsEqual<A,B>::val return true iff
// -- A and B is the same type.
// --
template<class A, class B>
struct IsEqual { enum {val = false }; };

template<class A>
struct IsEqual<A,A> { enum {val = true }; };

template<class A, class B>
struct IsNotEqual { enum {val = !IsEqual<A,B>::val }; };

#define ISEQUAL(A,B) IsEqual<A,B>::val
#define ISNOTEQUAL(A,B) IsNotEqual<A,B>::val



// ==========================================================================
// --
// -- Characteristic for representation of vector and matrix
// --
// Sporadic is a Dense but may have quite a few zero elements
// Therefore some algorithms might be specialized and might 
// take advantage of this.
class Sporadic{}; 
class Dense : public Sporadic {}; 
class Sparse{}; 
class Diagonal{}; typedef Diagonal Diag; 
class Toeplitz{}; 
class Hensel{}; 

template<class CLASS>
struct Sparsity_Trait {
  typedef Undefined  Sparsity_t;  // To be defined in specialized class
};



// ==========================================================================
// --
// -- Envelop class to specialize format output
// --
class DefaultFormat {};
template<class T, class Tag = DefaultFormat>
struct StructFormat {
  const T& value;
  StructFormat(const T& val) : value(val) {}
  operator const T& () { return value; }
};

template<class T, class Tag>
StructFormat<T,Tag> Formatted( const T& val, Tag xx )
{ return StructFormat<T,Tag>(val); }

#endif
