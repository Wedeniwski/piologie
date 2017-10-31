/////////////////////////////////
//
// Piologie V 1.3.3
// multi-precision arithmetic
// Integer
//
// (c) 1996-2001 HiPiLib
//               www.hipilib.de
//
// Sebastian Wedeniwski
// 10/01/2011
//

#ifndef _Include_Integer_H_
#define _Include_Integer_H_

#include "natural.h"

struct Integer_plus_tag {};
struct Integer_negate_tag {};
struct Integer_minus_tag {};
struct Integer_multiplies_tag {};
struct Integer_divides_tag {};
struct Integer_modulus_tag {};
struct Integer_square_root_tag {};
struct Integer_lshift_tag {};
struct Integer_rshift_tag {};
struct Integer_and_tag {};
struct Integer_or_tag {};
struct Integer_xor_tag {};
struct Integer_not_tag {};

class Integer : public Natural {
public:
  struct rep {
    const int   sgn;
    const       Natural& nat;
    const bool  bin;
    rep(const int a, const Natural& b, const bool c = false)
      : sgn(a), nat(b), bin(c) {}
  };
  
private:
  int sgn;

  Integer(const size_t, char);

  Natural&  abs();
  void      neg();
  void      neg(const Integer&);
  void      add(const Integer&, const Integer&);
  void      sub(const Integer&, const Integer&);
  void      mul(const Integer&, const Integer&);
  void      div(const Integer&, const Integer&);
  void      sqrt(const Integer&);
  
  void      lshift(const Integer&, const size_t);
  void      rshift(const Integer&, const size_t);
  void      bitwise_and(const Integer&, const Integer&);
  void      bitwise_or(const Integer&, const Integer&);
  void      bitwise_xor(const Integer&, const Integer&);
  void      bitwise_not(const Integer&);
  
  void      add(const Integer&, const Natural&);
  void      sub(const Integer&, const Natural&);
  void      sub(const Natural&, const Integer&);
  void      mul(const Integer&, const Natural&);
  void      div(const Integer&, const Natural&);
  void      div(const Natural&, const Integer&);
  
  void      add(const Integer&, const SignDigit);
  void      sub(const Integer&, const SignDigit);
  void      mul(const Integer&, const SignDigit);
  void      muladd(const Integer&, const SignDigit);
  void      mulsub(const Integer&, const SignDigit);
  

public:
  Integer(const SignDigit = 0);
  Integer(const Natural&);
  Integer(const Integer&);
#ifndef _Old_STD_
explicit
#endif
  Integer(const char*, const Digit = 10);
  ~Integer();

  inline Integer&   operator=(const Integer&);
  inline Integer&   operator=(const Natural&);
  Integer&          operator+=(const Integer&);
  Integer&          operator-=(const Integer&);
  Integer&          operator*=(const Integer&);
  Integer&          operator/=(const Integer&);
  Integer&          operator%=(const Integer&);
  Integer&          operator&=(const Integer&);
  Integer&          operator|=(const Integer&);
  Integer&          operator^=(const Integer&);
    
  Integer&          operator+=(const Natural&);
  Integer&          operator-=(const Natural&);
  Integer&          operator*=(const Natural&);
  Integer&          operator/=(const Natural&);
  Integer&          operator%=(const Natural&);
    
  Integer&          operator+=(const SignDigit);
  Integer&          operator-=(const SignDigit);
  Integer&          operator*=(const SignDigit);
  Integer&          operator/=(const SignDigit);
  Integer&          operator%=(const SignDigit);
  Integer&          operator>>=(const size_t);
  Integer&          operator<<=(const size_t);
  inline SignDigit  operator=(const SignDigit);
  Digit             operator&=(const Digit);
  Integer&          operator|=(const Digit);
    
  const Integer&    operator++();
  const Integer&    operator--();
  Integer           operator++(int);
  Integer           operator--(int);
  
  
  Integer(const binder_arguments<Integer, Integer, Integer_negate_tag>&);
  Integer& operator=(const binder_arguments<Integer, Integer, Integer_negate_tag>&);
  
  Integer(const binder_arguments<Integer, Integer, Integer_plus_tag>&);
  Integer& operator=(const binder_arguments<Integer, Integer, Integer_plus_tag>&);
  
  Integer(const binder_arguments<Integer, Integer, Integer_minus_tag>&);
  Integer& operator=(const binder_arguments<Integer, Integer, Integer_minus_tag>&);
  
  Integer(const binder_arguments<Integer, Integer, Integer_multiplies_tag>&);
  Integer& operator=(const binder_arguments<Integer, Integer, Integer_multiplies_tag>&);
  
  Integer(const binder_arguments<Integer, Integer, Integer_divides_tag>&);
  Integer& operator=(const binder_arguments<Integer, Integer, Integer_divides_tag>&);
  
  Integer(const binder_arguments<Integer, Integer, Integer_modulus_tag>&);
  Integer& operator=(const binder_arguments<Integer, Integer, Integer_modulus_tag>&);
  
  
  Integer(const binder_arguments<Integer, size_t, Integer_lshift_tag>&);
  Integer& operator=(const binder_arguments<Integer, size_t, Integer_lshift_tag>&);
  
  Integer(const binder_arguments<Integer, size_t, Integer_rshift_tag>&);
  Integer& operator=(const binder_arguments<Integer, size_t, Integer_rshift_tag>&);
  
  Integer(const binder_arguments<Integer, Integer, Integer_and_tag>&);
  Integer& operator=(const binder_arguments<Integer, Integer, Integer_and_tag>&);
  
  Integer(const binder_arguments<Integer, Integer, Integer_or_tag>&);
  Integer& operator=(const binder_arguments<Integer, Integer, Integer_or_tag>&);
  
  Integer(const binder_arguments<Integer, Integer, Integer_xor_tag>&);
  Integer& operator=(const binder_arguments<Integer, Integer, Integer_xor_tag>&);
  
  Integer(const binder_arguments<Integer, Integer, Integer_not_tag>&);
  Integer& operator=(const binder_arguments<Integer, Integer, Integer_not_tag>&);
  
  
  Integer(const binder_arguments<Integer, Natural, Integer_plus_tag>&);
  Integer& operator=(const binder_arguments<Integer, Natural, Integer_plus_tag>&);
  
  Integer(const binder_arguments<Integer, Natural, Integer_minus_tag>&);
  Integer& operator=(const binder_arguments<Integer, Natural, Integer_minus_tag>&);
  
  Integer(const binder_arguments<Natural, Integer, Integer_minus_tag>&);
  Integer& operator=(const binder_arguments<Natural, Integer, Integer_minus_tag>&);
  
  Integer(const binder_arguments<Integer, Natural, Integer_multiplies_tag>&);
  Integer& operator=(const binder_arguments<Integer, Natural, Integer_multiplies_tag>&);
  
  Integer(const binder_arguments<Integer, Natural, Integer_divides_tag>&);
  Integer& operator=(const binder_arguments<Integer, Natural, Integer_divides_tag>&);
  
  Integer(const binder_arguments<Natural, Integer, Integer_divides_tag>&);
  Integer& operator=(const binder_arguments<Natural, Integer, Integer_divides_tag>&);
  
  
  Integer(const binder_arguments<Integer, SignDigit, Integer_minus_tag>&);
  Integer& operator=(const binder_arguments<Integer, SignDigit, Integer_minus_tag>&);
  Integer(const binder_arguments<SignDigit, Integer, Integer_minus_tag>&);
  Integer& operator=(const binder_arguments<SignDigit, Integer, Integer_minus_tag>&);
    
  Integer(const binder_arguments<Integer, SignDigit, Integer_plus_tag>&);
  Integer& operator=(const binder_arguments<Integer, SignDigit, Integer_plus_tag>&);
    
  Integer(const binder_arguments<Integer, SignDigit, Integer_multiplies_tag>&);
  Integer& operator=(const binder_arguments<Integer, SignDigit, Integer_multiplies_tag>&);
  Integer& operator+=(const binder_arguments<Integer, SignDigit, Integer_multiplies_tag>&);
  Integer& operator-=(const binder_arguments<Integer, SignDigit, Integer_multiplies_tag>&);
  
  Integer(const binder_arguments<Integer, Integer, Integer_square_root_tag>&);
  Integer& operator=(const binder_arguments<Integer, Integer, Integer_square_root_tag>&);
  
  friend void       div(const Integer&, const Integer&, Integer&, Integer&);
  friend void       div(const Integer&, const SignDigit, Integer&, SignDigit&);
  friend void       swap(Integer&, Integer&);
  friend void       sqrt(const Integer&, Integer&, Integer&);
  inline friend int sign(const Integer&);
  friend Integer    pow(const Integer&, const SignDigit);
  friend Integer    pow(const Integer&, const Integer&);
  friend ISTREAM&   operator>>(ISTREAM&, Integer&);
  
  
  void        split(const size_t, Integer&, Integer&) const;
  void        setbit(const size_t);
  void        clearbit(const size_t);
  bool        testbit(const size_t) const;
  void        rand(const size_t);
  const char* atoI(const char*, const Digit = 10);
  bool        scan(ISTREAM&);
  
};


inline binder_arguments<Integer, Integer, Integer_negate_tag>
 operator-(const Integer&);
inline binder_arguments<Integer, Integer, Integer_plus_tag>
 operator+(const Integer&, const Integer&);
inline binder_arguments<Integer, Integer, Integer_minus_tag>
 operator-(const Integer&, const Integer&);
inline binder_arguments<Integer, Integer, Integer_multiplies_tag>
 operator*(const Integer&, const Integer&);
inline binder_arguments<Integer, Integer, Integer_divides_tag>
 operator/(const Integer&, const Integer&);
inline binder_arguments<Integer, Integer, Integer_modulus_tag>
 operator%(const Integer&, const Integer&);
inline binder_arguments<Integer, size_t, Integer_lshift_tag>
 operator<<(const Integer&, const size_t&);
inline binder_arguments<Integer, size_t, Integer_rshift_tag>
 operator>>(const Integer&, const size_t&);
inline binder_arguments<Integer, Integer, Integer_and_tag>
 operator&(const Integer&, const Integer&);
inline binder_arguments<Integer, Integer, Integer_or_tag>
 operator|(const Integer&, const Integer&);
inline binder_arguments<Integer, Integer, Integer_xor_tag>
 operator^(const Integer&, const Integer&);
inline binder_arguments<Integer, Integer, Integer_not_tag>
 operator~(const Integer&);
inline binder_arguments<Integer, Natural, Integer_plus_tag>
 operator+(const Natural&, const Integer&);
inline binder_arguments<Integer, Natural, Integer_plus_tag>
 operator+(const Integer&, const Natural&);
inline binder_arguments<Natural, Integer, Integer_minus_tag>
 operator-(const Natural&, const Integer&);
inline binder_arguments<Integer, Natural, Integer_minus_tag>
 operator-(const Integer&, const Natural&);
inline binder_arguments<Integer, Natural, Integer_multiplies_tag>
 operator*(const Natural&, const Integer&);
inline binder_arguments<Integer, Natural, Integer_multiplies_tag>
 operator*(const Integer&, const Natural&);
inline binder_arguments<Natural, Integer, Integer_divides_tag>
 operator/(const Natural&, const Integer&);
inline binder_arguments<Integer, Natural, Integer_divides_tag>
 operator/(const Integer&, const Natural&);
inline binder_arguments<Integer, SignDigit, Integer_plus_tag>
 operator+(const SignDigit&, const Integer&);
inline binder_arguments<Integer, SignDigit, Integer_plus_tag>
 operator+(const Integer&, const SignDigit&);
inline binder_arguments<Integer, SignDigit, Integer_multiplies_tag>
 operator*(const SignDigit&, const Integer&);
inline binder_arguments<Integer, SignDigit, Integer_multiplies_tag>
 operator*(const Integer&, const SignDigit&);

inline Integer    operator/(const Integer&, const SignDigit);
inline SignDigit  operator%(const Integer&, const SignDigit);
inline Digit      operator&(const Integer&, const Digit);
inline Integer    operator|(const Integer&, const Digit);


inline bool       operator==(const Integer&, const Integer&);
inline bool       operator!=(const Integer&, const Integer&);
inline bool       operator<(const Integer&, const Integer&);
inline bool       operator<=(const Integer&, const Integer&);
inline bool       operator>(const Integer&, const Integer&);
inline bool       operator>=(const Integer&, const Integer&);

inline bool       operator==(const Integer&, const SignDigit);
inline bool       operator!=(const Integer&, const SignDigit);
inline bool       operator<(const Integer&, const SignDigit);
inline bool       operator<=(const Integer&, const SignDigit);
inline bool       operator>(const Integer&, const SignDigit);
inline bool       operator>=(const Integer&, const SignDigit);

inline binder_arguments<Integer, Integer, Integer_square_root_tag>
 sqrt(const Integer&);


inline const Natural& abs(const Integer&);
inline Digit          log2(const Integer&);
Integer               root(const Integer&, const SignDigit);
inline int            units(Integer&);
void                  gcd(Integer, Integer, Integer&, Integer&, Integer&);

inline Integer        atoI(const char*, const Digit = 10);
char*                 Itoa(const Integer&, char*, const Digit = 10);

inline Integer::rep   print(const Integer&, bool = false);
inline OSTREAM&       operator<<(OSTREAM&, const Integer::rep&);
inline OSTREAM&       operator<<(OSTREAM&, const Integer&);

void                  div(const Integer&, const Integer&, Integer&, Integer&);
void                  div(const Integer&, const SignDigit, Integer&, SignDigit&);
void                  swap(Integer&, Integer&);
void                  sqrt(const Integer&, Integer&, Integer&);
inline int            sign(const Integer&);
Integer               pow(const Integer&, const SignDigit);
Integer               pow(const Integer&, const Integer&);
ISTREAM&              operator>>(ISTREAM&, Integer&);

///////////////////////// Inline-Implementation ////////////////////

inline const Natural& abs(const Integer& a)
// Algorithm:  c := abs(a)
// Input:      a in Integer.
// Output:     c in Natural such that c = |a| ||
{
  return (const Natural&)a;
}

inline Natural& Integer::abs()
// Algorithm:  c := a.abs()
// Input:      a in Integer.
// Output:     c in Natural such that c = |a| ||
{
  return (Natural&)*this;
}

inline int sign(const Integer& a)
// Algorithm:  c := sign(a)
// Input:      a in Integer.
// Output:     c in int such that if a = 0 then c = 0
//             else if a > 0 then c = 1 else c = -1 ||
{
  return a.sgn;
}

inline void Integer::neg()
// Algorithm:  a.neg()
// Input:      a in Integer.
// Output:     a in Integer such that a := -a ||
{
  sgn = -sgn;
}
inline void Integer::neg(const Integer& a)
// Algorithm:  b.neg(a)
// Input:      a in Integer.
// Output:     b in Integer such that b = -a ||
{
  *this = a; neg();
}

inline Integer::Integer(const size_t a, char b)
 : Natural(a, b)
// Algorithm:  c := Integer(a, b)
// Input:      b in char, a in size_t where a >= 1.
// Output:     c in Integer such that L(c) = R(c) = a;
//             map Integer(a, b) to Natural(a, b) ||
//
// Note:       This constructor don't fulfill the conditions for Integers.
//
// internal constructor without the initialization of the elements.
{
}

inline Integer::Integer(const SignDigit a)
 : Natural((a < 0)? -a : a)
// Algorithm:  c := Integer(a)
// Input:      a in SignDigit.
// Output:     c in Integer such that c = a ||
{
  if (a < 0) sgn = -1;
  else if (a > 0) sgn = 1;
  else sgn = 0;
}

inline Integer::Integer(const Natural& a)
 : Natural(a)
// Algorithm:  c := Integer(a)
// Input:      a in Natural.
// Output:     c in Integer such that c = a ||
{
  sgn = (a != 0);
}

inline Integer::Integer(const Integer& a)
: Natural(::abs(a))
// Algorithm:  c := Integer(a)
// Input:      a in Integer.
// Output:     c in Integer such that c = a ||
{
  sgn = a.sgn;
}

inline Integer::Integer(const char* a, const Digit b)
 : Natural(a + (*a == '-'), b)
// Algorithm:  c := Integer(a, b)
// Input:      a in String, b in Digit.
// Output:     c in Integer such that c = a ||
{
  if (*a == '-') sgn = -1;
  else sgn = (::abs(*this) != 0);
}

inline Integer::~Integer()
{
}

inline bool operator==(const Integer& a, const SignDigit b)
// Algorithm:  c := a == b
// Input:      a in Integer, b in SignDigit.
// Output:     c in bool such that if a = b then c = true else c = false ||
{
  return (b < 0)? (sign(a) < 0 && abs(a) == Digit(-b))
                : (sign(a) >= 0 && abs(a) == Digit(b));
}

inline bool operator!=(const Integer& a, const SignDigit b)
// Algorithm:  c := a != b
// Input:      a in Integer, b in SignDigit.
// Output:     c in bool such that if not a = b then c = true else c = false ||
{
  return (b < 0)? (sign(a) >= 0 || abs(a) != Digit(-b))
                : (sign(a) < 0 || abs(a) != Digit(b));
}

inline bool operator<(const Integer& a, const SignDigit b)
// Algorithm:  c := a < b
// Input:      a in Integer, b in SignDigit.
// Output:     c in bool such that if a < b then c = true else c = false ||
{
  return (b < 0)? (sign(a) < 0 && abs(a) > Digit(-b))
                : (sign(a) < 0 || abs(a) < Digit(b));
}

inline bool operator<=(const Integer& a, const SignDigit b)
// Algorithm:  c := a <= b
// Input:      a in Integer, b in SignDigit.
// Output:     c in bool such that if a <= b then c = true else c = false ||
{
  return (b < 0)? (sign(a) < 0 && abs(a) >= Digit(-b))
                : (sign(a) <= 0 || abs(a) <= Digit(b));
}

inline bool operator>(const Integer& a, const SignDigit b)
// Algorithm:  c := a > b
// Input:      a in Integer, b in SignDigit.
// Output:     c in bool such that if a > b then c = true else c = false ||
{
  return (b <= 0)? (sign(a) > 0 || abs(a) < Digit(-b))
                 : (sign(a) > 0 && abs(a) > Digit(b));
}

inline bool operator>=(const Integer& a, const SignDigit b)
// Algorithm:  c := a >= b
// Input:      a in Integer, b in SignDigit.
// Output:     c in bool such that if a >= b then c = true else c = false ||
{
  return (b <= 0)? (sign(a) >= 0 || abs(a) <= Digit(-b))
                 : (sign(a) > 0 && abs(a) >= Digit(b));
}

inline Integer& Integer::operator+=(const Integer& a)
// Algorithm:  c := c += a
// Input:      a,c in Integer.
// Output:     c in Integer such that c := c+a ||
{
  add(a, *this);
  return *this;
}

inline Integer& Integer::operator-=(const Integer& a)
// Algorithm:  c := c -= a
// Input:      a,c in Integer.
// Output:     c in Integer such that c := c-a ||
{
  sub(*this, a);
  return *this;
}

inline Integer& Integer::operator*=(const Integer& a)
// Algorithm:  c := c *= a
// Input:      a,c in Integer.
// Output:     c in Integer such that c := c*a ||
{
  sgn *= a.sgn;
  abs() *= ::abs(a);
  return *this;
}

inline Integer& Integer::operator/=(const Integer& a)
// Algorithm:  c := c /= a
// Input:      a,c in Integer where not a = 0.
// Output:     c in Integer such that c := sign(a*c)*[|c/a|] ||
{
  abs() /= ::abs(a);
  sgn = (::abs(*this) != 0)? sgn*a.sgn : 0;
  return *this;
}

inline Integer& Integer::operator%=(const Integer& a)
// Algorithm:  c := c %= a
// Input:      a,c in Integer where not a = 0.
// Output:     c in Integer such that c := c - sign(a*c)*[|c/a|]*a ||
{
  abs() %= ::abs(a);
  if (::abs(*this) == 0) sgn = 0;
  return *this;
}

inline Integer& Integer::operator&=(const Integer& a)
// Algorithm:  c := c &= a
// Input:      a,c in Integer.
// Output:     c in Integer such that c := c and a ||
{
  bitwise_and(*this, a);
  return *this;
}

inline Integer& Integer::operator|=(const Integer& a)
// Algorithm:  c := c |= a
// Input:      a,c in Integer.
// Output:     c in Integer such that c := c or a ||
{
  bitwise_or(*this, a);
  return *this;
}

inline Integer& Integer::operator^=(const Integer& a)
// Algorithm:  c := c ^= a
// Input:      a,c in Integer.
// Output:     c in Integer such that c := c xor a ||
{
  bitwise_xor(*this, a);
  return *this;
}

inline Integer& Integer::operator+=(const Natural& a)
// Algorithm:  c := c += a
// Input:      a in Natural, c in Integer.
// Output:     c in Integer such that c := c+a ||
{
  add(*this, a);
  return *this;
}

inline Integer& Integer::operator-=(const Natural& a)
// Algorithm:  c := c -= a
// Input:      a in Natural, c in Integer.
// Output:     c in Integer such that c := c-a ||
{
  sub(*this, a);
  return *this;
}

inline Integer& Integer::operator*=(const Natural& a)
// Algorithm:  c := c *= a
// Input:      a in Natural, c in Integer.
// Output:     c in Integer such that c := c*a ||
{
  if (a != 0) abs() *= a;
  else *this = 0;
  return *this;
}

inline Integer& Integer::operator/=(const Natural& a)
// Algorithm:  c := c /= a
// Input:      c in Integer, a in Natural where not a = 0.
// Output:     c in Integer such that c := sign(c)*[|c|/a] ||
{
  abs() /= a;
  if (::abs(*this) == 0) sgn = 0;
  return *this;
}

inline Integer& Integer::operator%=(const Natural& a)
// Algorithm:  c := c %= a
// Input:      c in Integer, a in Natural where not a = 0.
// Output:     c in Integer such that c := c - sign(c)*[|c|/a]*a ||
{
  abs() %= a;
  if (::abs(*this) == 0) sgn = 0;
  return *this;
}

inline void Integer::split(const size_t b, Integer& c, Integer& d) const
// Algorithm:  a.split(b, c, d)
// Input:      a,c,d in Integer, b in size_t where not c = d.
// Output:     c,d in Integer such that c = sign(a)*[|a|/2^(BETA*b)],
//             d = a - c*2^(BETA b) ||
{
  Natural::split(b, c.abs(), d.abs());
  c.sgn = (::abs(c) != 0)? sgn : 0;
  d.sgn = (::abs(d) != 0)? sgn : 0;
}

inline bool operator==(const Integer& a, const Integer& b)
// Algorithm:  c := a == b
// Input:      a,b in Integer.
// Output:     c in bool such that if a = b then c = true else c = false ||
{
  return (sign(a) == sign(b) && abs(a) == abs(b));
}

inline bool operator!=(const Integer& a, const Integer& b)
// Algorithm:  c := a != b
// Input:      a,b in Integer.
// Output:     c in bool such that if a = b then c = false else c = true ||
{
  return (sign(a) != sign(b) || abs(a) != abs(b));
}

inline bool operator<(const Integer& a, const Integer& b)
// Algorithm:  c := a < b
// Input:      a,b in Integer.
// Output:     c in bool such that if a < b then c = true else c = false ||
{
  return (sign(a) < sign(b) || sign(a) == sign(b)
   && ((sign(a) > 0)? abs(a) < abs(b) : abs(b) < abs(a)));
}

inline bool operator<=(const Integer& a, const Integer& b)
// Algorithm:  c := a <= b
// Input:      a,b in Integer.
// Output:     c in bool such that if a <= b then c = true else c = false ||
{
  return !(b < a);
}

inline bool operator>(const Integer& a, const Integer& b)
// Algorithm:  c := a > b
// Input:      a,b in Integer.
// Output:     c in bool such that if a > b then c = true else c = false ||
{
  return (b < a);
}

inline bool operator>=(const Integer& a, const Integer& b)
// Algorithm:  c := a >= b
// Input:      a,b in Integer.
// Output:     c in bool such that if a >= b then c = true else c = false ||
{
  return !(a < b);
}

inline Integer::Integer(const binder_arguments<Integer, Integer,
                                               Integer_negate_tag>& a)
 : Natural(' ', ' ')
{
  get_memory(a.x.length()+DELTA);
  neg(a.x);
}

inline Integer& Integer::operator=(const binder_arguments<Integer, Integer,
                                                          Integer_negate_tag>& a)
{
  if (this == &a.x) neg();
  else neg(a.x);
  return *this;
}

inline binder_arguments<Integer, Integer, Integer_negate_tag>
 operator-(const Integer& a)
// Algorithm:  c := -a
// Input:      a in Integer.
// Output:     c in Integer such that c = -a ||
{
  return binder_arguments<Integer, Integer, Integer_negate_tag>(a, a);
}

inline Integer::Integer(const binder_arguments<Integer, Integer,
                                               Integer_plus_tag>& a)
 : Natural(' ', ' ')
{
  get_memory(::max(a.x.length(), a.y.length())+DELTA);
  add(a.x, a.y);
}

inline Integer& Integer::operator=(const binder_arguments<Integer, Integer,
                                                          Integer_plus_tag>& a)
{
  add(a.x, a.y);
  return *this;
}

inline binder_arguments<Integer, Integer, Integer_plus_tag>
 operator+(const Integer& a, const Integer& b)
// Algorithm:  c := a+b
// Input:      a,b in Integer.
// Output:     c in Integer such that c = a+b ||
{
  return binder_arguments<Integer, Integer, Integer_plus_tag>(a, b);
}

inline Integer::Integer(const binder_arguments<Integer, Integer,
                                               Integer_minus_tag>& a)
 : Natural(' ', ' ')
{
  get_memory(::max(a.x.length(), a.y.length())+DELTA);
  sub(a.x, a.y);
}

inline Integer& Integer::operator=(const binder_arguments<Integer, Integer,
                                                          Integer_minus_tag>& a)
{
  sub(a.x, a.y);
  return *this;
}

inline binder_arguments<Integer, Integer, Integer_minus_tag>
 operator-(const Integer& a, const Integer& b)
// Algorithm:  c := a-b
// Input:      a,b in Integer.
// Output:     c in Integer such that c = a-b ||
{
  return binder_arguments<Integer, Integer, Integer_minus_tag>(a, b);
}

inline Integer::Integer(const binder_arguments<Integer, Natural,
                                               Integer_plus_tag>& a)
 : Natural(' ', ' ')
{
  get_memory(::max(a.x.length(), a.y.length())+DELTA);
  add(a.x, a.y);
}

inline Integer& Integer::operator=(const binder_arguments<Integer, Natural,
                                                          Integer_plus_tag>& a)
{
  add(a.x, a.y);
  return *this;
}

inline binder_arguments<Integer, Natural, Integer_plus_tag>
 operator+(const Integer& a, const Natural& b)
// Algorithm:  c := a+b
// Input:      a in Integer, b in Natural.
// Output:     c in Integer such that c = a+b ||
{
  return binder_arguments<Integer, Natural, Integer_plus_tag>(a, b);
}

inline binder_arguments<Integer, Natural, Integer_plus_tag>
 operator+(const Natural& a, const Integer& b)
// Algorithm:  c := a+b
// Input:      a in Natural, b in Integer.
// Output:     c in Integer such that c = a+b ||
{
  return binder_arguments<Integer, Natural, Integer_plus_tag>(b, a);
}

inline Integer::Integer(const binder_arguments<Integer, Natural,
                                               Integer_minus_tag>& a)
 : Natural(' ', ' ')
{
  get_memory(::max(a.x.length(), a.y.length())+DELTA);
  sub(a.x, a.y);
}

inline Integer& Integer::operator=(const binder_arguments<Integer, Natural,
                                                          Integer_minus_tag>& a)
{
  sub(a.x, a.y);
  return *this;
}

inline Integer::Integer(const binder_arguments<Natural, Integer,
                                               Integer_minus_tag>& a)
 : Natural(' ', ' ')
{
  get_memory(::max(a.x.length(), a.y.length())+DELTA);
  sub(a.x, a.y);
}

inline Integer& Integer::operator=(const binder_arguments<Natural, Integer,
                                                          Integer_minus_tag>& a)
{
  sub(a.x, a.y);
  return *this;
}

inline binder_arguments<Integer, Natural, Integer_minus_tag>
 operator-(const Integer& a, const Natural& b)
// Algorithm:  c := a-b
// Input:      a in Integer, b in Natural.
// Output:     c in Integer such that c = a-b ||
{
  return binder_arguments<Integer, Natural, Integer_minus_tag>(a, b);
}

inline binder_arguments<Natural, Integer, Integer_minus_tag>
 operator-(const Natural& a, const Integer& b)
// Algorithm:  c := a-b
// Input:      a in Natural, b in Integer.
// Output:     c in Integer such that c = a-b ||
{
  return binder_arguments<Natural, Integer, Integer_minus_tag>(a, b);
}

inline void Integer::mul(const Integer& a, const SignDigit b)
// Algorithm:  c.mul(a, b)
// Input:      a in Integer, b in SignDigit.
// Output:     c in Integer such that c = a*b ||
{
  if (b > 0) {
    abs() = ::abs(a) * Digit(b);
    sgn = a.sgn;
  } else if (b < 0) {
    abs() = ::abs(a) * Digit(-b);
    sgn = -a.sgn;
  } else *this = 0;
}

inline void Integer::mul(const Integer& a, const Natural& b)
// Algorithm:  c.mul(a, b)
// Input:      a in Integer, b in Natural.
// Output:     c in Integer such that c = a*b ||
{
  abs() = ::abs(a) * b;
  sgn = a.sgn;
}

inline void Integer::mul(const Integer& a, const Integer& b)
// Algorithm:  c.mul(a, b)
// Input:      a,b in Integer.
// Output:     c in Integer such that c = a*b ||
{
  abs() = ::abs(a) * ::abs(b);
  sgn = a.sgn*b.sgn;
}

inline Integer::Integer(const binder_arguments<Integer, Integer,
                                               Integer_multiplies_tag>& a)
 : Natural(' ', ' ')
{
  get_memory(a.x.length()+a.y.length()+DELTA);
  mul(a.x, a.y);
}

inline Integer& Integer::operator=(const binder_arguments<Integer, Integer,
                                                          Integer_multiplies_tag>& a)
{
  mul(a.x, a.y);
  return *this;
}

inline binder_arguments<Integer, Integer, Integer_multiplies_tag>
 operator*(const Integer& a, const Integer& b)
// Algorithm:  c := a*b
// Input:      a,b in Integer.
// Output:     c in Integer such that c = a*b ||
{
  return binder_arguments<Integer, Integer, Integer_multiplies_tag>(a, b);
}

inline Integer::Integer(const binder_arguments<Integer, Natural,
                                               Integer_multiplies_tag>& a)
 : Natural(' ', ' ')
{
  get_memory(a.x.length()+a.y.length()+DELTA);
  mul(a.x, a.y);
}

inline Integer& Integer::operator=(const binder_arguments<Integer, Natural,
                                                          Integer_multiplies_tag>& a)
{
  mul(a.x, a.y);
  return *this;
}

inline binder_arguments<Integer, Natural, Integer_multiplies_tag>
 operator*(const Integer& a, const Natural& b)
// Algorithm:  c := a*b
// Input:      a in Integer, b in Natural.
// Output:     c in Integer such that c = a*b ||
{
  return binder_arguments<Integer, Natural, Integer_multiplies_tag>(a, b);
}

inline binder_arguments<Integer, Natural, Integer_multiplies_tag>
 operator*(const Natural& a, const Integer& b)
// Algorithm:  c := a*b
// Input:      a in Natural, b in Integer.
// Output:     c in Integer such that c = a*b ||
{
  return binder_arguments<Integer, Natural, Integer_multiplies_tag>(b, a);
}

inline void Integer::bitwise_not(const Integer& a)
// Algorithm:  b.not(a)
// Input:      a in Integer.
// Output:     b in Integer such that b = not a = -a-1 ||
{
  *this = a; neg(); --*this;
}

inline void Integer::lshift(const Integer& a, const size_t b)
// Algorithm:  c.lshift(a, b)
// Input:      a in Integer, b in size_t.
// Output:     c in Integer such that c = a*2^b ||
{
  abs() = ::abs(a) << b;
  sgn = a.sgn;
}

inline void Integer::rshift(const Integer& a, const size_t b)
// Algorithm:  c.rshift(a, b)
// Input:      a in Integer, b in size_t.
// Output:     c in Integer such that c = sign(a)*[|a|/2^b] ||
{
  abs() = ::abs(a) >> b;
  if (::abs(*this) == 0) sgn = 0;
  else sgn = a.sgn;
}

inline Integer& Integer::operator*=(const SignDigit a)
// Algorithm:  c := c *= a
// Input:      a in SignDigit, c in Integer.
// Output:     c in Integer such that c := c*a ||
{
  if (a > 0) abs() *= a;
  else if (a < 0) { abs() *= -a; sgn = -sgn; }
  else *this = 0;
  return *this;
}

inline Integer& Integer::operator/=(const SignDigit a)
// Algorithm:  c := c /= a
// Input:      c in Integer, a in SignDigit where not a = 0.
// Output:     c in Integer such that c := sign(a*c)*[|c/a|] ||
{
  if (a >= 0)
    if (::abs(*this) < Digit(a)) *this = 0;
    else abs() /= a;
  else
    if (::abs(*this) < Digit(-a)) *this = 0;
    else { abs() /= -a; sgn = -sgn; }
  return *this;
}

inline Integer& Integer::operator%=(const SignDigit a)
// Algorithm:  c := c %= a
// Input:      c in Integer, a in SignDigit where not a = 0.
// Output:     c in Integer such that c := c - sign(a*c)*[|c/a|]*a ||
{
  abs() %= (a >= 0)? a : -a;
  if (::abs(*this) == 0) sgn = 0;
  return *this;
}

inline Integer& Integer::operator>>=(const size_t a)
// Algorithm:  c := c >>= a
// Input:      a in size_t, c in Integer.
// Output:     c in Integer such that c := sign(c)*[|c|/2^a] ||
{
  abs() >>= a;
  if (::abs(*this) == 0) sgn = 0;
  return *this;
}

inline Integer& Integer::operator<<=(const size_t a)
// Algorithm:  c := c <<= a
// Input:      a in size_t, c in Integer.
// Output:     c in Integer such that c := c*2^a ||
{
  abs() <<= a;
  return *this;
}

inline Integer::Integer(const binder_arguments<Integer, size_t,
                                               Integer_rshift_tag>& a)
 : Natural(' ', ' ')
{
  get_memory(a.x.length()+DELTA);
  rshift(a.x, a.y);
}

inline Integer& Integer::operator=(const binder_arguments<Integer, size_t,
                                                          Integer_rshift_tag>& a)
{
  rshift(a.x, a.y);
  return *this;
}

inline binder_arguments<Integer, size_t, Integer_rshift_tag>
 operator>>(const Integer& a, const size_t& b)
// Algorithm:  c := a >> b
// Input:      a in Integer, b in size_t.
// Output:     c in Integer such that c = sign(a)*[|a|/2^b] ||
{
  return binder_arguments<Integer, size_t, Integer_rshift_tag>(a, b);
}

inline Integer::Integer(const binder_arguments<Integer, size_t,
                                               Integer_lshift_tag>& a)
 : Natural(' ', ' ')
{
  get_memory(a.x.length()+a.y%BETA+DELTA);
  lshift(a.x, a.y);
}

inline Integer& Integer::operator=(const binder_arguments<Integer, size_t,
                                                          Integer_lshift_tag>& a)
{
  lshift(a.x, a.y);
  return *this;
}

inline binder_arguments<Integer, size_t, Integer_lshift_tag>
 operator<<(const Integer& a, const size_t& b)
// Algorithm:  c := a << b
// Input:      a in Integer, b in size_t.
// Output:     c in Integer such that c = a*2^b ||
{
  return binder_arguments<Integer, size_t, Integer_lshift_tag>(a, b);
}

inline Integer& Integer::operator=(const Integer& a)
// Algorithm:  c := c = a
// Input:      a,c in Integer.
// Output:     c in Integer such that c = a ||
{
  sgn = a.sgn;
  abs() = ::abs(a);
  return *this;
}

inline Integer& Integer::operator=(const Natural& a)
// Algorithm:  c := c = a
// Input:      a in Natural,c in Integer.
// Output:     c in Integer such that c = a ||
{
  sgn = (a > 0);
  abs() = a;
  return *this;
}

inline SignDigit Integer::operator=(const SignDigit a)
// Algorithm:  c := b = a
// Input:      a in SignDigit, b in Integer.
// Output:     b in Integer, c in SignDigit such that b = a, c = a ||
{
  if (a >= 0) { sgn = (a > 0); abs() = a; }
  else { sgn = -1; abs() = -a; }
  return a;
}

inline Integer& Integer::operator|=(const Digit a)
// Algorithm:  c := c |= a
// Input:      c in Integer, a in Digit.
// Output:     c in Integer such that c := c or a ||
{
  if (sgn >= 0) abs() |= a;
  else *this |= Integer(a);
  return *this;
}

inline Integer Integer::operator++(int)
// Algorithm:  c := a++
// Input:      a in Integer.
// Output:     a,c in Integer such that c := a, a := a+1 ||
{
  Integer a(*this);
  ++(*this);
  return a;
}

inline Integer Integer::operator--(int)
// Algorithm:  c := a--
// Input:      a in Integer.
// Output:     a,c in Integer such that c := a, a := a-1 ||
{
  Integer a(*this);
  --(*this);
  return a;
}

inline Integer::Integer(const binder_arguments<Integer, Integer,
                                               Integer_divides_tag>& a)
 : Natural(' ', ' ')
{
  get_memory(a.x.length()+DELTA);
  div(a.x, a.y);
}

inline Integer& Integer::operator=(const binder_arguments<Integer, Integer,
                                                          Integer_divides_tag>& a)
{
  div(a.x, a.y);
  return *this;
}

inline binder_arguments<Integer, Integer, Integer_divides_tag>
 operator/(const Integer& a, const Integer& b)
// Algorithm:  c := a/b
// Input:      a,b in Integer where not b = 0.
// Output:     c in Integer such that c = sign(a*b)*[|a/b|] ||
{
  return binder_arguments<Integer, Integer, Integer_divides_tag>(a, b);
}

inline Integer::Integer(const binder_arguments<Integer, Integer,
                                               Integer_modulus_tag>& a)
 : Natural(' ', ' ')
{
  Integer t(a.x.length()+DELTA, ' ');
  get_memory(a.y.length()+DELTA);
  ::div(a.x, a.y, t, *this);
}

inline Integer& Integer::operator=(const binder_arguments<Integer, Integer,
                                                          Integer_modulus_tag>& a)
{
  Integer t(a.x.length()+DELTA, ' ');
  ::div(a.x, a.y, t, *this);
  return *this;
}

inline binder_arguments<Integer, Integer, Integer_modulus_tag>
 operator%(const Integer& a, const Integer& b)
// Algorithm:  c := a%b
// Input:      a,b in Integer where not b = 0.
// Output:     c in Integer such that c = a - sign(a*b)*[|a/b|]*b ||
{
  return binder_arguments<Integer, Integer, Integer_modulus_tag>(a, b);
}

inline Integer::Integer(const binder_arguments<Integer, Natural,
                                               Integer_divides_tag>& a)
 : Natural(' ', ' ')
{
  get_memory(a.x.length()+DELTA);
  div(a.x, a.y);
}

inline Integer& Integer::operator=(const binder_arguments<Integer, Natural,
                                                          Integer_divides_tag>& a)
{
  div(a.x, a.y);
  return *this;
}

inline Integer::Integer(const binder_arguments<Natural, Integer,
                                               Integer_divides_tag>& a)
 : Natural(' ', ' ')
{
  get_memory(a.x.length()+DELTA);
  div(a.x, a.y);
}

inline Integer& Integer::operator=(const binder_arguments<Natural, Integer,
                                                          Integer_divides_tag>& a)
{
  div(a.x, a.y);
  return *this;
}

inline binder_arguments<Integer, Natural, Integer_divides_tag>
 operator/(const Integer& a, const Natural& b)
// Algorithm:  c := a/b
// Input:      a in Integer, b in Natural where not b = 0.
// Output:     c in Integer such that c = sign(a)*[|a|/b] ||
{
  return binder_arguments<Integer, Natural, Integer_divides_tag>(a, b);
}

inline binder_arguments<Natural, Integer, Integer_divides_tag>
 operator/(const Natural& a, const Integer& b)
// Algorithm:  c := a/b
// Input:      a in Natural, b in Integer where not b = 0.
// Output:     c in Integer such that c = sign(b)*[a/|b|] ||
{
  return binder_arguments<Natural, Integer, Integer_divides_tag>(a, b);
}

inline Integer::Integer(const binder_arguments<Integer, Integer,
                                               Integer_not_tag>& a)
 : Natural(' ', ' ')
{
  get_memory(a.x.length()+DELTA);
  bitwise_not(a.x);
}

inline Integer& Integer::operator=(const binder_arguments<Integer, Integer,
                                                          Integer_not_tag>& a)
{
  bitwise_not(a.x);
  return *this;
}

inline binder_arguments<Integer, Integer, Integer_not_tag>
 operator~(const Integer& a)
// Algorithm:  c := ~a
// Input:      a in Integer.
// Output:     c in Integer such that c = not a ||
{
  return binder_arguments<Integer, Integer, Integer_not_tag>(a, a);
}

inline Integer::Integer(const binder_arguments<Integer, Integer,
                                               Integer_and_tag>& a)
 : Natural(' ', ' ')
{
  get_memory(::min(a.x.length(), a.y.length())+DELTA);
  bitwise_and(a.x, a.y);
}

inline Integer& Integer::operator=(const binder_arguments<Integer, Integer,
                                                          Integer_and_tag>& a)
{
  bitwise_and(a.x, a.y);
  return *this;
}

inline binder_arguments<Integer, Integer, Integer_and_tag>
 operator&(const Integer& a, const Integer& b)
// Algorithm:  c := a & b
// Input:      a,b in Integer.
// Output:     c in Integer such that c = a and b ||
{
  return binder_arguments<Integer, Integer, Integer_and_tag>(a, b);
}

inline Integer::Integer(const binder_arguments<Integer, Integer,
                                               Integer_or_tag>& a)
 : Natural(' ', ' ')
{
  get_memory(::max(a.x.length(), a.y.length())+DELTA);
  bitwise_or(a.x, a.y);
}

inline Integer& Integer::operator=(const binder_arguments<Integer, Integer,
                                                          Integer_or_tag>& a)
{
  bitwise_or(a.x, a.y);
  return *this;
}

inline binder_arguments<Integer, Integer, Integer_or_tag>
 operator|(const Integer& a, const Integer& b)
// Algorithm:  c := a | b
// Input:      a,b in Integer.
// Output:     c in Integer such that c = a or b ||
{
  return binder_arguments<Integer, Integer, Integer_or_tag>(a, b);
}

inline Integer::Integer(const binder_arguments<Integer, Integer,
                                               Integer_xor_tag>& a)
 : Natural(' ', ' ')
{
  get_memory(::max(a.x.length(), a.y.length())+DELTA);
  bitwise_xor(a.x, a.y);
}

inline Integer& Integer::operator=(const binder_arguments<Integer, Integer,
                                                          Integer_xor_tag>& a)
{
  bitwise_xor(a.x, a.y);
  return *this;
}

inline binder_arguments<Integer, Integer, Integer_xor_tag>
 operator^(const Integer& a, const Integer& b)
// Algorithm:  c := a ^ b
// Input:      a,b in Integer.
// Output:     c in Integer such that c = a xor b ||
{
  return binder_arguments<Integer, Integer, Integer_xor_tag>(a, b);
}

inline Integer::Integer(const binder_arguments<Integer, SignDigit,
                                               Integer_plus_tag>& a)
 : Natural(' ', ' ')
{
  get_memory(a.x.length()+DELTA);
  add(a.x, a.y);
}

inline Integer& Integer::operator=(const binder_arguments<Integer, SignDigit,
                                                          Integer_plus_tag>& a)
{
  if (this == &a.x) return *this += a.y;
  else { add(a.x, a.y); return *this; }
}

inline binder_arguments<Integer, SignDigit, Integer_plus_tag>
 operator+(const Integer& a, const SignDigit& b)
// Algorithm:  c := a+b
// Input:      a in Integer, SignDigit b.
// Output:     c in Integer such that c = a+b ||
{
  return binder_arguments<Integer, SignDigit, Integer_plus_tag>(a, b);
}

inline binder_arguments<Integer, SignDigit, Integer_plus_tag>
 operator+(const SignDigit& a, const Integer& b)
// Algorithm:  c := a+b
// Input:      a in SignDigit, Integer b.
// Output:     c in Integer such that c = a+b ||
{
  return binder_arguments<Integer, SignDigit, Integer_plus_tag>(b, a);
}

inline Integer::Integer(const binder_arguments<Integer, SignDigit,
                                               Integer_minus_tag>& a)
 : Natural(' ', ' ')
{
  get_memory(a.x.length()+DELTA);
  sub(a.x, a.y);
}

inline Integer& Integer::operator=(const binder_arguments<Integer, SignDigit,
                                                          Integer_minus_tag>& a)
{
  if (this == &a.x) return *this -= a.y;
  else { sub(a.x, a.y); return *this; }
}

inline Integer::Integer(const binder_arguments<SignDigit, Integer,
                                               Integer_minus_tag>& a)
 : Natural(' ', ' ')
{
  get_memory(a.y.length()+DELTA);
  sub(a.y, a.x);
  neg();
}

inline Integer& Integer::operator=(const binder_arguments<SignDigit, Integer,
                                                          Integer_minus_tag>& a)
{
  if (this == &a.y) *this -= a.x;
  else sub(a.y, a.x);
  neg();
  return *this;
}

inline binder_arguments<Integer, SignDigit, Integer_minus_tag>
 operator-(const Integer& a, const SignDigit& b)
// Algorithm:  c := a-b
// Input:      a in Integer, SignDigit b.
// Output:     c in Integer such that c = a-b ||
{
  return binder_arguments<Integer, SignDigit, Integer_minus_tag>(a, b);
}

inline binder_arguments<SignDigit, Integer, Integer_minus_tag>
 operator-(const SignDigit& a, const Integer& b)
// Algorithm:  c := a-b
// Input:      a in SignDigit, Integer b.
// Output:     c in Integer such that c = a+b ||
{
  return binder_arguments<SignDigit, Integer, Integer_minus_tag>(a, b);
}

inline Integer::Integer(const binder_arguments<Integer, SignDigit,
                                               Integer_multiplies_tag>& a)
 : Natural(' ', ' ')
{
  get_memory(a.x.length()+DELTA);
  mul(a.x, a.y);
}

inline Integer& Integer::operator=(const binder_arguments<Integer, SignDigit,
                                                          Integer_multiplies_tag>& a)
{
  mul(a.x, a.y);
  return *this;
}

inline Integer& Integer::operator+=(const binder_arguments<Integer, SignDigit,
                                                           Integer_multiplies_tag>& a)
{
  muladd(a.x, a.y);
  return *this;
}

inline Integer& Integer::operator-=(const binder_arguments<Integer, SignDigit,
                                                           Integer_multiplies_tag>& a)
{
  mulsub(a.x, a.y);
  return *this;
}

inline binder_arguments<Integer, SignDigit, Integer_multiplies_tag>
 operator*(const Integer& a, const SignDigit& b)
// Algorithm:  c := a*b
// Input:      a in Integer, b in SignDigit.
// Output:     c in Integer such that c = a*b ||
{
  return binder_arguments<Integer, SignDigit, Integer_multiplies_tag>(a, b);
}

inline binder_arguments<Integer, SignDigit, Integer_multiplies_tag>
 operator*(const SignDigit& a, const Integer& b)
// Algorithm:  c := a*b
// Input:      a in SignDigit, b in Integer.
// Output:     c in Integer such that c = a*b ||
{
  return binder_arguments<Integer, SignDigit, Integer_multiplies_tag>(b, a);
}

inline Integer operator/(const Integer& a, const SignDigit b)
// Algorithm:  c := a/b
// Input:      a in Integer, b in SignDigit where not b = 0.
// Output:     c in Integer such that c = sign(a*b)*[|a/b|] ||
{
  return Integer(a) /= b;
}

inline SignDigit operator%(const Integer& a, const SignDigit b)
// Algorithm:  c := a%b
// Input:      a in Integer, b in SignDigit where not b = 0.
// Output:     c in SignDigit such that c = a - sign(a*b)*[|a/b|]*b ||
{
  SignDigit c = abs(a) % Digit((b >= 0)? b : -b);
  if (sign(a) == -1) c = -c;
  return c;
}

inline Digit operator&(const Integer& a, const Digit b)
// Algorithm:  c := a & b
// Input:      a in Integer, b in Digit.
// Output:     c in Integer such that c = a and b ||
{
  return (sign(a) >= 0)? abs(a)&b : ((~a.lowest())+1)&b;
}

inline Integer operator|(const Integer& a, const Digit b)
// Algorithm:  c := a | b
// Input:      a in Integer, b in Digit.
// Output:     c in Integer such that c = a or b ||
{
  return Integer(a) |= b;
}

inline Digit log2(const Integer& a)
// Algorithm:  b := log2(a)
// Input:      a in Integer.
// Output:     b in Digit
//             such that if |a| > 0 then b = [log2(|a|)] else b = 0 ||
{
  return log2(abs(a));
}

inline int units(Integer& a)
// Algorithm:  c := units(a)
// Input:      a in Integer.
// Output:     a in Integer, c in int such that a := |a|,
//             if a >= 0 then c = 1 else c = -1 ||
{
  if (sign(a) >= 0) return 1;
  else { a = -a; return -1; }
}

inline Integer::Integer(const binder_arguments<Integer, Integer,
                                               Integer_square_root_tag>& a)
 : Natural(' ', ' ')
{
  get_memory(a.x.length()/2+DELTA);
  sqrt(a.x);
}

inline Integer& Integer::operator=(const binder_arguments<Integer, Integer,
                                                          Integer_square_root_tag>& a)
{
  sqrt(a.x);
  return *this;
}

inline binder_arguments<Integer, Integer, Integer_square_root_tag>
 sqrt(const Integer& a)
// Algorithm:  b := sqrt(a)
// Input:      a in Integer where a >= 0.
// Output:     b in Integer such that b = [sqrt(a)] ||
{
  return binder_arguments<Integer, Integer, Integer_square_root_tag>(a, a);
}

inline Integer atoI(const char* a, const Digit b)
// Algorithm:  c := atoI(a, b)
// Input:      a in String, b in Digit where 2 <= b <= 36.
// Output:     c in Integer such that c = a ||
//
// Note:       conversion string to Integer; return 0 by conversion error.
{
  Integer result;
  result.atoI(a, b);
  return result;
}

inline Integer::rep print(const Integer& a, bool b)
// Algorithm:  o := o << print(a, b)
// Input:      o in ostream, a in Integer, b in bool.
// Output:     o in ostream ||
//
// Note:       puts internal representation of Integer a on an output stream.
{
  return Integer::rep(sign(a), abs(a), b);
}

inline OSTREAM& operator<<(OSTREAM& out, const Integer::rep& a)
// puts internal representation of Integer a on output stream.
{
  if (a.sgn == -1) out << '-';
  return out << print(a.nat, a.bin);
}

inline OSTREAM& operator<<(OSTREAM& out, const Integer& a)
// Algorithm:  o := o << a
// Input:      o in ostream, a in Integer.
// Output:     o in ostream ||
//
// Note:       puts Integer a on output stream.
{
  if (sign(a) < 0) {
    const int b = out.width();
    if (b > 0) { out.width(0); out << '-'; out.width(b-1); }
    else out << '-';
  }
  return out << abs(a);
}


#endif
