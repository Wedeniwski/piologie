/////////////////////////////////
//
// Piologie V 1.3.2
// multi-precision arithmetic
// Rational
//
// (c) 1996-2001 HiPiLib
//               www.hipilib.de
//
// Sebastian Wedeniwski
// 10/21/2001
//

#ifndef _Include_Rational_H_
#define _Include_Rational_H_

#include "integer.h"

struct Rational_plus_tag {};
struct Rational_negate_tag {};
struct Rational_minus_tag {};
struct Rational_multiplies_tag {};
struct Rational_divides_tag {};
struct Rational_modulus_tag {};
struct Rational_square_root_tag {};
struct Rational_lshift_tag {};
struct Rational_rshift_tag {};

//#ifdef _Old_STD_
//template <class I, class N>
//#else
//template <class I, class N = I>
//#endif
class Rational {
public:
  struct rep {
    const Integer& num;
    const Natural& den;
    const bool     bin;
    rep(const Integer& a, const Natural& b, const bool c = false)
      : num(a), den(b), bin(c) {}
  };
  

private:
  Integer num;
  Natural den;

  void add(const Rational&, const Rational&);
  void sub(const Rational&, const Rational&);
  void mul(const Rational&, const Rational&);
  void sqr(const Rational&);
  void div(const Rational&, const Rational&);
  void lshift(const Rational&, size_t);
  void rshift(const Rational&, size_t);

public:
  Rational(const SignDigit = 0);
  Rational(const Integer&);
  Rational(const Integer&, const Natural&);
  Rational(const Rational&);
  ~Rational();

  const Integer&   numerator() const;
  const Natural&   denominator() const;

  SignDigit        operator=(const SignDigit);
  Rational&        operator=(const Integer&);
  Rational&        operator=(const Rational&);
  Rational&        operator+=(const Rational&);
  Rational&        operator-=(const Rational&);
  Rational&        operator*=(const Rational&);
  Rational&        operator/=(const Rational&);
  Rational&        operator<<=(const size_t);
  Rational&        operator>>=(const size_t);
  
  const Rational&  operator++();
  const Rational&  operator--();
  Rational         operator++(int);
  Rational         operator--(int);


  Rational(const binder_arguments<Rational, Rational, Rational_negate_tag>&);
  Rational& operator=(const binder_arguments<Rational, Rational, Rational_negate_tag>&);

  Rational(const binder_arguments<Rational, Rational, Rational_plus_tag>&);
  Rational& operator=(const binder_arguments<Rational, Rational, Rational_plus_tag>&);

  Rational(const binder_arguments<Rational, Rational, Rational_minus_tag>&);
  Rational& operator=(const binder_arguments<Rational, Rational, Rational_minus_tag>&);

  Rational(const binder_arguments<Rational, Rational, Rational_multiplies_tag>&);
  Rational& operator=(const binder_arguments<Rational, Rational, Rational_multiplies_tag>&);

  Rational(const binder_arguments<Rational, Rational, Rational_divides_tag>&);
  Rational& operator=(const binder_arguments<Rational, Rational, Rational_divides_tag>&);

  Rational(const binder_arguments<Rational, size_t, Rational_lshift_tag>&);
  Rational& operator=(const binder_arguments<Rational, size_t, Rational_lshift_tag>&);

  Rational(const binder_arguments<Rational, size_t, Rational_rshift_tag>&);
  Rational& operator=(const binder_arguments<Rational, size_t, Rational_rshift_tag>&);

  void        rand(const size_t);
  inline void swap(Rational&);
  bool        scan(ISTREAM&);
  const char* atoR(const char*, const Digit = 10);

  friend ISTREAM&     operator>>(ISTREAM&, Rational&);
  friend Rational     inv(const Rational&);
};

inline binder_arguments<Rational, Rational, Rational_negate_tag>
 operator-(const Rational&);
inline binder_arguments<Rational, Rational, Rational_plus_tag>
 operator+(const Rational&, const Rational&);
inline binder_arguments<Rational, Rational, Rational_minus_tag>
 operator-(const Rational&, const Rational&);
inline binder_arguments<Rational, Rational, Rational_multiplies_tag>
 operator*(const Rational&, const Rational&);
inline binder_arguments<Rational, Rational, Rational_divides_tag>
 operator/(const Rational&, const Rational&);
inline binder_arguments<Rational, size_t, Rational_lshift_tag>
 operator<<(const Rational&, const size_t&);
inline binder_arguments<Rational, size_t, Rational_rshift_tag>
 operator>>(const Rational&, const size_t&);


inline bool operator==(const Rational&, const Rational&);
inline bool operator!=(const Rational&, const Rational&);
inline bool operator<(const Rational&, const Rational&);
inline bool operator<=(const Rational&, const Rational&);
inline bool operator>(const Rational&, const Rational&);
inline bool operator>=(const Rational&, const Rational&);

inline bool operator==(const Rational&, const SignDigit);
inline bool operator!=(const Rational&, const SignDigit);
inline bool operator<(const Rational&, const SignDigit);
inline bool operator<=(const Rational&, const SignDigit);
inline bool operator>(const Rational&, const SignDigit);
inline bool operator>=(const Rational&, const SignDigit);

inline void           swap(Rational&, Rational&);
inline const Integer& numerator(const Rational&);
inline const Natural& denominator(const Rational&);
inline Rational       abs(const Rational&);
inline int            sign(const Rational&);
inline Rational       atoR(const char*, const Digit = 10);
char*                 Rtoa(const Rational&, char*, const Digit = 10);
inline Rational::rep  print(const Rational&, bool = false);
ISTREAM&              operator>>(ISTREAM&, Rational&);
inline OSTREAM&       operator<<(OSTREAM&, const Rational::rep&);
inline OSTREAM&       operator<<(OSTREAM&, const Rational&);
Integer               ceil(const Rational&);
Integer               floor(const Rational&);
inline Integer        round(const Rational&);
Integer               trunc(const Rational&);
Rational              pow(const Rational&, const SignDigit);
Rational              pow(const Rational&, const Integer&);

ISTREAM&              operator>>(ISTREAM&, Rational&);
Rational              inv(const Rational&);

// Unfortunatelly not:
//typedef Rational<Integer, Natural> Rational;



///////////////////////// Inline-Implementation ////////////////////

inline void Rational::sqr(const Rational& a)
// Algorithm:  b.sqr(a)
// Input:      a in Rational.
// Output:     b in Rational such that b = a^2 ||
{
  num = a.num*a.num;
  den = a.den*a.den;
}

inline Rational::Rational(const SignDigit a)
 : num(a), den(1)
// Algorithm:  c := Rational(a)
// Input:      a in SignDigit.
// Output:     c in Rational such that c = a/1 ||
{
}

inline Rational::Rational(const Integer& a)
 : num(a), den(1)
// Algorithm:  c := Rational(a)
// Input:      a in Integer.
// Output:     c in Rational such that c = a/1 ||
{
}

inline Rational::Rational(const Rational& a)
 : num(a.num), den(a.den)
// Algorithm:  c := Rational(a)
// Input:      a in Rational.
// Output:     c in Rational such that c = a ||
{
}

inline Rational::~Rational()
{
}

inline const Integer& Rational::numerator() const
// Algorithm:  c := a.numerator()
// Input:      a in Rational.
// Output:     c in Integer such that c = a_1 where a_1/a_2 = a ||
{
  return num;
}

inline const Natural& Rational::denominator() const
// Algorithm:  c := a.denominator()
// Input:      a in Rational.
// Output:     c in Natural such that c = a_2 where a_1/a_2 = a ||
{
  return den;
}

inline SignDigit Rational::operator=(const SignDigit a)
// Algorithm:  c := b = a
// Input:      a in SignDigit, b in Rational.
// Output:     b in Rational, c in SignDigit such that b = a/1, c = a ||
{
  num = a; den = 1;
  return a;
}

inline Rational& Rational::operator=(const Integer& a)
// Algorithm:  c := c = a
// Input:      a in Integer, c in Rational.
// Output:     c in Rational such that c = a/1 ||
{
  num = a; den = 1;
  return *this;
}

inline Rational& Rational::operator=(const Rational& a)
// Algorithm:  c := c = a
// Input:      a,c in Rational.
// Output:     c in Rational such that c = a ||
{
  num = a.num; den = a.den;
  return *this;
}

inline Rational& Rational::operator+=(const Rational& a)
// Algorithm:  c := c += a
// Input:      a,c in Rational.
// Output:     c in Rational such that c := c+a ||
{
  add(*this, a);
  return *this;
}

inline Rational& Rational::operator-=(const Rational& a)
// Algorithm:  c := c -= a
// Input:      a,c in Rational.
// Output:     c in Rational such that c := c-a ||
{
  sub(*this, a);
  return *this;
}

inline Rational& Rational::operator*=(const Rational& a)
// Algorithm:  c := c *= a
// Input:      a,c in Rational.
// Output:     c in Rational such that c := c*a ||
{
  if (this == &a) sqr(*this);
  else mul(*this, a);
  return *this;
}

inline Rational& Rational::operator/=(const Rational& a)
// Algorithm:  c := c /= a
// Input:      a,c in Rational where not a = 0.
// Output:     c in Rational such that c := c/a ||
{
  div(*this, a);
  return *this;
}

inline Rational& Rational::operator<<=(const size_t a)
// Algorithm:  c := c <<= a
// Input:      a in size_t, c in Rational.
// Output:     c in Rational such that c := c * 2^a ||
{
  lshift(*this, a);
  return *this;
}

inline Rational& Rational::operator>>=(const size_t a)
// Algorithm:  c := c >>= a
// Input:      a in size_t, c in Rational.
// Output:     c in Rational such that c := c / 2^a ||
{
  rshift(*this, a);
  return *this;
}

inline const Rational& Rational::operator++()
// Algorithm:  c := ++a
// Input:      a in Rational.
// Output:     a,c in Rational such that a := a+1, c := a ||
{
  num += den;
  return *this;
}

inline const Rational& Rational::operator--()
// Algorithm:  c := --a
// Input:      a in Rational.
// Output:     a,c in Rational such that a := a-1, c := a ||
{
  num -= den;
  return *this;
}

inline Rational Rational::operator++(int)
// Algorithm:  c := a++
// Input:      a in Rational.
// Output:     a,c in Rational such that c := a, a := a+1 ||
{
  Rational a(*this);
  ++(*this);
  return a;
}

inline Rational Rational::operator--(int)
// Algorithm:  c := a--
// Input:      a in Rational.
// Output:     a,c in Rational such that c := a, a := a-1 ||
{
  Rational a(*this);
  --(*this);
  return a;
}

inline Rational::Rational(const binder_arguments<Rational, Rational,
                                                   Rational_negate_tag>& a)
 : num(-a.x.num), den(a.x.den)
{
}

inline Rational& Rational::operator=(const binder_arguments<Rational, Rational,
                                                            Rational_negate_tag>& a)
{
  num = -a.x.num;
  if (this != &a.x) den = a.x.den;
  return *this;
}

inline binder_arguments<Rational, Rational, Rational_negate_tag>
 operator-(const Rational& a)
// Algorithm:  c := -a
// Input:      a in Rational.
// Output:     c in Rational such that c = -a ||
{
  return binder_arguments<Rational, Rational, Rational_negate_tag>(a, a);
}

inline Rational::Rational(const binder_arguments<Rational, Rational,
                                                   Rational_plus_tag>& a)
{
  add(a.x, a.y);
}

inline Rational& Rational::operator=(const binder_arguments<Rational, Rational,
                                                            Rational_plus_tag>& a)
{
  add(a.x, a.y);
  return *this;
}

inline binder_arguments<Rational, Rational, Rational_plus_tag>
 operator+(const Rational& a, const Rational& b)
// Algorithm:  c := a+b
// Input:      a,b in Rational.
// Output:     c in Rational such that c = a+b ||
{
  return binder_arguments<Rational, Rational, Rational_plus_tag>(a, b);
}

inline Rational::Rational(const binder_arguments<Rational, Rational,
                                                   Rational_minus_tag>& a)
{
  sub(a.x, a.y);
}

inline Rational& Rational::operator=(const binder_arguments<Rational, Rational,
                                                            Rational_minus_tag>& a)
{
  sub(a.x, a.y);
  return *this;
}

inline binder_arguments<Rational, Rational, Rational_minus_tag>
 operator-(const Rational& a, const Rational& b)
// Algorithm:  c := a-b
// Input:      a,b in Rational.
// Output:     c in Rational such that c = a-b ||
{
  return binder_arguments<Rational, Rational, Rational_minus_tag>(a, b);
}

inline Rational::Rational(const binder_arguments<Rational, Rational,
                                                 Rational_multiplies_tag>& a)
{
  if (&a.x == &a.y) sqr(a.x);
  else mul(a.x, a.y);
}

inline Rational& Rational::operator=(const binder_arguments<Rational, Rational,
                                                            Rational_multiplies_tag>& a)
{
  if (&a.x == &a.y) sqr(a.x);
  else mul(a.x, a.y);
  return *this;
}

inline binder_arguments<Rational, Rational, Rational_multiplies_tag>
 operator*(const Rational& a, const Rational& b)
// Algorithm:  c := a*b
// Input:      a,b in Rational.
// Output:     c in Rational such that c = a*b ||
{
  return binder_arguments<Rational, Rational, Rational_multiplies_tag>(a, b);
}

inline Rational::Rational(const binder_arguments<Rational, Rational,
                                                 Rational_divides_tag>& a)
{
  div(a.x, a.y);
}

inline Rational& Rational::operator=(const binder_arguments<Rational, Rational,
                                                            Rational_divides_tag>& a)
{
  div(a.x, a.y);
  return *this;
}

inline binder_arguments<Rational, Rational, Rational_divides_tag>
 operator/(const Rational& a, const Rational& b)
// Algorithm:  c := a/b
// Input:      a,b in Rational.
// Output:     c in Rational such that c = a/b ||
{
  return binder_arguments<Rational, Rational, Rational_divides_tag>(a, b);
}

inline Rational::Rational(const binder_arguments<Rational, size_t,
                                                   Rational_lshift_tag>& a)
{
  lshift(a.x, a.y);
}

inline Rational& Rational::operator=(const binder_arguments<Rational, size_t,
                                                            Rational_lshift_tag>& a)
{
  lshift(a.x, a.y);
  return *this;
}

inline binder_arguments<Rational, size_t, Rational_lshift_tag>
 operator<<(const Rational& a, const size_t& b)
// Algorithm:  c := a << b
// Input:      a in Rational, b in size_t.
// Output:     c in Rational such that c = a * 2^b ||
{
  return binder_arguments<Rational, size_t, Rational_lshift_tag>(a, b);
}

inline Rational::Rational(const binder_arguments<Rational, size_t,
                                                 Rational_rshift_tag>& a)
{
  rshift(a.x, a.y);
}

inline Rational& Rational::operator=(const binder_arguments<Rational, size_t,
                                                            Rational_rshift_tag>& a)
{
  rshift(a.x, a.y);
  return *this;
}

inline binder_arguments<Rational, size_t, Rational_rshift_tag>
 operator>>(const Rational& a, const size_t& b)
// Algorithm:  c := a >> b
// Input:      a in Rational, b in size_t.
// Output:     c in Rational such that c = a / 2^b ||
{
  return binder_arguments<Rational, size_t, Rational_rshift_tag>(a, b);
}

inline void Rational::rand(const size_t n)
// Algorithm:  a.rand(n)
// Input:      n in size_t.
// Output:     a in Rational
//             such that |a_1| < 2^n, a_2 <= 2^n where a_1/a_2 = a (random number) ||
{
  num.rand(n);
  den.rand(n); ++den;
  const Natural t = gcd(abs(num), den);
  num /= t; den /= t;
}

inline void Rational::swap(Rational& a)
// Algorithm:  a.swap(b)
// Input:      a,b in Rational.
// Output:     a,b in Rational such that t := a, a := b, b := t
//             where t in Rational ||
{
  ::swap(num, a.num);
  ::swap(den, a.den);
}

inline void swap(Rational& a, Rational& b)
// Algorithm:  swap(a, b)
// Input:      a,b in Rational.
// Output:     a,b in Rational such that t := a, a := b, b := t
//             where t in Rational ||
{
  a.swap(b);
}

inline bool operator==(const Rational& a, const Rational& b)
// Algorithm:  c := a == b
// Input:      a,b in Rational.
// Output:     c in bool such that if a = b then c = true else c = false ||
{
  return (a.numerator() == b.numerator() && a.denominator() == b.denominator());
}


inline bool operator!=(const Rational& a, const Rational& b)
// Algorithm:  c := a != b
// Input:      a,b in Rational.
// Output:     c in bool such that if a = b then c = false else c = true ||
{
  return (a.numerator() != b.numerator() || a.denominator() != b.denominator());
}

inline bool operator<(const Rational& a, const Rational& b)
// Algorithm:  c := a < b
// Input:      a,b in Rational.
// Output:     c in bool such that if a < b then c = true else c = false ||
{
  return (a.numerator()*b.denominator() < a.denominator()*b.numerator());
}

inline bool operator<=(const Rational& a, const Rational& b)
// Algorithm:  c := a <= b
// Input:      a,b in Rational.
// Output:     c in bool such that if a <= b then c = true else c = false ||
{
  return (a.numerator()*b.denominator() <= a.denominator()*b.numerator());
}

inline bool operator>(const Rational& a, const Rational& b)
// Algorithm:  c := a > b
// Input:      a,b in Rational.
// Output:     c in bool such that if a > b then c = true else c = false ||
{
  return (a.numerator()*b.denominator() > a.denominator()*b.numerator());
}

inline bool operator>=(const Rational& a, const Rational& b)
// Algorithm:  c := a >= b
// Input:      a,b in Rational.
// Output:     c in bool such that if a >= b then c = true else c = false ||
{
  return (a.numerator()*b.denominator() >= a.denominator()*b.numerator());
}

inline bool operator==(const Rational& a, const SignDigit b)
// Algorithm:  c := a == b
// Input:      a in Rational, b in SignDigit.
// Output:     c in bool such that if a = b then c = true else c = false ||
{
  return (a.numerator() == b && a.denominator() == 1);
}

inline bool operator!=(const Rational& a, const SignDigit b)
// Algorithm:  c := a != b
// Input:      a in Rational, b in SignDigit.
// Output:     c in bool such that if a = b then c = false else c = true ||
{
  return (a.numerator() != b || a.denominator() != 1);
}

inline bool operator<(const Rational& a, const SignDigit b)
// Algorithm:  c := a < b
// Input:      a in Rational, b in SignDigit.
// Output:     c in bool such that if a < b then c = true else c = false ||
{
  return (b < 0)? (a.numerator() > a.denominator()*Digit(-b))
                : (a.numerator() < a.denominator()*Digit(b));
}

inline bool operator<=(const Rational& a, const SignDigit b)
// Algorithm:  c := a <= b
// Input:      a in Rational, b in SignDigit.
// Output:     c in bool such that if a <= b then c = true else c = false ||
{
  return (b < 0)? (a.numerator() >= a.denominator()*Digit(-b))
                : (a.numerator() <= a.denominator()*Digit(b));
}

inline bool operator>(const Rational& a, const SignDigit b)
// Algorithm:  c := a > b
// Input:      a in Rational, b in SignDigit.
// Output:     c in bool such that if a > b then c = true else c = false ||
{
  return (b < 0)? (a.numerator() < a.denominator()*Digit(-b))
                : (a.numerator() > a.denominator()*Digit(b));
}

inline bool operator>=(const Rational& a, const SignDigit b)
// Algorithm:  c := a >= b
// Input:      a in Rational, b in SignDigit.
// Output:     c in bool such that if a >= b then c = true else c = false ||
{
  return (b < 0)? (a.numerator() <= a.denominator()*Digit(-b))
                : (a.numerator() >= a.denominator()*Digit(b));
}

inline const Integer& numerator(const Rational& a)
// Algorithm:  c := numerator(a)
// Input:      a in Rational.
// Output:     c in Integer such that c = a_1 where a_1/a_2 = a ||
{
  return a.numerator();
}

inline const Natural& denominator(const Rational& a)
// Algorithm:  c := denominator(a)
// Input:      a in Rational.
// Output:     c in Natural such that c = a_2 where a_1/a_2 = a ||
{
  return a.denominator();
}

inline Rational abs(const Rational& a)
// Algorithm:  c := abs(a)
// Input:      a in Rational.
// Output:     c in Rational such that c = |a| ||
{
  return Rational(abs(a.numerator()), a.denominator());
}

inline int sign(const Rational& a)
// Algorithm:  c := sign(a)
// Input:      a in Rational.
// Output:     c in int such that if a = 0 then c = 0
//             else if a > 0 then c = 1 else c = -1 ||
{
  return sign(a.numerator());
}

inline Rational atoR(const char* a, const Digit b)
// Algorithm:  c := atoR(a, b)
// Input:      a in String, b in Digit where 2 <= b <= 36.
// Output:     c in Rational such that c = a ||
//
// Note:       conversion string to Rational; return 0 by conversion error.
{
  Rational result;
  result.atoR(a, b);
  return result;
}

inline Rational::rep print(const Rational& a, bool b)
// Algorithm:  o := o << print(a, b)
// Input:      o in ostream, a in Rational, b in bool.
// Output:     o in ostream ||
//
// Note:       puts internal representation of Rational a on an output stream.
{
  return Rational::rep(a.numerator(), a.denominator(), b);
}

inline OSTREAM& operator<<(OSTREAM& out, const Rational::rep& a)
// puts internal representation of Rational a on output stream.
{
  return out << print(a.num, a.bin) << '/' << print(a.den, a.bin);
}

inline OSTREAM& operator<<(OSTREAM& out, const Rational& a)
// Algorithm:  o := o << a
// Input:      o in ostream, a in Rational.
// Output:     o in ostream ||
//
// Note:       puts Rational a on output stream.
{
  return out << a.numerator() << '/' << a.denominator();
}

inline Integer round(const Rational& a)
// Algorithm:  c := round(a)
// Input:      a in Rational.
// Output:     c in Integer such that c = floor(a+1/2) ||
{
  return floor(a + Rational(1, 2));
}


#endif
