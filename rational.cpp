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

#include "rational.h"
#include <string.h>


#include <cstdio>
using namespace std;


/////////////////// Rational Arithmetic ////////////////////////

void Rational::add(const Rational& a, const Rational& b)
// Algorithm:  c.add(a, b)
// Input:      a,b in Rational.
// Output:     c in Rational such that c = a + b ||
{
  if (sign(a.num) == 0) { num = b.num; den = b.den; }
  else if (sign(b.num) == 0) { num = a.num; den = a.den; }
  else {
    Natural g = gcd(a.den, b.den);
    if (g == 1) {
      const Integer h = a.den * b.num;
      num = a.num * b.den;
      num += h;
      den = a.den * b.den;
    } else {
      Natural s = b.den / g;
      Integer t = a.num * s;
      s = a.den / g;
      t += b.num * s;
      g = gcd(abs(t), g);
      num = t / g;
      g = b.den / g;
      den = s * g;
    }
  }
}

void Rational::sub(const Rational& a, const Rational& b)
// Algorithm:  c.sub(a, b)
// Input:      a,b in Rational.
// Output:     c in Rational such that c = a - b ||
{
  if (sign(a.num) == 0) { num = -b.num; den = b.den; }
  else if (sign(b.num) == 0) { num = a.num; den = a.den; }
  else {
    Natural g = gcd(a.den, b.den);
    if (g == 1) {
      const Integer h = a.den * b.num;
      num = a.num * b.den;
      num -= h;
      den = a.den * b.den;
    } else {
      Natural s = b.den / g;
      Integer t = a.num * s;
      s = a.den / g;
      t -= b.num * s;
      g = gcd(abs(t), g);
      num = t / g;
      g = b.den / g;
      den = s * g;
    }
  }
}

void Rational::mul(const Rational& a, const Rational& b)
// Algorithm:  c.mul(a, b)
// Input:      a,b in Rational.
// Output:     c in Rational such that c = a * b ||
{
  Natural g = gcd(abs(a.num), b.den);
  Natural h = gcd(a.den, abs(b.num));
  if (g == 1) {
    if (h == 1) {
      num = a.num * b.num;
      den = a.den * b.den;
    } else {
      Integer t = b.num / h;
      g = a.den / h;
      num = a.num * t;
      den = g * b.den;
    }
  } else if (h == 1) {
    Integer t = a.num / g;
    h = b.den / g;
    num = t * b.num;
    den = a.den * h;
  } else {
    Integer s = a.num / g;
    Integer t = b.num / h;
    h = a.den / h;
    g = b.den / g;
    num = s * t;
    den = h * g;
  }
}

void Rational::div(const Rational& a, const Rational& b)
// Algorithm:  c.div(a, b)
// Input:      a,b in Rational where not b = 0.
// Output:     c in Rational such that c = a / b ||
{
  if (b.num == 0) b.num.errmsg(4, "(div)");
  Natural g = gcd(abs(a.num), abs(b.num));
  Natural h = gcd(a.den, b.den);
  if (g == 1) {
    if (h == 1) {
      if (this == &b) {
        const bool i = (sign(num) == -1);
        g = a.den * abs(num);
        num = a.num * den;
        den = g;
        if (i) num = -num;
      } else {
        num = a.num * b.den;
        den = a.den * abs(b.num);
        if (sign(b.num) == -1) num = -num;
      }
    } else {
      Natural t = b.den / h;
      g = a.den / h;
      den = g * abs(b.num);
      if (sign(b.num) == -1) {
        num = a.num * t;
        num = -num;
      } else num = a.num * t;
    }
  } else if (h == 1) {
    Integer t = a.num / g;
    Integer s = b.num / g;
    num = t * b.den;
    den = a.den * abs(s);
    if (sign(s) == -1) num = -num;
  } else {
    Integer s = a.num / g;
    Natural t = b.den / h;
    Integer r = b.num / g;
    num = s * t;
    t = a.den / h;
    den = t * abs(r);
    if (sign(r) == -1) num = -num;
  }
}

void Rational::lshift(const Rational& a, size_t b)
// Algorithm:  c.lshift(a, b)
// Input:      a in Rational, b in size_t.
// Output:     c in Rational such that c = a * 2^b ||
{
  size_t i = 0;
  const Natural& c = a.den;
  while (i < b && !c.testbit(i)) ++i;
  den = a.den >> i;
  num = a.num << (b-i);
}

void Rational::rshift(const Rational& a, size_t b)
// Algorithm:  c.rshift(a, b)
// Input:      a in Rational, b in size_t.
// Output:     c in Rational such that c = a / 2^b ||
{
  size_t i = 0;
  const Natural& c = abs(a.num);
  while (i < b && !c.testbit(i)) ++i;
  num = a.num >> i;
  den = a.den << (b-i);
}

Rational::Rational(const Integer& a, const Natural& b)
 : num(a), den(b)
// Algorithm:  c := Rational(a, b)
// Input:      a in Integer, b in Natural where not b = 0.
// Output:     c in Rational such that c = a/b ||
{
  if (b == 0) den.errmsg(4, "(constructor)");
  const Natural t = gcd(abs(a), b);
  num /= t; den /= t;
}

bool Rational::scan(ISTREAM& in)
// Algorithm:  b := a.scan(i)
// Input:      a in Rational, i in istream.
// Output:     a in Rational, i in istream, b in bool ||
//
// Note:       gets Rational a as an internal representation from input stream
//             if b is true.
{
  if (!num.scan(in)) return false;
  char c = 0;
  if (in.get(c) && c != '/') { in.putback(c); return false; }
  return den.scan(in);
}

ISTREAM& operator>>(ISTREAM& in, Rational& a)
// Algorithm:  i := i >> a
// Input:      i in istream.
// Output:     i in istream, a in Rational ||
//
// Note:       gets Rational a from input stream.
{
  in >> a.num;
  char ch = 0;
  if (in.get(ch) && ch != '/') { in.putback(ch); return in; }
  in >> a.den;
  if (a.den == 0) a.den = 1;
  Natural b = gcd(abs(a.num), a.den);
  a.num /= b; a.den /= b;
  return in;
}

char* Rtoa(const Rational& a, char* b, const Digit c)
// Algorithm:  c := Rtoa(a, c, b)
// Input:      a in Rational, b in Digit, c in String
//             where 2 <= b <= 36, sizeof(c) > BETA*L(a)/log2(b).
// Output:     c in String such that c = a ||
//
// Note:       conversion Rational to string.
{
  Itoa(a.numerator(), b, c);
  char* d = b+strlen(b);
  *d = '/';
  Ntoa(a.denominator(), ++d, c);
  return b;
}

const char* Rational::atoR(const char* a, const Digit b)
// Algorithm:  c := d.atoR(a, b)
// Input:      d in Rational, a in String, b in Digit where 2 <= b <= 36.
// Output:     d in Rational, c in String such that d = a ||
//
// Note: Returns a pointer to the first occurrence of a non-digit character
//       in a after the character '/'.
{
  a = num.atoI(a, b);
  if (*a == '/') a = den.atoN(++a, b);
  if (den == 0) den = 1;
  const Natural c = gcd(abs(num), den);
  num /= c; den /= c;
  return a;
}


Rational inv(const Rational& a)
// Algorithm:  b := inv(a)
// Input:      a in Rational where not a = 0.
// Output:     b in Rational such that b = a^(-1) ||
{
  if (a.num == 0) a.num.errmsg(4, "(inv)");
  Rational b(a.den, abs(a.num));
  if (sign(a.num) == -1) b.num = -b.num;
  return b;
}

Integer ceil(const Rational& a)
// Algorithm:  c := ceil(a)
// Input:      a in Rational.
// Output:     c in Integer such that c = min{x in Integer | x >= a} ||
{
  if (a == 0) return SignDigit(0);
  Natural q,r;
  div(abs(a.numerator()), a.denominator(), q, r);
  if (sign(a) == -1) return -Integer(q);
  r <<= 1;
  if (r >= q) ++q;
  return q;
}

Integer floor(const Rational& a)
// Algorithm:  c := floor(a)
// Input:      a in Rational.
// Output:     c in Integer such that c = max{x in Integer | x <= a} ||
{
  Natural q,r;
  div(abs(a.numerator()), a.denominator(), q, r);
  Integer t = q;
  if (sign(a) == -1) {
    t = -t;
    if (r != 0) --t;
  }
  return t;
}

Integer trunc(const Rational& a)
// Algorithm:  c := trunc(a)
// Input:      a in Rational.
// Output:     c in Integer such that c = sign(a)*[|a|] ||
{
  Natural q,r;
  div(abs(a.numerator()), a.denominator(), q, r);
  if (sign(a) == -1) return -Integer(q);
  return q;
}

Rational pow(const Rational& a, const SignDigit b)
// Algorithm:  c := pow(a, b)
// Input:      a in Rational, b in SignDigit.
// Output:     c in Rational such that c = a^b ||
{
  if (b >= 0) return Rational(pow(a.numerator(), b), pow(a.denominator(), Digit(b)));
  else {
    Rational c(pow(a.denominator(), Digit(-b)), pow(abs(a.numerator()), Digit(-b)));
    if ((b&1) && sign(a) == -1) c = -c;
    return c;
  }
}

Rational pow(const Rational& a, const Integer& b)
// Algorithm:  c := pow(a, b)
// Input:      a in Rational, b in Integer.
// Output:     c in Rational such that c = a^b ||
{
  if (sign(b) >= 0) return Rational(pow(a.numerator(), b), pow(a.denominator(), abs(b)));
  else {
    const Natural& d = abs(b);
    Rational c(pow(a.denominator(), d), pow(abs(a.numerator()), d));
    if ((b&1) && sign(a) == -1) c = -c;
    return c;
  }
}

