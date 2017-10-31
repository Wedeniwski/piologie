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

#include "integer.h"


#ifdef _Piologie_Debug_
# define CONDITION(expr)     assert(expr)

# define NATURALCONDITION(a)                            \
  if ((a).root) {                                       \
      CONDITION(*(a).root == 0 &&                       \
                *((a).p-1) == 0 &&                      \
                (a).root < (a).p &&                     \
                (a).size > 0 &&                         \
                ((a).size == 1 || *(a).p != 0));        \
      const Digit* _pE = (a).p;                         \
      for (Digit* _pA = (a).root; _pA != _pE; ++_pA)    \
        CONDITION(*_pA == 0);                           \
  }

# define NATURAL_FOR_CHECK(a, b)                        \
  Natural (a) = (b);

# define INTEGERCONDITION(a)                            \
  CONDITION((a).highest() == 0 && sign(a) == 0 ||       \
            (a).highest() != 0 &&                       \
            (sign(a) == 1 || sign(a) == -1));

# define INTEGER_FOR_CHECK(a, b)                        \
  Integer (a) = (b);

#else
# define CONDITION(__ignore)        ((void)0)
# define NATURALCONDITION(__ignore) ((void)0)
# define INTEGERCONDITION(__ignore) ((void)0)
# define NATURAL_FOR_CHECK(__ignore1, __ignore2) ((void)0)
# define INTEGER_FOR_CHECK(__ignore1, __ignore2) ((void)0)
#endif



/////////////////// Integer Arithmetic /////////////////////////

void Integer::add(const Integer& a, const Integer& b)
// Algorithm:  c.add(a, b)
// Input:      a,b in Integer.
// Output:     c in Integer such that c = a+b ||
{
  INTEGERCONDITION(a);
  INTEGERCONDITION(b);

  const int sA = a.sgn;
  const int sB = b.sgn;
  if (sA == 0) *this = b;
  else if (sB == 0) *this = a;
  else if (sA == sB) {
    abs() = ::abs(a) + ::abs(b);
    sgn = sA;
  } else sgn = sA*::abs(::abs(a), ::abs(b), abs());
  
  INTEGERCONDITION(*this);
}

void Integer::add(const Integer& a, const Natural& b)
// Algorithm:  c.add(a, b)
// Input:      a in Integer, b in Natural.
// Output:     c in Integer such that c = a+b ||
{
  INTEGERCONDITION(a);

  const int sA = a.sgn;
  if (sA == 0) *this = b;
  else if (sA == 1) {
    abs() = ::abs(a) + b;
    sgn = 1;
  } else sgn = -::abs(::abs(a), b, abs());
  
  INTEGERCONDITION(*this);
}

void Integer::sub(const Integer& a, const Integer& b)
// Algorithm:  c.sub(a, b)
// Input:      a,b in Integer.
// Output:     c in Integer such that c = a-b ||
{
  INTEGERCONDITION(a);
  INTEGERCONDITION(b);

  const int sA = a.sgn;
  const int sB = b.sgn;
  if (sA == 0) neg(b);
  else if (sB == 0) *this = a;
  else if (sA != sB) {
    abs() = ::abs(a) + ::abs(b);
    sgn = sA;
  } else sgn = sA*::abs(::abs(a), ::abs(b), abs());
  
  INTEGERCONDITION(*this);
}

void Integer::sub(const Integer& a, const Natural& b)
// Algorithm:  c.sub(a, b)
// Input:      a in Integer, b in Natural.
// Output:     c in Integer such that c = a-b ||
{
  INTEGERCONDITION(a);

  const int sA = a.sgn;
  if (sA == 0) { *this = b; neg(); }
  else if (sA == -1) {
    abs() = ::abs(a) + b;
    sgn = -1;
  } else sgn = ::abs(::abs(a), b, abs());
  
  INTEGERCONDITION(*this);
}

void Integer::sub(const Natural& a, const Integer& b)
// Algorithm:  c.sub(a, b)
// Input:      a in Natural, b in Integer.
// Output:     c in Integer such that c = a-b ||
{
  INTEGERCONDITION(b);

  const int sB = b.sgn;
  if (sB == 0) *this = a;
  else if (sB == -1) {
    abs() = a + ::abs(b);
    sgn = 1;
  } else sgn = ::abs(a, ::abs(b), abs());
  
  INTEGERCONDITION(*this);
}

void Integer::div(const Integer& a, const Integer& b)
// Algorithm:  c.div(a, b)
// Input:      a,b,c in Integer where not b = 0;
// Output:     c in Integer such that c = sign(a*b)*[|a/b|] ||
{
  INTEGERCONDITION(a);
  INTEGERCONDITION(b);

  Natural t;
  Natural::div(::abs(a), ::abs(b), t);
  sgn = (::abs(*this) != 0)? a.sgn*b.sgn : 0;

  INTEGERCONDITION(*this);
}

void Integer::sqrt(const Integer& a)
// Algorithm:  b.sqrt(a)
// Input:      a,b in Integer where a >= 0.
// Output:     b in Integer such that b = [sqrt(a)] ||
{
  INTEGERCONDITION(a);

  const int sA = a.sgn;
  if (sA == -1) errmsg(6, "(sqrt)");
  Natural::sqrt(::abs(a));
  sgn = sA;

  INTEGERCONDITION(*this);
}

Integer& Integer::operator+=(const SignDigit a)
// Algorithm:  c := c += a
// Input:      a SignDigit, c in Integer.
// Output:     c in Integer such that c := c+a ||
{
  INTEGERCONDITION(*this);

  if (a < 0) return *this -= -a;
  const int sT = sgn;
  if (sT > 0) abs() += a;
  else if (sT == 0) *this = a;
  else if (::abs(*this) > Digit(a)) abs() -= a;
  else *this = a - highest();

  INTEGERCONDITION(*this);

  return *this;
}

Integer& Integer::operator-=(const SignDigit a)
// Algorithm:  c := c -= a
// Input:      a SignDigit, c in Integer.
// Output:     c in Integer such that c := c-a ||
{
  INTEGERCONDITION(*this);

  if (a < 0) return *this += -a;
  const int sT = sgn;
  if (sT < 0) abs() += a;
  else if (sT == 0) neg(a);
  else if (::abs(*this) > Digit(a)) abs() -= a;
  else {
    *this = a - highest();
    neg();
  }

  INTEGERCONDITION(*this);

  return *this;
}

Digit Integer::operator&=(const Digit b)
// Algorithm:  c := a &= b
// Input:      a in Integer, b in Digit.
// Output:     a in Integer, c in Digit such that a := a and b, c = a ||
{
  INTEGERCONDITION(*this);

  Digit c;
  if (sgn >= 0) {
    c = abs() &= b;
    sgn = (c != 0);
  } else *this = c = ((~lowest())+1)&b;

  INTEGERCONDITION(*this);

  return c;
}

const Integer& Integer::operator++()
// Algorithm:  c := ++a
// Input:      a in Integer.
// Output:     a,c in Integer such that a := a+1, c := a ||
{
  INTEGERCONDITION(*this);

  const int sT = sgn;
  if (sT == 0) *this = 1;
  else if (sT > 0) ++abs();
  else if (::abs(*this) == 1) *this = 0;
  else --abs();

  INTEGERCONDITION(*this);

  return *this;
}

const Integer& Integer::operator--()
// Algorithm:  c := --a
// Input:      a in Integer.
// Output:     a,c in Integer such that a := a-1, c := a ||
{
  INTEGERCONDITION(*this);

  const int sT = sgn;
  if (sT == 0) { abs() = 1; sgn = -1; }
  else if (sT < 0) ++abs();
  else if (::abs(*this) == 1) *this = 0;
  else --abs();

  INTEGERCONDITION(*this);

  return *this;
}

void Integer::add(const Integer& a, const SignDigit b)
// Algorithm:  c.add(a, b)
// Input:      a in Integer, b in SignDigit where not &a = &c.
// Output:     c in Integer such that c = a+b ||
{
  INTEGERCONDITION(a);
  CONDITION(&a != this);

  if (b < 0) sub(a, -b);
  else {
    const int sA = a.sgn;
    sgn = sA;
    if (sA > 0) Natural::add(::abs(a), b);
    else if (sA == 0) *this = b; 
    else if (::abs(a) > Digit(b)) Natural::sub(::abs(a), b);
    else *this = b - a.highest();
  }

  INTEGERCONDITION(*this);
}

void Integer::sub(const Integer& a, const SignDigit b)
// Algorithm:  c.sub(a, b)
// Input:      a,c in Integer, b in SignDigit where not &a = &c.
// Output:     c in Integer such that c = a-b ||
{
  INTEGERCONDITION(a);
  CONDITION(&a != this);

  if (b < 0) add(a, -b);
  else {
    const int sA = a.sgn;
    sgn = sA;
    if (sA < 0) Natural::add(::abs(a), b);
    else if (sA == 0) neg(b);
    else if (::abs(a) > Digit(b)) Natural::sub(::abs(a), b);
    else {
      *this = b - a.highest();
      neg();
    }
  }

  INTEGERCONDITION(*this);
}

void Integer::muladd(const Integer& a, const SignDigit b)
// Algorithm:  c.muladd(a, b)
// Input:      a,c in Integer, b in SignDigit.
// Output:     c in Integer such that c := c + a*b ||
{
  INTEGERCONDITION(*this);
  INTEGERCONDITION(a);

  if (b < 0) mulsub(a, -b);
  else {
    const int sT = sgn;
    if (sT == 0) mul(a, b);
    else if (sT == a.sgn) Natural::muladd(::abs(a), b);
    else if (a.length()+1 < length()) Natural::mulsub(::abs(a), b);
    else {
      Integer c = a*b;
      *this += c;
    }
  }

  INTEGERCONDITION(*this);
}

void Integer::mulsub(const Integer& a, const SignDigit b)
// Algorithm:  c.mulsub(a, b)
// Input:      a,c in Integer, b in SignDigit.
// Output:     c in Integer such that c := c - a*b ||
{
  INTEGERCONDITION(*this);
  INTEGERCONDITION(a);

  if (b < 0) muladd(a, -b);
  else {
    const int sT = sgn;
    if (sT == 0) { mul(a, b); neg(); }
    else if (sT != a.sgn) Natural::muladd(::abs(a), b);
    else if (a.length()+1 < length()) Natural::mulsub(::abs(a), b);
    else {
      Integer c = a*b;
      *this -= c;
    }
  }

  INTEGERCONDITION(*this);
}

void div(const Integer& a, const Integer& b, Integer& c, Integer& d)
// Algorithm:  div(a, b, c, d)
// Input:      a,b,c,d in Integer where not b = 0, not c = d;
// Output:     c,d in Integer such that c = sign(a*b)*[|a/b|], d = a-c*b ||
{
  INTEGERCONDITION(a);
  INTEGERCONDITION(b);
  INTEGER_FOR_CHECK(_a, a);
  INTEGER_FOR_CHECK(_b, b);

  div(abs(a), abs(b), c.abs(), d.abs());
  const int sA = a.sgn;
  const int sB = b.sgn;
  c.sgn = (abs(c) != 0)? sA*sB : 0;
  d.sgn = (abs(d) != 0)? sA : 0;

  CONDITION(d+c*_b == _a);
  INTEGERCONDITION(c);
  INTEGERCONDITION(d);
}

void Integer::div(const Integer& a, const Natural& b)
// Algorithm:  c.div(a, b)
// Input:      a,c in Integer, b in Natural where not b = 0;
// Output:     c in Integer such that c = sign(a)*[|a|/b] ||
{
  INTEGERCONDITION(a);

  Natural t;
  Natural::div(::abs(a), b, t);
  sgn = (::abs(*this) != 0)? a.sgn : 0;

  INTEGERCONDITION(*this);
}

void Integer::div(const Natural& a, const Integer& b)
// Algorithm:  c.div(a, b)
// Input:      a in Natural, b,c in Integer where not b = 0;
// Output:     c in Integer such that c = sign(b)*[a/|b|] ||
{
  INTEGERCONDITION(b);

  Natural t;
  Natural::div(a, ::abs(b), t);
  sgn = (::abs(*this) != 0)? b.sgn : 0;

  INTEGERCONDITION(*this);
}

void div(const Integer& a, const SignDigit b, Integer& c, SignDigit& d)
// Algorithm:  div(a, b, c, d)
// Input:      a in Integer, b in SignDigit where not b = 0.
// Output:     c in Integer, d in SignDigit
//             such that c = sign(a*b)*[|a/b|], d = a - c*b ||
{
  INTEGERCONDITION(a);

  Digit r;
  if (b >= 0) {
    div(abs(a), Digit(b), c.abs(), r);
    c.sgn = (abs(c) != 0)? a.sgn : 0;
  } else {
    div(abs(a), Digit(-b), c.abs(), r);
    c.sgn = (abs(c) != 0)? -a.sgn : 0;
  }
  d = (sign(a) == -1)? -r : r;

  INTEGERCONDITION(c);
}

void swap(Integer& a, Integer& b)
// Algorithm:  swap(a, b)
// Input:      a,b in Integer.
// Output:     a,b in Integer such that t := a, a := b, b := t
//             where t in Integer ||
{
  INTEGERCONDITION(a);
  INTEGERCONDITION(b);

  swap(a.abs(), b.abs());
  int t = a.sgn; a.sgn = b.sgn; b.sgn = t;

  INTEGERCONDITION(a);
  INTEGERCONDITION(b);
}

void sqrt(const Integer& a, Integer& b, Integer& c)
// Algorithm:  sqrt(a, b, c)
// Input:      a in Integer where a >= 0.
// Output:     b,c in Integer such that b = [sqrt(a)], c = a - b^2 ||
{
  INTEGERCONDITION(a);
  INTEGER_FOR_CHECK(_a, a);

  const int sA = a.sgn;
  if (sA == -1) a.errmsg(6, "(sqrt)");
  else if (sA == 0) b = c = 0;
  else {
    sqrt(abs(a), b.abs(), c.abs());
    b.sgn = 1;
    c.sgn = (abs(c) != 0);
  }

  CONDITION(c+b*b == _a);
  INTEGERCONDITION(b);
  INTEGERCONDITION(c);
}

const char* Integer::atoI(const char* a, const Digit b)
// Algorithm:  c := d.atoI(a, b)
// Input:      d in Integer, a in String, b in Digit where 2 <= b <= 36.
// Output:     d in Integer, c in String such that d = a ||
//
// Note: Returns a pointer to the first occurrence of a non-digit character
//       in a.
{
  const bool d = (*a == '-');
  a = atoN(a + d, b);
  sgn = (::abs(*this) != 0);
  if (d) neg();

  INTEGERCONDITION(*this);

  return a;
}

char* Itoa(const Integer& a, char* b, const Digit c)
// Algorithm:  c := Itoa(a, c, b)
// Input:      a in Integer, b in Digit, c in String
//             where 2 <= b <= 36, sizeof(c) > BETA*L(a)/log2(b).
// Output:     c in String such that c = a ||
//
// Note:       conversion Integer to string.
{
  INTEGERCONDITION(a);

  if (sign(a) >= 0) Ntoa(abs(a), b, c);
  else { *b = '-'; Ntoa(abs(a), b+1, c); }
  return b;
}

bool Integer::scan(ISTREAM& in)
// Algorithm:  b := a.scan(i)
// Input:      a in Integer, i in istream.
// Output:     a in Integer, i in istream, b in bool ||
//
// Note:       gets Integer a as an internal representation from input stream
//             if b is true.
{
  if (!in.good()) return false;
  char c = 0;
  if (in.get(c) && c != '-') in.putback(c);
  const bool b = Natural::scan(in);
  sgn = (::abs(*this) != 0);
  if (c == '-') neg();
  return b;
}

ISTREAM& operator>>(ISTREAM& in, Integer& a)
// Algorithm:  i := i >> a
// Input:      i in istream.
// Output:     i in istream, a in Integer ||
//
// Note:       gets Integer a from input stream.
{
  INTEGERCONDITION(a);

  if (!in.good()) return in;
  char ch = 0;
  if (in.get(ch) && ch != '-') in.putback(ch);
  in >> a.abs();
  a.sgn = (abs(a) != 0);
  if (ch == '-') a.neg();

  INTEGERCONDITION(a);

  return in;
}

void Integer::bitwise_and(const Integer& a, const Integer& b)
// Algorithm:  c.bitwise_and(a, b)
// Input:      a,b in Integer.
// Output:     c in Integer such that c = a and b ||
{
  INTEGERCONDITION(a);
  INTEGERCONDITION(b);

  const int sA = a.sgn;
  const int sB = b.sgn;
  if (sA >= 0) {
    if (sB >= 0) {
      abs() = ::abs(a) & ::abs(b);
      sgn = (::abs(*this) != 0);
    } else {
      const Natural c = ::abs(a);
      Natural::bitwise_not(::abs(b));
      ++abs();
      const size_t sC = c.length();
      const size_t sT = length();
      abs() &= c;
      if (sT < sC) { copy(c, sC-sT); sgn = 1; }
      else sgn = (::abs(*this) != 0);
    }
  } else if (sB < 0) {
    Natural c = b;
    --c;
    *this = a;
    --abs(); abs() |= c; ++abs();
  } else {
    const Natural c = ::abs(b);
    Natural::bitwise_not(::abs(a));
    ++abs();
    const size_t sC = c.length();
    const size_t sT = length();
    abs() &= c;
    if (sT < sC) { copy(c, sC-sT); sgn = 1; }
    else sgn = (::abs(*this) != 0);
  }

  INTEGERCONDITION(*this);
}

void Integer::bitwise_or(const Integer& a, const Integer& b)
// Algorithm:  c.bitwise_or(a, b)
// Input:      a,b in Integer.
// Output:     c in Integer such that c = a or b ||
{
  INTEGERCONDITION(a);
  INTEGERCONDITION(b);

  const int sA = a.sgn;
  const int sB = b.sgn;
  if (sA >= 0) {
    if (sB >= 0) {
      abs() = ::abs(a) | ::abs(b);
      sgn = (::abs(*this) != 0);
    } else {
      Natural c = b;
      --c;
      Natural::bitwise_not(::abs(a));
      const size_t sC = c.length();
      const size_t sT = length();
      abs() &= c;
      if (sT < sC) { copy(c, sC-sT); sgn = -1; }
      else { const int b = (::abs(*this) != 0); sgn = -b; }
      --(*this);
    }
  } else if (sB < 0) {
    Natural c = b;
    --c;
    *this = a;
    --abs(); abs() &= c; ++abs();
  } else {
    Natural c = a;
    --c;
    Natural::bitwise_not(::abs(b));
    const size_t sC = c.length();
    const size_t sT = length();
    abs() &= c;
    if (sT < sC) { copy(c, sC-sT); sgn = -1; }
    else { const int b = (::abs(*this) != 0); sgn = -b; }
    --(*this);
  }

  INTEGERCONDITION(*this);
}

void Integer::bitwise_xor(const Integer& a, const Integer& b)
// Algorithm:  c.bitwise_xor(a, b)
// Input:      a,b in Integer.
// Output:     c in Integer such that c = a xor b ||
{
  INTEGERCONDITION(a);
  INTEGERCONDITION(b);

  if (a.sgn >= 0 && b.sgn >= 0) {
    abs() = ::abs(a) ^ ::abs(b);
    sgn = (::abs(*this) != 0);
  } else {
    Integer s,t;
    s.bitwise_not(a); t.bitwise_not(b);
    s &= b; t &= a;
    bitwise_or(s, t);
  }

  INTEGERCONDITION(*this);
}

void Integer::setbit(const size_t a)
// Algorithm:  c.setbit(a)
// Input:      a in size_t, c in Integer.
// Output:     c in Integer such that c := c or 2^a ||
{
  INTEGERCONDITION(*this);
  INTEGER_FOR_CHECK(_t, *this);

  if (sgn >= 0) Natural::setbit(a);
  else {
    const size_t b = a/BETA;
    if (b <= length()) {
      Digit* pT;
      const size_t c = trailing_zeros(pT);
      if (c < b) Natural::clearbit(a);
      else if (c == b) *pT = ((*pT-1) & ~(Digit(1)<<(a%BETA)))+1;
      else {
        pT += c-b;                    // not very good!
        *pT = ~(Digit(1) << (a%BETA))+1;
        dec(pT);
      }
    }
  }

  CONDITION(*this == (_t | Integer(1) << a));
  INTEGERCONDITION(*this);
}

void Integer::clearbit(const size_t a)
// Algorithm:  c.clearbit(a)
// Input:      a in size_t, c in Integer.
// Output:     c in Integer such that c := c and not(2^a) ||
{
  INTEGERCONDITION(*this);
  INTEGER_FOR_CHECK(_t, *this);

  if (sgn >= 0) Natural::clearbit(a);
  else {
    const size_t b = a/BETA;
    Digit* pT;
    const size_t c = trailing_zeros(pT);
    if (c < b) Natural::setbit(a);
    else if (c == b) *pT = ((*pT-1) | (Digit(1)<<(a%BETA)))+1;
  }

  CONDITION(*this == (_t & ~(Integer(1) << a)));           
  INTEGERCONDITION(*this);
}

bool Integer::testbit(const size_t a) const
// Algorithm:  c := b.testbit(a)
// Input:      a in size_t, b in Integer.
// Output:     c in bool such that if b and 2^a then c = true else c = false ||
{
  INTEGERCONDITION(*this);

  if (sgn >= 0) return Natural::testbit(a);
  const size_t b = a/BETA;
  Digit* pT;
  const size_t c = trailing_zeros(pT);
  if (c < b) return !Natural::testbit(a);
  else if (c > b) return false;
  else return (((*pT-1) & (Digit(1)<<(a%BETA))) == 0);
}

void Integer::rand(const size_t n)
// Algorithm:  a.rand(n)
// Input:      n in size_t.
// Output:     a in Integer such that |a| < 2^n (random number) ||
{
  INTEGERCONDITION(*this);

  Natural::rand(n);
  if (::abs(*this) == 0) sgn = 0;
  else sgn = (::rand()&1)? 1 : -1;

  INTEGERCONDITION(*this);
}

Integer pow(const Integer& a, const SignDigit b)
// Algorithm:  c := pow(a, b)
// Input:      a in Integer, b in SignDigit.
// Output:     c in Integer such that c = a^b ||
{
  INTEGERCONDITION(a);

  if (b < 0 && abs(a) != 1) return SignDigit(0);
  Integer c(pow(abs(a), Digit(b)));
  c.sgn = (abs(c) != 0);
  if (a.sgn == -1 && (b&1)) c.neg();

  INTEGERCONDITION(c);

  return c;
}

Integer pow(const Integer& a, const Integer& b)
// Algorithm:  c := pow(a, b)
// Input:      a,b in Integer.
// Output:     c in Integer such that c = sign(a)^b*[|a|^b] ||
{
  INTEGERCONDITION(a);
  INTEGERCONDITION(b);

  if (b < 0 && abs(a) != 1) return SignDigit(0);
  Integer c;
  c.abs() = pow(abs(a), abs(b));
  c.sgn = (abs(c) != 0);
  if (a.sgn == -1 && (b&1)) c.neg();

  INTEGERCONDITION(c);

  return c;
}

Integer root(const Integer& a, const SignDigit b)
// Algorithm:  c := root(a, b)
// Input:      a in Integer, b in SignDigit where 2|b or a >= 0.
// Output:     c in Integer such that c = sign(a)^b*[|a|^(1/b)] ||
{
  INTEGERCONDITION(a);

  if (b < 0 && abs(a) != 1) return SignDigit(0);
  const int sA = sign(a);
  if ((b&1) && sA < 0) return -Integer(root(abs(a), Digit(b)));
  if (sA < 0) a.errmsg(6, "(root)");
  return root(abs(a), Digit(b));
}

void gcd(Integer a, Integer b, Integer& x, Integer& y, Integer& z)
// Algorithm:  gcd(a, b, x, y, z)
// Input:      a,b in Integer.
// Output:     x,y,z in Integer such that z = a*x + b*y = gcd(a, b) ||
{
  if (&x == &y || &x == &z || &y == &z) z.errmsg(5, "(gcd)");
  x = 0; y = 1;
  if (b == 0) z = a;
  else {
    const int as = sign(a);
    const int bs = sign(b);
    if (as == -1) a = -a;
    if (bs == -1) b = -b;
    Integer u = 1;
    Integer v = SignDigit(0);
    Integer t;

    while (true) {
      div(a, b, z, a);
      if (a == 0) { z = b; break; }
      t = x*z; u -= t;
      t = y*z; v -= t;

      div(b, a, z, b);
      if (b == 0) { z = a; x = u; y = v; break; }
      t = u*z; x -= t;
      t = v*z; y -= t;
    }
    if (as == -1) x = -x;
    if (bs == -1) y = -y;
  }
}

