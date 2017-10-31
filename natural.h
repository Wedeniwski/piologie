/////////////////////////////////
//
// Piologie V 1.3.4
// multi-precision arithmetic
// Natural
//
// (c) 1996-2001 HiPiLib
//               www.hipilib.de
//
// Sebastian Wedeniwski
// 10/21/2001
//

#ifndef _Include_Natural_H_
#define _Include_Natural_H_

#include "digit.h"

#ifdef setbit       // some "math.h" libraries have a defined setbit
# undef setbit
#endif

struct Natural_plus_tag {};
struct Natural_minus_tag {};
struct Natural_multiplies_tag {};
struct Natural_divides_tag {};
struct Natural_modulus_tag {};
struct Natural_square_root_tag {};
struct Natural_lshift_tag {};
struct Natural_rshift_tag {};
struct Natural_and_tag {};
struct Natural_or_tag {};
struct Natural_xor_tag {};
struct Natural_not_tag {};


class Natural : public NumberBase {
public:
  struct rep {
    const size_t size;
    const Digit* p;
    const bool   bin;
    rep(const size_t a, const Digit* b, const bool c = false)
      : size(a), p(b), bin(c) {}
  };
  
private:
  size_t  size;
  Digit*  root;
  Digit*  p;
  
  static size_t NaturalSize;
  static size_t NaturalSizeOld;
  

  void          inc(Digit*);
  void          add_with_inc(const Digit*, Digit*, const Digit*);
  void          add_with_inc(const Digit*, Digit*, const Digit*, const Digit*);
  bool          add_no_inc(const Digit*, Digit*, const Digit*) const;
  bool          add_no_inc(const Digit*, Digit*, const Digit*, const Digit*) const;
  void          sub(Digit*, Digit*, const Digit*);
  void          sub(const Digit*, Digit*, const Digit*, const Digit*);
  bool          sub_no_dec(const Digit*, Digit*, const Digit*) const;
  int           abs(Digit*, const Digit*, const Digit*, size_t) const;

  void          sqr(const Digit*, Digit*, size_t) const;
  Digit         mul(const Digit*, const Digit*, Digit*, const Digit) const;
  Digit*        muladd(const Digit*, const Digit*, Digit*, const Digit) const;
  void          mul(const Digit*, const Digit*, const size_t, Digit*) const;
  void          cmul(const Digit*, size_t, const Digit*, const size_t, Digit*) const;
  void          mul(const Digit*, size_t, const Digit*, const size_t, Digit*) const;
  Digit         mod_div(const Digit);

protected:
  void          dec(Digit*);      // bad! because Integer::setbit
  void          get_memory(const size_t);

  Natural(char, char);
  Natural(const size_t, char);
  Natural(const Digit, size_t);
  Natural(size_t, const Natural&);
  Natural(const Natural&, size_t);

  const Digit*  first() const;
  Digit*        last() const;
  size_t        rootsize() const;
  void          normalize();
  size_t        trailing_zeros(Digit*&) const;
  Digit*        setsize(const size_t);
  Natural&      copy(const Natural&, const size_t);
  void          enlarge(const size_t);
  void          fast_rshift(const size_t);
  void          fast_append(const size_t);

  int           compare(const Natural&) const;
  void          add(const Natural&, const Natural&);
  void          sub(const Natural&, const Natural&);
  void          mul(const Natural&, const Natural&);
  void          sqr(const Natural&);
  void          div(const Natural&, Natural, Natural&);
  void          sqrt(const Natural&);
  void          newton_sqrt(Natural);

  void          lshift(const Natural&, size_t);
  void          rshift(const Natural&, size_t);

  void          bitwise_and(const Natural&, const Natural&);
  void          bitwise_or(const Natural&, const Natural&);
  void          bitwise_xor(const Natural&, const Natural&);
  void          bitwise_not(const Natural&);

  void          add(const Natural&, const Digit);
  void          sub(const Natural&, const Digit);
  void          mul(const Natural&, const Digit);
  void          muladd(const Natural&, const Digit);
  void          mulsub(const Natural&, const Digit);


public:
  static size_t NumberOfDecimals(const size_t);
  static size_t NumberOfDigits(const size_t);
  static void   RestoreSize();

  Natural(const Digit = 0);                   // default constructor
  Natural(const Natural&);                    // copy constructor
  #ifndef _Old_STD_
  explicit
  #endif
  Natural(const char*, const Digit = 10);     // constructor for string conversion
  ~Natural();                                 // destructor
  
  inline Natural& operator=(const Natural&);
  Natural&        operator+=(const Natural&);
  Natural&        operator-=(const Natural&);
  Natural&        operator*=(const Natural&);
  Natural&        operator/=(const Natural&);
  Natural&        operator%=(const Natural&);
  Natural&        operator&=(const Natural&);
  Natural&        operator|=(const Natural&);
  Natural&        operator^=(const Natural&);

  Digit           operator=(const Digit);
  Natural&        operator+=(const Digit);
  Natural&        operator-=(const Digit);
  Natural&        operator*=(const Digit);
  Natural&        operator/=(const Digit);
  Digit           operator%=(const Digit);
  Digit           operator&=(const Digit);
  Natural&        operator|=(const Digit);
  Natural&        operator^=(const Digit);
  Natural&        operator>>=(size_t);
  Natural&        operator<<=(size_t);
  void            lmove(size_t);
  void            rmove(size_t);

  Natural&        operator++();
  Natural&        operator--();
  const Natural   operator++(int);
  const Natural   operator--(int);


  Natural(const binder_arguments<Natural, Natural, Natural_plus_tag>&);
  Natural& operator=(const binder_arguments<Natural, Natural, Natural_plus_tag>&);

  Natural(const binder_arguments<Natural, Natural, Natural_minus_tag>&);
  Natural& operator=(const binder_arguments<Natural, Natural, Natural_minus_tag>&);

  Natural(const binder_arguments<Natural, Natural, Natural_multiplies_tag>&);
  Natural& operator=(const binder_arguments<Natural, Natural, Natural_multiplies_tag>&);

  Natural(const binder_arguments<Natural, Natural, Natural_divides_tag>&);
  Natural& operator=(const binder_arguments<Natural, Natural, Natural_divides_tag>&);

  Natural(const binder_arguments<Natural, Natural, Natural_modulus_tag>&);
  Natural& operator=(const binder_arguments<Natural, Natural, Natural_modulus_tag>&);


  Natural(const binder_arguments<Natural, size_t, Natural_lshift_tag>&);
  Natural& operator=(const binder_arguments<Natural, size_t, Natural_lshift_tag>&);

  Natural(const binder_arguments<Natural, size_t, Natural_rshift_tag>&);
  Natural& operator=(const binder_arguments<Natural, size_t, Natural_rshift_tag>&);

  Natural(const binder_arguments<Natural, Natural, Natural_and_tag>&);
  Natural& operator=(const binder_arguments<Natural, Natural, Natural_and_tag>&);

  Natural(const binder_arguments<Natural, Natural, Natural_or_tag>&);
  Natural& operator=(const binder_arguments<Natural, Natural, Natural_or_tag>&);

  Natural(const binder_arguments<Natural, Natural, Natural_xor_tag>&);
  Natural& operator=(const binder_arguments<Natural, Natural, Natural_xor_tag>&);

  Natural(const binder_arguments<Natural, Natural, Natural_not_tag>&);
  Natural& operator=(const binder_arguments<Natural, Natural, Natural_not_tag>&);


  Natural(const binder_arguments<Natural, Digit, Natural_minus_tag>&);
  Natural& operator=(const binder_arguments<Natural, Digit, Natural_minus_tag>&);

  Natural(const binder_arguments<Natural, Digit, Natural_plus_tag>&);
  Natural& operator=(const binder_arguments<Natural, Digit, Natural_plus_tag>&);

  Natural(const binder_arguments<Natural, Digit, Natural_multiplies_tag>&);
  Natural& operator=(const binder_arguments<Natural, Digit, Natural_multiplies_tag>&);
  Natural& operator+=(const binder_arguments<Natural, Digit, Natural_multiplies_tag>&);
  Natural& operator-=(const binder_arguments<Natural, Digit, Natural_multiplies_tag>&);


  friend Digit        operator%(const Natural&, const Digit);

  inline friend bool  operator==(const Natural&, const Natural&);
  inline friend bool  operator!=(const Natural&, const Natural&);
  inline friend bool  operator<(const Natural&, const Natural&);
  inline friend bool  operator<=(const Natural&, const Natural&);
  inline friend bool  operator>(const Natural&, const Natural&);
  inline friend bool  operator>=(const Natural&, const Natural&);

  inline friend bool  operator==(const Natural&, const Digit);
  inline friend bool  operator!=(const Natural&, const Digit);
  inline friend bool  operator<(const Natural&, const Digit);
  inline friend bool  operator<=(const Natural&, const Digit);
  inline friend bool  operator>(const Natural&, const Digit);
  inline friend bool  operator>=(const Natural&, const Digit);


  Natural(const binder_arguments<Natural, Natural, Natural_square_root_tag>&);
  Natural& operator=(const binder_arguments<Natural, Natural, Natural_square_root_tag>&);


  friend void         div(const Natural&, const Natural&, Natural&, Natural&);
  friend void         div(const Natural&, const Digit, Natural&, Digit&);
  friend void         sqrt(const Natural&, Natural&, Natural&);
  friend void         swap(Natural&, Natural&);
  inline friend bool  sign(const Natural&);


  void    split(const size_t, Natural&, Natural&) const;
  Digit   highest() const;
  Digit   lowest() const;
  size_t  length() const;
  void    setbit(const Digit);
  void    clearbit(const Digit);
  bool    testbit(const Digit) const;
  bool    odd() const;
  bool    even() const;
  void    rand(size_t);
  bool    scan(ISTREAM&);

  const char* atoN(const char*, const Digit = 10);

  inline friend rep print(const Natural&, bool);
  inline friend rep print(const Natural&);

  friend Natural    abs(const Natural&, const Natural&);
  friend int        abs(const Natural&, const Natural&, Natural&);
  friend Natural    gcd(Natural, Natural);

  friend class   FFT;
};

inline binder_arguments<Natural, Natural, Natural_plus_tag>
 operator+(const Natural&, const Natural&);
inline binder_arguments<Natural, Natural, Natural_minus_tag>
 operator-(const Natural&, const Natural&);
inline binder_arguments<Natural, Natural, Natural_multiplies_tag>
 operator*(const Natural&, const Natural&);
inline binder_arguments<Natural, Natural, Natural_divides_tag>
 operator/(const Natural&, const Natural&);
inline binder_arguments<Natural, Natural, Natural_modulus_tag>
 operator%(const Natural&, const Natural&);
inline binder_arguments<Natural, size_t, Natural_lshift_tag>
 operator<<(const Natural&, const size_t&);
inline binder_arguments<Natural, size_t, Natural_rshift_tag>
 operator>>(const Natural&, const size_t&);
inline binder_arguments<Natural, Natural, Natural_and_tag>
 operator&(const Natural&, const Natural&);
inline binder_arguments<Natural, Natural, Natural_or_tag>
 operator|(const Natural&, const Natural&);
inline binder_arguments<Natural, Natural, Natural_xor_tag>
 operator^(const Natural&, const Natural&);
inline binder_arguments<Natural, Natural, Natural_not_tag>
 operator~(const Natural&);
inline binder_arguments<Natural, Digit, Natural_minus_tag>
 operator-(const Natural&, const Digit&);
inline binder_arguments<Natural, Digit, Natural_plus_tag>
 operator+(const Natural&, const Digit&);
inline binder_arguments<Natural, Digit, Natural_plus_tag>
 operator+(const Digit&, const Natural&);
inline binder_arguments<Natural, Digit, Natural_multiplies_tag>
 operator*(const Natural&, const Digit&);
inline binder_arguments<Natural, Digit, Natural_multiplies_tag>
 operator*(const Digit&, const Natural&);
inline Natural  operator/(const Natural&, const Digit);
inline Digit    operator&(const Natural&, const Digit);
inline Natural  operator|(const Natural&, const Digit);
inline binder_arguments<Natural, Natural, Natural_square_root_tag>
 sqrt(const Natural&);

char*           Ntoa(const Natural&, char*, const Digit = 10);
inline Natural  atoN(const char*, const Digit = 10);
OSTREAM&        operator<<(OSTREAM&, const Natural&);
OSTREAM&        operator<<(OSTREAM&, const Natural::rep&);
ISTREAM&        operator>>(ISTREAM&, Natural&);

inline Digit  log2(const Natural&);
Natural       lcm(const Natural&, const Natural&);
Natural       pow(const Natural&, Natural);
Natural       pow(const Natural&, Digit);
Natural       root(const Natural&, const Digit);

Digit         operator%(const Natural&, const Digit);
void          div(const Natural&, const Natural&, Natural&, Natural&);
void          div(const Natural&, const Digit, Natural&, Digit&);
void          sqrt(const Natural&, Natural&, Natural&);
void          swap(Natural&, Natural&);
Natural       abs(const Natural&, const Natural&);
int           abs(const Natural&, const Natural&, Natural&);
Natural       gcd(Natural, Natural);


///////////////////////// Inline-Implementation ////////////////////

inline Natural::Natural(char, char)
{
}

inline bool operator==(const Natural& a, const Natural& b)
// Algorithm:  c := a == b
// Input:      a,b in Natural.
// Output:     c in bool such that if a = b then c = true else c = false ||
{
  return (a.compare(b) == 0);
}

inline bool operator!=(const Natural& a, const Natural& b)
// Algorithm:  c := a != b
// Input:      a,b in Natural.
// Output:     c in bool such that if a = b then c = false else c = true ||
{
  return (a.compare(b) != 0);
}

inline bool operator<(const Natural& a, const Natural& b)
// Algorithm:  c := a < b
// Input:      a,b in Natural.
// Output:     c in bool such that if a < b then c = true else c = false ||
{
  return (a.compare(b) < 0);
}

inline bool operator<=(const Natural& a, const Natural& b)
// Algorithm:  c := a <= b
// Input:      a,b in Natural.
// Output:     c in bool such that if a <= b then c = true else c = false ||
{
  return (a.compare(b) <= 0);
}

inline bool operator>(const Natural& a, const Natural& b)
// Algorithm:  c := a > b
// Input:      a,b in Natural.
// Output:     c in bool such that if a > b then c = true else c = false ||
{
  return (a.compare(b) > 0);
}

inline bool operator>=(const Natural& a, const Natural& b)
// Algorithm:  c := a >= b
// Input:      a,b in Natural.
// Output:     c in bool such that if a >= b then c = true else c = false ||
{
  return (a.compare(b) >= 0);
}

inline bool operator==(const Natural& a, const Digit b)
// Algorithm:  c := a == b
// Input:      a in Natural, b in Digit.
// Output:     c in bool such that if a = b then c = true else c = false ||
{
  return (a.size == 1 && *a.p == b);
}

inline bool operator!=(const Natural& a, const Digit b)
// Algorithm:  c := a != b
// Input:      a in Natural, b in Digit.
// Output:     c in bool such that if not a = b then c = true else c = false ||
{
  return (a.size != 1 || *a.p != b);
}

inline bool operator<(const Natural& a, const Digit b)
// Algorithm:  c := a < b
// Input:      a in Natural, b in Digit.
// Output:     c in bool such that if a < b then c = true else c = false ||
{
  return (a.size == 1 && *a.p < b);
}

inline bool operator<=(const Natural& a, const Digit b)
// Algorithm:  c := a <= b
// Input:      a in Natural, b in Digit.
// Output:     c in bool such that if a <= b then c = true else c = false ||
{
  return (a.size == 1 && *a.p <= b);
}

inline bool operator>(const Natural& a, const Digit b)
// Algorithm:  c := a > b
// Input:      a in Natural, b in Digit.
// Output:     c in bool such that if a > b then c = true else c = false ||
{
  return (a.size > 1 || *a.p > b);
}

inline bool operator>=(const Natural& a, const Digit b)
// Algorithm:  c := a >= b
// Input:      a in Natural, b in Digit.
// Output:     c in bool such that if a >= b then c = true else c = false ||
{
  return (a.size > 1 || *a.p >= b);
}

inline Natural::Natural(const size_t a, char)
// Algorithm:  c := Natural(a, b)
// Input:      b in char, a in size_t where a >= 1.
// Output:     c in Natural such that L(c) = R(c) = a ||
//
// Note:       This constructor don't fulfill the conditions for Naturals.
//
// internal constructor without the initialization of the elements.
{
  get_memory(a);
}

inline bool sign(const Natural& a)
// Algorithm:  c := sign(a)
// Input:      a in Natural.
// Output:     c in bool such that if a = 0 then c = false else c = true ||
{
  return (*a.p != 0);
}
inline Digit Natural::highest() const
// Algorithm:  c := a.highest()
// Input:      a in Natural.
// Output:     c in Digit such that c = [a/2^(BETA*(L(a)-1))] ||
{
  return *p;
}

inline Digit Natural::lowest() const
// Algorithm:  c := a.lowest()
// Input:      a in Natural.
// Output:     c in Digit such that c = a and GAMMA ||
{
  return p[size-1];
}

inline size_t Natural::length() const
// Algorithm:  c := a.length()
// Input:      a in Natural.
// Output:     c in size_t such that c = L(a) ||
{
  return size;
}

inline bool Natural::odd() const
// Algorithm:  c := a.odd()
// Input:      a in Natural.
// Output:     c in bool such that if 2|a then c = false else c = true ||
{
  return ((lowest() & 1) != 0);
}

inline bool Natural::even() const
// Algorithm:  c := a.even()
// Input:      a in Natural.
// Output:     c in bool such that if 2|a then c = true else c = false ||
{
  return ((lowest() & 1) == 0);
}

inline const Digit* Natural::first() const
// Algorithm:  c := a.first()
// Input:      a in Natural.
// Output:     c in [a.p, a.p+L(a)[ such that c = a.p ||
{
  return p;
}

inline Digit* Natural::last() const
// Algorithm:  c := a.last()
// Input:      a in Natural.
// Output:     c in [a.p, a.p+L(a)[ such that c = a.p+L(a)-1 ||
{
  return p+size-1;
}
inline size_t Natural::rootsize() const
// Algorithm:  c := a.rootsize()
// Input:      a in Natural.
// Output:     c in size_t such that c = R(a) ||
{
  return (p+size) - root;
}

inline void Natural::normalize()
// Algorithm:  a.normalize()
// Input:      a in Natural.
// Output:     a in Natural such that not (a.p) = (0) or a.size = 1 ||
{
  Digit* pT = p;
  if (*pT == 0) {
    size_t sT = size;
    if (sT > 2) {
      do { ++pT; --sT; } while (*pT == 0 && sT > 1);
      p = pT; size = sT;
    } else if (sT == 2) { p = ++pT; size = --sT; }
  }
}

inline void Natural::fast_rshift(const size_t a)
// Algorithm:  b.fast_rshift(a)
// Input:      b in Natural, a in size_t where L(b) > a.
// Output:     b in Natural such that b := [b/2^(BETA*a)] ||
{
  size -= a;
}

inline void Natural::fast_append(const size_t a)
// Algorithm:  b.fast_append(a)
// Input:      b in Natural, a in size_t.
// Output:     b in Natural such that L(b) := L(b)+a ||
{
  size += a;
}

inline void Natural::sub(Digit* pT, Digit* pDif, const Digit* pSub)
// Subtraktion *pDif -= *pSub until pDif = pT,
// pDif-pT <= pSub.rootsize().
{
  if (sub_no_dec(pT, pDif, pSub)) dec(pT);
}

inline Digit Natural::mod_div(const Digit b)
// Algorithm:  c := a.mod_div(b)
// Input:      a in Natural, b in Digit where not b = 0.
// Output:     a in Natural, c in Digit such that c = a - [a/b]*b, a := [a/b] ||
{
  Digit c;
  ::div(*this, b, *this, c);
  return c;
}

inline size_t Natural::NumberOfDecimals(const size_t sz)
{
  NaturalSizeOld = NaturalSize;
  return NaturalSize =
         min(size_t(sz/(BETA*0.301029995664))+1, size_t(GAMMA/BETA));
}

inline size_t Natural::NumberOfDigits(const size_t sz)
{
  NaturalSizeOld = NaturalSize;
  return NaturalSize = min(max(sz, size_t(1)), size_t(GAMMA/BETA));
}

inline void Natural::RestoreSize()
{
  NaturalSize = NaturalSizeOld;
}

inline Natural& Natural::operator*=(const Natural& a)
// Algorithm:  c := c *= a
// Input:      a,c in Natural.
// Output:     c in Natural such that c := c*a ||
{
  if (this == &a) sqr(*this);
  else mul(*this, a);
  return *this;
}

inline Natural& Natural::operator/=(const Natural& a)
// Algorithm:  c := c /= a
// Input:      a,c in Natural where not a = 0.
// Output:     c in Natural such that c := [c/a] ||
{
  Natural t(size+DELTA, ' ');
  div(*this, a, t);
  return *this;
}

inline Natural& Natural::operator%=(const Natural& a)
// Algorithm:  c := c %= a
// Input:      a,c in Natural where not a = 0.
// Output:     c in Natural such that c := c - [c/a]*a ||
{
  Natural t(size+DELTA, ' ');
  ::div(*this, a, t, *this);
  return *this;
}

inline Natural& Natural::operator/=(const Digit a)
// Algorithm:  c := c /= a
// Input:      c in Natural, a in Digit where not a = 0.
// Output:     c in Natural such that c := [c/a] ||
{
  mod_div(a);
  return *this;
}

inline Digit Natural::operator%=(const Digit b)
// Algorithm:  c := a %= b
// Input:      a in Natural, b in Digit where not b = 0.
// Output:     c in Digit, a in Natural such that c = a - [a/b]*b, a = c ||
{
  Digit c = (*this)%b;
  *this = c;
  return c;
}

inline Natural::Natural(const binder_arguments<Natural, size_t,
                                               Natural_rshift_tag>& a)
{
  get_memory(a.x.size+DELTA);
  rshift(a.x, a.y);
}

inline Natural& Natural::operator=(const binder_arguments<Natural, size_t,
                                                          Natural_rshift_tag>& a)
{
  if (this == &a.x) return *this >>= a.y;
  else { rshift(a.x, a.y); return *this; }
}

inline binder_arguments<Natural, size_t, Natural_rshift_tag>
 operator>>(const Natural& a, const size_t& b)
// Algorithm:  c := a >> b
// Input:      a in Natural, b in size_t.
// Output:     c in Natural such that c = [a/2^b] ||
{
  return binder_arguments<Natural, size_t, Natural_rshift_tag>(a, b);
}

inline Natural::Natural(const binder_arguments<Natural, size_t,
                                               Natural_lshift_tag>& a)
{
  get_memory(a.x.size+a.y%BETA+DELTA);
  lshift(a.x, a.y);
}

inline Natural& Natural::operator=(const binder_arguments<Natural, size_t,
                                                          Natural_lshift_tag>& a)
{
  if (this == &a.x) return *this <<= a.y;
  else { lshift(a.x, a.y); return *this; }
}

inline binder_arguments<Natural, size_t, Natural_lshift_tag>
 operator<<(const Natural& a, const size_t& b)
// Algorithm:  c := a << b
// Input:      a in Natural, b in size_t.
// Output:     c in Natural such that c = a*2^b ||
{
  return binder_arguments<Natural, size_t, Natural_lshift_tag>(a, b);
}

inline Natural& Natural::operator=(const Natural& a)
// Algorithm:  c := c = a
// Input:      a,c in Natural.
// Output:     c in Natural such that c = a ||
{
  return copy(a, a.size);
}

inline Natural& Natural::operator|=(const Digit a)
// Algorithm:  c := c |= a
// Input:      c in Natural, a in Digit.
// Output:     c in Natural such that c := c or a ||
{
  p[size-1] |= a;
  return *this;
}

inline Digit Natural::operator&=(const Digit b)
// Algorithm:  c := a &= b
// Input:      a in Natural, b in Digit.
// Output:     a in Natural, c in Digit such that a := a and b, c = a ||
{
  return *this = lowest() & b;
}

inline Natural& Natural::operator^=(const Digit a)
// Algorithm:  c := c ^= a
// Input:      c in Natural, a in Digit.
// Output:     c in Natural such that c := c xor a ||
{
  p[size-1] ^= a;
  return *this;
}

inline Natural& Natural::operator++()
// Algorithm:  c := ++a
// Input:      a in Natural.
// Output:     a,c in Natural such that a := a+1, c := a ||
{
  inc(p+size);
  return *this;
}

inline Natural& Natural::operator--()
// Algorithm:  c := --a
// Input:      a in Natural.
// Output:     a,c in Natural such that a := a-1, c := a ||
{
  dec(p+size);
  return *this;
}

inline const Natural Natural::operator++(int)
// Algorithm:  c := a++
// Input:      a in Natural.
// Output:     a,c in Natural such that c := a, a := a+1 ||
{
  const Natural a(*this);
  inc(p+size);
  return a;
}

inline const Natural Natural::operator--(int)
// Algorithm:  c := a--
// Input:      a in Natural.
// Output:     a,c in Natural such that c := a, a := a-1 ||
{
  const Natural a(*this);
  dec(p+size);
  return a;
}

inline Natural::Natural(const binder_arguments<Natural, Natural,
                                               Natural_not_tag>& a)
{
  get_memory(a.x.size+DELTA);
  bitwise_not(a.x);
}

inline Natural& Natural::operator=(const binder_arguments<Natural, Natural,
                                                          Natural_not_tag>& a)
{
  bitwise_not(a.x);
  return *this;
}

inline binder_arguments<Natural, Natural, Natural_not_tag>
 operator~(const Natural& a)
// Algorithm:  c := ~a
// Input:      a in Natural.
// Output:     c in Natural such that c = not a ||
{
  return binder_arguments<Natural, Natural, Natural_not_tag>(a, a);
}

inline Natural::Natural(const binder_arguments<Natural, Natural,
                                               Natural_plus_tag>& a)
{
  get_memory(max(a.x.size, a.y.size)+DELTA);
  add(a.x, a.y);
}

inline Natural& Natural::operator=(const binder_arguments<Natural, Natural,
                                                          Natural_plus_tag>& a)
{
  if (this == &a.x) return *this += a.y;
  else if (this == &a.y) return *this += a.x;
  else { add(a.x, a.y); return *this; }
}

inline binder_arguments<Natural, Natural, Natural_plus_tag>
 operator+(const Natural& a, const Natural& b)
// Algorithm:  c := a+b
// Input:      a,b in Natural.
// Output:     c in Natural such that c = a+b ||
{
  return binder_arguments<Natural, Natural, Natural_plus_tag>(a, b);
}

inline Natural::Natural(const binder_arguments<Natural, Natural,
                                               Natural_minus_tag>& a)
{
  get_memory(max(a.x.size, a.y.size)+DELTA);
  sub(a.x, a.y);
}

inline Natural& Natural::operator=(const binder_arguments<Natural, Natural,
                                                          Natural_minus_tag>& a)
{
  sub(a.x, a.y);
  return *this;
}

inline binder_arguments<Natural, Natural, Natural_minus_tag>
 operator-(const Natural& a, const Natural& b)
// Algorithm:  c := a-b
// Input:      a,b in Natural.
// Output:     c in Natural such that c = a-b ||
{
  return binder_arguments<Natural, Natural, Natural_minus_tag>(a, b);
}

inline Natural::Natural(const binder_arguments<Natural, Natural,
                                               Natural_multiplies_tag>& a)
{
  if (&a.x == &a.y) {
    get_memory(2*a.x.size+DELTA);
    sqr(a.x);
  } else {
    get_memory(a.x.size+a.y.size+DELTA);
    mul(a.x, a.y);
  }
}

inline Natural& Natural::operator=(const binder_arguments<Natural, Natural,
                                                          Natural_multiplies_tag>& a)
{
  if (&a.x == &a.y) sqr(a.x);
  else mul(a.x, a.y);
  return *this;
}

inline binder_arguments<Natural, Natural, Natural_multiplies_tag>
 operator*(const Natural& a, const Natural& b)
// Algorithm:  c := a*b
// Input:      a,b in Natural.
// Output:     c in Natural such that c = a*b ||
{
  return binder_arguments<Natural, Natural, Natural_multiplies_tag>(a, b);
}

inline Natural::Natural(const binder_arguments<Natural, Natural,
                                               Natural_divides_tag>& a)
{
  get_memory(a.x.size+DELTA);
  Natural t(a.y.size+DELTA, ' ');
  div(a.x, a.y, t);
}

inline Natural& Natural::operator=(const binder_arguments<Natural, Natural,
                                                          Natural_divides_tag>& a)
{
  Natural t(a.y.size+DELTA, ' ');
  div(a.x, a.y, t);
  return *this;
}

inline binder_arguments<Natural, Natural, Natural_divides_tag>
 operator/(const Natural& a, const Natural& b)
// Algorithm:  c := a/b
// Input:      a,b in Natural where not b = 0.
// Output:     c in Natural such that c = [a/b] ||
{
  return binder_arguments<Natural, Natural, Natural_divides_tag>(a, b);
}

inline Natural::Natural(const binder_arguments<Natural, Natural,
                                               Natural_modulus_tag>& a)
{
  Natural t(a.x.size+DELTA, ' ');
  get_memory(a.y.size+DELTA);
  ::div(a.x, a.y, t, *this);
}

inline Natural& Natural::operator=(const binder_arguments<Natural, Natural,
                                                          Natural_modulus_tag>& a)
{
  Natural t(a.x.size+DELTA, ' ');
  ::div(a.x, a.y, t, *this);
  return *this;
}

inline binder_arguments<Natural, Natural, Natural_modulus_tag>
 operator%(const Natural& a, const Natural& b)
// Algorithm:  c := a%b
// Input:      a,b in Natural where not b = 0.
// Output:     c in Natural such that c = a - [a/b]*b ||
{
  return binder_arguments<Natural, Natural, Natural_modulus_tag>(a, b);
}

inline Natural::Natural(const binder_arguments<Natural, Natural,
                                               Natural_and_tag>& a)
{
  get_memory(min(a.x.size, a.y.size)+DELTA);
  bitwise_and(a.x, a.y);
}

inline Natural& Natural::operator=(const binder_arguments<Natural, Natural,
                                                          Natural_and_tag>& a)
{
  if (this == &a.x) return *this &= a.y;
  else if (this == &a.y) return *this &= a.x;
  else { bitwise_and(a.x, a.y); return *this; }
}

inline binder_arguments<Natural, Natural, Natural_and_tag>
 operator&(const Natural& a, const Natural& b)
// Algorithm:  c := a & b
// Input:      a,b in Natural.
// Output:     c in Natural such that c = a and b ||
{
  return binder_arguments<Natural, Natural, Natural_and_tag>(a, b);
}

inline Natural::Natural(const binder_arguments<Natural, Natural,
                                               Natural_or_tag>& a)
{
  get_memory(max(a.x.size, a.y.size)+DELTA);
  bitwise_or(a.x, a.y);
}

inline Natural& Natural::operator=(const binder_arguments<Natural, Natural,
                                                          Natural_or_tag>& a)
{
  if (this == &a.x) return *this |= a.y;
  else if (this == &a.y) return *this |= a.x;
  else { bitwise_or(a.x, a.y); return *this; }
}

inline binder_arguments<Natural, Natural, Natural_or_tag>
 operator|(const Natural& a, const Natural& b)
// Algorithm:  c := a | b
// Input:      a,b in Natural.
// Output:     c in Natural such that c = a or b ||
{
  return binder_arguments<Natural, Natural, Natural_or_tag>(a, b);
}

inline Natural::Natural(const binder_arguments<Natural, Natural,
                                               Natural_xor_tag>& a)
{
  get_memory(max(a.x.size, a.y.size)+DELTA);
  bitwise_xor(a.x, a.y);
}

inline Natural& Natural::operator=(const binder_arguments<Natural, Natural,
                                                          Natural_xor_tag>& a)
{
  if (this == &a.x) return *this ^= a.y;
  else if (this == &a.y) return *this ^= a.x;
  else { bitwise_xor(a.x, a.y); return *this; }
}

inline binder_arguments<Natural, Natural, Natural_xor_tag>
 operator^(const Natural& a, const Natural& b)
// Algorithm:  c := a ^ b
// Input:      a,b in Natural.
// Output:     c in Natural such that c = a xor b ||
{
  return binder_arguments<Natural, Natural, Natural_xor_tag>(a, b);
}

inline Natural::Natural(const binder_arguments<Natural, Digit,
                                               Natural_minus_tag>& a)
{
  get_memory(a.x.size+DELTA);
  sub(a.x, a.y);
}

inline Natural& Natural::operator=(const binder_arguments<Natural, Digit,
                                                          Natural_minus_tag>& a)
{
  if (this == &a.x) return *this -= a.y;
  else { sub(a.x, a.y); return *this; }
}

inline binder_arguments<Natural, Digit, Natural_minus_tag>
 operator-(const Natural& a, const Digit& b)
// Algorithm:  c := a - b
// Input:      a in Natural, b in Digit.
// Output:     c in Natural such that c = a-b ||
{
  return binder_arguments<Natural, Digit, Natural_minus_tag>(a, b);
}

inline Natural::Natural(const binder_arguments<Natural, Digit,
                                               Natural_plus_tag>& a)
{
  get_memory(a.x.size+DELTA);
  add(a.x, a.y);
}

inline Natural& Natural::operator=(const binder_arguments<Natural, Digit,
                                                          Natural_plus_tag>& a)
{
  if (this == &a.x) return *this += a.y;
  else { add(a.x, a.y); return *this; }
}

inline binder_arguments<Natural, Digit, Natural_plus_tag>
 operator+(const Natural& a, const Digit& b)
// Algorithm:  c := a + b
// Input:      a in Natural, b in Digit.
// Output:     c in Natural such that c = a+b ||
{
  return binder_arguments<Natural, Digit, Natural_plus_tag>(a, b);
}

inline binder_arguments<Natural, Digit, Natural_plus_tag>
 operator+(const Digit& a, const Natural& b)
// Algorithm:  c := a + b
// Input:      a in Digit, b in Natural.
// Output:     c in Natural such that c = a+b ||
{
  return binder_arguments<Natural, Digit, Natural_plus_tag>(b, a);
}

inline Natural::Natural(const binder_arguments<Natural, Digit,
                                               Natural_multiplies_tag>& a)
{
  get_memory(a.x.size+DELTA);
  mul(a.x, a.y);
}

inline Natural& Natural::operator=(const binder_arguments<Natural, Digit,
                                                          Natural_multiplies_tag>& a)
{
  if (this == &a.x) return *this *= a.y;
  else { mul(a.x, a.y); return *this; }
}

inline Natural& Natural::operator+=(const binder_arguments<Natural, Digit,
                                                           Natural_multiplies_tag>& a)
{
  muladd(a.x, a.y);
  return *this;
}

inline Natural& Natural::operator-=(const binder_arguments<Natural, Digit,
                                                           Natural_multiplies_tag>& a)
{
  mulsub(a.x, a.y);
  return *this;
}

inline binder_arguments<Natural, Digit, Natural_multiplies_tag>
 operator*(const Natural& a, const Digit& b)
// Algorithm:  c := a*b
// Input:      a in Natural, b in Digit.
// Output:     c in Natural such that c = a*b ||
{
  return binder_arguments<Natural, Digit, Natural_multiplies_tag>(a, b);
}

inline binder_arguments<Natural, Digit, Natural_multiplies_tag>
 operator*(const Digit& a, const Natural& b)
// Algorithm:  c := a*b
// Input:      a in Digit, b in Natural.
// Output:     c in Natural such that c = a*b ||
{
  return binder_arguments<Natural, Digit, Natural_multiplies_tag>(b, a);
}

inline Natural operator/(const Natural& a, const Digit b)
// Algorithm:  c := a/b
// Input:      a in Natural, b in Digit where not b = 0.
// Output:     c in Natural such that c = [a/b] ||
{
  const size_t sA = a.length();
  if (sA == 1) return Natural(a.highest() / b);
  else return Natural(a) /= b;
}

inline Digit operator&(const Natural& a, const Digit b)
// Algorithm:  c := a & b
// Input:      a in Natural, b in Digit.
// Output:     c in Digit such that c = a and b ||
{
  return a.lowest() & b;
}

inline Natural operator|(const Natural& a, const Digit b)
// Algorithm:  c := a | b
// Input:      a in Natural, b in Digit.
// Output:     c in Natural such that c = a or b ||
{
  return Natural(a) |= b;
}

inline Digit log2(const Natural& a)
// Algorithm:  b := log2(a)
// Input:      a in Natural.
// Output:     b in Digit
//             such that if a > 0 then b = [log2(a)] else b = 0 ||
{
  return log2(a.highest()) + (a.length()-1)*BETA;
}

inline Natural atoN(const char* a, const Digit b)
// Algorithm:  c := atoN(a, b)
// Input:      a in String, b in Digit where 2 <= b <= 36.
// Output:     c in Natural such that c = a ||
//
// Note:       conversion string to Natural; return 0 by conversion error.
{
  Natural result;
  result.atoN(a, b);
  return result;
}

inline Natural::rep print(const Natural& a, bool b)
// Algorithm:  o := o << print(a, b)
// Input:      o in ostream, a in Natural, b in bool.
// Output:     o in ostream ||
//
// Note:       puts internal representation of Natural a on output stream.
{
  return Natural::rep(a.size, a.p, b);
}

inline Natural::rep print(const Natural& a)
// Algorithm:  o := o << print(a)
// Input:      o in ostream, a in Natural.
// Output:     o in ostream ||
//
// Note:       puts internal representation of Natural a on output stream.
{
  return Natural::rep(a.size, a.p, false);
}

inline Natural::Natural(const binder_arguments<Natural, Natural,
                                               Natural_square_root_tag>& a)
{
  get_memory(a.x.size/2+DELTA);
  sqrt(a.x);
}

inline Natural& Natural::operator=(const binder_arguments<Natural, Natural,
                                                          Natural_square_root_tag>& a)
{
  sqrt(a.x);
  return *this;
}

inline binder_arguments<Natural, Natural, Natural_square_root_tag>
 sqrt(const Natural& a)
// Algorithm:  b := sqrt(a)
// Input:      a in Natural.
// Output:     b in Natural such that b = [sqrt(a)] ||
{
  return binder_arguments<Natural, Natural, Natural_square_root_tag>(a, a);
}

#endif
