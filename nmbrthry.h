/////////////////////////////////
//
// Piologie V 1.3.2
// multi-precision arithmetic
// Number Theory Package
//
// (c) 1996-2001 HiPiLib
//               www.hipilib.de
//
// Sebastian Wedeniwski
// 10/21/2001
//

#ifndef _Include_NumberTheory_H_
#define _Include_NumberTheory_H_

#include "natural.h"
#ifndef _Old_STD_
# include "modulo.h"
#endif

Natural factorial(const Digit a);
Natural binomial(const Digit a, const Digit b);
Natural fibonacci(Digit n);
void    pell(const Digit x, Natural& u, Natural& v);

#ifndef _Old_STD_
bool    MillerRabin(unsigned int i, const Natural& n);
bool    isprime(const Natural& n);
Natural euler(Natural a);
#endif


class Primes : public NumberBase {
private:
  static unsigned char* primes;
  unsigned char* primPos;
  Digit          primNumber;
  unsigned char  idx;

public:
  Primes();
  ~Primes();

  void   destroy();
  Digit  firstPrime();
  bool   nextPrime(Digit&);
  Digit  lastPrime() const;
  size_t numberOfPrimes(Digit) const;
};

inline Primes::~Primes()
{
}

inline Digit Primes::firstPrime()
// Algorithm:  c := n.firstPrime()
// Input:      n in Primes.
// Output:     c in Digit such that c = 2 ||
{
  primPos = primes;
  idx = 1;
  return primNumber = 2;
}

#ifndef _Old_STD_
template<class T>
bool spsp(const T& base, const T& n)
// Algorithm:  c := spsp(b, n)
// Input:      b,n in T where n == 1 (mod 2).
// Output:     c in bool
//             such that if b^((n-1)/2^k) == 1 (mod n) or b^((n-1)/2^l) == -1 (mod n) 
//             for 0 <= l < k, where 2^k|n-1 and not 2^(k+1)|n-1,
//             then c = true else c = false ||
{
  const T m = n-1;
  T l = m;
  int k = 0;
  do { ++k; l >>= 1; } while ((l&1) == 0);
  l = pow(base, l, n);
  if (l != 1 && l != m) {
    int j = k-1;
    if (j == 0) return false;
    do {
      l *= l; l %= n;
      if (l == 1) return false;
    } while (l != m && --j);
    if (j == 0) return false;
  }
  return true;
}

template <class Container>
void factoring(Natural a, Container& p)
// Algorithm:  factoring(a, c)
// Input:      a in Natural.
// Output:     c is a container over Natural with output-iterator
//             such that if a >= 1 then a = prod_{c.begin() <= i < c.end()} *i
//             else c.begin() = c.end(),
//             for all i in [c.begin(), c.end()[, j in [i, c.end()[ : *i <= *j ||
{
  p.erase(p.begin(), p.end());
  if (a == 0) return;
  while (a.even()) {                      // less factors
    p.push_back(2);
    a >>= 1;
  }
  if (a == 1) return;
  Natural t;
  Digit i;
  Primes primes;
  Digit prim = primes.firstPrime();            // didn't need 2
  while (primes.nextPrime(prim)) {
    while (true) {
      div(a, prim, t, i);
      if (i) break;
      p.push_back(prim);
      if (t == 1) return;
      swap(a, t);
    }
    if (t < prim) { p.push_back(a); return; }
  }
  
  if (isprime(a)) { p.push_back(a); return; }   // greater
  Natural s = prim;
  const Digit n[8] = { 4, 2, 4, 2, 4, 6, 2, 6 };
  i = 0;
  switch (s % 30) {
    case 1:  ++i;
    case 29: ++i;
    case 23: ++i;
    case 19: ++i;
    case 17: ++i;
    case 13: ++i;
    case 11: ++i;
  }
  Natural q,r;
  t = root(a, 3);
  for (s += n[i]; s.length() == 1 && s <= t; s += n[i&7]) {
    if (a%s.highest() == 0) {
      p.push_back(s);
      a /= s.highest();
      if (isprime(a)) { p.push_back(a); return; }
      t = root(a, 3);
    }
    ++i;
  }
  while (s <= t) {
    div(a, s, q, r);
    if (r == 0) {
      p.push_back(s);
      if (isprime(q)) { p.push_back(q); return; }
      swap(a, q);
      t = root(a, 3);
    }
    ++i;
    s += n[i&7];
  }
  
  Natural w;                                        // large factors
  sqrt(a, s, w);
  if (w == 0) { p.push_back(s); p.push_back(s); return; }
  s = root(a, 6);
  q = a << 2;
  t *= a; t <<= 2;

  Natural x,y,z,d = 4;
  Natural e = 2;
  r = s >> 2;
  const char c[32] = { 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
                       1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 };
  Natural k = q;
  for (; k <= t; k += q) {
    if (e == 0) {
      d += 4;
      if (d > s) break;
      r = d >> 1;
      div(s, d, r, x);
    } else --e;
    sqrt(k, x, y);
    z = x << 1; ++x; ++z; z -= y;
    while (true) {
      if (c[z.lowest() & 31]) {
        sqrt(z, y, w);
        if (w == 0) {
          y = gcd(x+y, a);
          a /= y;
          if (y < a) { p.push_back(y); p.push_back(a); }
          else { p.push_back(a); p.push_back(y); }
          return;
        }
      }
      if (r == 0) break;
      --r; z = x << 1; ++x; ++z;
    }
  } 
  while (k <= t) {
    sqrt(k, x, z);
    y = x << 1; ++x; ++y; y -= z;
    if (c[y.lowest() & 31]) {
      sqrt(y, y, w);
      if (w == 0) {
        y = gcd(x+y, a);
        a /= y;
        if (y < a) { p.push_back(y); p.push_back(a); }
        else { p.push_back(a); p.push_back(y); }
        return;
      }
    } 
    k += q;
  }
  p.push_back(a);
}

template<class T>
int jacobi(T a, T b)
// Algorithm:  c := jacobi(a, b)
// Input:      a,b in T where a,b > 0 and not 2 | b.
// Output:     c in {-1,0,1} such that c = (a|b) ||
{
  if (a == 1) return 1;
  if (a == 2) {
    const int c = b&7;
    return (c == 1 || c == 7)? 1 : -1;
  }
  int h,c = 1;
  while (a != 0) {
    for (h = 0; (a&1) == 0; ++h) a >>= 1;
    int d = b&7;
    if (h&1 && d != 1 && d != 7) c = -c;
    if (a < b) {
      if ((d&3) == 3 && (a&3) == 3) c = -c;
      d = a&7;
      swap(a, b);
    }
    a -= b; a >>= 1;
    if (d == 3 || d == 5) c = -c;
  }
  return (b == 1)? c : 0;
}

#endif

#endif
