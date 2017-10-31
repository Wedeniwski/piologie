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

#include "nmbrthry.h"

#ifndef _Old_STD_
# include <list>
#endif

#ifndef _Old_STD_
using namespace std;
# define NOTHROW_NEW  new(nothrow)
#else
# define NOTHROW_NEW  new
#endif


const size_t     SieveSize = (size_t(1) << ((BETA > 32)? 16 : BETA/2)) - 1;

unsigned char*  Primes::primes     = 0;

Primes::Primes()
{
  if (primes) { firstPrime(); return; }
  // Generating prime numbers
  const Digit n = 8*SieveSize;
  unsigned char* p = primes = NOTHROW_NEW unsigned char[SieveSize];
  if (!primes) errmsg(2, "(constructor)");
  firstPrime();
  const unsigned char* pE = p+SieveSize;
  do *p++ = (unsigned char)~0; while (p != pE);
  p -= SieveSize;
  
  Digit c = GAMMA;
  Digit k = 1;
  Digit q = 2*(sqrt(Digit(24*SieveSize))/3);
  
  for (Digit i = 3; i <= q; i += 2) {
    k += 3;
    Digit j = c += 2*k;
    k += 3;
  
    while (j < n) {
      p[j/8] &= (unsigned char)~(1 << (j%8)); j += i;
      if (j >= n) break;
      p[j/8] &= (unsigned char)~(1 << (j%8)); j += k;
    }
  
    i += 2;
    if (i > q) break;
    j = c += ++k; ++k;
  
    while (j < n) {
      p[j/8] &= (unsigned char)~(1 << (j%8)); j += k;
      if (j >= n) break;
      p[j/8] &= (unsigned char)~(1 << (j%8)); j += i;
    }
  }
  
}

void Primes::destroy()
{
  delete[] primes;
  primPos = primes = 0;
}

bool Primes::nextPrime(Digit& a)
// Algorithm:  c := n.nextPrime(a)
// Input:      n in Primes.
// Output:     a in Digit, c in bool
//             such that if n.primNumber <= a.lastPrime()
//             then c = true else c = false ||
{
  if (primNumber == 2) { a = 3; primNumber = 5; }
  else {
    Digit p;
    do {
      if (primPos - primes == SieveSize) return false;
      p = *primPos & idx;
      if (p) a = primNumber;
      if (idx == 1 << (CHAR_BIT-1)) { idx = 1; ++primPos; primNumber += 4; }
      else {
        if (idx == 2 || idx == 8 || idx == 32) primNumber += 2;
        primNumber += 2;
        idx *= 2;
      }
    } while (p == 0);
  }
  return true;
}

Digit Primes::lastPrime() const
// Algorithm:  c := n.lastPrime()
// Input:      n in Primes.
// Output:     c in Digit such that c is a prime ||
{
  static Digit lprim = 0;
  if (lprim) return lprim;
  Digit i = SieveSize;
  while (primes[--i] == 0);
  const Digit p = log2(Digit(primes[i]));
  i *= 24; i += 3*p + 5;
  return lprim = i - (p&1);
}

size_t Primes::numberOfPrimes(Digit n) const
// Algorithm:  c := n.numberOfPrimes(a)
// Input:      n in Primes, a in Digit.
// Output:     c in Digit
//             such that if a <= n.lastPrime() then c = pi(a) else c = 0 ||
{
  if (n < 5) return size_t(n/2 + (n == 3));
  n -= 5;
  if (n >= 24*SieveSize) return 0;
  const unsigned char* q = primes + n/24;
  n %= 24;
  size_t i = 2;
  for (unsigned char j = char(1 << (n/3 + ((n%6) == 2))); j >= 1; j /= 2)
    if (*q & j) ++i;

  for (unsigned char* p = primes; p < q; ++p) {
    for (unsigned char j = 1; j; j *= 2)
      if (*p & j) ++i;
  }
  return i;
}

Natural factorial(const Digit a)
// Algorithm:  b := factorial(a)
// Input:      a in Digit.
// Output:     b in Natural such that b = a! ||
{
  if (a <= 1) return 1;
  Digit d = a >> 1;
  Digit i = a;
  while (d > 1) { d >>= 1; i += a; }

  Natural::NumberOfDigits(size_t(i)/BETA+1);
  Natural c(1);
  Natural::RestoreSize();
  Primes primes;
  if (a > primes.lastPrime()) {
    d = 1;
    for (i = 2; i <= a; ++i) {
      Digit d1,d2;
      c.digitmul(i, d, d1, d2);
      if (d1) { c *= d; d = i; }
      else d = d2;
    }
    return c *= d;
  }
  const size_t t = primes.numberOfPrimes(a) - 1;
  Digit* e = NOTHROW_NEW Digit[t];
  if (!e) c.errmsg(2, "(factorial)");
  Digit* p = NOTHROW_NEW Digit[t];
  if (!p) c.errmsg(2, "(factorial)");
  primes.firstPrime();
  for (i = 0; i < t; ++i) {
    primes.nextPrime(p[i]);
    e[i] = d = a/p[i];
    while (d >= p[i]) e[i] += d /= p[i];
  }
  for (Digit j = Digit(1) << log2(e[0]);; j >>= 1) {
    d = 1;
    for (i = 0; i < t; ++i)
      if (e[i]&j) {
        Digit d1,d2;
        c.digitmul(d, p[i], d1, d2);
        if (d1) { c *= d; d = p[i]; }
        else d = d2;
      }
    c *= d;
    if (j == 1) break;
    c = c*c;
  }
  d = a;
  for (i = a; i; i >>= 1) d -= i&1;

  delete[] e;
  delete[] p;

  return c <<= size_t(d);
}

Natural binomial(const Digit a, const Digit b)
// Algorithm:  c := binomial(a, b)
// Input:      a,b in Digit.
// Output:     c in Natural such that c = a!/((a-b)!b!) ||
{
  if (b > a) return Digit(0);
  if (a <= 1) return 1;
  Natural c(1);
  Primes primes;
  if (a > primes.lastPrime()) {
    for (Digit i = 1; i <= b; ++i) {
      c *= a-i+1;
      c /= i;
    }
    return c;
  }
  const size_t t = primes.numberOfPrimes(a) - 1;
  Digit* e = NOTHROW_NEW Digit[t];
  if (!e) c.errmsg(2, "(binomial)");
  Digit* p = NOTHROW_NEW Digit[t];
  if (!p) c.errmsg(2, "(binomial)");
  primes.firstPrime();
  Digit d,i,j = 0;
  for (i = 0; i < t; ++i) {
    primes.nextPrime(p[i]);
    e[i] = d = a/p[i];
    while (d >= p[i]) e[i] += d /= p[i];
    e[i] -= d = b/p[i];
    while (d >= p[i]) e[i] -= d /= p[i];
    e[i] -= d = (a-b)/p[i];
    while (d >= p[i]) e[i] -= d /= p[i];
    if (e[i] > j) j = e[i];
  }
  for (j = Digit(1) << log2(j);; j >>= 1) {
    d = 1;
    for (i = 0; i < t; ++i)
      if (e[i]&j) {
        Digit d1,d2;
        c.digitmul(d, p[i], d1, d2);
        if (d1) { c *= d; d = p[i]; }
        else d = d2;
      }
    c *= d;
    if (j == 1) break;
    c = c*c;
  }
  d = 0;
  for (i = b; i; i >>= 1) d += i&1;
  for (i = a-b; i; i >>= 1) d += i&1;
  for (i = a; i; i >>= 1) d -= i&1;

  delete[] e;
  delete[] p;

  return c <<= size_t(d);
}

Natural fibonacci(Digit n)
// Algorithm:  c := fibonacci(n)
// Input:      n in Digit.
// Output:     c in Natural such that c = F_n
//             where F_0 = 0, F_1 = 1, F_k = F_(k-1)+F_(k-2) for k >= 2 ||
{
  if (n <= 1) return n;
  Natural::NumberOfDecimals(size_t((n*209)/1000));
  Natural k,i,j(1);
  
  if (n >= 50) {
    Natural t;
    Natural::RestoreSize();
    Digit a = Digit(1) << log2(n);
    do {
      k = j*j; t = i*i;
      j = k+t;                    // (4.2)
      --k; i = t-k;
      i <<= 1; i += t;            // (4.4)
      while (n&a) {
        j += i; swap(i, j);
        a >>= 1;
        if (a == 0) return i;
        k = j*j; t = i*i;
        j = k+t;                  // (4.2)
        i = t-k;
        --i; i <<= 1; i += t;     // (4.4)
      }
    } while (a >>= 1);
    return i;
  } else {
    Natural::RestoreSize();
    while (--n) {
      k = i+j;
      if (--n == 0) return k;
      i = j+k;
      if (--n == 0) return i;
      j = k+i;
    }
  }
  return j;
}

void pell(const Digit x, Natural& u, Natural& v)
// Algorithm:  pell(a, b, c)
// Input:      a in Digit where not sqrt(a) = [sqrt(a)].
// Output:     b,c in Natural such that b^2 + a * c^2 = 1 ||
{
  short sign = -1;
  Digit q;
  if (x == 5) {
    u = 1; v = 1; q = 4;
  } else {
    Digit w;
    sqrt(x, w, q);
    if (q == 0) { u = v = 0; return; }
    u = w; v = 1;
    Digit q0 = 1;
    Digit m = w;
    Natural u0 = 1;
    Natural v0 = Digit(0);
    while (q != 4 && q != 1) {
      const Digit a = (m+w)/q;
      const Digit m1 = a*q - m;
      const Digit q1 = (m-m1)*a + q0;
      Natural u1 = a*u;
      u1 += u0;
      Natural v1 = a*v;
      v1 += v0;
      m = m1; q0 = q; q = q1;
      u0 = u; u = u1; v0 = v; v = v1;
      sign = -sign;
    }
  }
  if (q == 4) {
    Natural u2 = u*u;
    if (sign < 0) {
      ++u2;
      Natural t = u2 >> 1;
      v *= t;
      u2 += 2;
    } else {
      --u2;
      Natural t = u2 >> 1;
      v *= t;
      u2 -= 2;
    }
    u2 >>= 1;
    u *= u2;
  }
  if (sign < 0) {
    v *= u; v <<= 1;
    u *= u; u <<= 1; ++u;
  }
}

#ifndef _Old_STD_
bool MillerRabin(unsigned int i, const Natural& n)
// Algorithm:  c := MillerRabin(i, n)
// Input:      i in unsigned int, n in Natural where n == 1 (mod 2).
// Output:     c in bool
//             such that if for all primes b less than or equal the ith prime
//             b^((n-1)/2^k) == 1 (mod n) or b^((n-1)/2^l) == -1 (mod n) 
//             for 0 <= l < k, where 2^k|n-1 and not 2^(k+1)|n-1,
//             then c = true else c = false ||
{
  Primes p;
  Digit a = p.firstPrime();
  Natural b;
  do {
    b = a;
    if (!spsp(b, n)) return false;
  } while (p.nextPrime(a) && --i);
  return true;
}

bool isprime(const Natural& n)
// Algorithm:  c := isprime(a)
// Input:      a in Natural.
// Output:     c in bool such that
//             if a is a prime element then c = true else c = false ||
{
  static const bool prim[32] = { 0, 0, 1, 1, 0, 1, 0, 1, 0, 0,
                                 0, 1, 0, 1, 0, 0, 0, 1, 0, 1,
                                 0, 0, 0, 1, 0, 0, 0, 0, 0, 1,
                                 0, 1 };
  static const bool mod[30] = { 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1,
                                1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0 };
  const size_t sT = n.length();
  const Digit z = n.highest();
  if (sT == 1 && z < 32) return prim[z];
  if (mod[n%30]) return false;
  Primes primes;
  if (sT == 1) {
    const Digit y = primes.lastPrime();
    if (z == y) return true;
    else if (z < y) {
      Digit q = primes.firstPrime();
      while (primes.nextPrime(q) && q < z);
      if (q == z) return true;
      else if (q > z) return false;
    }
  }

  static const Natural y0("4DQKT65", 32);
  static const Natural y1("T3AAA400", 32);
  static const Natural y2("1UKFJVRHR", 32);
  static const Natural y3("3543JG76V", 32);
  static const Natural y4("9MKD9B5I61", 32);

  if (n < Digit(9080191)) return (spsp(Natural(31), n) && spsp(Natural(73), n));
  else if (n < y0) return (spsp(Natural(2), n) && spsp(Natural(7), n)
                           && spsp(Natural(61), n));
  else if (n < y1) return (spsp(Natural(2), n) && spsp(Natural(13), n)
                           && spsp(Natural(23), n) && spsp(Natural(1662803), n));
  else if (n < y2) return MillerRabin(5, n);
  else if (n < y3) return MillerRabin(6, n);
  else if (n < y4) return MillerRabin(7, n);

  list<Natural> p;
  Natural b = n-1;
  factoring(b, p);
  Digit k = primes.firstPrime();
  while (true) {
    if (!spsp(Natural(k), n)) return false;
    list<Natural>::iterator i = p.begin();
    while (true) {
      const Natural c = b / *i;
      const Digit k2 = k;
      while (pow(Natural(k), c, n) == 1) {
        if (!primes.nextPrime(k)) k += 2;
      }
      if (k != k2) break;

      list<Natural>::const_iterator j = i;
      do
        if (++i == p.end()) return true;
      while (*j == *i);
    }
  }
}

Natural euler(Natural a)
// Algorithm:  c := euler(a)
// Input:      a in Natural.
// Output:     c in Natural such that c = phi(a) ||
{
  list<Natural> p;
  factoring(a, p);
  Natural c,d;

  list<Natural>::iterator i = p.begin();
  while (i != p.end()) {
    div(a, *i, c, d); a -= c;
    list<Natural>::const_iterator j = i;
    while (++i != p.end() && *j == *i);
  }
  return a;
}

#endif
