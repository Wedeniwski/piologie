/////////////////////////////////
//
// Piologie V 1.3.2
// multi-precision arithmetic
// Calculation of Constants
//
// (c) 1996-2001 HiPiLib
//               www.hipilib.de
//
// Sebastian Wedeniwski
// 10/21/2001
//

#include "pi.h"

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


static void output(OSTREAM& out, const Natural& a, const Natural* c,
                   const size_t m, size_t n)
{
  Natural q,r;
  div(a, *c, q, r);
  if (--n == 0) {
    out.width(m);
    out.fill('0');
    out << q;
    out.width(m);
    out.fill('0');
    out << r;
  } else {
    output(out, q, c+1, m, n);
    output(out, r, c+1, m, n);
  }
}

void Fixed::output(OSTREAM& out)
// Algorithm:  a.output(o)
// Input:      o in ostream, a in Fixed.
// Output:     o in ostream ||
//
// Note:       puts Fixed a on output stream.
{
  size_t n = decimals();
  size_t steps = 1;

  if (n >= 80000000) steps = 4096;
  else if (n >= 40000000) steps = 2048;
  else if (n >= 20000000) steps = 1024;
  else if (n >= 8000000) steps = 512;
  else if (n >= 4000000) steps = 256;
  else if (n >= 2000000) steps = 128;
  else if (n >= 800000) steps = 64;
  else if (n >= 400000) steps = 32;
  else if (n >= 200000) steps = 16;
  else if (n >= 100000) steps = 8;
  else if (n >= 50000) steps = 4;
  else if (n >= 10000) steps = 2;

  const size_t block_size = 100*steps;
  const size_t d          = log2(steps)+1;
  size_t m = precision();

  Natural* c = new Natural[d];
  if (c == 0) errmsg(2, "(output)");
  c[d-1] = pow(Natural(ALPHA), Digit(block_size/steps));
  for (int i = d-2; i >= 0; --i) c[i] = c[i+1]*c[i+1];

  const size_t sB = length();
  enlarge((sB > m)? 3 : 3+m-sB);
  if (sB > m) {
    fast_rshift(m);
    out << value();
    if (n) out << '.';
    value() = 0; fast_append(m); normalize();
  } else out << "0.";
  while (n >= ALPHA_WIDTH*(block_size/steps)) {
    int j = 0;
    do {
      if (n >= ((ALPHA_WIDTH*block_size) >> j)) {
        value() *= c[j];
        out.width((ALPHA_WIDTH*block_size) >> j);
        break;
      }
      ++j;
    } while (j < d);
    out.fill('0');
    if (length() > m) {
      fast_rshift(m);
      if (j == d-1) out << value();
      else ::output(out, value(), c+j+1, ALPHA_WIDTH*(block_size/steps), d-j-1);
      value() = 0; fast_append(m); normalize();
      if (value() != 0) {
        Digit* pT;
        const size_t sz = trailing_zeros(pT);
        fast_rshift(sz); m -= sz;
      }
    } else out << '0';
    n -= (block_size*ALPHA_WIDTH) >> j;
  }
  Digit* pE = last();
  Digit* pT = pE-m;
  while (n) {
    value() *= ALPHA;
    if (n <= ALPHA_WIDTH) {
       out.width(n);
       out.fill('0');
       Digit d = highest();
       for (n = ALPHA_WIDTH-n; n; --n) d /= 10;
       out << d;
       break;
    }
    out.width(ALPHA_WIDTH);
    out.fill('0');
    if (length() > m) {
      out << *pT;
      *pT = 0;
      if (*pE == 0) { --pE; fast_rshift(1); --m; }
      normalize();
    } else out << '0';
    n -= ALPHA_WIDTH;
  }
  delete[] c;
}

void Pi::stoermer(const size_t n, Natural& pi)
// Algorithm:  p.stoermer(n, a)
// Input:      p in Pi, a in Natural, n in size_t
//             where n < 3(GAMMA_LOW+1)/BETA, BETA >= 32.
// Output:     a in Natural such that |a-pi*2^(BETA*n)| < 1 ||
{
  const size_t sz = BETA*n;
  Natural t;
  Natural a = Digit(8*24);
  Natural b = Digit(57*8);
  Natural c = Digit(239*4);
  pi        = 0;
  Digit   k = 0;

  a <<= sz; b <<= sz; c <<= sz;
  do {
    a >>= 12;
    t = a * (63*k+191);
    if (b != 0) {
      b /= Digit(57*57*57*57);
      t += b * Digit(3248*k+9746);
      if (c != 0) {
        c /= Digit(239UL*239*239*239);
        t += c * Digit(57120*k+171362);
      }
    }
    pi += t /= (k+1)*(k+3);
    k += 4;
  } while (t != 0);
}

void Pi::schoenhage(const size_t n, Natural& pi)
// Algorithm:  p.schoenhage(n, a)
// Input:      p in Pi, a in Natural, n in size_t.
// Output:     a in Natural such that |a-pi*2^(BETA*n)| < 1 ||
{
  const size_t m = size_t(log2(BETA*n+BETA)-1.68025927);
  const size_t sz = BETA*n;
 
  Integer s,t,a(1),A(1),b(1),B(1),c(1);
  a <<= sz; A <<= sz; B <<= sz-1; c <<= sz-2;
  for (size_t k = 0; k <= m; ++k) {
    t = A+B; t >>= 2;
    s = B << sz; b = sqrt(s);
    a += b; a >>= 1;
    A = a*a; A >>= sz;
    B = A-t; B <<= 1;
    t = A-B; t <<= k; c -= t;
  }
  A <<= sz;
  div(A, c, t, a);
  pi = abs(t);
}

void Pi::chudnovsky(Digit n, const Digit m, Integer& p, Integer& q, Integer& t) const
// Algorithm:  p.chudnovsky(n, m, p, q, t)
// Input:      p in Pi, p,q,t in Integer, n,m in Digit where n < 2^BETA/6, n < m.
// Output:     p,q,t in Integer ||
{
  CONDITION(n < m && n < ((size_t(1) << (BETA-1))/3));

  const SignDigit A = SignDigit(13591409);
  const SignDigit B = SignDigit(545140134);
  const SignDigit C = SignDigit(640320);
  switch (m-n) {
    case 1: {
      if (n <= GAMMA_LOW/8) p = (6*n-5)*(2*n-1);
      else { p = 6*n-5; p *= 2*n-1; }
      p *= 6*n-1; p = -p;
      if (n <= GAMMA_LOW/2) q = n*n;
      else { q = n; q *= n; }
      q *= n; q *= 230*C; q *= 116*C;
      t = B; t *= n; t += A; t *= p;
      
      break;
    }
    case 2: {
      Integer s1,s2;
      if (n <= GAMMA_LOW/8) {
        p = (6*n-5)*(2*n-1);
        s1 = p *= 6*n-1; s1 = -s1;
        p *= 6*n+5; p *= (2*n+1)*(6*n+1);
      } else {
        p = 6*n-5; p *= 2*n-1;
        s1 = p *= 6*n-1; s1 = -s1;
        p *= 6*n+5; p *= 2*n+1; p *= 6*n+1;
      }
      if (n < GAMMA_LOW/2) {
        const size_t k = n+1;
        q = k*k; q *= k; q *= 230*C;
        s2 = q *= 116*C; q *= n*n;
      } else {
        const size_t k = n+1;
        q = k; q *= k; q *= k; q *= 230*C;
        s2 = q *= 116*C; q *= n; q *= n;
      }
      q *= n; q *= 230*C; q *= 116*C;
      t = B; t *= n; t += A;
      const Integer s3 = t;
      t += B;
      const Integer s4 = t*p;
      t = s1*s2*s3;
      t += s4;
      
      break;
    }
    case 3: {
      Integer s1,s2;
      if (n <= GAMMA_LOW/12) {
        p = (6*n+1)*(2*n-1); p *= (6*n-1)*(6*n-5);
        s1 = p *= (2*n+1)*(6*n+5);
        p *= (6*n+7)*(2*n+3);
      } else {
        p = 6*n-5; p *= 2*n-1; p *= 6*n-1;
        p *= 6*n+5; p *= 2*n+1;
        s1 = p *= 6*n+1;
        p *= 6*n+7; p *= 2*n+3;
      }
      p *= 6*n+11; p = -p;
      if (n < GAMMA_LOW/2-1) {
        const size_t k = n+1;
        q = k*k; q *= k; q *= 230*C; q *= 116*C;
        const size_t k2 = n+2;
        q *= k2*k2; q *= k2; q *= 230*C;
        s2 = q *= 116*C; q *= n*n;
      } else {
        const size_t k = n+1;
        q = k; q *= k; q *= k; q *= 230*C; q *= 116*C;
        const size_t k2 = n+2;
        q *= k2; q *= k2; q *= k2; q *= 230*C;
        s2 = q *= 116*C; q *= n; q *= n;
      }
      q *= n; q *= 230*C; q *= 116*C;
      t = B; t *= n; t += A;
      const Integer s3 = t;
      t += B;
      const Integer s4 = t;
      t += B;
      const Integer s5 = t*p;
      t = s2*s3;
      if (n <= GAMMA_LOW/8) t *= (6*n-5)*(2*n-1);
      else { t *= 6*n-5; t *= 2*n-1; }
      t *= 6*n-1; t = -t;
      s2 = s1*s4;
      n += 2;
      if (n <= GAMMA_LOW/2) s2 *= n*n;
      else { s2 *= n; s2 *= n; }
      s2 *= n; s2 *= 230*C; s2 *= 116*C;
      t += s2; t += s5;
      
      break;
    }
    default: {
      Integer p1,q1,t1;
      Integer p2,q2,t2;

      const size_t l = (n+m)/2;
      chudnovsky(n, l, p1, q1, t1);
      chudnovsky(l, m, p2, q2, t2);
      t = t1*q2; t += t1 = t2*p1;
      p = p1*p2; q = q1*q2;
    }
  }
}
void Pi::chudnovsky2(const Digit n, const Digit m, Integer& q, Integer& t) const
// Algorithm:  p.chudnovsky(n, m, q, t)
// Input:      p in Pi, q,t in Integer, n,m in Digit where n < m.
// Output:     q,t in Integer ||
{
  CONDITION(n < m);

  Integer p1,q1,t1;
  Integer q2,t2;

  const Digit l = (n+m)/2;
  chudnovsky(n, l, p1, q1, t1);
  if (m-l <= 3) {
    Integer p2;
    chudnovsky(l, m, p2, q2, t2);
  } else chudnovsky2(l, m, q2, t2);
  t = t1*q2; t += t1 = t2*p1;
  q = q1*q2;
}

void Pi::sqrt_series(Digit n, const Digit m, Integer& p, Integer& q, Integer& t) const
// Algorithm:  z.series(n, m, p, q, t)
// Input:      z in Zeta3, p,q,t in Integer, n,m in Digit where n < 2^BETA/64-32, n < m.
// Output:     p,q,t in Integer ||
{
  CONDITION(n < m && n < GAMMA/2);

  const static Natural u2("2JSQ4J31843NN94KLN88HF3G1EA6DH2NQ24PCG707KEGHQ1", 32);
  switch (m-n) {
    case 1: {
      p = 2*n-3;
      q = u2; q *= 2*n;
      t = p; t = -t;
      break;
    }
    case 2: {
      t = q = u2;
      q *= 2*n; t *= 2*n+2; q *= t;
      t += 2*n-1; t *= 2*n-3; t = -t;
      if (n < GAMMA/8) p = (2*n-3)*(2*n-1);
      else { p = 2*n-3; p *= 2*n-1; }
      break;
    }
    case 3: {
      if (n < GAMMA/8) p = (2*n-1)*(2*n+1);
      else { p = 2*n-1; p *= 2*n+1; }
      p *= 2*n-3;
      t = q = u2;
      q *= 2*n+2; t *= 2*n+4; q *= t;
      t += 2*n+1; t *= 2*n-1; t = -t;
      t -= q; t *= 2*n-3;
      q *= 2*n; q *= u2;
      break;
    }
    default: {
      Integer p1,q1,p2,q2;
      Integer t1,t2;

      const Digit l = (n+m)/2;
      sqrt_series(n, l, p1, q1, t1);
      sqrt_series(l, m, p2, q2, t2);
      t = t1*q2; t += t1 = t2*p1;
      p = p1*p2; q = q1*q2;
    }
  }
}
void Pi::sqrt_series2(const Digit n, const Digit m, Integer& q, Integer& t) const
// Algorithm:  z.series2(n, m, q, t)
// Input:      z in Zeta3, q,t in Integer, n,m in Digit where n < m.
// Output:     q,t in Integer ||
{
  CONDITION(n < m);

  Integer p1,q1,q2;
  Integer t1,t2;

  const Digit l = (n+m)/2;
  sqrt_series(n, l, p1, q1, t1);
  if (m-l <= 3) {
    Integer p2;
    sqrt_series(l, m, p2, q2, t2);
  } else sqrt_series2(l, m, q2, t2);
  t = t1*q2; t += t1 = t2*p1;
  q = q1*q2;
}

Pi::Pi(const size_t n)
// Algorithm:  a := Pi(n)
// Input:      n in size_t.
// Output:     a in Pi such that |a-pi*10^n| < 1 ||
 : pi(n)
{
  const size_t sz = pi.precision();

  if (BETA < 32) schoenhage(sz, pi.value());
  else if (sz < (3*GAMMA_LOW+3)/BETA && n < 10000)
    stoermer(sz, pi.value());
  else if (n >= Digit(100000000) || n >= ((size_t(1) << (BETA-1))/3))
    schoenhage(sz, pi.value());
  else {
    Integer p1,q1,t1;
    Integer q2,t2;

    // sqrt(C):
    const static Natural u("1JQ2P3OQIFH25E09DDDA1KD1", 32);
    const static Natural v("229Q8MPR3DNV36SDE9P7HD", 32);
    Digit k = 2 + Digit(n/69.65551373);   // 2*log10(u), k >= 5

    Digit l = (k+1)/2;
    sqrt_series(2, l, p1, q1, t1);
    q2 = u; q2 *= u; q2 <<= 1; t1 -= q1; q1 *= q2; t1 += q1;
    sqrt_series2(l, k, q2, t2);
    Integer t = t1*q2;
    t += t1 = t2*p1; t *= u;
    q1 *= v;
    Integer t3 = q1*q2;


    // Pi:
    const SignDigit A = SignDigit(13591409);
    const SignDigit C = SignDigit(640320);
    k = size_t(n/14.1643);
    l = (k+1)/2;
    chudnovsky(1, l, p1, q1, t1);
    chudnovsky2(l, k, q2, t2);
    Integer t4 = t1*q2;
    t4 += t1 = t2*p1;
    p1 = q1*q2;
    t4 += p1*A;
    p1 *= SignDigit(C/12);


    t4 *= t3;
    t *= p1;
    t.lmove(sz);
    pi.value() = abs(t) / abs(t4);

    /*
    t2 = C; t2.lmove(sz*2);
    t1 = sqrt(t2);
    t2 = p1*t1;
    pi.value() = abs(t2) / abs(t4);
    */
  }
}

#ifdef _Old_STD_
#include <math.h>
#else
#include <cmath>
#endif
#ifdef log2
# undef log2
#endif
void Sqrt::series(Digit n, const Digit m, Natural& p, Natural& q, Integer& t) const
// Algorithm:  z.series(n, m, p, q, t)
// Input:      z in Sqrt, p,q,t in Natural, n,m in Digit where n < GAMMA/2, n < m.
// Output:     p,q,t in Natural ||
{
  CONDITION(n < m && n < GAMMA/2);

  switch (m-n) {
    case 1: {
      p = 2*n-3;
      if (d == 0) q = u*u;
      else { q = d; q *= d; }
      q *= 2*n;
      t = p; t = -t;
      break;
    }
    case 2: {
      if (d == 0) t = q = u*u;
      else { q = d; t = q *= d; }
      q *= 2*n; t *= 2*n+2; q *= t;
      t += 2*n-1; t *= 2*n-3; t = -t;
      if (n < GAMMA/8) p = (2*n-3)*(2*n-1);
      else { p = 2*n-3; p *= 2*n-1; }
      break;
    }
    case 3: {
      if (n < GAMMA/8) p = (2*n-1)*(2*n+1);
      else { p = 2*n-1; p *= 2*n+1; }
      p *= 2*n-3;
      if (d == 0) t = q = u*u;
      else { q = d; t = q *= d; }
      q *= 2*n+2; t *= 2*n+4; q *= t;
      t += 2*n+1; t *= 2*n-1; t = -t;
      t -= q; t *= 2*n-3;
      if (d == 0) q *= u*u;
      else { q *= d; q *= d; }
      q *= 2*n;
      break;
    }
    default: {
      Natural p1,q1,p2,q2;
      Integer t1,t2;

      const Digit l = (n+m)/2;
      series(n, l, p1, q1, t1);
      series(l, m, p2, q2, t2);
      t = t1*q2; t += t1 = t2*p1;
      p = p1*p2; q = q1*q2;
    }
  }
}
void Sqrt::series2(const Digit n, const Digit m, Natural& q, Integer& t) const
// Algorithm:  z.series2(n, m, q, t)
// Input:      z in Sqrt, q,t in Natural, n,m in Digit where n < m.
// Output:     q,t in Natural ||
{
  CONDITION(n < m);

  Natural p1,q1,q2;
  Integer t1,t2;

  const Digit l = (n+m)/2;
  series(n, l, p1, q1, t1);
  if (m-l <= 3) {
    Natural p2;
    series(l, m, p2, q2, t2);
  } else series2(l, m, q2, t2);
  t = t1*q2; t += t1 = t2*p1;
  q = q1*q2;
}

Sqrt::Sqrt(const Digit a, const size_t n)
// Algorithm:  b := Sqrt(a, n)
// Input:      a in Digit, n in size_t.
// Output:     b in Sqrt such that |b-sqrt(a)*10^n| < 1 ||
 : sqrt(n), d(0)
{
  const size_t sz = sqrt.precision();
  Digit b,c;
  ::sqrt(a, b, c);
  if (c == 0) {
    sqrt.value() = b;
    sqrt.value().lmove(sz);
  } else {
    Natural p1,q1,q2;
    Integer t1,t2;

    pell(a, u, v);
    if (u.length() == 1) {
      Natural c = u << 1;
      Natural u0 = c*u;
      Natural v0 = c*v;
      --u0;
      if (u0.length() == 1) {
        Natural u1,v1;
        while (true) {
          u1 = c*u0; u1 -= u;
          v1 = c*v0; v1 -= v;
          u = u0; v = v0;
          if (u1.length() > 1) break;
          u0 = u1; v0 = v1;
        }
      }
    }
    if (u.length() == 1) d = u&GAMMA;

    Digit k = 2 + n/Digit((d == 0)? (2*log10(double(GAMMA))) : (2*log10(double(d))));
    if (k < 5) k = 5;
    const Digit l = (k+1)/2;

    series(2, l, p1, q1, t1);
    if (d == 0) q2 = u*u;
    else { q2 = d; q2 *= d; }
    q2 <<= 1; t1 -= q1; q1 *= q2; t1 += q1;
    series2(l, k, q2, t2);

    Integer t = t1*q2;
    t += t1 = t2*p1; t *= u;
    t2 = q1*q2; t2 *= v;
    t.lmove(sz);
    sqrt.value() = abs(t) / abs(t2);
  }
}

void Zeta3::series(Digit n, const Digit m, Integer& p, Integer& q, Integer& t) const
// Algorithm:  z.series(n, m, p, q, t)
// Input:      z in Zeta3, p,q,t in Integer, n,m in Digit where n < 2^BETA/64-32, n < m.
// Output:     p,q,t in Integer ||
{
  CONDITION(n < m && n < ((size_t(1) << (BETA-1))/64)-32);

  if (n+1 == m) {
    if (n <= GAMMA_LOW/2) { p = n*n; p *= n*n; }
    else { p = n; p = p*p; p = p*p; }
    const Digit k = 2*n-1;
    if (n <= GAMMA_LOW/8) { p *= k*n; p *= k*k; }
    else { p *= n; p *= k; p *= k; p *= k; }
    p = -p;
    t = 126392; t *= n;
    t += 412708; t *= n;
    t += 531578; t *= n;
    t += 336367; t *= n;
    t += 104000; t *= n;
    t += 12463; t *= p;

    const Digit k2 = 4*n+1;
    const Digit k3 = k2+2;
    if (k3 <= GAMMA_LOW/2) { q = k2*k2; q *= k2*k3; q *= k3*k3; }
    else {
      q = k2; q *= k2; q *= k2;
      q *= k3; q *= k3; q *= k3;
    }
    q *= 72*n+24; q *= 3*n+2;
  } else {
    Integer p1,q1,t1;
    Integer p2,q2,t2;

    const Digit l = (n+m)/2;
    series(n, l, p1, q1, t1);
    series(l, m, p2, q2, t2);
    t = t1*q2; t += t1 = t2*p1;
    p = p1*p2; q = q1*q2;
  }
}
void Zeta3::series2(const Digit n, const Digit m, Integer& q, Integer& t) const
// Algorithm:  z.series2(n, m, q, t)
// Input:      z in Zeta3, q,t in Integer, n,m in Digit where n < m.
// Output:     q,t in Integer ||
{
  CONDITION(n < m);

  Integer p1,q1,t1;
  Integer q2,t2;

  const Digit l = (n+m)/2;
  series(n, l, p1, q1, t1);
  if (m-l <= 1) {
    Integer p2;
    series(l, m, p2, q2, t2);
  } else series2(l, m, q2, t2);
  t = t1*q2; t += t1 = t2*p1;
  q = q1*q2;
}

void Zeta3::linear(const size_t n, Natural& z) const
// Algorithm:  b.linear(n, a)
// Input:      b in Zeta3, a in Natural, n in size_t where b.decimals()/574 <= GAMMA.
// Output:     a in Natural such that |a-Zeta(3)*2^(BETA*n)| < 1 ||
{
  CONDITION((zeta.decimals()/1.4)/410 <= GAMMA);

  Natural a = 1;
  Natural c = 77;
  Natural t;

  a <<= n*BETA+2;
  z = 0;
  Digit k = 1;
  do {
    t = a * c;
    if (k&1) z += t;
    else z -= t;

    c += 410*k+45;
    if (k <= GAMMA_LOW) { a *= k*k; a *= k*k; }
    else { a *= k; a *= k; a *= k; a *= k; }
    a *= k;
    const Digit i = 2*k+1;
    if (i <= GAMMA_LOW) { a /= i*i; a /= i*i; }
    else { a /= i; a /= i; a /= i; a /= i; }
    a /= 32*i;
    ++k;
  } while (a != 0);
  z >>= 8;
}

Zeta3::Zeta3(const size_t n)
// Algorithm:  a := Zeta3(n)
// Input:      n in size_t.
// Output:     a in Zeta3 such that |a-zeta(3)*10^n| < 1 ||
 : zeta(n)
{
  const size_t sz = zeta.precision();

  if (n <= 4000 && n/574 <= GAMMA) linear(sz, zeta.value());
  else {
    Integer p1,q1,t1;
    Integer q2,t2;

    const size_t k = (n >= 100000)? size_t(n/5.03) : n/5;
    const size_t l = (k+1)/2;
    series(1, l, p1, q1, t1);
    series2(l, k, q2, t2);
    Integer t = t1*q2;
    t += t1 = t2*p1;
    t2 = q1*q2;
    t += t2*SignDigit(12463); t2 *= SignDigit(10368);
    size_t i = 0;
    if (t2 != 0)
      while (!abs(t2).testbit(i)) ++i;
    t2 >>= i; t.lmove(sz);
    zeta.value() = abs(t) / abs(t2);
    zeta.value() >>= i;
  }
}

void Exp1::linear(const size_t n, Natural& z) const
// Algorithm:  b.linear(n, a)
// Input:      b in Exp1, a in Natural, n in size_t.
// Output:     a in Natural such that |a-Exp(1)*2^(BETA*n)| < 1 ||
{
  Natural a = 1;
  a <<= n*BETA;
  z = a << 1;
  Digit k = 2;
  do {
    a /= k;
    z += a;
    ++k;
  } while (a != 0);
}

void Exp1::series(Digit n, const Digit m, Natural& q, Natural& t) const
// Algorithm:  e.series(n, m, q, t)
// Input:      e in Exp1, q,t in Natural, n,m in Digit where n < m.
// Output:     q,t in Natural ||
{
  CONDITION(n < m);

  switch (m-n) {
    case 1: {
      q = (n == 0)? 1 : n;
      t = 1;
      break;
    }
    case 2: {
      const Digit l = (n+m)/2;
      if (n == 0) q = l;
      else if (l < GAMMA_LOW) q = n*l;
      else { q = n; q *= l; }
      t = l+1;
      break;
    }
    case 3: {
      const Digit l = (n+m)/2;
      const Digit l2 = (l+m)/2;
      if (l2 < GAMMA_LOW) {
        const Digit x = l*l2;
        t = x; q = x;
      } else {
        t = l; t *= l2; q = t;
      }
      if (n != 0) q *= n;
      t += l2+1;
      break;
    }
    default: {
      Natural q1,t1;
      Natural q2,t2;

      const Digit l = (n+m)/2;
      series(n, l, q1, t1);
      series(l, m, q2, t2);
      t = t1*q2; t += t2;
      q = q1*q2;
      break;
    }
  }
}

Exp1::Exp1(const size_t n)
// Algorithm:  a := Exp1(n)
// Input:      n in size_t.
// Output:     a in Exp1 such that |a-exp(1)*10^n| < 1 ||
 : exp(n)
{
  const size_t sz = exp.precision();

  if (n <= 4000) linear(sz, exp.value());
  else {
    Natural t;
    // n < (m + 1/2)*log10(double(m)) - 11/25*m
    Digit m = Digit(n/4.239174208864);
    Digit m2 = m/2;
    Digit d = Digit(log10(double(m)) * (m+0.5) - (m/25.0)*11);
    while (d <= n) {
      m *= 2;
      d = Digit(log10(double(m)) * (m+0.5) - (m/25.0)*11);
    }
    Digit d2 = Digit(log10(double(m2)) * (m2+0.5) - (m2/25.0)*11);
    while (d2 > n) {
      m2 /= 2;
      d2 = Digit(log10(double(m2)) * (m2+0.5) - (m2/25.0)*11);
    }
    while (m-m2 > 1) {
      Digit m3 = (m+m2)/2;
      Digit d3 = Digit(log10(double(m3)) * (m3+0.5) - (m3/25.0)*11);
      if (d3 > n) m = m3;
      else m2 = m3;
    }
    series(0, m, ((Natural&)*this), t);
    size_t i = 0;
    if (((Natural&)*this) != 0) {
      for (const Digit* p = last(); *p == 0 && i < BETA*sz; --p) i += BETA;
      while (!testbit(i) && i < BETA*sz) ++i;
    }
    ((Natural&)*this) >>= i; t <<= BETA*sz-i;
    exp.value() = t / ((Natural&)*this);
  }
}

Natural Ln::atanh_inv_linear(const Digit a, const size_t n) const
// Algorithm:  c := b.atanh_inv_linear(a, n)
// Input:      a in Digit, b in Ln, n in size_t.
// Output:     c in Natural such that |c-atanh(1/a)*2^(BETA*n)| < 1 ||
{
  Natural c = Digit(0);
  Natural s,t = a;
  t <<= BETA*n;
  Digit k = 1;
  if (a < GAMMA_LOW) {
    const Digit b = a*a;
    do {
      t /= b;
      c += s = t / k;
      k += 2;
    } while (s != 0);
  } else {
    do {
      t /= a; t /= a;
      c += s = t / k;
      k += 2;
    } while (s != 0);
  }
  return c;
}

void Ln::atanh_inv_series(Digit n, const Digit m,
                          Natural& b, Natural& q, Natural& t) const
// Algorithm:  l.atanh_inv_series(n, m, b, q, t)
// Input:      l in Ln, b,q,t in Natural, n,m in Digit where n < m.
// Output:     b,q,t in Natural ||
{
  CONDITION(n < m);

  if (n+1 == m) {
    b = 2*n+1;
    t = 1;
    if (n == 0) q = x;
    else if (x < GAMMA_LOW) q = x*x;
    else { q = x; q *= q; }
  } else {
    Natural b1,q1,t1;
    Natural b2,q2,t2;
    const Digit l = (n+m)/2;
    atanh_inv_series(n, l, b1, q1, t1);
    atanh_inv_series(l, m, b2, q2, t2);
    q = q1*q2; b = b1*b2;
    t = t2*b1; b1 = b2*t1; t += b2 = b1*q2;
  }
}

Natural Ln::atanh_inv_series(const Digit a)
// Algorithm:  c := b.atanh_inv_linear(a)
// Input:      a in Digit, b in Ln.
// Output:     c in Natural such that |c-atanh(1/a)*2^(BETA*b.precision())| < 1 ||
{
  const size_t n = ln.precision();
  Digit m = Digit(ln.decimals()/log10(double(a)))+1;
  x = a;
  Natural b,q,t;
  atanh_inv_series(0, m, b, q, t);
  size_t i = 0;
  if (q != 0)
    while (!q.testbit(i)) ++i;
  size_t j = 0;
  if (b != 0)
    while (!b.testbit(j)) ++j;
  q >>= i; b >>= j; t <<= BETA*n-i-j;
  q *= b;
  b = t / q;
  return b;
}

Natural Ln::ln2(const size_t n)
// Algorithm:  c := b.ln2(n)
// Input:      b in Ln, n in size_t.
// Output:     c in Natural such that |c-ln(2)*2^(BETA*n)| < 1 ||
{
  if (n <= 400) {
    Natural ln = 144*atanh_inv_linear(255, n);
    ln += 54*atanh_inv_linear(449, n);
    ln += 24*atanh_inv_linear(4801, n);
    ln += 20*atanh_inv_linear(31751, n);
    return ln += 82*atanh_inv_linear(32257, n);
  }
  Natural ln = 144*atanh_inv_series(255);
  ln += 54*atanh_inv_series(449);
  ln += 24*atanh_inv_series(4801);
  ln += 20*atanh_inv_series(31751);
  return ln += 82*atanh_inv_series(32257);
}

void Ln::linear(const size_t n, Natural& ln) const
// Algorithm:  a.linear(n, b)
// Input:      a in Ln, b in Natural, n in size_t.
// Output:     b in Natural such that |b-ln(a.u/2^a.v)*2^(BETA*n)| < 1 ||
{
  Natural s,t = u;
  t <<= BETA*n-v;
  Digit k = 2;
  ln = t;
  do {
    t *= u; t >>= v;
    s = t / k;
    if (k&1) ln += s;
    else ln -= s;
    ++k;
  } while (s != 0);
}

void Ln::series(Digit n, const Digit m,
                  Natural& b, Integer& p, Digit& q, Integer& t) const
// Algorithm:  l.series(n, m, b, p, q, t)
// Input:      l in Ln, b,q in Natural, p,t in Integer, n,m in Digit where n < m.
// Output:     b,q in Natural, p,t in Integer ||
{
  CONDITION(n < m);

  if (n+1 == m) {
    b = n+1;
    p = u; q = v;
    if (n != 0) p = -p;
    t = p;
  } else {
    Natural b1,b2;
    Digit   q1,q2;
    Integer p1,t1,p2,t2;

    const Digit l = (n+m)/2;
    series(n, l, b1, p1, q1, t1);
    series(l, m, b2, p2, q2, t2);
    p = p1*p2; q = q1+q2; b = b1*b2;
    p2 = t2*b1; t = p2*p1;
    b1 = b2 << q2; t += p2 = b1*t1;
  }
}

Ln::Ln(const Digit a, const size_t n)
// Algorithm:  b := Ln(a, n)
// Input:      a in Digit, n in size_t.
// Output:     b in Ln such that |b-ln(a)*10^n| < 1 ||
 : ln(n)
{
  const size_t sz = ln.precision();
  switch (a) {
    case 0: ln.value().errmsg(0, "value undefined!"); break;
    case 1: ln.value() = 0; break;
    case 2: ln.value() = ln2(sz); break;
    default: {
      v = log2(a);
      u = a - (Digit(1) << v);
      if (sz < 4000) linear(sz, ln.value());
      else {
        const Digit m = Digit(n/log10((Digit(1) << v)/double(u)))+1;
        Integer p,t;
        Digit q;
        series(0, m, (Natural&)*this, p, q, t);
        size_t i = 0;
        if (((Natural&)*this) != 0) {
          for (const Digit* p = last(); *p == 0; --p) i += BETA;
          while (!testbit(i)) ++i;
        }
        ((Natural&)*this) >>= i;
        const size_t j = BETA*sz;
        if (i >= j) t >>= i-j+q;
        else {
          i = j-i;
          if (i >= q) t <<= i-q;
          else t >>= q-i;
        }
        ln.value() = ::abs(t) / ((Natural&)*this);
      }
      ln.value() += v*ln2(sz);
      break;
    }
  }
}

void EulerGamma::linear(const size_t n, Natural& euler)
// Algorithm:  a.linear(n, b)
// Input:      a in EulerGamma, b in Natural, n in size_t.
// Output:     b in Natural such that |b-gamma*2^(BETA*n)| < 1 ||
{
  const Digit a = Digit(0.17328679514*n*BETA)+1;  // ln(2)/4
  const Digit m = Digit(3.591121477*a);

  Natural f = 1;
  Natural g = Digit(0);
  f.lmove(n);
  Natural s = f;
  Natural s2 = Digit(0);
  if (a < GAMMA_LOW) {
    const Digit b = a*a;
    if (m < GAMMA_LOW) {
      for (Digit k = 1; k < m; ++k) {
        f *= b; g *= b;
        f /= k*k; g /= k; g += f; g /= k;
        s += f; s2 += g;
      }
    } else {
      Digit k;
      for (k = 1; k < GAMMA_LOW; ++k) {
        f *= b; g *= b;
        f /= k*k; g /= k; g += f; g /= k;
        s += f; s2 += g;
      }
      while (k < m) {
        f *= b; g *= b;
        f /= k; f /= k;
        g /= k; g += f; g /= k;
        s += f; s2 += g;
        ++k;
      }
    }
  } else {
    if (m < GAMMA_LOW) {
      for (Digit k = 1; k < m; ++k) {
        f *= a; f *= a; g *= a; g *= a;
        f /= k*k; g /= k; g += f; g /= k;
        s += f; s2 += g;
      }
    } else {
      Digit k;
      for (k = 1; k < GAMMA_LOW; ++k) {
        f *= a; f *= a; g *= a; g *= a;
        f /= k*k; g /= k; g += f; g /= k;
        s += f; s2 += g;
      }
      while (k < m) {
        f *= a; f *= a; g *= a; g *= a;
        f /= k; f /= k;
        g /= k; g += f; g /= k;
        s += f; s2 += g;
        ++k;
      }
    }
  }
  s2.lmove(n);
  euler = s2/s;
  Ln l(a, this->euler.decimals());
  euler -= l.value();
}

void EulerGamma::series(Digit n, const Digit m,
                          Natural& p, Natural& q, Natural& t,
                          Natural& c, Natural& d, Natural& v) const
// Algorithm:  e.series(n, m, p, q, t, c, d, v)
// Input:      e in EulerGamma, p,q,t,c,d,v in Natural, n,m in Digit where n < m.
// Output:     p,q,t,c,d,v in Natural ||
{
  CONDITION(n < m);

  switch (m-n) {
    case 1: {
      const Digit l = n+1;
      if (l < GAMMA_LOW) q = l*l;
      else { q = l; q *= l; }
      v = t = p = ((Natural&)*this); c = 1; d = l;
      break;
    }
    case 2: {
      const Digit l1 = n+1;
      const Digit l2 = n+2;
      const Digit l3 = l1+l2;
      if (l2 < GAMMA_LOW) { d = q = l1*l2; t = ((Natural&)*this)*(l2*l2); }
      else { q = l1; d = q *= l2; t = ((Natural&)*this)*l2; t *= l2; }
      q *= q; t += p = ((Natural&)*this)*((Natural&)*this);
      if (l3 >= l1) c = l3;
      else { c = l1; c += l2; }
      v = l2*t; v += l1*p;
      break;
    }
    case 3: {
      const Digit l1 = n+1;
      const Digit l2 = n+2;
      const Digit l3 = n+3;
      const Digit l4 = l1+l2;
      Natural s = ((Natural&)*this)*((Natural&)*this);
      v = l1*s; p = s*((Natural&)*this);
      if (l3 < GAMMA_LOW) {
        d = q = l1*l2; q *= l3;
        s += ((Natural&)*this)*(l2*l2); c = l4*l3;
        v += s*l2;
        const Digit l5 = l3*l3;
        t = l5*s; v *= l5; v += l4*p;
      } else {
        q = l1; d = q *= l2; q *= l3;
        s += (((Natural&)*this)*l2)*l2;
        v += s*l2;
        const Digit l5 = l3*l3;
        t = l5*s; v *= l5;
        if (l4 >= l1) { c = l4; v += l4*p; }
        else { c = l1; c += l2; v += c*p; }
        c *= l3;
      }
      q *= q; c += d; t += p;
      v *= l3; v += d*p; d *= l3;
      break;
    }
    default: {
      Natural p1,q1,t1,c1,d1,v1;
      Natural p2,q2,t2,c2,d2,v2;

      const Digit l = (n+m)/2;
      series(n, l, p1, q1, t1, c1, d1, v1);
      series(l, m, p2, q2, t2, c2, d2, v2);
      p = p1*p2; q = q1*q2; p2 = p1*t2;
      t = t1*q2; t += p2;
      c = c1*d2; d = d1*d2;
      p2 *= c1; q1 = c2*d1;
      c += q1;
      q1 = q2*v1; q1 += p2;
      v = d2*q1; q1 = d1*p1; p2 = q1*v2; v += p2;
      break;
    }
  }
}

void EulerGamma::series2(Digit n, const Digit m,
                         Natural& q, Natural& t, Natural& d, Natural& v) const
// Algorithm:  e.series2(n, m, q, t, d, v)
// Input:      e in EulerGamma, q,t,d,v in Natural, n,m in Digit where n < m.
// Output:     q,t,d,v in Natural ||
{
  CONDITION(n < m);

  switch (m-n) {
    case 1: {
      const Digit l = n+1;
      if (l < GAMMA_LOW) q = l*l;
      else { q = l; q *= l; }
      v = t = ((Natural&)*this); d = l;
      break;
    }
    case 2: {
      const Natural p = ((Natural&)*this)*((Natural&)*this);
      const Digit l1 = n+1;
      const Digit l2 = n+2;
      if (l2 < GAMMA_LOW) { d = q = l1*l2; t = ((Natural&)*this)*(l2*l2); }
      else { q = l1; d = q *= l2; t = ((Natural&)*this)*l2; t *= l2; }
      q *= q; t += p;
      v = l2*t; v += l1*p;
      break;
    }
    case 3: {
      const Digit l1 = n+1;
      const Digit l2 = n+2;
      const Digit l3 = n+3;
      const Digit l4 = l1+l2;
      Natural s = ((Natural&)*this)*((Natural&)*this);
      Natural p = s*((Natural&)*this);
      v = l1*s;
      if (l3 < GAMMA_LOW) {
        d = q = l1*l2; q *= l3;
        s += ((Natural&)*this)*(l2*l2);
        v += s*l2;
        const Digit l5 = l3*l3;
        t = l5*s; v *= l5; v += l4*p;
      } else {
        q = l1; d = q *= l2; q *= l3;
        s += (((Natural&)*this)*l2)*l2;
        v += s*l2;
        const Digit l5 = l3*l3;
        t = l5*s; v *= l5;
        if (l4 >= l1) v += l4*p;
        else { s = l1; s += l2; v += s*p; }
      }
      q *= q; t += p;
      v *= l3; v += d*p; d *= l3;
      break;
    }
    default: {
      Natural p1,q1,t1,c1,d1,v1;
      Natural q2,t2,d2,v2;

      const Digit l = (n+m)/2;
      series(n, l, p1, q1, t1, c1, d1, v1);
      series2(l, m, q2, t2, d2, v2);
      q = q1*q2; d = d1*d2;
      t = t1*q2; v = p1*t2; t += v; c1 *= v;
      q1 = q2*v1; q1 += c1;
      v = d2*q1; d2 = d1*p1; q1 = d2*v2; v += q1;
      break;
    }
  }
}

EulerGamma::EulerGamma(const size_t n)
// Algorithm:  a := EulerGamma(n)
// Input:      n in size_t.
// Output:     a in EulerGamma such that |a-gamma*10^n| < 1 ||
 : euler(n)
{
  const size_t sz = euler.precision();

  if (n <= 10000) linear(sz, euler.value());
  else {
    const Digit a = Digit(0.17328679514*sz*BETA)+1;  // ln(2)/4
    Digit m = Digit(3.591121477*a);
    ((Natural&)*this) = a; ((Natural&)*this) *= a;

    Natural* q = new Natural(); // To reduce memory consumption
    Natural* t = new Natural();
    Natural* d = new Natural();
    Natural* v = new Natural();
    series2(0, m, *q, *t, *d, *v);
    *q += *t;
    delete t;
    // ToDo: optimize by 2^k|*this and 2^k|*v
    ((Natural&)*this) = (*q)*(*d);
    delete d;
    delete q;
    size_t i = 0;
    if (((Natural&)*this) != 0) {
      for (const Digit* p = last(); *p == 0 && i < BETA*sz; --p) i += BETA;
      while (!testbit(i) && i < BETA*sz) ++i;
    }
    ((Natural&)*this) >>= i; *v <<= BETA*sz-i;

    euler.value() = *v / ((Natural&)*this);
    delete v;
    Ln ln(a, euler.decimals());
    euler.value() -= ln.value();
  }
}
