/////////////////////////////////
//
// Piologie V 1.3.3
// multi-precision arithmetic
// Natural
//
// (c) 1996-2001 HiPiLib
//
// Sebastian Wedeniwski
// 10/01/2011
// div bug fixed by Grant Atoyan
//

#include "natural.h"
#include <string.h>


const size_t KARATSUBA_MUL_MARK = 8;
const size_t KARATSUBA_SQR_MARK = 2*KARATSUBA_MUL_MARK;

#ifdef _DigitAsm_
const size_t FFT_SQR_MARK = 8125;
const size_t FFT_MUL_MARK = 5000;
#else
const size_t FFT_SQR_MARK = 11875;
const size_t FFT_MUL_MARK = 16875;
#endif

const size_t NEWTON_DIV_MARK1 = 500;
const size_t NEWTON_DIV_MARK2 = 2500;       // >= 2

const size_t NEWTON_SQRT_MARK = 500;       // >= 2

const size_t NTOA_MARK = 17000;

const size_t ATON_MARK = 20000;


#define FASTER_BINSQRT      // binary square root without shifting
                            // need BETA/2 times more memory!


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

#define FILL_ZERO(a, b)                                 \
  do *a++ = 0; while (a < b)
//do *a++ = 0; while (a != b)

#define FILL_DELTA(a)                                   \
  a[0] = a[1] = a[2] = a[3] = a[4] = a[5] = a[6] = a[7] = 0;

#define COPY(a, b, c, d)                                \
  do *a++ = *b++; while (c != d);

#define COPY_BACKWARD(a, b, c, d)                       \
  do *--a = *--b; while (c != d);

#define MOVE(a, b, c, d)                                \
  do { *a++ = *b; *b++ = 0; } while (c != d);

#define MOVE_BACKWARD(a, b, c, d)                       \
  do { *--a = *--b; *b = 0; } while (c != d);

#ifndef _Old_STD_
using namespace std;
# define NOTHROW_NEW  new(nothrow)
#else
# define NOTHROW_NEW  new
#endif



class FFT : private NumberBase {
private:
  Natural             t[3];
  FFT*                factor;
  size_t              shift;
  const Natural&      arg;
  static const Digit* moduli;
  static const Digit* primroots;
  static Digit        m;
  static size_t       n1,n2;
  static Digit*       omega;
  static Digit*       omega2;
  static size_t*      order;

  static const Digit* init_moduli();
  static const Digit* init_primroots();
  static size_t       max_size();
  static void         setmodulo(const Digit);

  void  digitmulmod(const Digit, const Digit, const Digit, Digit&) const;
  Digit pow(Digit, Digit, const Digit) const;
  Digit digitinv(Digit, const Digit) const;
  void  init_omega(const Digit);
  void  innerfft(const Digit*, Digit*, Digit*, const size_t, const Digit) const;
  void  innerfft(const Digit*, Digit*, Digit*, const size_t, const Digit,
                 const Digit) const;
  void  innerfftinv(const Digit*, Digit*, Digit*, const size_t, const Digit,
                    const Digit) const;
  void  innerfft3(Digit*, const size_t, const Digit) const;
  void  innerfftinv3(Digit*, const size_t, const Digit) const;
  void  fft(Digit*, const Digit*, const size_t) const;
  void  multiply_matrix(Digit*, const Digit, const Digit) const;
  void  fftinv(Digit*, const Digit*, const size_t) const;
  void  five_step(const Digit*, Digit*, const Digit) const;
  void  five_step(const Digit*, Digit*, const Digit, const Digit) const;
  void  chinese_remainder();

  void  square(const size_t);
  void  multiply(const size_t);
  void  result(Natural&) const;

  void  base_conversion(Natural&, const Natural&, const size_t);

  FFT(const Natural&, const bool = true, FFT* = 0);
  ~FFT();

  size_t size() const;
  void sqr(Natural&);
  void mul(Natural&);

  friend class Natural;
};


inline size_t FFT::size() const
// Algorithm:  a := f.size()
// Input:      f in FFT.
// Output:     a in size_t such that a = L(f.t_0) ||
{
  return t[0].size;
}

inline void FFT::setmodulo(const Digit a)
// Algorithm:  b.setmodulo(a)
// Input:      b in FFT, a in Digit where a in b.moduli.
// Output:     b in FFT such that b.m = a ||
{
  CONDITION(a == moduli[0] || a == moduli[1] || a == moduli[2]);

  m = a;
}

#if defined(_DigitAsm_) && _M_IX86 >= 300 && defined(_MSC_VER)
static Digit __declspec(naked) __fastcall _digitmulmod(const Digit a, const Digit b,
                                                       const Digit m)
{
  __asm {
    mov eax,ecx
    mul edx
    div dword ptr [esp+4]
    mov eax,edx
    ret 4
  }
}

inline void FFT::digitmulmod(const Digit a, const Digit b, const Digit m, Digit& c) const
// Algorithm:  f.digitmulmod(a, b, m, c)
// Input:      f in FFT, a,b,m in Digit where m > 0, a,b <= m.
// Output:     c in Digit such that c = a*b - f.m*[(a*b)/f.m] ||
{
  c = _digitmulmod(a, b, m);
}
#elif defined(_DigitAsm_) && (defined (__i386__) || defined (__i486__)) && defined (__GNUC__)
inline void FFT::digitmulmod(const Digit a, const Digit b, const Digit m, Digit& c) const
{
  int dummy;
  __asm__ ("mull %3; divl %4"
    : "=&d" (c), "=a" (dummy)
    : "%a" (a), "rm" (b), "rm" (m));
}
#else
inline void FFT::digitmulmod(const Digit a, const Digit b, const Digit M, Digit& c) const
// Algorithm:  f.digitmulmod(a, b, m, c)
// Input:      f in FFT, a,b,m in Digit where m > 0, a,b <= m, BETA in {32, 64},
//             m in f.moduli.
// Output:     c in Digit such that c = a*b - m*[(a*b)/m] ||
{
  CONDITION(BETA == 32 || BETA == 64);

  Digit x,y;
  digitmul(a, b, x, y);

#if defined(_DigitAsm_)
  digitdiv(x, y, M, x, c);
#elif defined(_Unknown_Apogee_Bug_)
  digitmod(x, y, M, c);
#else
  if (BETA == 32) {
    const Digit M1 = (Digit(1) << (BETA-1)) - (Digit(1) << (BETA-4))
                                            - (Digit(1) << (BETA-6)) + 1;
    const Digit M2 = (Digit(1) << (BETA-1)) - (Digit(1) << (BETA-5)) + 1;
    const Digit M3 = (Digit(1) << (BETA-1)) - (Digit(1) << (BETA-7)) + 1;

    CONDITION(M == M1 || M == M2 || M == M3);

    if (M == M1) {
      do {
        const Digit r = x << (BETA-5);
        const Digit s = x << (BETA-3);
        const Digit u = y - (x << 1);
        x >>= BETA-29; x += x >> 2;
        if (y < u) --x;
        y = u+s; if (y < s) ++x;
        y += r; if (y < r) ++x;
      } while (x);
      if (y >= M1) {
        y -= M1;
        if (y >= M1) {
          y -= M1;
          if (y >= M1) y -= M1;
        }
      }
      c = y;
    } else if (M == M2) {
      do {
        const Digit s = x << (BETA-4);
        const Digit u = y - (x << 1);
        x >>= BETA-28; if (y < u) --x;
        y = u+s; if (y < s) ++x;
      } while (x);
      if (y >= M2) {
        y -= M2;
        if (y >= M2) {
          y -= M2;
          if (y >= M2) y -= M2;
        }
      }
      c = y;
    } else { // M = M3
      do {
        const Digit s = x << (BETA-6);
        const Digit u = y - (x << 1);
        x >>= BETA-26; if (y < u) --x;
        y = u+s; if (y < s) ++x;
      } while (x);
      if (y >= M3) {
        y -= M3;
        if (y >= M3) {
          y -= M3;
          if (y >= M3) y -= M3;
        }
      }
      c = y;
    }
    
  } else {  // BETA == 64
    const Digit M1 = Digit(1 - (Digit(1) << (BETA-24)));
    const Digit M2 = Digit(1 - (Digit(1) << (BETA-30)));
    const Digit M3 = Digit(1 - (Digit(1) << (BETA/2)));

    CONDITION(M == M1 || M == M2 || M == M3);

    if (M == M1) {
      Digit s = x << (BETA-24);
      Digit u = y - x;
      x >>= 24; if (y < u) --x;
      y = u+s; if (y < s) ++x;
    
      s = x << (BETA-24); u = y - x;
      x >>= 24; if (y < u) --x;
      y = u+s; if (y < s) ++x;
    
      s = x << (BETA-24); u = y - x;
      x >>= 24; if (y < u) --x;
      y = u+s; if (y < s) ++x;
      if (x || y >= M1) y -= M1;
      c = y;
    } else if (M == M2) {
      Digit s = x << (BETA-30);
      Digit u = y - x;
      x >>= 30; if (y < u) --x;
      y = u+s; if (y < s) ++x;
    
      s = x << (BETA-30); u = y - x;
      x >>= 30; if (y < u) --x;
      y = u+s; if (y < s) ++x;
    
      s = x << (BETA-30); u = y - x;
      x >>= 30; if (y < u) --x;
      y = u+s; if (y < s) ++x;
      if (x || y >= M2) y -= M2;
      c = y;
    } else {  // M == M3
      Digit s = x << BETA/2;
      Digit u = y - x;
      x >>= BETA/2; if (y < u) --x;
      y = u+s; if (y < s) ++x;
    
      s = x << BETA/2; u = y - x;
      x >>= BETA/2; if (y < u) --x;
      y = u+s; if (y < s) ++x;
      if (x || y >= M3) y -= M3;
      c = y;
    }
    
  }
#endif
}

#endif


/////////////////// Natural arithmetic /////////////////////////

size_t Natural::NaturalSize    = DELTA;
size_t Natural::NaturalSizeOld = DELTA;

void Natural::get_memory(const size_t a)
// Algorithm:  c.get_memory(a)
// Input:      a in size_t where a >= 1.
// Output:     c in Natural such that L(c) = R(c) = a ||
//
// internal allocation without the initialization of the elements.
{
  CONDITION(a >= 1);

  Digit* pT = NOTHROW_NEW Digit[size = a];
  p = root = pT;
  if (!pT) errmsg(2, "(get_memory)");
}

Natural::Natural(const Digit a, size_t b)
// Algorithm:  c := Natural(a, b)
// Input:      a in Digit, b in size_t where b >= 2.
// Output:     c in Natural such that c = a*2^{BETA*(b-1)}, R(c) = b+DELTA ||
//
// Note:       This constructor don't fulfill the conditions for Naturals
//             e.g. a = 0.
{
  CONDITION(b >= 2);

  Digit* pT = NOTHROW_NEW Digit[b+DELTA];
  root = pT; size = b;
  if (!pT) errmsg(2, "(internal constructor)");
  FILL_DELTA(pT);
  p = pT += DELTA;
  const Digit* pE = pT+b;
  *pT++ = a;
  FILL_ZERO(pT, pE);
}
Natural::Natural(size_t a, const Natural& b)
// Algorithm:  c := Natural(a, b)
// Input:      a in size_t, b in Natural where a >= 1.
// Output:     c in Natural such that c = b, R(c) = L(b)+a ||
{
  NATURALCONDITION(b);
  CONDITION(a >= 1);

  const size_t sB = b.size;
  Digit* pT = NOTHROW_NEW Digit[sB+a];
  root = pT; size = sB;
  if (!pT) errmsg(2, "(internal constructor)");
  const Digit* pE = pT+a;
  FILL_ZERO(pT, pE);
  pE += sB; p = pT;
  const Digit* pB = b.p;
  COPY(pT, pB, pT, pE);

  NATURALCONDITION(*this);
}
Natural::Natural(const Natural& a, size_t b)
// Algorithm:  c := Natural(a, b)
// Input:      a in Natural, b in size_t where b >= 1.
// Output:     c in Natural such that c = a*2^{BETA*b},
//             R(c) = L(a)+b+DELTA ||
//
// Note:       This constructor don't fulfill the conditions for Naturals (a=0).
{
  NATURALCONDITION(a);
  CONDITION(b >= 1);

  const size_t sA = a.size;
  const size_t sT = sA+b;
  Digit* pT = NOTHROW_NEW Digit[sT+DELTA];
  root = pT; size = sT;
  if (!pT) errmsg(2, "(internal constructor)");
  FILL_DELTA(pT);
  p = pT += DELTA;
  const Digit* pE = pT+sA;
  const Digit* pA = a.p;
  COPY(pT, pA, pT, pE);
  pE += b;
  FILL_ZERO(pT, pE);
}

size_t Natural::trailing_zeros(Digit*& a) const
// Algorithm:  c := b.trailing_zeros(a)
// Input:      b in Natural where not b = 0.
// Output:     a in [b.p, b.p+b.size[, c in size_t
//             such that not (a) = (0), [a+1, a+c] = 0^c
//             and a+c = b.p+b.size-1 ||
{
  NATURALCONDITION(*this);
  CONDITION(*this != 0);

  Digit* pT = p+size-1;
  if (*pT) { a = pT; return 0; }
  else {
    size_t c = 1;
    while (*--pT == 0) ++c;
    a = pT;
    return c;
  }
}

Digit* Natural::setsize(const size_t b)
// Algorithm:  r := a.setsize(b)
// Input:      a in Natural, b in size_t where b >= 1.
// Output:     a in Natural, r in [a.root, a.p+L(a)[
//             such that R(a) >= b+DELTA, L(a) = b, r = a.p ||
//
// Note:       a is not normalize!
{
  CONDITION(b >= 1);

  Digit* pT = p;
  const size_t sT = size;
  Digit* rT = root;
  if (pT+sT < rT+b+DELTA) {                // (rootsize() < b+DELTA)?
    delete[] rT;
    root = pT = NOTHROW_NEW Digit[b+DELTA];
    if (!pT) errmsg(2, "(setsize)");
    FILL_DELTA(pT);
    pT += DELTA;
  } else if (sT > b) {
    const Digit* pE = pT+sT-b;
    FILL_ZERO(pT, pE);
  } else pT -= b-sT;
  p = pT; size = b;
  return pT;
}

Natural& Natural::copy(const Natural& a, const size_t b)
// Algorithm:  c.copy(a, b)
// Input:      a in Natural, b in size_t where 1 <= b <= L(a).
// Output:     c in Natural such that
//             c = [a/2^(BETA*(L(a)-b))]*2^(BETA*(L(a)-b))
//                 + c mod 2^(BETA*(L(a)-b)) ||
{
  NATURALCONDITION(a);
  CONDITION(b >= 1 && b <= a.size);

  const size_t sA = a.size;
  const size_t sT = size;
  Digit* rT = root;
  Digit* pT = p;
  if (sA >= sT+size_t(pT-rT)) {        // (a.size >= rootsize())?
    Digit* pX = NOTHROW_NEW Digit[sA+DELTA];
    if (!pX) errmsg(2, "(copy)");
    root = pX;
    FILL_DELTA(pX);
    pX += DELTA;
    const size_t sz = sA-b;
    if (sz) {
      pX += b; pT += sT;
      const Digit* pE = pT;
      pT -= sz;
      COPY(pX, pT, pT, pE);
      pX -= sA;
    }
    pT = pX;
    delete[] rT;
  } else if (sT > sA) {
    const Digit* pE = pT+sT-sA;
    FILL_ZERO(pT, pE);
  } else pT -= sA-sT;
  p = pT; size = sA;
  const Digit* pA = a.p;
  const Digit* pE = pA+b;
  COPY(pT, pA, pA, pE);

  NATURALCONDITION(*this);

  return *this;
}

void Natural::enlarge(const size_t b)
// Algorithm:  a.enlarge(b)
// Input:      a in Natural, b in size_t where b >= 1.
// Output:     a in Natural such that R(a) := b+R(a) ||
{
  CONDITION(b >= 1);

  const size_t sA = p-root+b;
  const size_t sT = size;
  Digit* pT = NOTHROW_NEW Digit[sA+sT];
  if (!pT) errmsg(2, "(enlarge)");
  Digit* pA = pT+sA;
  FILL_ZERO(pT, pA);
  pA = p; p = pT;
  const Digit* pE = pA+sT;
  COPY(pT, pA, pA, pE);
  delete[] root;
  root = pT-sT-sA;
}

void Natural::inc(Digit* b)
// Algorithm:  a.inc(b)
// Input:      a in Natural, b in [a.p, a.p+L(a)].
// Output:     a in Natural such that a := a + 2^(BETA*(a.p+L(a)-b)) ||
{
  NATURALCONDITION(*this);
  CONDITION(b >= p && b <= p+size);

  while (++(*--b) == 0);
  if (b < p) {
    ++size; p = b;
    if (b == root) enlarge(DELTA);
  }

  NATURALCONDITION(*this);
}

void Natural::dec(Digit* b)
// Algorithm:  a.dec(b)
// Input:      a in Natural, b in [a.p, a.p+L(a)].
// Output:     a in Natural such that a := a - 2^(BETA*(a.p+L(a)-b)) ||
{
  NATURALCONDITION(*this);
  CONDITION(b >= p && b <= p+size);

  Digit c;
  Digit* pT = p;
  while (b != pT) {
    c = --(*--b);
    if (c != GAMMA) {
      if (c == 0 && b == pT) {
        size_t sT = size;
        if (sT > 2) {       // is possible, e.g. after sub!
          do { ++pT; --sT; } while (*pT == 0 && sT > 1);
          p = pT; size = sT;
        } else if (sT == 2) { p = ++pT; size = --sT; }
      }

      NATURALCONDITION(*this);

      return;
    }
  }
  errmsg(3, "(dec)");
}

void Natural::add_with_inc(const Digit* pT, Digit* pSum, const Digit* pSmd)
// Algorithm:  a.add_with_inc(r, s, t)
//             Let b in Natural.
// Input:      a in Natural where R(b) > L(a) >= L(b) or R(a) > L(b) >= L(a),
//             r,s in [a.p, a.p+L(a)] where r < s,
//             t in [b.p, b.p+L(b)] where t-(s-r) also in [b.p, b.p+L(b)].
// Output:     a in Natural such that [r, s[ := [r, s[ + [t-(s-r), t[ ||
//
// Note:       R(b) > L(a) >= L(b) or R(a) > L(b) >= L(a)!
{
  CONDITION(pT < pSum && pT >= p && pT <= p+size);
  CONDITION(pSum >= p && pSum <= p+size);

  Digit c,d;
  do {
    c = *--pSmd;                          // non carry addition
    d = *--pSum;
    *pSum = d += c;
    if (c > d)
      do {
        c = *--pSmd;
        d = *--pSum;
        *pSum = d += c+1;                 // addition with carry
      } while (c >= d);
  } while (pSum > pT);
  pT = p;
  if (pSum < pT) { p = pSum; ++size; }
}

void Natural::add_with_inc(const Digit* pT, Digit* pSum,
                         const Digit* pSmd1, const Digit* pSmd2)
// Algorithm:  a.add_with_inc(r, s, u, v)
//             Let b,c in Natural.
// Input:      a in Natural where L(a) >= L(b) and L(a) >= L(c)
//             and (R(b) > L(c) >= L(b) or R(c) > L(b) >= L(c)),
//             r,s in [a.p, a.p+L(a)] where r < s,
//             u in [b.p, b.p+L(b)] where u-(s-r) also in [b.p, b.p+L(b)].
//             v in [c.p, c.p+L(c)] where v-(s-r) also in [c.p, c.p+L(c)].
// Output:     a in Natural such that [r, s[ := [u-(s-r), u[ + [v-(s-r), v[ ||
//
// Note:       R(b) > L(c) >= L(b) or R(c) > L(b) >= L(c)!
{
  CONDITION(pT < pSum && pT >= p && pT <= p+size);
  CONDITION(pSum >= p && pSum <= p+size);

  Digit c,d;
  do {
    c = *--pSmd1;                 // non carry addition
    d = *--pSmd2;
    *--pSum = c += d;
    if (d > c)
      do {
        c = *--pSmd1;
        d = *--pSmd2;
        *--pSum = c += d+1;
      } while (d >= c);
  } while (pSum > pT);
  pT = p;
  if (pSum < pT) { p = pSum; ++size; }
}

bool Natural::add_no_inc(const Digit* pT, Digit* pSum, const Digit* pSmd) const
// Algorithm:  d := x.add_no_inc(r, s, t)
//             Let a,b in Natural.
// Input:      x in Natural,
//             r,s in [a.root, a.p+L(a)] where r < s,
//             t in [b.root, b.p+L(b)] where t-(s-r) also in [b.root, b.p+L(b)].
// Output:     d in bool such that [d] x [r, s[ := [r, s[ + [t-(s-r), t[ ||
{
  CONDITION(pT < pSum);

  Digit c,d;
  do {
    c = *--pSmd;                          // non carry addition
    d = *--pSum;
    *pSum = d += c;
    if (c > d)
      do {
        if (pSum == pT) return true;
        c = *--pSmd;
        d = *--pSum;
        *pSum = d += c+1;                 // addition with carry
      } while (c >= d);
  } while (pSum != pT);
  return false;
}

bool Natural::add_no_inc(const Digit* pT, Digit* pSum,
                         const Digit* pSmd1, const Digit* pSmd2) const
// Algorithm:  d := x.add_no_inc(r, s, u, v)
//             Let a,b,c in Natural.
// Input:      x in Natural,
//             r,s in [a.root, a.p+L(a)] where r < s,
//             u in [b.root, b.p+L(b)] where u-(s-r) also in [b.root, b.p+L(b)].
//             v in [c.root, c.p+L(c)] where v-(s-r) also in [c.root, c.p+L(c)].
// Output:     d in bool such that [d] x [r, s[ := [u-(s-r), u[ + [v-(s-r), v[ ||
{
  CONDITION(pT < pSum);

  do {
    Digit c = *--pSmd1;                 // non carry addition
    Digit d = *--pSmd2;
    *--pSum = c += d;
    if (d > c)
      do {
        if (pSum == pT) return true;    // addition with carry
        c = *--pSmd1;
        d = *--pSmd2;
        *--pSum = c += d+1;
      } while (d >= c);
  } while (pSum != pT);
  return false;
}

void Natural::sub(const Digit* pT, Digit* pDif,
                  const Digit* pMin, const Digit* pSub)
// Algorithm:  a.sub(r, s, u, v)
//             Let b,c in Natural.
// Input:      a in Natural,
//             r,s in [a.root, a.p+L(a)] where r < s,
//             u in [b.root, b.p+L(b)] where u-(s-r) also in [b.root, b.p+L(b)].
//             v in [c.root, c.p+L(c)] where v-(s-r) also in [c.root, c.p+L(c)].
// Output:     [r, s[ := [u-(s-r), u[ - [v-(s-r), v[ ||
{
  CONDITION(pT < pDif && pT >= p && pT <= p+size);
  CONDITION(pDif >= p && pDif <= p+size);

  Digit c,d;
  do {
    c = *--pSub;
    d = *--pMin;
    if (d >= c) *--pDif = d - c;
    else {
      *--pDif = d - c;
      do {
        if (pT == pDif) { dec(pDif); return; }
        c = ~(*--pSub);
        d = *--pMin;
        *--pDif = d += c;
      } while (d >= c);
    }
  } while (pT != pDif);
}

bool Natural::sub_no_dec(const Digit* pT, Digit* pDif, const Digit* pSub) const
// Algorithm:  d := x.sub_no_inc(r, s, t)
//             Let a,b in Natural.
// Input:      x in Natural,
//             r,s in [a.p a.p+L(a)] where r < s,
//             t in [b.p b.p+L(b)] where t-(s-r) also in [b.p b.p+L(b)].
// Output:     d in bool such that [d] x [r, s[ := [r, s[ - [t-(s-r), t[ ||
{
  CONDITION(pT < pDif);

  Digit c,d;
  do {
    c = *--pSub;
    d = *--pDif;
    if (d >= c) *pDif = d -= c;
    else {
      *pDif = d -= c;
      do {
        if (pT == pDif) return true;
        c = ~(*--pSub);
        d = *--pDif;
        *pDif = d += c;
      } while (d >= c);
    }
  } while (pT != pDif);
  return false;
}

inline void subpos(const Digit* pT, Digit* pDif, const Digit* pSub)
// Algorithm:  subpos(r, s, t)
//             Let a,b in Natural.
// Input:      r,s in [a.p, a.p+L(a)] where r < s,
//             t in [b.p, b.p+L(b)] where t-(s-r) also in [b.p, b.p+L(b)],
//             where R(a) > L(b), [r, s[ >= [t-(s-r), t[.
// Output:     [r, s[ := [r, s[ - [t-(s-r), t[ ||
//
// Note:       R(a) > L(b)!
{
  CONDITION(pT < pDif);

  Digit c,d;
  do {
    c = *--pSub;
    d = *--pDif;
    if (d >= c) *pDif = d -= c;
    else {
      *pDif = d -= c;
      do {
        c = ~(*--pSub); 
        d = *--pDif;
        *pDif = d += c;
      } while (d >= c);
    }
  } while (pT < pDif);
}

inline void subpos(const Digit* pT, Digit* pDif,
                   const Digit* pMin, const Digit* pSub)
// Algorithm:  subpos(r, s, u, v)
//             Let b,c in Natural.
// Input:      a in Natural,
//             r,s in [a.root, a.p+L(a)] where r < s,
//             u in [b.root, b.p+L(b)] where u-(s-r) also in [b.root, b.p+L(b)].
//             v in [c.root, c.p+L(c)] where v-(s-r) also in [c.root, c.p+L(c)].
//             where R(a) > L(b), [u-(s-r), u[ >= [v-(s-r), v[.
// Output:     [r, s[ := [u-(s-r), u[ - [v-(s-r), v[ ||
//
// Note:       R(a) > L(b)!
{
  CONDITION(pT < pDif);

  Digit c,d;
  do {
    c = *--pSub;
    d = *--pMin;
    if (d >= c) *--pDif = d - c;
    else {
      *--pDif = d - c;
      do {
        c = ~(*--pSub);
        d = *--pMin;
        *--pDif = d += c;
      } while (d >= c);
    }
  } while (pT < pDif);
}

int Natural::abs(Digit* pC, const Digit* pA, const Digit* pB, size_t n) const
// Algorithm:  d := x.abs(r, s, t, n)
//             Let a,b,c in Natural.
// Input:      x in Natural,
//             n in size_t where n > 0,
//             r,r+n in [a.p, a.p+L(a)], s,s+n in [b.p, b.p+L(b)],
//             t,t+n in [c.p, c.p+L(c)].
// Output:     d in int such that d = sign([s, s+n[ - [t, t+n[),
//             [r, r+n[ = |[s, s+n[ - [t, t+n[| ||
{
  CONDITION(n > 0);

  do {
    const Digit a = *pA;
    const Digit b = *pB;
    if (a > b) {
      subpos(pC, pC+n, pA+n, pB+n);
      return 1;
    } else if (a < b) {
      subpos(pC, pC+n, pB+n, pA+n);
      return -1;
    }
    *pC = 0; ++pA; ++pB; ++pC;
  } while (--n);
  return 0;
}

#if defined(_DigitAsm_) && _M_IX86 >= 300 && defined(_MSC_VER)
Digit Natural::mul(const Digit* pE, const Digit* pA, Digit* pT,
                   const Digit c) const
{
  Digit r;
  __asm {
    xor ebx,ebx
    mov edi,pA
    mov esi,pT
    mov ecx,pE
L1: sub edi,4
    sub esi,4
    mov eax,[edi]
    mul c
    add eax,ebx
    adc edx,0
    mov [esi],eax
    mov ebx,edx
    cmp edi,ecx
    jne L1
    mov [esi-4],ebx
    mov r,ebx
  }
  return r;
}

Digit* Natural::muladd(const Digit* pE, const Digit* pA, Digit* pT,
                       const Digit c) const
{
  Digit* r;
  __asm {
    xor ebx,ebx
    mov edi,pA
    mov esi,pT
    mov ecx,pE
L1: sub edi,4
    sub esi,4
    mov eax,[edi]
    mul c
    add eax,ebx
    adc edx,0
    add eax,[esi]
    adc edx,0
    mov [esi],eax
    mov ebx,edx
    cmp edi,ecx
    jne L1
    sub esi,4
    mov eax,[esi]
    add eax,ebx
    mov [esi],eax
    jnc L3
L2: sub esi,4
    mov eax,[esi]
    inc eax
    mov [esi],eax
    jz  L2
L3: mov r,esi
  }
  return r;
}
#else
Digit Natural::mul(const Digit* pE, const Digit* pA, Digit* pT,
                   const Digit c) const
// Algorithm:  d.mul(r, s, t, c)
//             Let a,b in Natural where R(b) > L(a)+1.
// Input:      d in Natural,
//             r,s in [a.p, a.p+L(a)] where r < s,
//             t in [b.p, b.p+L(b)] where t-s+r-1 also in [b.root+1, b.p+L(b)].
// Output:     [t-(s-r)-1, t[ = [r, s[ * c ||
//
// Note:       t-s+r-1 in [b.root+1, b.p+L(b)].
{
  CONDITION(pE < pA);

  Digit x,y,z = 0;
  do {
    digitmul(*--pA, c, x, y);
    y += z;
    z = x + (y < z);
    *--pT = y;
  } while (pA != pE);
  *--pT = z;
  return z;
}

Digit* Natural::muladd(const Digit* pE, const Digit* pA, Digit* pT,
                       const Digit c) const
// Algorithm:  b.muladd(r, s, t, c)
//             Let a in Natural.
// Input:      b in Natural where R(b) > L(a)+2,
//             r,s in [a.p, a.p+L(a)] where r < s,
//             t in [b.p, b.p+L(b)] where t-s+r-1 also in [b.root+2, b.p+L(b)].
// Output:     [t-(s-r)-1, t[ := [t-(s-r)-1, t[ + [r, s[ * c ||
//
// Note:       t-s+r-1 in [b.root+2, b.p+L(b)].
{
  CONDITION(pE < pA);

  Digit x,y,z = 0;
  do {
    digitmul(*--pA, c, x, y);
    y += z;
    z = x + (y < z);
    x = *--pT;
    y += x;
    z += (y < x);
    *pT = y;
  } while (pA != pE);
  x = *--pT;
  *pT = x += z;
  if (x < z)
    while (++(*--pT) == 0);
  return pT;
}

#endif
void Natural::sqr(const Digit* pA, Digit* pB, size_t n) const
// Algorithm:  x.sqr(r, s, n)
//             Let a,b in Natural where not a.p = b.p.
// Input:      x in Natural,
//             n in size_t where n > 0,
//             r,r+n in [a.p, a.p+L(a)],
//             s,s+2*n in [b.p, b.p+L(b)].
// Output:     [s, s+2*n[ = ([r, r+n[)^2 ||
{
  CONDITION(n > 0);

  if (n == 1) {
    const Digit x = *pA;
    digitmul(x, x, *pB, pB[1]);
    return;
  } else if (n == 2) {
    digitsqr(pA[0], pA[1], pB);
    return;
  }
  const Digit* pE = pA + (n&(GAMMA-1));
  do {
    digitsqr(*pA, pA[1], pB);
    pA += 2; pB += 4;
  } while (pA != pE);
  pE -= n;
  if (n&1) {
    digitmul(*pA, *pA, *pB, pB[1]);
    ++pE; ++pB; --n;
    const Digit c = *pA;
    Digit x,y;
    do {
      digitmul(*--pA, c, x, y);
      Digit* pC = --pB;
      Digit z = *pB;
      *pC = z += y;
      if (z < y) {
        z = *--pC;
        *pC = z += x+1;
        if (z <= x)
          while (++(*--pC) == 0);
      } else {
        z = *--pC;
        *pC = z += x;
        if (z < x)
          while (++(*--pC) == 0);
      }
      pC = pB; z = *pB;
      *pC = z += y;
      if (z < y) {
        z = *--pC;
        *pC = z += x+1;
        if (z <= x)
          while (++(*--pC) == 0);
      } else {
        z = *--pC;
        *pC = z += x;
        if (z < x)
          while (++(*--pC) == 0);
      }
    } while (pA != pE);
    pB += n-1; pA += n;
  }
  pA -= 2; --pB; n -= 2;
  while (n) {
    const Digit c[2] = { *pA, pA[1] };
    Digit x[4];
    do {
      pA -= 2;
      digitmul(*pA, pA[1], c[0], c[1], x);
      Digit y = x[0] >> (BETA-1);
      pB -= 2;
      x[0] <<= 1; x[0] |= x[1] >> (BETA-1); x[1] <<= 1; x[1] |= x[2] >> (BETA-1);
      x[2] <<= 1; x[2] |= x[3] >> (BETA-1); x[3] <<= 1;
      Digit* pC = pB;
      *pB += x[3];
      if (*pB < x[3]) {
        *--pC += x[2]+1;
        if (*pC <= x[2]) {
          *--pC += x[1]+1;
          if (*pC <= x[1]) {
            *--pC += x[0]+1;
            if (*pC <= x[0]) ++y;
          } else {
            *--pC += x[0];
            if (*pC < x[0]) ++y;
          }
        } else {
          *--pC += x[1];
          if (*pC < x[1]) {
            *--pC += x[0]+1;
            if (*pC <= x[0]) ++y;
          } else {
            *--pC += x[0];
            if (*pC < x[0]) ++y;
          }
        }
      } else {
        *--pC += x[2];
        if (*pC < x[2]) {
          *--pC += x[1]+1;
          if (*pC <= x[1]) {
            *--pC += x[0]+1;
            if (*pC <= x[0]) ++y;
          } else {
            *--pC += x[0];
            if (*pC < x[0]) ++y;
          }
        } else {
          *--pC += x[1];
          if (*pC < x[1]) {
            *--pC += x[0]+1;
            if (*pC <= x[0]) ++y;
          } else {
            *--pC += x[0];
            if (*pC < x[0]) ++y;
          }
        }
      }
      *--pC += y;
      if (*pC < y)
        while (++(*--pC) == 0);
    } while (pA != pE);
    n -= 2; pB += n-2; pA += n;
  }
}

void Natural::cmul(const Digit* pA, size_t sA, const Digit* pB, const size_t sB,
                   Digit* pC) const
// Algorithm:  x.cmul(r, n, s, m, t)
//             Let a,b,c in Natural where c.p not in {a.p, b.p}.
// Input:      x in Natural, n,m in size_t where n >= m > 0,
//             r,r+n in [a.p, a.p+L(a)],
//             s,s+m in [b.p, b.p+L(b)],
//             t,t+n+m in [c.p, c.p+L(c)].
// Output:     [t, t+n+m[ = [r, r+n[ * [s, s+m[ ||
{
  CONDITION(sA >= sB && sB > 0);

  const Digit* pE = pC+sB;
  FILL_ZERO(pC, pE);
  const Digit* pF = pB;
  pE = pA; pA += sA; pB += sB; pC += sA;
  mul(pE, pA, pC, *--pB);
  while (pB > pF) muladd(pE, pA, --pC, *--pB);
}

void Natural::mul(const Digit* pA, const Digit* pB, const size_t sz, Digit* pC) const
// Algorithm:  c.mul(r, s, n, t)
//             Let a,b in Natural.
// Input:      c in Natural where not c.p in {a.p, b.p},
//             n in size_t where n > KARATSUBA_MUL_MARK,
//             r,r+n in [a.root, a.p+L(a)],
//             s,s+n in [b.root, b.p+L(b)],
//             t,t+2*n in [c.root, c.p+L(c)].
// Output:     [t, t+2*n[ = [r, r+n[ * [s, s+n[ ||
{
  struct stack {
    const Digit* a;
    const Digit* b;
    Digit* c;
    char   cr,d;
  };
      
  size_t n = sz;
  size_t i = 0;
  Digit* h = NOTHROW_NEW Digit[8*n];
  if (!h) errmsg(2, "(mul)");
  while (n >>= 1) ++i;
  stack* st = NOTHROW_NEW stack[i];
  for (n = 0; n < i; ++n) st[n].d = 0;
  n = sz;
  --i;
  
  do {
    if (n <= KARATSUBA_MUL_MARK) {
      cmul(pA, n, pB, n, pC);
      h -= 4*n;
      pA = st[++i].a; pB = st[i].b; pC = st[i].c;
    }
    switch (++st[i].d) {
      case 1:
        st[i].cr = 0;
        if (n&1) {
          ++st[i].cr;
          digitmul(*pA, *pB, *pC, pC[1]);
          ++pA; ++pB; pC += 2;
        }
        n >>= 1; st[i].a = pA; st[i].b = pB; st[i].c = pC;
        h += 4*n; --i;
        break;
      case 2:
        pA += n; pB += n; pC += 2*n;
        h += 4*n; --i;
        break;
      case 3:
        if (abs(h, pA, pA+n, n)*abs(h+n, pB, pB+n, n) == -1) ++st[i].d;
        pA = h; pB = h+n; pC = h+2*n;
        h += 4*n; --i;
        break;
      case 4:
        pC += n;
        {
          Digit k = (sub_no_dec(h+2*n, h+4*n, pC+n) == true)
                  + (sub_no_dec(h+2*n, h+4*n, pC+3*n) == true);   // k >= 1!
          k -= (sub_no_dec(pC, pC+2*n, h+4*n) == true);
          Digit* pPos = pC;
          *--pPos += k;
          if (*pPos < k)
            while (++(*--pPos) == 0);
        }
        st[i].d = 0;
        if (st[i].cr) {
          muladd(pB, pB+2*n, pC+n, *(pA-1));
          muladd(pA, pA+2*n, pC+n, *(pB-1));
          n <<= 1; ++n;
        } else n <<= 1;
        if (n < sz) { pA = st[++i].a; pB = st[i].b; pC = st[i].c; h -= 4*n; }
        break;
      case 5:
        pC += n;
        {
          Digit k = (add_no_inc(h+2*n, h+4*n, pC+n) == true)
                  + (add_no_inc(h+2*n, h+4*n, pC+3*n) == true);
          k += (add_no_inc(pC, pC+2*n, h+4*n) == true);
          Digit* pPos = pC;
          *--pPos += k;
          if (*pPos < k)
            while (++(*--pPos) == 0);
        }
        st[i].d = 0;
        if (st[i].cr) {
          muladd(pB, pB+2*n, pC+n, *(pA-1));
          muladd(pA, pA+2*n, pC+n, *(pB-1));
          n <<= 1; ++n;
        } else n <<= 1;
        if (n < sz) { pA = st[++i].a; pB = st[i].b; pC = st[i].c; h -= 4*n; }
        break;
    }
  } while (n < sz);
  delete[] h;
  delete[] st;
}

void Natural::mul(const Digit* pA, size_t sA, const Digit* pB, const size_t sB,
                  Digit* pC) const
// Algorithm:  c.mul(r, n, s, m, t)
//             Let a,b in Natural.
// Input:      c in Natural where not c.p in {a.p, b.p},
//             n,m in size_t where n >= m > 0,
//             r,r+n in [a.root, a.p+L(a)],
//             s,s+m in [b.root, b.p+L(b)],
//             t,t+n+m in [c.root, c.p+L(c)].
// Output:     [t, t+n+m[ = [r, r+n[ * [s, s+m[ ||
{
  CONDITION(sA >= sB && sB > 0);

  if (sB <= KARATSUBA_MUL_MARK) cmul(pA, sA, pB, sB, pC);
  else {
    sA -= sB;
    if (sA) {
      pA += sA;
      const Digit* pE = pC+sA;
      FILL_ZERO(pC, pE);
    }
    mul(pA, pB, sB, pC);
    Digit* pT = NOTHROW_NEW Digit[sA+sB];
    if (!pT) errmsg(2, "(mul)");
    while (sA >= sB) {
      pA -= sB; pC -= sB; sA -= sB;
      mul(pA, pB, sB, pT);
      if (add_no_inc(pC, pC+2*sB, pT+2*sB))
        for (Digit* pPos = pC; ++(*--pPos) == 0;);
    }
    if (sA) {
      pA -= sA; pC -= sA;
      mul(pB, sB, pA, sA, pT);
      sA += sB;
      add_no_inc(pC, pC+sA, pT+sA);
    }
    delete[] pT;
  }
}

Natural::Natural(const Digit a)
// Algorithm:  c := Natural(a)
// Input:      a in Digit.
// Output:     c in Natural such that c = a ||
{
  const size_t sT = NaturalSize+DELTA;
  Digit* pT = NOTHROW_NEW Digit[sT];
  root = pT; size = 1;
  if (!pT) errmsg(2, "(default constructor)");
  Digit* pE = pT+sT;
  *--pE = a; p = pE;
  FILL_ZERO(pT, pE);

  NATURALCONDITION(*this);
}

Natural::Natural(const Natural& a)
// Algorithm:  c := Natural(a)
// Input:      a in Natural.
// Output:     c in Natural such that c = a ||
{
  NATURALCONDITION(a);

  const size_t sA = a.size;
  Digit* pT = NOTHROW_NEW Digit[sA+DELTA];
  root = pT; size = sA;
  if (!pT) errmsg(2, "(copy constructor)");
  FILL_DELTA(pT);
  p = pT += DELTA;
  const Digit* pA = a.p;
  const Digit* pE = pA+sA;
  COPY(pT, pA, pA, pE);

  NATURALCONDITION(*this);
}

Natural::Natural(const char* a, const Digit b)
// Algorithm:  c := Natural(a, b)
// Input:      a in String, b in Digit.
// Output:     c in Natural such that c = a ||
{
  Digit* pT = NOTHROW_NEW Digit[DELTA];
  p = root = pT; size = DELTA;
  if (!pT) errmsg(2, "(constructor, string conversion)");
  atoN(a, b);
}

Natural::~Natural()
{
//  NATURALCONDITION(*this);

  delete[] root;
}

Natural& Natural::operator+=(const Natural& a)
// Algorithm:  c := c += a
// Input:      a,c in Natural.
// Output:     c in Natural such that c := c+a ||
{
  NATURALCONDITION(*this);
  NATURALCONDITION(a);

  const size_t sT = size;
  const size_t sA = a.size;
  const Digit* pA = a.p;
  Digit* rT = root;
  Digit* pT = p;
  if (sA == sT) {
    if (pT-rT == 1) {
      Digit* pC = NOTHROW_NEW Digit[sT+DELTA];
      root = pC;
      if (!pC) errmsg(2, "(operator+=)");
      FILL_DELTA(pC);
      p = pC += DELTA;
      add_with_inc(pC, pC+sT, pT+sT, pA+sT);
      delete[] rT;
    } else add_with_inc(pT, pT+sT, pA+sT);
  } else if (sA < sT) {
    if (pT-rT == 1) {
      Digit* pC = NOTHROW_NEW Digit[sT+DELTA];
      root = pC;
      if (!pC) errmsg(2, "(operator+=)");
      FILL_DELTA(pC);
      p = pC += DELTA;
      const Digit* pE = pT+sT-sA;
      COPY(pC, pT, pT, pE);
      if (pA+sA > a.root+sT) add_with_inc(pC, pC+sA, pT+sA, pA+sA);
      else if (add_no_inc(pC, pC+sA, pT+sA, pA+sA)) inc(pC);
      delete[] rT;
    } else {
      pT += sT;
      if (add_no_inc(pT-sA, pT, pA+sA)) inc(pT-sA);
    }
  } else {      // sA > sT
    if (pT+sT <= rT+sA+1) {                // (rootsize() <= a.size+1)?
      Digit* pC = NOTHROW_NEW Digit[sA+DELTA];
      root = pC;
      if (!pC) errmsg(2, "(operator+=)");
      FILL_DELTA(pC);
      p = pC += DELTA; size = sA;
      const Digit* pE = pA+sA-sT;
      COPY(pC, pA, pA, pE);
      if (add_no_inc(pC, pC+sT, pT+sT, pA+sT)) inc(pC);
      delete[] rT;
    } else {
      Digit* pE = pT+sT-sA;
      p = pE; size = sA;
      COPY(pE, pA, pE, pT);
      if (add_no_inc(pE, pE+sT, pA+sT)) inc(pE);    // because COPY
    }
  }

  NATURALCONDITION(*this);

  return *this;
}

Natural& Natural::operator-=(const Natural& a)
// Algorithm:  c := c -= a
// Input:      a,c in Natural.
// Output:     c in Natural such that c := c-a ||
{
  NATURALCONDITION(*this);
  NATURALCONDITION(a);

  const size_t sT = size;
  const size_t sA = a.size;
  const Digit* pA = a.p;
  Digit* pT = p;
  if (sA <= sT) sub(pT+sT-sA, pT+sT, pA+sA);
  else errmsg(3, "(sub)");
  normalize();

  NATURALCONDITION(*this);

  return *this;
}

Natural& Natural::operator&=(const Natural& a)
// Algorithm:  c := c &= a
// Input:      a,c in Natural.
// Output:     c in Natural such that c := c and a ||
{
  NATURALCONDITION(*this);

  const Digit* pA = a.p;
  const size_t sA = a.size;
  const Digit* pE = pA+sA;
  size_t sT = size;
  Digit* pT = p;
  if (sT > sA) {
    const Digit* pE = pT+sT-sA;
    FILL_ZERO(pT, pE);
    p = pT; size = sA;
  } else pA += sA-sT;
  do *pT++ &= *pA++; while (pA != pE);
  normalize();

  NATURALCONDITION(*this);

  return *this;
}

Natural& Natural::operator|=(const Natural& a)
// Algorithm:  c := c |= a
// Input:      a,c in Natural.
// Output:     c in Natural such that c := c or a ||
{
  NATURALCONDITION(*this);

  const Digit* pA = a.p;
  size_t sA = a.size;
  Digit* rT = root;
  size_t sT = size;
  Digit* pT = p;
  if (sA > sT)
    if (rT+sA >= pT+sT) {           // a.size >= rootsize()
      Digit* pC = NOTHROW_NEW Digit[sA+DELTA];
      root = pC; size = sA;
      if (!pC) errmsg(2, "(operator|=)");
      FILL_DELTA(pC);
      p = pC += DELTA;
      const Digit* pE = pA+sA-sT;
      COPY(pC, pA, pA, pE);
      pE += sT;
      do *pC++ = *pT++ | *pA++; while (pA != pE);
      delete[] rT;
    } else {
      size = sA; sA -= sT; p = pT -= sA;
      const Digit* pE = pA+sA;
      COPY(pT, pA, pA, pE);
      pE += sT;
      do *pT++ |= *pA++; while (pA != pE);
    }
  else {
    const Digit* pE = pA+sA;
    pT += sT-sA;
    do *pT++ |= *pA++; while (pA != pE);
  }

  NATURALCONDITION(*this);

  return *this;
}

Natural& Natural::operator^=(const Natural& a)
// Algorithm:  c := c ^= a
// Input:      a,c in Natural.
// Output:     c in Natural such that c := c xor a ||
{
  NATURALCONDITION(*this);

  const Digit* pA = a.p;
  size_t sA = a.size;
  Digit* rT = root;
  size_t sT = size;
  Digit* pT = p;
  if (sA > sT)
    if (rT+sA >= pT+sT) {       // a.size >= rootsize()
      Digit* pC = NOTHROW_NEW Digit[sA+DELTA];
      root = pC; size = sA;
      if (!pC) errmsg(2, "(operator^=)");
      FILL_DELTA(pC);
      p = pC += DELTA;
      const Digit* pE = pA+sA-sT;
      COPY(pC, pA, pA, pE);
      pE += sT;
      do *pC++ = *pT++ ^ *pA++; while (pA != pE);
      delete[] rT;
    } else {
      size = sA; sA -= sT; p = pT -= sA;
      const Digit* pE = pA+sA;
      COPY(pT, pA, pA, pE);
      pE += sT;
      do *pT++ ^= *pA++; while (pA != pE);
    }
  else if (sA < sT) {
    const Digit* pE = pA+sA;
    pT += sT-sA;
    do *pT++ ^= *pA++; while (pA != pE);
  } else {
    const Digit* pE = pA+sA;
    do *pT++ ^= *pA++; while (pA != pE);
    normalize();
  }

  NATURALCONDITION(*this);

  return *this;
}

int Natural::compare(const Natural& b) const
// Algorithm:  c := a.compare(b)
// Input:      a,b in Natural.
// Output:     c in int such that c = sign(a - b) ||
{
  NATURALCONDITION(*this);
  NATURALCONDITION(b);

  const size_t sT = size;
  const size_t sB = b.size;
  if (sT > sB) return 1;
  else if (sT < sB) return -1;

  Digit* pT = p;
  Digit* pB = b.p;
  for (const Digit* pE = pT+sT; *pT == *pB; ++pB)
    if (++pT == pE) return 0;
  return (*pT < *pB)? -1 : 1;
}

void swap(Natural& a, Natural& b)
// Algorithm:  swap(a, b)
// Input:      a,b in Natural.
// Output:     a,b in Natural such that t := a, a := b, b := t
//             where t in Natural ||
{
  NATURALCONDITION(a);
  NATURALCONDITION(b);

  Digit* pA = a.p; a.p = b.p; b.p = pA;
  pA = a.root; a.root = b.root; b.root = pA;
  size_t sz = a.size; a.size = b.size; b.size = sz;

  NATURALCONDITION(a);
  NATURALCONDITION(b);
}

void Natural::split(const size_t n, Natural& a, Natural& b) const
// Algorithm:  c.split(n, a, b)
// Input:      a,b,c in Natural, n in size_t where not a = b.
// Output:     a,b in Natural such that a = [c/2^(BETA*n)], b = c - a*2^(BETA*n) ||
{
  NATURALCONDITION(*this);

  Digit* pA = a.p;
  Digit* pB = b.p;
  size_t sT = size;
  if (pA == pB) a.errmsg(5, "(split)");
  if (sT <= n) { b = *this; a = 0; return; }
  if (n == 0) { a = *this; b = 0; return; }

  const Digit* pT = p;
  if (pT == pB) {
    sT -= n;
    pA = a.setsize(sT);
    const Digit* pE = pB+sT;
    MOVE(pA, pB, pB, pE);
    b.p = pB; b.size = n;
  } else if (pT == pA) {
    pB = b.setsize(n);
    pA += sT;
    Digit* pE = pA-n;
    COPY(pB, pE, pE, pA);
    pE -= n; sT -= n;
    MOVE_BACKWARD(pA, pE, pE, pT);
    a.p = pA; a.size = sT;
    FILL_ZERO(pE, pA);
  } else {
    sT -= n;
    pA = a.setsize(sT);
    const Digit* pE = pT+sT;
    COPY(pA, pT, pT, pE);
    pB = b.setsize(n);
    pE += n;
    COPY(pB, pT, pT, pE);
  }
  b.normalize();

  NATURALCONDITION(a);
  NATURALCONDITION(b);
}

void Natural::add(const Natural& a, const Natural& b)
// Algorithm:  c.add(a, b)
// Input:      a,b in Natural where not a.p = c.p and not b.p = c.p.
// Output:     c in Natural such that c = a+b ||
{
  NATURALCONDITION(a);
  NATURALCONDITION(b);
  CONDITION(a.p != p && b.p != p);

  const size_t sA = a.size;
  const size_t sB = b.size;
  const Digit* pA = a.p;
  const Digit* pB = b.p;
  if (sA == sB) {
    Digit* pT = setsize(sA);
    add_with_inc(pT, pT+sA, pA+sA, pB+sA);
  } else if (sA > sB) {
    Digit* pT = setsize(sA);
    const Digit* pE = pA+sA-sB;
    COPY(pT, pA, pA, pE);
    if (pB+sB > b.root+sA) add_with_inc(pT, pT+sB, pA+sB, pB+sB);
    else if (add_no_inc(pT, pT+sB, pA+sB, pB+sB)) inc(pT);
  } else {
    Digit* pT = setsize(sB);
    const Digit* pE = pB+sB-sA;
    COPY(pT, pB, pB, pE);
    if (pA+sA > a.root+sB) add_with_inc(pT, pT+sA, pA+sA, pB+sA);
    else if (add_no_inc(pT, pT+sA, pA+sA, pB+sA)) inc(pT);
  }

  NATURALCONDITION(*this);
}

void Natural::sub(const Natural& a, const Natural& b)
// Algorithm:  c.sub(a, b)
// Input:      a,b in Natural.
// Output:     c in Natural such that c = a-b ||
{
  NATURALCONDITION(a);
  NATURALCONDITION(b);

  size_t sA = a.size;
  const size_t sB = b.size;
  const Digit* pA = a.p;
  const Digit* pB = b.p;
  Digit* pT = p;
  if (pA == pT) *this -= b;
  else if (pB == pT) {
    if (sB == sA) {
      if (abs(pT, pA, pT, sB) == -1) errmsg(3, "(sub)");
      normalize();
    } else if (sB < sA) {
      Digit* rT = root;
      if (pT+sB <= rT+sA) {
        pT = NOTHROW_NEW Digit[sA+DELTA];
        root = pT; size = sA;
        if (!pT) errmsg(2, "(sub)");
        FILL_DELTA(pT);
        p = pT += DELTA;
        const Digit* pE = pA+sA-sB;
        COPY(pT, pA, pA, pE);
        sub(pT, pT+sB, pA+sB, pB+sB);
        normalize();
        delete[] rT;
      } else {
        size = sA;
        sA -= sB;
        p = pT -= sA;
        const Digit* pE = pA+sA;
        COPY(pT, pA, pA, pE);
        sub(pT, pT+sB, pA+sB, pT+sB);
        normalize();
      }
    } else errmsg(3, "(sub)");
  } else if (sA == sB) {
    pT = setsize(sA);
    if (abs(pT, pA, pB, sA) == -1) errmsg(3, "(sub)");
    normalize();
  } else if (sA > sB) {
    pT = setsize(sA);
    const Digit* pE = pA+sA-sB;
    COPY(pT, pA, pA, pE);
    sub(pT, pT+sB, pA+sB, b.p+sB);
    normalize();
  } else errmsg(3, "(sub)");

  NATURALCONDITION(*this);
}

void Natural::mul(const Natural& a, const Natural& b)
// Algorithm:  c.mul(a, b)
// Input:      a,b in Natural.
// Output:     c in Natural such that c = a*b ||
{
  NATURALCONDITION(a);
  NATURALCONDITION(b);
  NATURAL_FOR_CHECK(_a, a);
  NATURAL_FOR_CHECK(_b, b);

  const size_t sA = a.size;
  const size_t sB = b.size;
  const Digit* pA = a.p;
  const Digit* pB = b.p;

  if (sA <= 2 && sB <= 2) {
    if (sA == 1) {
      const Digit x = *pA;
      const Digit y = *pB;
      if (sB == 1) {
        Digit* pT = setsize(sB+1);
        digitmul(x, y, pT[0], pT[1]);
      } else {
        const Digit z = pB[1];
        Digit* pT = setsize(sB+1);
        digitmul(x, y, z, pT);
      }
    } else {
      const Digit x0 = pA[0];
      const Digit x1 = pA[1];
      const Digit y0 = pB[0];
      if (sB == 1) {
        Digit* pT = setsize(sB+2);
        digitmul(y0, x0, x1, pT);
      } else {
        const Digit y1 = pB[1];
        Digit* pT = setsize(sB+2);
        digitmul(x0, x1, y0, y1, pT);
      }
    }
    normalize();
  } else if (sA >= FFT_MUL_MARK && sB >= FFT_MUL_MARK
             && max(sA, sB) <= FFT::max_size()) {
    FFT f(a, false);
    FFT g(b, true, &f);
    g.mul(*this);
  } else {
    Digit* pT = p;
    if (sA < sB)
      if (pA == pT) {
        Natural t(a);
        pT = setsize(sA+sB);
        mul(pB, sB, t.p, sA, pT);
      } else if (pB == pT) {
        Natural t(b);
        pT = setsize(sA+sB);
        mul(t.p, sB, pA, sA, pT);
      } else {
        pT = setsize(sA+sB);
        mul(b.p, sB, pA, sA, pT);
      }
    else
      if (pA == pT) {
        Natural t(a);
        pT = setsize(sA+sB);
        mul(t.p, sA, pB, sB, pT);
      } else if (pB == pT) {
        Natural t(b);
        pT = setsize(sA+sB);
        mul(pA, sA, t.p, sB, pT);
      } else {
        pT = setsize(sA+sB);
        mul(pA, sA, pB, sB, pT);
      }
    normalize();
  }

  NATURAL_FOR_CHECK(_c, _a+_b);
  CONDITION((*this*2) == (_c*_c-_a*_a-_b*_b));
  NATURALCONDITION(*this);
}

void Natural::sqr(const Natural& a)
// Algorithm:  b.sqr(a)
// Input:      a,b in Natural.
// Output:     b in Natural such that b = a^2 ||
{
  NATURALCONDITION(a);

  size_t n = a.length();
  if (n == 1) {
    Digit x,y,z = *a.p;
    a.digitmul(z, z, x, y);
    setsize(2);
    if (x) { *p = x; p[1] = y; }
    else { *p = 0; *++p = y; --size; }
    return;
  }
  if (n <= KARATSUBA_SQR_MARK) {
    if (a.p != p) {
      setsize(2*n);
      sqr(a.p, p, n);
    } else {
      const Natural c(a);
      setsize(2*n);
      sqr(c.p, p, n);
    }
    normalize();
    return;
  }
  if (n >= FFT_SQR_MARK && n <= FFT::max_size()) {
    FFT f(a); f.sqr(*this);
    return;
  }
  
  struct stack {
    const Digit* a;
    Digit* b;
    char   c,d;
  };
    
  const Natural c(a);
  setsize(2*n);
  size_t i = 0;
  Digit* h = NOTHROW_NEW Digit[6*n];
  if (!h) errmsg(2, "(sqr)");
  while (n >>= 1) ++i;
  stack* st = NOTHROW_NEW stack[i];
  for (n = 0; n < i; ++n) st[n].d = 0;
  const size_t sC = c.size;
  n = sC;
  --i;
  
  const Digit* pA = c.p;
  Digit* pB = p;
  do {
    if (n <= KARATSUBA_SQR_MARK) {
      sqr(pA, pB, n);
      h -= 3*n;
      pA = st[++i].a; pB = st[i].b;
    }
    switch (++st[i].d) {
      case 1:
        st[i].c = 0;
        if (n&1) {
          ++st[i].c;
          digitmul(*pA, *pA, *pB, pB[1]);
          ++pA; pB += 2;
        }
        n >>= 1; st[i].a = pA; st[i].b = pB;
        h += 3*n; --i;
        break;
      case 2:
        pA += n; pB += 2*n;
        h += 3*n; --i;
        break;
      case 3:
        a.abs(h, pA, pA+n, n);
        pA = h; pB = h+n;
        h += 3*n; --i;
        break;
      case 4:
        pB += n;
        {
          int i = sub_no_dec(h+n, h+3*n, pB+n) + sub_no_dec(h+n, h+3*n, pB+3*n);
          i -= sub_no_dec(pB, pB+2*n, h+3*n);
          if (i == 1)
            for (Digit* pPos = pB; ++(*--pPos) == 0;);
          else if (i == 2) {
            Digit* pPos = pB;
            *--pPos += 2;
            if (*pPos < 2)
              while (++(*--pPos) == 0);
          }
        }
        st[i].d = 0;
        if (st[i].c) {
          const Digit d = *(pA-1);
          muladd(pA, pA+2*n, pB+n, d);
          muladd(pA, pA+2*n, pB+n, d);
          n <<= 1; ++n;
        } else n <<= 1;
        if (n < sC) { pA = st[++i].a; pB = st[i].b; h -= 3*n; }
        break;
    }
  } while (n < sC);
  delete[] h;
  delete[] st;
  normalize();

  NATURALCONDITION(*this);
}

void div(const Natural& a, const Natural& b, Natural& q, Natural& r)
// Algorithm:  div(a, b, q, r)
// Input:      a,b,q,r in Natural where not b = 0, not q = r.
// Output:     q,r in Natural such that q = [a/b], r = a - q*b ||
{
  NATURALCONDITION(a);
  NATURALCONDITION(b);
  NATURAL_FOR_CHECK(_a, a);
  NATURAL_FOR_CHECK(_b, b);

  if (b == 0) b.errmsg(4, "(div)");
  if (q.p == r.p) r.errmsg(5, "(div)");
  switch (b.compare(a)) {
    case 0: r = 0; q = 1; return;
    case 1: r = a; q = 0; return;
  }
  const size_t sA = a.size;
  const size_t sB = b.size;
  if (sB == 1) {
    const Digit c = *b.p;     // important! (q = b)
    q = a; r = q.mod_div(c);
    return;
  } else if (sB >= NEWTON_DIV_MARK2 && sA >= sB+NEWTON_DIV_MARK1) {
    Natural* a2 = 0;
    Natural* b2 = 0;
    if (&q == &a) {
      a2 = new Natural(a);
      if (!a2) a.errmsg(2, "(div)");
    }
    if (&q == &b) {
      b2 = new Natural(b);
      if (!b2) b.errmsg(2, "(div)");
    }
    if (&r == &a) {
      a2 = new Natural(a);
      if (!a2) a.errmsg(2, "(div)");
    }
    if (a2)
      if (b2) q.div(*a2, *b2, r);
      else q.div(*a2, b, r);
      else if (b2) q.div(a, *b2, r);
      else q.div(a, b, r);
    delete b2;
    delete a2;
    return;
  }

  Digit m = 1;
  int m2 = 0;
  Natural c(3, b);    // because subpos
  r = a;
  Digit d = *c.p;
  if (d != GAMMA) {
    m = GAMMA/(d+1);
    c *= m; r *= m;
    d = *c.p;
  }
  if (d == GAMMA) { m2 = 1; c *= GAMMA; r *= GAMMA; d = *c.p; }
  if (d <= GAMMA/2) {
    m2 += 2; c <<= 1; r <<= 1; d = *c.p;
    if (d == GAMMA) { m2 += 4; c *= GAMMA; r *= GAMMA; d = *c.p; }
  }
  CONDITION(d > GAMMA/2);
  const Digit d2 = d+1;
  const size_t sQ = r.size-c.size+(*r.p >= d);
  Digit* pQ = q.setsize(sQ);
  const Digit* pE = pQ+sQ;
  FILL_ZERO(pQ, pE);
  pQ -= sQ;
  size_t sC = c.size;
  const Digit* pC = c.p+sC;
  do {
    Digit* pR = r.p;
    const Digit k = *pR;
    const size_t sR = r.size;
    r.size = sC;
    bool e = (k < d || k == d && r < c);
    if (e) { // sR > sC is always fulfilled if r < c
      Digit s,t;
      q.digitdiv(*pR, pR[1], d2, s, t);
      r.size = sC+1;
      if (s) r.mulsub(c, s);
      else r.normalize();
      while (r >= c) {
        subpos(r.p, r.p + r.size, pC);
        r.normalize();
        ++s;
      }
      *pQ = s;
    } else {
      ++(*pQ);
      subpos(pR, pR+sC, pC);
      r.normalize();
    }
    const Digit* pF = r.p;
    r.size = pR+sR-pF; // restore the size back
    const size_t l = pE-pQ;
    size_t m = (pF-pR) - (*pF >= d);
    if (!e || m == 0) ++m;
    pQ += (l <= m)? l : m;
  } while (pQ != pE);
  q.normalize(); r.normalize();
  if (r == c) {
    r = 0; ++q;
  } else {
    r /= m;
    if (m2&1) r /= GAMMA;
    if (m2&2) r >>= 1;
    if (m2&4) r /= GAMMA;
  }
  
  CONDITION(r < _b);
  CONDITION(r+q*_b == _a);
  NATURALCONDITION(q);
  NATURALCONDITION(r);
}

void Natural::div(const Natural& a, Natural b, Natural& r)
// Algorithm:  c.div(a, b, r)
// Input:      a,b,r in Natural where not b = 0.
// Output:     c,r in Natural such that c = [a/b], r = a - c*b ||
{
  NATURALCONDITION(a);
  NATURALCONDITION(b);
  NATURAL_FOR_CHECK(_b, b);

  const size_t sA = a.size;
  const size_t sB = b.size;
  if (sB == 1) {
    const Digit x = *b.p;
    if (x == 0) b.errmsg(4, "(div)");
    *this = a / x;
    return;
  } else if (sB < NEWTON_DIV_MARK2 || sA < sB+NEWTON_DIV_MARK1) {
    Natural t(sA+DELTA, ' ');
    ::div(a, b, *this, t);
    return;
  }

  CONDITION(b.length() >= 2 && a.length() > b.length());

  Natural* a2 = 0;
  if (this == &a) {
    a2 = new Natural(a); // if this == &a
    if (!a2) a.errmsg(2, "(div)");
  }
  const size_t l = BETA-1-size_t(log2(b.highest()));
  b <<= l;
  Digit q,q2 = b.highest();
  if (q2 == ~(GAMMA/2)) q = GAMMA;
  else digitdiv(~(GAMMA/2), 0, q2, q, q2);
  Natural y = q;
  size_t m;
  size_t* s = quad_convergence_sizes(a.size-sB+3, m);
  size_t sY = s[--m];
  do {
    const size_t k = s[--m];
    *this = y*y;
    if (size > k) this->rmove(size-k);
    if (sY < sB) { b.size = sY+1; *this *= b; b.size = sB; }
    else *this *= b;
    *this >>= (size-k)*BETA - 1;
    y <<= (k-sY)*BETA + 1; y -= *this;
    sY = k;
  } while (m);
  sY += 2;
  *this = y*y; this->rmove(size-sY-2);
  *this *= b; *this >>= (size-sY)*BETA - 1;
  y <<= 1+2*BETA; y -= *this;
  if (a2) *this = y * (*a2);
  else *this = y * a;
  *this >>= (sB+sY)*BETA-l-1;
  delete[] s;
  b >>= l;
  y = (*this) * b;
  if (a2) r = (*a2) - y;
  else r = a - y;
  if (r == b) {
    r = 0; ++(*this);
  }
  delete a2;

  CONDITION(r < _b);
  NATURALCONDITION(*this);
}

void Natural::bitwise_and(const Natural& a, const Natural& b)
// Algorithm:  c.bitwise_and(a, b)
// Input:      a,b in Natural where not a.p = c.p and not b.p = c.p.
// Output:     c in Natural such that c = a and b ||
{
  NATURALCONDITION(a);
  NATURALCONDITION(b);
  CONDITION(a.p != p && b.p != p);

  const Digit* pA = a.p;
  const size_t sA = a.size;
  const Digit* pB = b.p;
  const size_t sB = b.size;
  size_t sT;
  if (sA >= sB) { pA += sA-sB; sT = sB; }
  else { pB += sB-sA; sT = sA; }
  Digit d;
  do d = *pA++ & *pB++; while (d == 0 && --sT);
  if (sT) {
    Digit* pT = setsize(sT);
    *pT = d;
    const Digit* pE = pT+sT;
    while (++pT != pE) *pT = *pA++ & *pB++;
  } else *this = 0;

  NATURALCONDITION(*this);
}

void Natural::bitwise_or(const Natural& a, const Natural& b)
// Algorithm:  c.bitwise_or(a, b)
// Input:      a,b in Natural where not a.p = c.p and not b.p = c.p.
// Output:     c in Natural such that c = a or b ||
{
  NATURALCONDITION(a);
  NATURALCONDITION(b);
  CONDITION(a.p != p && b.p != p);

  const Digit* pA = a.p;
  const size_t sA = a.size;
  const Digit* pB = b.p;
  const size_t sB = b.size;
  if (sA > sB) {
    Digit* pT = setsize(sA);
    const Digit* pE = pA+sA-sB;
    COPY(pT, pA, pA, pE);
    pE += sB;
    do *pT++ = *pA++ | *pB++; while (pA != pE);
  } else if (sA < sB) {
    Digit* pT = setsize(sB);
    const Digit* pE = pB+sB-sA;
    COPY(pT, pB, pB, pE);
    pE += sA;
    do *pT++ = *pA++ | *pB++; while (pB != pE);
  } else {
    Digit* pT = setsize(sA);
    const Digit* pE = pA+sA;
    do *pT++ = *pA++ | *pB++; while (pA != pE);
  }

  NATURALCONDITION(*this);
}

void Natural::bitwise_xor(const Natural& a, const Natural& b)
// Algorithm:  c.bitwise_xor(a, b)
// Input:      a,b in Natural where not a.p = c.p and not b.p = c.p.
// Output:     c in Natural such that c = a xor b ||
{
  NATURALCONDITION(a);
  NATURALCONDITION(b);
  CONDITION(a.p != p && b.p != p);

  const Digit* pA = a.p;
  const size_t sA = a.size;
  const Digit* pB = b.p;
  const size_t sB = b.size;
  if (sA > sB) {
    Digit* pT = setsize(sA);
    const Digit* pE = pA+sA-sB;
    COPY(pT, pA, pA, pE);
    pE += sB;
    do *pT++ = *pA++ ^ *pB++; while (pA != pE);
  } else if (sA < sB) {
    Digit* pT = setsize(sB);
    const Digit* pE = pB+sB-sA;
    COPY(pT, pB, pB, pE);
    pE += sA;
    do *pT++ = *pA++ ^ *pB++; while (pB != pE);
  } else {
    Digit* pT = setsize(sA);
    const Digit* pE = pA+sA;
    do *pT++ = *pA++ ^ *pB++; while (pA != pE);
    normalize();
  }

  NATURALCONDITION(*this);
}

void Natural::bitwise_not(const Natural& a)
// Algorithm:  b.bitwise_not(a)
// Input:      a in Natural.
// Output:     b in Natural such that b = not a ||
{
  NATURALCONDITION(a);

  const size_t sA = a.size;
  const Digit* pA = a.p;
  const Digit* pE = pA+sA;
  Digit* pT = setsize(sA);
  do *pT++ = ~*pA++; while (pA != pE);
  normalize();

  NATURALCONDITION(*this);
}

void Natural::lshift(const Natural& a, size_t b)
// Algorithm:  c.lshift(a, b)
// Input:      a in Natural, b in size_t where not a.p = c.p.
// Output:     c in Natural such that c = a*2^b ||
{
  NATURALCONDITION(a);
  CONDITION(a.p != p);

  const Digit* pA = a.p;
  Digit d = *pA;
  if (d == 0) *this = 0;
  else {
    const size_t sA = a.size;
    const size_t b2 = b/BETA;
    size_t sT = sA+b2+1;
    Digit* pT = setsize(sT);
    if (b2) {
      const Digit* pE = pT+sT;
      pT += sA+1;
      FILL_ZERO(pT, pE);
      pT -= sT;
    }
    b %= BETA;
    if (b) {
      const Digit b2 = BETA-b;
      const Digit* pE = pA+sA;
      Digit e = d >> b2;
      *pT++ = e;
      if (e == 0) { size = --sT; p = pT; }
      while (++pA != pE) {
        e = *pA;
        *pT++ = (d << b) | (e >> b2);
        if (++pA == pE) { d = e; break; }
        d = *pA;
        *pT++ = (e << b) | (d >> b2);
      }
      *pT = d << b;
    } else {
      *pT = 0;
      size = --sT; p = ++pT;
      const Digit* pE = pA+sA;
      COPY(pT, pA, pA, pE);
    }
  }

  NATURALCONDITION(*this);
}

void Natural::lmove(size_t b)
// Algorithm:  a.lmove(b)
// Input:      a in Natural, b in size_t.
// Output:     a in Natural such that a := a*2^(BETA*b) ||
{
  NATURALCONDITION(*this);

  Digit* pT = p;
  size_t sT = size;
  Digit c = *pT;
  if (c) {
    ++b;
    Digit* rT = root;
    size_t sC = sT+b;
    Digit* pC = pT-b;
    if (rT+b >= pT) {
      root = pC = NOTHROW_NEW Digit[sC+DELTA];
      if (!pC) errmsg(2, "(lmove)");
      FILL_DELTA(pC);
      pC += DELTA;
    }
    *pC = 0;
    size = --sC; p = ++pC;
    const Digit* pE = pT+sT;
    COPY(pC, pT, pT, pE);
    if (--b) {
      const Digit* pE = pC+b;
      FILL_ZERO(pC, pE);
    }
    if (root != rT) delete[] rT;
  }

  NATURALCONDITION(*this);
}

void Natural::rshift(const Natural& a, size_t b)
// Algorithm:  c.rshift(a, b)
// Input:      a in Natural, b in size_t where not a.p = c.p.
// Output:     c in Natural such that c = [a/2^b] ||
{
  NATURALCONDITION(a);
  CONDITION(a.p != p);

  const Digit* pA = a.p;
  const size_t sA = a.size;
  const size_t b2 = b/BETA;
  if (b2 >= sA) *this = 0;
  else {
    b %= BETA;
    size_t sT = sA-b2;
    Digit* pT = setsize(sT);
    if (b) {
      pT += sT;
      const size_t b2 = BETA-b;
      const Digit* pE = pA+sT;
      Digit e,d = *--pE;
      while (pE != pA) {
        e = *--pE;
        *--pT = (e << b2) | (d >> b);
        if (pE == pA) { d = e; break; }
        d = *--pE;
        *--pT = (d << b2) | (e >> b);
      }
      *--pT = d >>= b;
      if (d == 0 && sT > 1) { p = ++pT; size = --sT; }
    } else {
      const Digit* pE = pA+sT;
      COPY(pT, pA, pA, pE);
    }
  }

  NATURALCONDITION(*this);
}

void Natural::rmove(size_t a)
// Algorithm:  b.rmove(a)
// Input:      a in size_t, b in Natural.
// Output:     b in Natural such that b := [b/2^(BETA*a)] ||
{
  NATURALCONDITION(*this);

  size_t sT = size;
  if (a >= sT) *this = 0;
  else {
    size_t sA = sT-a;
    Digit* pE = p;
    Digit* pT = pE+sT;
    if (a) {
      Digit* pA = pT-a;
      COPY_BACKWARD(pT, pA, pA, pE);
      p = pT; size = sA;
      FILL_ZERO(pA, pT);
    }
  }

  NATURALCONDITION(*this);
}

void Natural::add(const Natural& a, const Digit b)
// Algorithm:  c.add(a, b)
// Input:      a,c in Natural, b in Digit where not a.p = c.p.
// Output:     c in Natural such that c = a+b ||
{
  NATURALCONDITION(a);
  CONDITION(a.p != p);
  
  const Digit* pA = a.p;
  const size_t sA = a.size;
  if (sA == 1) {
    Digit* pT = setsize(1);
    Digit x = *pA;
    *pT = x += b;
    if (x < b) { *--pT = 1; p = pT; size = 2; }
  } else {
    Digit* pT = setsize(sA);
    const Digit* pE = pT+sA-1;
    COPY(pT, pA, pT, pE);
    Digit x = *pA;
    *pT = x += b;
    if (x < b) {
      while (++(*--pT) == 0);
      const Digit* pC = p;
      if (pT < pC) { p = pT; size = sA+1; }
    }
  }

  NATURALCONDITION(*this);
}

void Natural::sub(const Natural& a, const Digit b)
// Algorithm:  c.sub(a, b)
// Input:      a,c in Natural, b in Digit where not a.p = c.p.
// Output:     c in Natural such that c = a-b ||
{
  NATURALCONDITION(a);
  CONDITION(a.p != p);
  
  const Digit* pA = a.p;
  const size_t sA = a.size;
  if (sA == 1) {
    Digit* pT = setsize(1);
    Digit x = *pA;
    if (x < b) errmsg(3, "(sub)");
    *pT = x-b;
  } else {
    Digit* pT = setsize(sA);
    Digit* pE = pT+sA-1;
    COPY(pT, pA, pT, pE);
    Digit x = *pA;
    if (x < b) {
      Digit y;
      do y = --(*--pE); while (y == GAMMA);
      if (y == 0 && pE == p) { p = pE+1; size = sA-1; }
    }
    *pT = x-b;
  }

  NATURALCONDITION(*this);
}

void Natural::mul(const Natural& a, const Digit b)
// Algorithm:  c.mul(a, b)
// Input:      a in Natural, b in Digit where not a.p = c.p.
// Output:     c in Natural such that c = a*b ||
{
  NATURALCONDITION(a);
  CONDITION(a.p != p);

  const Digit* pA = a.p;
  const size_t sA = a.size;
  const size_t sT = sA+1;
  Digit* pT = setsize(sT);
  if (mul(pA, pA+sA, pT+sT, b) == 0) normalize();

  NATURALCONDITION(*this);
}

#if defined(_DigitAsm_) && _M_IX86 >= 300 && defined(_MSC_VER)
static void x_div(const Digit* pE, const Digit* pA, Digit* pC, const Digit b, Digit& d)
{
  __asm {
    xor edx,edx
    mov edi,pA
    mov esi,pC
    mov ebx,b
    mov ecx,pE
L1: mov eax,[edi]
    div ebx
    mov [esi],eax
    add edi,4
    add esi,4
    cmp edi,ecx
    jne L1
    mov eax,d
    mov [eax],edx
  }
}
#endif

void div(const Natural& a, const Digit b, Natural& c, Digit& d)
// Algorithm:  div(a, b, q, r)
// Input:      a in Natural, b in Digit where not b = 0.
// Output:     q in Natural, r in Digit such that q = [a/b], r = a - q*b ||
{
  NATURALCONDITION(a);

  if (b == 0) a.errmsg(4, "(div)");

  const size_t sA = a.size;
  const Digit* pA = a.p;
  Digit* pC = c.p;
  if (pA != pC) pC = c.setsize(sA);
  if (sA == 1) {
    Digit x = *pA;
    Digit y = x/b;
    *pC = y; d = x - y*b;
  } else {
    const Digit* pE = pA+sA;
#if defined(_DigitAsm_) && _M_IX86 >= 300 && defined(_MSC_VER)
    x_div(pE, pA, pC, b, d);
#elif defined(_DigitAsm_)
    Digit e = 0;
    do {
      a.digitdiv(e, *pA, b, *pC, e); 
      ++pC;
    } while (++pA != pE);
    d = e;
#else
    if (b > GAMMA/2) {
      Digit e = 0;
      do {
        a.digitdiv(e, *pA, b, *pC, e); 
        ++pC;
      } while (++pA != pE);
      d = e;
    } else {
      const Digit n = log2(b)+1;
      const Digit n2 = BETA-n;
      const Digit b2 = b << n2;
      Digit x,y = *pA++;
      Digit z = y >> n;
      do {
        x = *pA;
        a.digitdiv(z, (y << n2) | (x >> n), b2, *pC, z);
        ++pC;
        if (++pA == pE) { y = x; break; }
        y = *pA;
        a.digitdiv(z, (x << n2) | (y >> n), b2, *pC, z);
        ++pC;
      } while (++pA != pE);
      a.digitdiv(z, y << n2, b2, *pC, z);
      d = z >> n2;
    }
#endif
    c.normalize();
  }

  NATURALCONDITION(c);
}

Natural& Natural::operator+=(const Digit a)
// Algorithm:  c := c += a
// Input:      a in Digit, c in Natural.
// Output:     c in Natural such that c := c+a ||
{
  NATURALCONDITION(*this);
  
  Digit* pT = p+size-1;
  Digit b = *pT;
  *pT = b += a;
  if (b < a) inc(pT);

  NATURALCONDITION(*this);

  return *this;
}

Natural& Natural::operator-=(const Digit a)
// Algorithm:  c := c -= a
// Input:      a in Digit, c in Natural.
// Output:     c in Natural such that c := c-a ||
{
  NATURALCONDITION(*this);
  
  Digit* pT = p+size;
  Digit b = *--pT;
  if (b < a) dec(pT);
  *pT = b-a;             // no normalize important!

  NATURALCONDITION(*this);

  return *this;
}

Natural& Natural::operator*=(const Digit a)
// Algorithm:  c := c *= a
// Input:      a in Digit, c in Natural.
// Output:     c in Natural such that c := c*a ||
{
  NATURALCONDITION(*this);

  Digit* pT = p;
  size_t sT = size;
  if (sT == 1) {
    const Digit x = *pT;
    pT = setsize(2);
    digitmul(x, a, pT[0], pT[1]);
    normalize();
  } else if (sT == 2) {
    const Digit x0 = pT[0];
    const Digit x1 = pT[1];
    pT = setsize(3);
    digitmul(a, x0, x1, pT);
    normalize();
  } else {
    const Digit* rT = root;
    if (pT-rT == 1) {
      enlarge(DELTA);
      pT = p;
    }
    Digit* pE = pT+sT;
    if (mul(pT, pE, pE, a)) {
      p = --pT; size = ++sT;
    } else normalize();
  }

  NATURALCONDITION(*this);

  return *this;
}

Natural& Natural::operator>>=(size_t a)
// Algorithm:  c := c >>= a
// Input:      a in size_t, c in Natural.
// Output:     c in Natural such that c := [c/2^a] ||
{
  NATURALCONDITION(*this);

  const size_t a2 = a/BETA;
  size_t sT = size;
  if (a2 >= sT) *this = 0;
  else {
    a %= BETA;
    size_t sA = sT-a2;
    Digit* pE = p;
    Digit* pT = pE+sT;
    if (a) {
      const Digit* pA = pE+sA;
      const size_t a2 = BETA-a;
      Digit e,d = *--pA;
      while (pA != pE) {
        e = *--pA;
        *--pT = (e << a2) | (d >> a);
        if (pA == pE) { d = e; break; }
        d = *--pA;
        *--pT = (d << a2) | (e >> a);
      }
      *--pT = d >>= a;
      if (pE != pT) FILL_ZERO(pE, pT);
      if (d == 0 && sA > 1) { ++pT; --sA; }
      p = pT; size = sA;
    } else if (a2) {
      Digit* pA = pT-a2;
      COPY_BACKWARD(pT, pA, pA, pE);
      p = pT; size = sA;
      FILL_ZERO(pA, pT);
    }
  }

  NATURALCONDITION(*this);

  return *this;
}

Natural& Natural::operator<<=(size_t a)
// Algorithm:  c := c <<= a
// Input:      a in size_t, c in Natural.
// Output:     c in Natural such that c := c*2^a ||
{
  NATURALCONDITION(*this);

  if (a >> (CHAR_BIT*sizeof(size_t)-1)) errmsg(2, "(operator<<=)");
  Digit* pT = p;
  size_t sT = size;
  Digit c = *pT;
  if (c) {
    size_t b = a/BETA+1;
    a %= BETA;
    Digit* rT = root;
    size_t sC = sT+b;
    Digit* pC = pT-b;
    if (rT+b >= pT) {
      root = pC = NOTHROW_NEW Digit[sC+DELTA];
      if (!pC) errmsg(2, "(operator<<=)");
      FILL_DELTA(pC);
      pC += DELTA;
    }
    p = pC; size = sC;
    if (a) {
      const Digit a2 = BETA-a;
      const Digit* pE = pT+sT;
      Digit d = c >> a2;
      *pC++ = d;
      if (d == 0) { size = --sC; p = pC; }
      while (++pT != pE) {
        d = *pT;
        *pC++ = (c << a) | (d >> a2);
        if (++pT == pE) { c = d; break; }
        c = *pT;
        *pC++ = (d << a) | (c >> a2);
      }
      *pC++ = c << a;
    } else {
      *pC = 0;
      size = --sC; p = ++pC;
      const Digit* pE = pT+sT;
      COPY(pC, pT, pT, pE);
    }
    if (--b) {
      const Digit* pE = pC+b;
      FILL_ZERO(pC, pE);
    }
    if (root != rT) delete[] rT;
  }

  NATURALCONDITION(*this);

  return *this;
}

Digit Natural::operator=(const Digit a)
// Algorithm:  c := b = a
// Input:      a in Digit, b in Natural.
// Output:     b in Natural, c in Digit such that b = a, c = a ||
{
  Digit* pT = p;
  size_t sT = size;
  if (sT == 1) *pT = a;
  else {
    Digit* pE = pT+sT;
    *--pE = a; size = 1; p = pE;
    FILL_ZERO(pT, pE);
  }

  NATURALCONDITION(*this);

  return a;
}

#if defined(_DigitAsm_) && _M_IX86 >= 300 && defined(_MSC_VER)
static Digit x_mod(const Digit* pE, const Digit* pA, const Digit b)
{
  Digit r;
  __asm {
    xor edx,edx
    mov edi,pA
    mov ebx,b
    mov ecx,pE
L1: mov eax,[edi]
    div ebx
    add edi,4
    cmp edi,ecx
    jne L1
    mov r,edx
  }
  return r;
}
#endif

Digit operator%(const Natural& a, const Digit b)
// Algorithm:  c := a%b
// Input:      a in Natural, b in Digit where not b = 0.
// Output:     c in Digit such that c = a - [a/b]*b ||
{
  NATURALCONDITION(a);

  if (b == 0) a.errmsg(4, "(operator%)");

  const size_t sA = a.size;
  const Digit* pA = a.p;
  if (sA == 1) return *pA % b;
  else {
    const Digit* pE = pA+sA;
#if defined(_DigitAsm_) && _M_IX86 >= 300 && defined(_MSC_VER)
    return x_mod(pE, pA, b);
#elif defined(_DigitAsm_)
    Digit d,c = 0;
    do a.digitdiv(c, *pA, b, d, c); while (++pA != pE);
    return c;
#else
    if (b > GAMMA/2) {
      Digit d,c = 0;
      do a.digitdiv(c, *pA, b, d, c); while (++pA != pE);
      return c;
    } else {
      const Digit n = log2(b)+1;
      const Digit n2 = BETA-n;
      const Digit b2 = b << n2;
      Digit x,y = *pA++;
      Digit d,z = y >> n;
      do {
        x = *pA;
        a.digitdiv(z, (y << n2) | (x >> n), b2, d, z);
        if (++pA == pE) { y = x; break; }
        y = *pA;
        a.digitdiv(z, (x << n2) | (y >> n), b2, d, z);
      } while (++pA != pE);
      a.digitdiv(z, y << n2, b2, d, z);
      return z >> n2;
    }
#endif
  }
}

void Natural::muladd(const Natural& a, const Digit b)
// Algorithm:  c.muladd(a, b)
// Input:      a,c in Natural, b in Digit.
// Output:     c in Natural such that c := c + a*b ||
{
  NATURALCONDITION(*this);
  NATURALCONDITION(a);

  const size_t sA = a.size;
  const size_t sT = size;
  const Digit* rT = root;
  Digit* pT = p;
  if (rT+sA+2 >= pT+sT) {       // (a.size+2 >= rootsize())?
    enlarge(sA-sT+2+DELTA);
    pT = p;
  }
  Digit* pA = a.p;
  pT += sT;
  Digit* pC = muladd(pA, pA+sA, pT, b);
  const size_t sz = pT-pC;
  if (sz > sT) { size = sz; p = pC; }
  normalize();

  NATURALCONDITION(*this);
}

void Natural::mulsub(const Natural& a, const Digit b)
// Algorithm:  c.mulsub(a, b)
// Input:      a,c in Natural, b in Digit.
// Output:     c in Natural such that c := c - a*b ||
{
  NATURALCONDITION(*this);
  NATURALCONDITION(a);

  const size_t sT = size;
  const size_t sA = a.size;
  if (sA > sT) errmsg(3, "(mulsub)");
  const Digit* pA = a.p;
  const Digit* pE = pA+sA;
  Digit* pT = p+sT;
  Digit x,y,z = 0;

  do {
    digitmul(*--pE, b, x, y);
    y += z;
    z = (y < z) + x;
    x = *--pT;
    z += (x < y);
    *pT = x-y;
  } while (pA != pE);
  if (z) {
    const Digit* pE = p;
    if (pT == pE) errmsg(3, "(mulsub)");
    Digit x = *--pT;
    if (x < z) dec(pT);
    *pT = x -= z;
    if (x == 0 && pT == pE) {       // normalize()
      size_t sT = size;
      if (sT > 2) {
        do { ++pT; --sT; } while (*pT == 0 && sT > 1);
        p = pT; size = sT;
      } else if (sT == 2) { p = ++pT; size = --sT; }
    }
  }

  NATURALCONDITION(*this);
}

inline void sqrtsub(const Digit* pT, Digit* pDif, const Digit* pSub)
// Algorithm:  sqrtsub(r, s, t)
//             Let a,b in Natural.
// Input:      r,s in [a.root, a.p+L(a)] where r < s,
//             t in [b.root, b.p+L(b)] where t-(s-r) also in [b.root, b.p+L(b)],
//             where R(a) > L(b), [r, s] >= [t-(s-r), t].
// Output:     [r, s] := [r, s] - [t-(s-r), t] ||
{
  do {
    Digit c = *pDif;
    Digit d = *pSub;
    *pDif -= d;
    if (c < d)
      do {
        c = *--pDif;
        d = *--pSub;
        *pDif -= d+1;
      } while (c <= d);
    --pDif;
  } while (pT <= --pSub);
}

void Natural::sqrt(const Natural& a)
// Algorithm:  b.sqrt(a)
// Input:      a in Natural.
// Output:     b in Natural such that b = [sqrt(a)] ||
{
  NATURALCONDITION(a);

  const size_t sA = a.size;
  if (sA == 1) *this = ::sqrt(*a.p);
  else if (sA == 2) *this = ::sqrt(*a.p, a.p[1]);

#ifdef FASTER_BINSQRT
  else if (sA <= 2*BETA) {
    Natural c;
    ::sqrt(a, *this, c);
  } else if (sA <= NEWTON_SQRT_MARK) {
    Natural x(sA+2+DELTA, ' ');
    Digit* y[BETA];
    Digit  i;
    --x.size;
    const size_t SHIFT_H = BETA - size_t(log2(*a.p))-1;
    const size_t SHIFT   = (SHIFT_H&1)? SHIFT_H-1 : SHIFT_H;
    x = a << SHIFT;
    size_t j = x.size + 2;
    for (i = 0; i < BETA; ++i) {
      Digit* pY = y[i] = NOTHROW_NEW Digit[j];
      const Digit* pE = pY+j;
      FILL_ZERO(pY, pE);
      ++y[i];                       // because subpos!
    }
    size_t sX = x.size;
    size_t sY = sX;
    Digit* pX = x.p;
    pX[sX] = 1;

    i = Digit(1) << (BETA-2); j = 0;
    size_t k = 0;
    do {
      y[j][k] |= i;
      if (sX == sY) {
        Digit* pA = pX;
        Digit* pB = y[j];
        Digit x2,y2;
        do { x2 = *pA++; y2 = *pB++; } while (x2 == y2);
        if (x2 > y2) {
          sqrtsub(y[j], pX+k, y[j]+k);
          while (*pX == 0) { ++pX; --sX; }
          size_t j2 = j;
          size_t k2 = k;
          do {
            if (++j2 == BETA) { j2 = 0; --k2; }
            y[j2][k2] |= i; i >>= 1;
            if (i == 0) { i = Digit(1) << (BETA-1); ++k2; }
          } while (j2 != j);
        }
      } else if (sX > sY) {
        if (j) sqrtsub(y[j], pX+k, y[j]+k);
        else sqrtsub(y[j], pX+k+1, y[j]+k);
        while (*pX == 0) { ++pX; --sX; }
        size_t j2 = j;
        size_t k2 = k;
        do {
          if (++j2 == BETA) { j2 = 0; --k2; }
          y[j2][k2] |= i; i >>= 1;
          if (i == 0) { i = Digit(1) << (BETA-1); ++k2; }
        } while (j2 != j);
      }
      y[j][k] &= ~i;

      if (++j == BETA) {
        i = Digit(1) << (BETA-2);
        j = 0; --sY;
      } else if (i < 4) { i = Digit(1) << (BETA-2); ++k; }
      else i >>= 2;
    } while ((x.size&1) && (j == BETA/2-1 || k <= x.size/2)
      || (x.size&1) == 0 && (j || k < x.size/2));

    sX = x.size;
    if ((sX&1) == 0) ++sY;
    pX = setsize(sY);
    Digit* pA = (sX&1)? y[BETA/2-1] : y[BETA-1];
    const Digit* pE = pA+sY;
    COPY(pX, pA, pA, pE);
    *this >>= SHIFT/2 + 1;
    // b = *this, a*2^SHIFT = (b/2)^2 + x,
    // but a != (b/2^(SHIFT/2+1))^2 + x/2^SHIFT
    // because b is not minimal!
    for (i = 0; i < BETA; ++i) delete[] --y[i];
    x.normalize();  // just to be on the safe side
#else
  } else if (sA <= NEWTON_SQRT_MARK) {
    Natural c;
    ::sqrt(a, *this, c);
#endif
  } else newton_sqrt(a);

  NATURALCONDITION(*this);
}

void Natural::newton_sqrt(Natural a)
// Algorithm:  b.newton_sqrt(a)
// Input:      a in Natural where L(a) >= 2.
// Output:     b in Natural such that b = [sqrt(a)] ||
{
  NATURALCONDITION(a);
  CONDITION(a.size >= 2);

  const size_t sA = a.size;
  const size_t l = BETA-1-size_t(log2(a.highest()));
  a <<= l & (~size_t(0)-1);
  const Digit d = ::sqrt(a.p[0], a.p[1]);
  Digit q,r;
  if (d == Digit(~(GAMMA/2))) q = GAMMA;
  else a.digitdiv(~(GAMMA/2), 0, d, q, r);
  Natural t,b = q;
  size_t m;
  size_t* s = quad_convergence_sizes(sA/2+1, m);
  size_t sB = s[--m];
  do {
    *this = b*b;
    if (sB < sA) a.size = sB+1;
    t = *this * a; a.size = sA;
    const size_t k = s[--m];
    t >>= (t.size-k)*BETA - 1;
    *this = t * b;
    *this >>= (size-k)*BETA - 1;
    b.lmove(k-sB); t = b;
    b <<= 1; b += t; b -= *this; b >>= 1;
    sB = k;
  } while (m);
  *this = b*b; t = *this * a; t >>= (t.size-sB)*BETA - 1;
  *this = t * b; *this >>= (size-sB)*BETA - 1;
  t = b; b <<= 1; b += t; b -= *this;
  *this = b * a;
  *this >>= sA*(BETA/2) + sB*BETA + l/2;
  delete[] s;

  NATURALCONDITION(*this);
}

void sqrt(const Natural& b, Natural& c, Natural& d)
// Algorithm:  sqrt(b, c, d)
// Input:      b in Natural.
// Output:     c,d in Natural such that c = [sqrt(b)], d = b - c^2 ||
{
  NATURALCONDITION(b);
  NATURAL_FOR_CHECK(_b, b);

  const size_t sB = b.size;
  if (b == 0) { c = d = 0; return; }
  else if (sB == 1) {
    Digit c2,d2;
    sqrt(*b.p, c2, d2);
    c = c2; d = d2;
    return;
  } else if (sB == 2) {
    Digit d1,d2;
    c = sqrt(*b.p, b.p[1], d1, d2);
    if (d1) {
      Digit* pD = d.setsize(2);
      *pD = d1; pD[1] = d2;
    } else d = d2;
    return;
  }
  
  static Digit ax[BETA] = { Digit(1) << (BETA-2),
                            Digit(1) << (BETA-3) | Digit(1) << (BETA-4),
                            0 };
  if (ax[2] == 0) {
    size_t i = 2;
    do ax[i] = ax[i-1] >> 2; while (++i < BETA/2);
    ax[i] = Digit(1) << (BETA-1);
    while (++i < BETA) ax[i] = ax[i-1] >> 2;
  }
  
  Natural a(b, 1);
  size_t sT = a.size;
  Natural y((Digit)0, sT);      // a.size >= 2!
  size_t sY       = --sT;
  Digit* pY       = y.p;
  Digit* pT       = a.p;
  Digit* pC       = pY;
  Digit* pD       = pT;
  const Digit* pE = pY+sY;
  size_t s        = 0;
  pT[sT] = 1;    // because compare
  
  do {
    *pC ^= ax[s];
    
    if (sT == sY) {    // (a >= y)?
      Digit* pA = pT;
      Digit* pB = pY;
      Digit x,y;
      do { x = *pA++; y = *pB++; } while (x == y);
      if (x > y) {
        sqrtsub(pY, pD, pC);
        while (*pT == 0) { ++pT; --sT; }
        *pC |= ax[s+BETA/2];
      }
    } else if (sT > sY) {
      sqrtsub(pY, pD, pC);
      while (*pT == 0) { ++pT; --sT; }
      *pC |= ax[s+BETA/2];
    }
    
    Digit* pA = pC;
    Digit x = *pA;
    x >>= 1;
    while (pA != pY) {        // y >>= 1
      Digit y = *--pA;
      pA[1] = x | ((y&1) << (BETA-1));
      x = y >> 1;
    }
    *pA = x;
    if (x == 0) { ++pY; --sY; }
    
    if (++s == BETA/2) { s = 0; ++pC; ++pD; }
    
  } while (pC < pE);
  y.p = pY; y.size = sY;
  if (sT) { a.size = sT; a.p = pT; }
  else { a.size = 1; a.p = --pT; }
  c = y; d = a;

  CONDITION(d+c*c == _b);
  NATURALCONDITION(c);
  NATURALCONDITION(d);
}

Natural root(const Natural& a, const Digit n)
// Algorithm:  b := root(a, n)
// Input:      a in Natural, n in Digit.
// Output:     b in Natural such that b = [a^(1/n)] ||
{
  if (n == 0) return 1;
  else if (n == 1) return a;
  else if (n == 2) return sqrt(a);
  const Digit k = log2(a)+1;
  const Digit m = n-1;
  const Digit l = k/n;
  Natural b(4);
  if (l >= 2) { b = (4*(k - l*n))/n + 5; b <<= size_t(l)-2; }
  Natural c,q,r;
  while (true) {
    c = pow(b, m);
    div(a, c, q, r);
    if (b == q) return b;
    else if (b < q) {
      c = pow(++b, n);
      if (a < c) --b;
      return b;
    }
    b *= m;
    b += q;
    b /= n;
  }
}

void Natural::setbit(const Digit a)
// Algorithm:  c.setbit(a)
// Input:      a in Digit, c in Natural.
// Output:     c in Natural such that c := c or 2^a ||
{
  NATURALCONDITION(*this);

  const size_t  b = size_t(a/BETA) + 1;
  const Digit   c = Digit(1) << (a%BETA);
  const size_t sT = size;
  const Digit* rT = root;
  Digit* pT = p;
  if (b <= sT) pT[sT-b] |= c;
  else if (rT+b < pT+sT) {                   // b < rootsize()
    p = pT -= b-sT; size = b; *pT |= c;
  } else {
    enlarge(DELTA+b-sT);
    pT = p;
    p = pT -= b-sT; size = b; *pT |= c;
  }

  NATURALCONDITION(*this);
}

void Natural::clearbit(const Digit a)
// Algorithm:  c.clearbit(a)
// Input:      a in Digit, c in Natural.
// Output:     c in Natural such that c := c and not(2^a) ||
{
  NATURALCONDITION(*this);

  const size_t b = size_t(a/BETA) + 1;
  const Digit c = Digit(1) << (a%BETA);
  size_t sT = size;
  Digit* pT = p;
  if (b == sT) {
    Digit d = *pT;
    *pT = d &= ~c;
    if (d == 0) {                       // normalize
      if (sT > 2) {
        do { ++pT; --sT; } while (*pT == 0 && sT > 1);
        p = pT; size = sT;
      } else if (sT == 2) { p = ++pT; size = --sT; }
    }
  } else if (b < sT) pT[sT-b] &= ~c;

  NATURALCONDITION(*this);
}

bool Natural::testbit(const Digit a) const
// Algorithm:  c := b.testbit(a)
// Input:      a in Digit, b in Natural.
// Output:     c in bool such that if b and 2^a then c = true else c = false ||
{
  NATURALCONDITION(*this);

  const size_t b = size_t(a/BETA) + 1;
  const size_t sT = size;
  if (b > sT) return false;
  return ((p[sT-b] & (Digit(1) << (a%BETA))) != 0);
}

// Correct ?!
#ifndef RAND_MAX
# define RAND_MAX ((unsigned(1) << (CHAR_BIT*sizeof(int)-1)) - 1)
#endif

void Natural::rand(size_t n)
// Algorithm:  a.rand(n)
// Input:      n in size_t.
// Output:     a in Natural such that a < 2^n (random number) ||
{
  NATURALCONDITION(*this);

  const size_t sT = max(n/BETA + (n%BETA != 0), size_t(1));
  const size_t s  = min(1+size_t(log2(Digit(RAND_MAX))), BETA);
  Digit* pT = setsize(sT);
  pT += sT;
  size_t i = 0;
  Digit k = 0;
  while (n >= s) {
    const Digit j = (Digit)::rand();
    k |= j << i;
    const size_t l = BETA-i;
    if (s >= l) {
      *--pT = k;
      k = j >> l;
      i -= BETA;
    }
    i += s; n -= s;
  }
  if (n) {
    const Digit k2 = Digit(1) << n;
    const Digit j = Digit(::rand()) & (k2-1);
    *--pT = k | (j << i);
    const size_t l = BETA-i;
    if (s > l) *--pT = j >> l;
  } else *--pT = k;
  normalize();

  NATURALCONDITION(*this);
}

const char* Natural::atoN(const char* a, const Digit b)
// Algorithm:  c := d.atoN(a, b)
// Input:      d in Natural, a in String, b in Digit where 2 <= b <= 36.
// Output:     d in Natural, c in String such that d = a ||
//
// Note: Returns a pointer to the first occurrence of a non-digit character
//       in a.
{
  CONDITION(b >= 2 && b <= 36);

  *this = 0;
  if (b < 2 || b > 36) return a;
  const char* a2 = a;
  if (b < 10) {
    while (true) {
      const char d = *a2;
      if (d == 0 || !isdigit(d) || d-'0' >= b) break;
      ++a2;
    }
  } else {
    while (true) {
      const char d = *a2;
      if (d == 0 || !isdigit(d) && toupper(d)-'A'+10 >= b) break;
      ++a2;
    }
  }
  if ((b & (~b+1)) == b) {
    const Digit j = log2(b);
    Digit l = (a2-a)*j;
    while (a != a2) {
      const Digit d = Digit((isdigit(*a))? *a - '0' : toupper(*a)-'A'+10);
      Digit i = Digit(1) << j;
      do {
        i >>= 1; --l;
        if (d&i) setbit(l);
      } while (i != 1);
      ++a;
    }
    return a;
  }

  Digit  b2 = b;
  const Digit w = GAMMA/b;
  size_t j = 1;
  while (b2 <= w) { b2 *= b; ++j; }
  size_t i = 0;
  size_t l = (a2-a)/j;
  Digit d = 0;
  if (l < ATON_MARK) {
    while (a != a2) {
      d *= b; d += Digit((isdigit(*a))? *a - '0' : toupper(*a)-'A'+10);
      if (++i == j) {
        *this *= b2; *this += d;
        d = i = 0;
      }
      ++a;
    }
  } else {
    Digit* v = new Digit[l];
    Digit* v2 = v;
    while (a != a2) {
      d *= b; d += Digit((isdigit(*a))? *a - '0' : toupper(*a)-'A'+10);
      if (++i == j) {
        *v2 = d; ++v2;
        d = i = 0;
      }
      ++a;
    }
    
    CONDITION(l > 1);
    size_t l2 = l >> 1;
    Natural* c = new Natural[l2];
    size_t k2 = l2;
    size_t k = l-1;
    while (true) {
      if (--k2 == 0) {
        if (k == 2) { c[0] = v[k-2]; c[0] *= b2; c[0] += v[k-1]; }
        else c[0] = v[k-1];
        c[0] *= b2; c[0] += v[k];
        break;
      }
      c[k2] = v[k-1]; c[k2] *= b2; c[k2] += v[k];
      k -= 2;
    }
    delete[] v;
    Natural b3 = b2;
    while (l > 3) {
      l = l2;
      l2 >>= 1;
      Natural* c2 = new Natural[l2];
      b3 *= b3;
      k2 = l2;
      k = l-1;
      while (true) {
        if (--k2 == 0) {
          if (k == 2) { c2[0] = c[k-2]; c2[0] *= b3; c2[0] += c[k-1]; }
          else c2[0] = c[k-1];
          c2[0] *= b3; c2[0] += c[k];
          break;
        }
        c2[k2] = c[k-1]; c2[k2] *= b3; c2[k2] += c[k];
        k -= 2;
      }
      delete[] c;
      c = c2;
    }
    *this = c[0];
    delete[] c;
    
  }
  if (i) {
    b2 = 1;
    do b2 *= b; while (--i);
    *this *= b2; *this += d;
  }

  NATURALCONDITION(*this);

  return a;
}

static void Ntoa2(const Natural& a, char* c, const Digit b, size_t width, bool active)
// Algorithm:  c := Ntoa2(a, c, b)
// Input:      a in Natural, b in Digit, c in String
//             where b = 2^x for x in N, sizeof(c) > BETA*L(a)/x.
// Output:     c in String such that c = a ||
//
// Note:       conversion Natural to string.
{
  const Digit j = log2(b);
  Digit k = log2(a)+1;
  Digit k2 = (k/j)*j;
  char* str = c;
  if (active && width) {
    char* str2 = str + width-k/j - (k != k2);
    do *str = '0'; while (++str != str2);
  }
  Digit d;
  if (k != k2) {
    d = 0;
    do {
      d <<= 1; d |= Digit(a.testbit(--k));
    } while (k != k2);
    *str++ = (d <= 9)? char(d + '0') : char(d + 'A' - 10);
  }
  while (k > 0) {
    d = 0;
    k2 = k-j;
    do {
      d <<= 1; d |= Digit(a.testbit(--k));
    } while (k != k2);
    *str++ = (d <= 9)? char(d + '0') : char(d + 'A' - 10);
  }
  *str = 0;
}
static void Ntoa(Natural a, char* c, const Digit b, size_t width, bool active)
// Algorithm:  c := Ntoa(a, c, b)
// Input:      a in Natural, b in Digit, c in String
//             where 2 <= b <= 36, sizeof(c) > BETA*L(a)/log2(b).
// Output:     c in String such that c = a ||
//
// Note:       conversion Natural to string.
{
  char* str = c;
  Digit  b2 = b;
  const Digit w = GAMMA/b;
  size_t j = 1;
  while (b2 <= w) { b2 *= b; ++j; }
  const size_t s = 1 + a.length()*BETA/size_t(log2(b2));
  if (s < NTOA_MARK) {
    Digit d;
    while (true) {
      div(a, b2, a, d);
      if (a == 0) break;
      for (size_t i = 0; i < j; ++i) {
        const char c = char(d%b);
        *str++ = (c <= 9)? char(c + '0') : char(c + 'A' - 10);
        d /= b;
      }
      width -= j;
    }
    while (d > 0) {
      const char c = char(d%b);
      *str++ = (c <= 9)? char(c + '0') : char(c + 'A' - 10);
      d /= b;
      --width;
    }
    if (active && width) {
      char* str2 = str+width;
      do *str = '0'; while (++str != str2);
    }

    *str = 0;
    for (char* str2 = c; str2 < --str; ++str2) {
      const char c = *str2;
      *str2 = *str;
      *str = c;
    }
  } else {
    Natural c = pow(Natural(b2), s/2);
    Natural q,r;
    div(a, c, q, r);

    const size_t sz = j*(s/2);
    Ntoa(q, str, b, (sz < width)? width-sz : 0, active);
    Ntoa(r, str+strlen(str), b, sz, true);
  }
}
char* Ntoa(const Natural& a, char* c, const Digit b)
// Algorithm:  c := Ntoa(a, c, b)
// Input:      a in Natural, b in Digit, c in String
//             where 2 <= b <= 36, sizeof(c) > BETA*L(a)/log2(b).
// Output:     c in String such that c = a ||
//
// Note:       conversion Natural to string.
{
  CONDITION(b >= 2 && b <= 36);

  *c = 0;
  if (b < 2 || b > 36) return c;
  if ((b & (~b+1)) == b) Ntoa2(a, c, b, 0, false);
  else Ntoa(a, c, b, 0, false);
  return c;
}

static void output_stream(OSTREAM& out, const Natural& a)
{
  const size_t s = 1 + a.length()*BETA/size_t(log2(ALPHA));
  const size_t t = out.width();
  if (s < 17000) {
    Digit* pC = NOTHROW_NEW Digit[s];
    if (!pC) a.errmsg(2, "(stream operator<<)");
    Digit* pE = pC + s;
    Digit* pOut = pE;
    size_t sz = 0;
    Natural b = a;
    while (b.length() != 1) {
      div(b, ALPHA, b, *--pOut);
      sz += ALPHA_WIDTH;
    }
    out.width((sz < t)? t-sz : 0);
    out << b.highest();
    while (pOut != pE) {
      out.width(ALPHA_WIDTH);
      out.fill('0');
      out << *pOut;
      ++pOut;
    }
    out.width(0);
    delete[] pC;
    return;
  }

  Natural c = pow(Natural(ALPHA), s/2);
  Natural q,r;
  div(a, c, q, r);

  const size_t sz = ALPHA_WIDTH*(s/2);
  out.width((sz < t)? t-sz : 0);
  output_stream(out, q);
  out.width(sz);
  out.fill('0');
  output_stream(out, r);
  out.width(0);
}

OSTREAM& operator<<(OSTREAM& out, const Natural& a)
// Algorithm:  o := o << a
// Input:      o in ostream, a in Natural.
// Output:     o in ostream ||
//
// Note:       Only decimal output is supported.
{
  output_stream(out, a);
  out.fill(' ');
  return out;
}

ISTREAM& operator>>(ISTREAM& in, Natural& a)
// Algorithm:  i := i >> a
// Input:      i in istream.
// Output:     i in istream, a in Natural ||
//
// Note:       Only decimal input is supported.
{
  if (!in.good()) return in;
  Digit  d = 0;
  size_t i = 0;
  char ch  = 0;
  a = 0;
  while (!in.eof() && in.get(ch) && isdigit(ch)) {
    d *= 10; d += ch - '0';
    if (++i == ALPHA_WIDTH) {
      a *= ALPHA; a += d;
      d = i = 0;
    }
  }
  if (i) { a *= pow10(i); a += d; }
  if (in.good() && ch != '\n') in.putback(ch);
  return in;
}

OSTREAM& operator<<(OSTREAM& out, const Natural::rep& a)
// Note: puts internal representation of Natural a on an output stream.
{
  const size_t sA = a.size;
  const Digit* pA = a.p;
  const Digit* pE = pA+sA;
  out << CHAR_BIT << '*' << sizeof(Digit) << '*' << sA;
  if (a.bin) {
    out << 'B';
    do {
      const Digit d = *--pE;
      for (size_t i = 0; i < sizeof(Digit); ++i) out << char((d >> (CHAR_BIT*i))&0xff);
    } while (pA != pE);
  } else {
    out << '(' << *--pE;
    while (pA != pE) out << ',' << *--pE;
  }
  return out << ')';
}

bool Natural::scan(ISTREAM& in)
// Algorithm:  b := a.scan(i)
// Input:      a in Natural, i in istream.
// Output:     a in Natural, i in istream, b in bool ||
//
// Note:       gets Natural a as an internal representation from input stream
//             if b is true.
{
  if (!in.good()) return false;
  char c = 0;
  size_t l = 0;
  if (!(in >> l) || l != CHAR_BIT || !in.get(c) || c != '*') return false;
  l = 0;
  if (!(in >> l) || l < 1 || (l%sizeof(Digit) && sizeof(Digit)%l)
      || !in.get(c) || c != '*')  return false;
  size_t sT = 0;
  if (!(in >> sT) || sT < 1) return false;
  sT *= l;
  if (sT % sizeof(Digit)) sT += sizeof(Digit);
  sT /= sizeof(Digit);
  Digit* pT = setsize(sT);
  Digit* pE = pT+sT;
  if (in.get(c)) {
    if (c == 'B') {
      if (l <= sizeof(Digit)) {
        do {
          Digit d = 0;
          for (size_t i = 0; i < sizeof(Digit); ++i) {
            if (!in.get(c)) return false;
            d |= Digit((unsigned char)c) << (CHAR_BIT*i);
          }
          *--pE = d;
        } while (--sT);
        if (!in.get(c) || c != ')') ++pE;
      } else {
        Natural x,y;
        do {
          size_t i;
          x = 0;
          for (i = 0; i < sizeof(Digit); ++i) {
            if (!in.get(c)) return false;
            y = Digit((unsigned char)c); y <<= CHAR_BIT*i;
            x |= y;
          }
          for (i = 0; i < l; i += sizeof(Digit)) {
            *--pE = x & GAMMA; x >>= BETA;
          }
        } while (--sT);
        if (!in.get(c) || c != ')') ++pE;
      }
    } else if (c == '(') {
      if (l == sizeof(Digit)) {
        do
          if (!(in >> *--pE)) break;
        while (in.get(c) && c == ',' && --sT);
        if (c != ')') ++pE;
      } else if (l < sizeof(Digit)) {
        do {
          Digit x = 0;
          for (size_t i = 0; i < sizeof(Digit); i += l) {
            Digit y;
            in >> y;
            x |= y << (CHAR_BIT*i);
            if (!in.get(c) || c != ',') break;
          }
          *--pE = x;
        } while (c == ',' && --sT);
      } else {
        Natural x;
        do {
          in >> x;
          for (size_t i = 0; i < l; i += sizeof(Digit)) {
            *--pE = x & GAMMA; x >>= BETA;
          }
        } while (in.get(c) && c == ',' && --sT);
      }
    }
  }
  normalize();

  NATURALCONDITION(*this);

  return (pT == pE && in.good() && c == ')');
}

Natural abs(const Natural& a, const Natural& b)
// Algorithm:  c := abs(a, b)
// Input:      a,b in Natural.
// Output:     c in Natural such that c = |a-b| ||
{
  NATURALCONDITION(a);
  NATURALCONDITION(b);

  const size_t sA = a.size;
  const size_t sB = b.size;
  Natural c(max(sA, sB)+DELTA, ' ');

  if (sA > sB) c = a-b;
  else if (sA < sB) c = b-a;
  else {
    const Digit* pA = a.p;
    const Digit* pB = b.p;
    Digit* pC = c.p;
    if (pA != pC && pB != pC) pC = c.setsize(sA);
    c.abs(pC, pA, pB, sA);
    c.normalize();
  }

  NATURALCONDITION(c);

  return c;
}
int abs(const Natural& a, const Natural& b, Natural& c)
// Algorithm:  d := abs(a, b, c)
// Input:      a,b in Natural.
// Output:     d in int, c in Natural such that c = |a-b|, d*c = a-b ||
{
  NATURALCONDITION(a);
  NATURALCONDITION(b);

  const size_t sA = a.size;
  const size_t sB = b.size;

  if (sA > sB) { c = a-b; return 1; }
  else if (sA < sB) { c = b-a; return -1; }
  const Digit* pA = a.p;
  const Digit* pB = b.p;
  Digit* pC = c.p;
  if (pA != pC && pB != pC) pC = c.setsize(sA);
  const int i = c.abs(pC, pA, pB, sA);
  c.normalize();

  NATURALCONDITION(c);

  return i;
}

Natural gcd(Natural a, Natural b)
// Algorithm:  c := gcd(a, b)
// Input:      a,b in Natural.
// Output:     c in Natural such that c = max{x in Natural : x|a, x|b} ||
{
  size_t sA = a.size;
  size_t sB = b.size;
  // small sizes
  if (sA == 1) {
    Digit x = *a.p;
    if (sB == 2) {
      Digit y;
      a.gcd(b.p[0], b.p[1], 0, x, x, y);
      if (x == 0) return y;
      Natural t(x, size_t(2));
      t.p[1] = y;
      return t;
    }
    if (x == 0) return b;
    if (sB > 2) b %= x;
    return gcd(x, *b.p);
  } else if (sB == 1) {
    Digit x = *b.p;
    if (sA == 2) {
      Digit y;
      a.gcd(a.p[0], a.p[1], 0, x, x, y);
      if (x == 0) return y;
      Natural t(x, size_t(2));
      t.p[1] = y;
      return t;
    }
    if (x == 0) return a;
    a %= x;
    return gcd(x, *a.p);
  } else if (sA == 2 && sB == 2) {
    Digit x,y;
    a.gcd(a.p[0], a.p[1], b.p[0], b.p[1], x, y);
    if (x == 0) return y;
    Natural t(x, size_t(2));
    t.p[1] = y;
    return t;
  }
  if (a == b) return a;
  // extract factor 2
  size_t i = 0;
  Digit* pA = a.p+sA;
  if (*--pA == 0) {
    if (sA == 1) return b;
    do i += BETA; while (*--pA == 0);
  }
  size_t j = 0;
  Digit* pB = b.p+sB;
  if (*--pB == 0) {
    if (sB == 1) return a;
    do j += BETA; while (*--pB == 0);
  }
  // count trailing zeros
  Digit x = *pA;
  i += log2(x&-x);
  if (i) a >>= i;
  x = *pB; j += log2(x&-x);
  if (j) b >>= j;
  if (i > j) i = j;
  sA = a.size; sB = b.size;
  // small sizes
  if (sA == 1) {
    x = *a.p;
    if (sB == 2) {
      Digit y;
      a.gcd(b.p[0], b.p[1], 0, x, x, y);
      if (x == 0) {
        if (log2(y)+i < BETA) return y << i;
        Natural t = y;
        return t <<= i;
      }
      Natural t(x, size_t(2));
      t.p[1] = y;
      return t <<= i;
    } else if (sB > 2) b %= x;
    x = gcd(x, *b.p);
    if (log2(x)+i < BETA) return x << i;
    Natural t = x;
    return t <<= i;
  } else if (sB == 1) {
    x = *b.p;
    if (sA == 2) {
      Digit y;
      a.gcd(a.p[0], a.p[1], 0, x, x, y);
      if (x == 0) {
        if (log2(y)+i < BETA) return y << i;
        Natural t = y;
        return t <<= i;
      }
      Natural t(x, size_t(2));
      t.p[1] = y;
      return t;
    }
    a %= x;
    x = gcd(x, *a.p);
    if (log2(x)+i < BETA) return x << i;
    Natural t = x;
    return t <<= i;
  } else if (sA == 2 && sB == 2) {
    Digit y;
    a.gcd(a.p[0], a.p[1], b.p[0], b.p[1], x, y);
    if (x == 0) {
      if (log2(y)+i < BETA) return y << i;
      Natural t = y;
      return t <<= i;
    }
    Natural t(x, size_t(2));
    t.p[1] = y;
    return t <<= i;
  }
  
  // k-ary gcd algorithm
  Natural t1,t2;
  const size_t K_ARY_BETA = (BETA > 32)? 16 : BETA/2;
  const Digit K_ARY = Digit(1) << K_ARY_BETA;
  const Digit K_BASE = K_ARY-1;
  // first trial division phase
  Natural n = 1;
  pA = NumberBase::k_ary_gcd_trial_division;
  sA = NumberBase::k_ary_gcd_trial_division_size;
  for (j = 0; j < sA; ++j) {
    x = pA[j];
    while (true) {
      Digit y1,y2;
      div(a, x, t1, y1);
      if (y1 > 0) y1 = gcd(y1, x);
      if (y1 == 1) break;
      div(b, x, t2, y2);
      if (y2 > 0) y2 = gcd(y2, x);
      if (y2 == 1) break;
      if (y1 == y2) {
        if (y1 == 0) { n *= x; swap(a, t1); swap(b, t2); }
        else { x = y1; n *= x; a /= x; b /= x; }
      } else {
        x = gcd(y1, y2);
        if (x == 1) break;
        n *= x; a /= x; b /= x;
      }
    }
  }

  // Euclidean algorithm for numbers with a large difference
  if (a < b) swap(a, b);
  while (b != 0) {
    sA = a.size;
    sB = b.size;
    if (sA <= sB+1) break;
    if (sA == 2) {
      Digit x,y;
      if (sB == 2) a.gcd(a.p[0], a.p[1], b.p[0], b.p[1], x, y);
      else a.gcd(a.p[0], a.p[1], 0, b.p[0], x, y);
      b = 0;
      if (x == 0) a = y;
      else { a = x; a.lmove(1); a |= y; }
      break;
    } else if (sA == 1) { a = gcd(a.p[0], b.p[0]); b = 0; break; }
    div(a, b, t1, a);
    swap(a, b);
  }
  // main loop
  while (a != 0 && b != 0) {
    sA = a.size;
    sB = b.size;
    x = a&K_BASE;
    if (x == 0) a >>= K_ARY_BETA;
    else if ((x&1) == 0) a >>= log2(x&-x);
    else {
      Digit y = b&K_BASE;
      if (y == 0) b >>= K_ARY_BETA;
      else if ((y&1) == 0) b >>= log2(y&-y);
      else {
        y = (x*NumberBase::k_ary_gcd_inverse[y >> 1])&K_BASE;
        x = NumberBase::k_ary_gcd_linear_id[y];
        if (NumberBase::k_ary_gcd_linear_id_neg[y]) {
          y = (x*y)&K_BASE;
          if (a >= b) {
            a *= x; t1 = b*y; a = abs(a, t1) >> K_ARY_BETA;
          } else {
            b *= y; t1 = a*x; b = abs(b, t1) >> K_ARY_BETA;
          }
        } else {
          y = K_ARY - ((x*y)&K_BASE);
          if (a >= b) { a *= x; a += y*b; a >>= K_ARY_BETA; }
          else { b *= y; b += x*a; b >>= K_ARY_BETA; }
        }
      }
    }
  }
  // trial division phase
  if (a == 0) swap(a, b);
  sA = NumberBase::k_ary_gcd_trial_division_size;
  for (j = 0; j < sA; ++j) {
    x = pA[j];
    while (true) {
      const Digit y = gcd(a%x, x);
      if (y == 1) break;
      a /= y;
    }
  }
  if (n > 1) a *= n;
  return a << i;
}

Natural lcm(const Natural& a, const Natural& b)
// Algorithm:  c := lcm(a, b)
// Input:      a,b in Natural.
// Output:     c in Natural such that c = min{x in Natural : a|x, b|x} ||
{
  if (a == 0) return b;
  if (b == 0) return a;

  Natural q,r;
  div(a, gcd(a, b), q, r);
  return q*b;
}

Natural pow(const Natural& a, Digit b)
// Algorithm:  c := pow(a, b)
// Input:      a in Natural, b in Digit.
// Output:     c in Natural such that c = a^b ||
{
  if (b == 0) return 1;
  Digit c = Digit(1) << log2(b);

  Natural d(a);
  while (b != c) {
    b &= ~c; c >>= 1;
    d *= d;
    if (b & c) d *= a;
  }
  while (b >>= 1) d *= d;
  return d;
}

Natural pow(const Natural& a, Natural b)
// Algorithm:  c := pow(a, b)
// Input:      a,b in Natural.
// Output:     c in Natural such that c = a^b ||
{
  if (b == 1) return a;
  else if (b > 1) {
    Natural c = 1;
    Natural d = a;
    do {
      if (b.odd()) c *= d;
      b >>= 1;
      d *= d;
    } while (b > 1);
    return c *= d;
  } else return 1;
}



///////////////////////// FFT class ////////////////////////////

const Digit*  FFT::moduli     = FFT::init_moduli();
const Digit*  FFT::primroots  = FFT::init_primroots();
Digit         FFT::m          = 0;
size_t        FFT::n1         = 0;
size_t        FFT::n2         = 0;
Digit*        FFT::omega      = 0;
Digit*        FFT::omega2     = 0;
size_t*       FFT::order      = 0;

const Digit* FFT::init_moduli()
// Algorithm:  m := f.init_moduli()
// Input:      f in FFT where BETA in {32, 64}.
// Output:     m in Digit^3 such that m[0],m[1],m[2] prim
//             and m[0] < m[1] < m[2] ||
{
  if (BETA == 32) {
#ifdef _DigitAsm_
    static Digit M[3] = { Digit(2717908993U), Digit(2868903937U), Digit(3221225473U) };
#else
    static Digit M[3] = { (Digit(1) << (BETA-1)) - (Digit(1) << (BETA-4))
                                                 - (Digit(1) << (BETA-6)) + 1,
                          (Digit(1) << (BETA-1)) - (Digit(1) << (BETA-5)) + 1,
                          (Digit(1) << (BETA-1)) - (Digit(1) << (BETA-7)) + 1 };
#endif
    return M;
  } else if (BETA == 64) {
    static Digit M[3] = { Digit(1 - (Digit(1) << (BETA-24))),
                          Digit(1 - (Digit(1) << (BETA-30))),
                          Digit(1 - (Digit(1) << (BETA-32))) };
    return M;
  } else {
    static Digit M[3] = { 0, 0, 0 };
    return M;
  }
}

const Digit* FFT::init_primroots()
// Algorithm:  p := f.init_primroots()
// Input:      f in FFT where BETA in {32, 64}.
// Output:     p in Digit^3 such that p[j]^(f.n1*f.n2) = 1 (mod f.moduli[j])
//             and not p[j]^k = 1 (mod f.moduli[j]) for 0 < k < f.n1*f.n2
//             and 0 <= j <= 2 ||
{
  if (BETA == 32) {
#ifdef _DigitAsm_
    static Digit p[3] = { 5, 35, 5 };
#else
    static Digit p[3] = { 13, 31, 5 };
#endif
    return p;
  } else if (BETA == 64) {
    static Digit p[3] = { 19, 10, 7 };
    return p;
  } else {
    static Digit p[3] = { 0, 0, 0 };
    return p;
  }
}

size_t FFT::max_size()
// Algorithm:  d := f.max_size()
// Input:      f in FFT.
// Output:     d in size_t such that  k*d+1 = f.moduli[0] for k in N, k > 0 ||
{
  CONDITION(moduli[0] < moduli[1] && moduli[1] < moduli[2]);

  static size_t d = 0;
  if (d) return d;

  const Digit a = log2(moduli[0]);
  const Digit b = log2(-(moduli[2] << (BETA-1-a))) + a+2-BETA;
  const size_t c = 3*(1 << (b-log2(BETA)-1));
  return d = size_t(c*a);
}

Digit FFT::pow(Digit a, Digit b, const Digit m) const
// Algorithm:  c := f.pow(a, b, m)
// Input:      f in FFT, a,b,m in Digit.
// Output:     c in Digit such that c = a^b (mod m) ||
{
  if (b == 1) return a;
  else if (b > 1) {
    Digit c = 1;
    do {
      if (b&1) digitmulmod(a, c, m, c);
      b >>= 1;
      digitmulmod(a, a, m, a);
    } while (b > 1);
    digitmulmod(a, c, m, c);
    return c;
  }
  return 1;
}

Digit FFT::digitinv(Digit a, const Digit m) const
// Algorithm:  c := f.digitinv(a, m)
// Input:      f in FFT, a,m in Digit.
// Output:     c in Digit such that c = a^{-1} (mod m) ||
{
  Digit b = m;
  Digit q = 1;
  Digit p = 0;

  while (true) {
    Digit s = a/b;
    a -= s*b;
    if (a == 0) return m - p;
    q += s*p;
    s = b/a;
    b -= s*a;
    if (b == 0) return q;
    p += s*q;
  }
}

void FFT::init_omega(const Digit c)
// Algorithm:  b.init_omega(a)
// Input:      b in FFT, a in Digit ||
//
// Note:       initialize all necessary powers of c in b.omega and b.omega2.
{
  size_t sT = t[0].size;
  if (sT%3 == 0) sT /= 3;
  const size_t n = n2;
  const Digit  M = m;
  const Digit  s = pow(c, sT/n, M);
  if (n == n1) {
    Digit* a = omega;
    const Digit* e = a+n-1;
    Digit x = s;
    *a = 1; *++a = s;
    do {
      digitmulmod(x, s, M, x);
      *++a = x;
    } while (a != e);
  } else {
    Digit* a = omega;
    Digit* b = omega2;
    const Digit* e = b+n-1;
    Digit x = s;
    *a = *b = 1; *++b = s;
    do {
      digitmulmod(x, s, M, x);
      *++a = *++b = x;
      digitmulmod(x, s, M, x);
      *++b = x;
    } while (b != e);
  }
}

static void transpose_sqr(Digit* a, size_t n, size_t m)
{
  CONDITION(n <= m);

  size_t nn = m-n;
  const Digit* c = a+n;
  const Digit* d = a+n*m-nn;
  do {
    Digit* b = a+m;
    while (++a != c) { swap(*a, *b); b += m; }
    a += ++nn; c += m;
  } while (c != d);
}

static void transpose(Digit* a, size_t n, size_t m)
{
  if (n == m) transpose_sqr(a, n, n);
  else if (m == 2*n) {
    transpose_sqr(a, n, m);
    transpose_sqr(a+n, n, m);
    CONDITION(n >= 2);

    Digit* s = new Digit[n];
    bool*  r = new bool[m];
    bool*  t = r+m;
    FILL_ZERO(r, t);
    r -= m;

    --m;
    size_t j = 1;
    do {
      size_t i = j;
      Digit* p = a+n*i;
      Digit* t = p+n;
      COPY(s, p, p, t);
      while (true) {
        i = (i < n) ? 2*i : 2*(i-n)+1;
        if (i == j) break;
        r[i] = true;
        t = p-n;
        Digit* q = a+n*i;
        COPY(t, q, t, p);
        i = (i < n) ? 2*i : 2*(i-n)+1;
        if (i == j) { p = q; break; }
        r[i] = true;
        t = q-n;
        p = a+n*i;
        COPY(t, p, t, q);
      }
      t = p-n; s -= n;
      COPY(t, s, t, p);
      s -= n;
      while (r[++j]);
    } while (j < m);

    delete[] r;
    delete[] s;
  } else {
    CONDITION(n == 2*m && m >= 2);

    Digit* s = new Digit[m];
    bool*  r = new bool[n];
    bool*  t = r+n;
    FILL_ZERO(r, t);
    r -= n;

    --n;
    size_t j = 1;
    do {
      size_t i = j;
      Digit* p = a+m*i;
      Digit* t = p+m;
      COPY(s, p, p, t);
      while (true) {
        i = (i&1)? i/2+m : i/2;
        if (i == j) break;
        r[i] = true;
        t = p-m;
        Digit* q = a+m*i;
        COPY(t, q, t, p);
        i = (i&1)? i/2+m : i/2;
        if (i == j) { p = q; break; }
        r[i] = true;
        t = q-m;
        p = a+m*i;
        COPY(t, p, t, q);
      }
      t = p-m; s -= m;
      COPY(t, s, t, p);
      s -= m;
      while (r[++j]);
    } while (j < n);

    delete[] r;
    delete[] s;

    ++n;
    transpose_sqr(a, m, n);
    transpose_sqr(a+m, m, n);
  }
}

static void brevorder(const size_t a, Digit* b, const size_t* c)
// Algorithm:  brevorder(a, b, c)
// Input:      a in size_t, b,c in Digit^a where a > 2.
// Output:     b = [b_0, b_a[ in Digit^a such that
//             t := b_1, b_1 := b_(a/2), b_(a/2) := t,
//             t := b_i, b_i := b_j, b_j := t  if i < j
//             with j = c_i and 2 <= i < a ||
{
  CONDITION(a > 2);

  Digit i = 2;
  swap(b[1], b[a/2]);
  do {
    const Digit j = c[i];
    if (i < j) swap(b[i], b[j]);
  } while (++i < a);
}

void FFT::innerfft(const Digit* pE, Digit* pA, Digit* pB,
                   const size_t i, const Digit M) const
{
  do {
    Digit a = *pA;
    Digit b = *pB;
    Digit c = a+b;
    if (M-a <= b) c -= M;
    if (a < b) a += M;
    *pA = c; *pB = a-b;
    pA += i; pB += i;
  } while (pA < pE);
}

void FFT::innerfft(const Digit* pE, Digit* pA, Digit* pB,
                   const size_t i, const Digit M, const Digit o) const
{
  do {
    Digit a = *pA;
    Digit b = *pB;
    Digit c = a+b;
    if (M-a <= b) c -= M;
    if (a < b) a += M;
    a -= b;
    digitmulmod(o, a, M, a);
    *pA = c; *pB = a;
    pA += i; pB += i;
  } while (pA < pE);
}

void FFT::innerfftinv(const Digit* pE, Digit* pA, Digit* pB,
                     const size_t i, const Digit M, const Digit o) const
{
  do {
    Digit a = *pA;
    Digit b = *pB;
    digitmulmod(o, b, M, b);
    Digit c = a+b;
    if (M-a <= b) c -= M;
    if (a < b) a += M;
    *pA = c; *pB = a-b;
    pA += i; pB += i;
  } while (pA < pE);
}

void FFT::innerfft3(Digit* pT, const size_t sT, const Digit w) const
{
  const Digit M = m;
  Digit ww;
  digitmulmod(w, w, M, ww);
  Digit w1 = digitinv(2, M);
  Digit w2 = pow(w, sT, M);
  if (M-w2 <= w1) w2 -= M;
  w2 += w1;
  digitmulmod(3, w1, M, w1);
  w1 = M - w1;            // w1 != 0

  Digit o  = 1;
  Digit oo = 1;
  const Digit* pE = pT+sT;
  do {
    Digit x = pT[0];
    Digit y = pT[sT];
    Digit z = pT[2*sT];

    Digit s = y+z;
    if (M-y <= z) s -= M;
    if (y < z) z -= M;
    z = y - z;
    if (M-x <= s) x -= M;
    x += s;
    digitmulmod(s, w1, M, s);
    digitmulmod(z, w2, M, z);
    if (M-s <= x) s -= M;
    s += x;
    y = s + z;
    if (M-s <= z) y -= M;
    if (s < z) z -= M;
    z = s - z;
    digitmulmod(y, o, M, y);
    digitmulmod(z, oo, M, z);

    pT[0] = x;
    pT[sT] = y;
    pT[2*sT] = z;

    digitmulmod(o, w, M, o);
    digitmulmod(oo, ww, M, oo);
  } while (++pT != pE);
  //static Digit* x = pT-sT;       // to avoid a Microsoft Visual C++ 6.0 "/O2" bug
}

void FFT::innerfftinv3(Digit* pT, const size_t sT, const Digit w) const
{
  const Digit M = m;
  Digit ww;
  digitmulmod(w, w, M, ww);
  Digit w1 = digitinv(2, M);
  Digit w2 = pow(w, sT, M);
  if (M-w2 <= w1) w2 -= M;
  w2 += w1;
  digitmulmod(3, w1, M, w1);
  w1 = M - w1;            // w1 != 0

  Digit o  = 1;
  Digit oo = 1;
  const Digit* pE = pT+sT;
  do {
    Digit x = pT[0];
    Digit y = pT[sT];
    Digit z = pT[2*sT];

    digitmulmod(y, o, M, y);
    digitmulmod(z, oo, M, z);
    Digit s = y+z;
    if (M-y <= z) s -= M;
    if (y < z) z -= M;
    z = y - z;
    if (M-x <= s) x -= M;
    x += s;
    digitmulmod(s, w1, M, s);
    digitmulmod(z, w2, M, z);
    if (M-s <= x) s -= M;
    s += x;
    y = s + z;
    if (M-s <= z) y -= M;
    if (s < z) z -= M;
    z = s - z;

    pT[0] = x;
    pT[sT] = y;
    pT[2*sT] = z;

    digitmulmod(o, w, M, o);
    digitmulmod(oo, ww, M, oo);
  } while (++pT != pE);
  //static Digit* x = pT-sT;       // to avoid a Microsoft Visual C++ 6.0 "/O2" bug
}

#if defined(_DigitAsm_) && _M_IX86 >= 300 && defined(_MSC_VER)
void FFT::fft(Digit* a, const Digit* b, const size_t c) const
{
  Digit o;
  Digit* s;
  Digit* pE;
  Digit M = m;
  __asm {
    mov ebx,c   // j
    mov eax,1   // i
    mov ecx,ebx
    shr ebx,1
    cmp ecx,1
    je  L10
L1: mov edi,a
    mov ecx,c
    mov edx,M
    lea ecx,[edi+ecx*4]
    push eax
    push ebp
    shl ebx,1
    mov ebp,ecx

    // innerfft(a+c, a, a+j, 2*j, m);
L2: mov eax,[edi]         // a
    mov esi,[edi+ebx*2]   // b
    mov ecx,edx
    sub ecx,eax
    cmp ecx,esi
    lea ecx,[eax+esi]
    ja  L3
    sub ecx,edx
L3: sub eax,esi
    jae L4
    add eax,edx
L4: mov [edi],ecx
    mov [edi+ebx*2],eax
    lea edi,[edi+ebx*4]
    cmp edi,ebp
    jb  L2

    shr ebx,1
    pop ebp
    pop eax
    mov edx,b
    lea edx,[edx+eax*4]
    mov s,edx
    mov ecx,1   // k
    cmp ebx,1
    jbe L9
L5: mov esi,a
    mov edx,c
    lea edi,[esi+ecx*4]
    lea edx,[esi+edx*4]
    mov pE,edx
    push eax
    push ecx
    mov edx,s
    mov ecx,[edx]
    mov o,ecx
    shl ebx,1

    // innerfft(a+c, a+k, a+k+j, 2*j, m, *s);
L6: mov eax,[edi]
    mov esi,[edi+ebx*2]
    mov edx,M
    sub edx,eax
    cmp edx,esi
    mov edx,M
    lea ecx,[eax+esi]
    ja  L7
    sub ecx,edx
L7: sub eax,esi
    jae L8
    add eax,edx
L8: mul o
    div M
    mov [edi],ecx
    mov [edi+ebx*2],edx
    mov ecx,pE
    lea edi,[edi+ebx*4]
    cmp edi,ecx
    jb  L6

    shr ebx,1
    pop ecx
    pop eax
    mov edx,s
    lea edx,[edx+eax*4]
    mov s,edx
    inc ecx
    cmp ecx,ebx
    jne L5
L9: shr ebx,1
    shl eax,1
    mov edx,c
    cmp eax,edx
    jne L1
L10:
  }
}

void FFT::fftinv(Digit* a, const Digit* b, const size_t c) const
{
  Digit o;
  Digit* s;
  Digit* pE;
  Digit M = m;
  __asm {
    mov ebx,c   // j
    mov eax,1   // i
    mov ecx,ebx
    shr ebx,1
    cmp ecx,1
    je  L10
L1: mov edi,a
    mov ecx,c
    mov edx,M
    lea ecx,[edi+ecx*4]
    push ebx
    push ebp
    mov ebp,ecx
    shl eax,1

    // innerfft(a+c, a, a+i, 2*i, m);
L2: mov esi,[edi]
    mov ebx,[edi+2*eax]
    mov ecx,edx
    sub ecx,ebx
    cmp ecx,esi
    lea ecx,[ebx+esi]
    ja  L3
    sub ecx,edx
L3: sub esi,ebx
    jae L4
    add esi,edx
L4: mov [edi],ecx
    mov [edi+2*eax],esi
    lea edi,[edi+eax*4]
    cmp edi,ebp
    jb  L2

    shr eax,1
    pop ebp
    pop ebx
    mov edx,b
    lea edx,[edx+ebx*4]
    mov s,edx
    mov ecx,1   // k
    cmp eax,1
    jbe L9
L5: mov esi,a
    mov edx,c
    lea edi,[esi+ecx*4]
    lea edx,[esi+edx*4]
    mov pE,edx
    push ebx
    push ecx
    mov edx,s
    mov ecx,[edx]
    mov o,ecx
    lea esi,[eax*2]

    // innerfftinv(a+c, a+k, a+k+i, 2*i, m, *s);
    mov ecx,M
L6: mov eax,[edi+esi*2]
    mov ebx,[edi]
    mul o
    div ecx
    mov eax,ecx
    sub eax,ebx
    cmp eax,edx
    lea eax,[ebx+edx]
    ja  L7
    sub eax,ecx
L7: sub ebx,edx
    jae L8
    add ebx,ecx
L8: mov [edi],eax
    mov [edi+esi*2],ebx
    mov eax,pE
    lea edi,[edi+esi*4]
    cmp edi,eax
    jb  L6

    pop ecx
    pop ebx
    mov eax,esi
    shr eax,1
    mov edx,s
    lea edx,[edx+ebx*4]
    mov s,edx
    inc ecx
    cmp ecx,eax
    jne L5
L9: shr ebx,1
    shl eax,1
    mov edx,c
    cmp eax,edx
    jne L1
L10:
  }
}

void FFT::multiply_matrix(Digit* pA, const Digit w, const Digit iw) const
{
  size_t n1_ = n1;
  size_t n2_ = n2;
  Digit M = m;
  Digit u;
  __asm {
    mov u,1
    mov esi,pA
    mov eax,n1_
    mov ecx,n2_
    imul eax,ecx
    lea edi,[esi+eax*4]
L1: mov edx,n2_
    lea ecx,[esi+edx*4]
    mov ebx,iw
L2: mov eax,[esi]
    mul ebx
    div M
    mov [esi],edx
    mov eax,ebx
    mul u
    div M
    mov ebx,edx
    add esi,4
    cmp esi,ecx
    jne L2
    mov eax,u
    mul w
    div M
    mov u,edx
    cmp esi,edi
    jne L1
  }
}
#else
void FFT::fft(Digit* a, const Digit* b, const size_t c) const
{
  CONDITION((c & (~c+1)) == c);

  for (size_t i = 1, j = c/2; i != c; i *= 2, j /= 2) {
    innerfft(a+c, a, a+j, 2*j, m);

    const Digit* s = b+i;
    for (size_t k = 1; k < j; ++k) {
      innerfft(a+c, a+k, a+k+j, 2*j, m, *s);
      s += i;
    }
  }
}

void FFT::fftinv(Digit* a, const Digit* b, const size_t c) const
{
  CONDITION((c & (~c+1)) == c);

  for (size_t i = 1, j = c/2; i != c; i *= 2, j /= 2) {
    innerfft(a+c, a, a+i, 2*i, m);

    const Digit* s = b+j;
    for (size_t k = 1; k < i; ++k) {
      innerfftinv(a+c, a+k, a+k+i, 2*i, m, *s);
      s += j;
    }
  }
}

void FFT::multiply_matrix(Digit* pA, const Digit w, const Digit iw) const
// Algorithm:  c.multiply_matrix(a, w, v)
// Input:      c in FFT, v,w in Digit,
//             a = ([a_0, a_c.n2[,...,[a_((c.n1-1)*c.n2), a_(c.n1*c.n2)[)
//             in Digit^(c.n1 x c.n2).
// Output:     a in Digit^(c.n1*c.n2)
//             such that a := (a_ij * v * w^(i*j))_(0 <= i,j < c.n1*c.n2) ||
{
  const Digit M = m;
  Digit u = 1;
  const size_t k = n2;
  const Digit* pE = pA+n1*k;
  do {
    const Digit* pF = pA+k;
    Digit v = iw;
    do {
      digitmulmod(*pA, v, M, *pA);
      digitmulmod(v, u, M, v);
    } while (++pA != pF);
    digitmulmod(u, w, M, u);
  } while (pA != pE);
}

#endif
void FFT::five_step(const Digit* pE, Digit* pA, const Digit w) const
// Algorithm:  c.five_step(e, a, w)
// Input:      c in FFT, w in Digit, a in Digit^(c.n1*c.n2),
//             e in [a_0, a_(c.n1*c.n2)] where e = a_(c.n1*c.n2).
// Output:     a in Digit^(c.n1*c.n2) such that a := F_(c.n1*c.n2)(a) ||
{
  const size_t nn = n1*n2;
  do {
    transpose(pA, n1, n2);          // first step

    const Digit* pF = pA+nn;        // second step
    do {
      fft(pA, omega, n1);
      brevorder(n1, pA, order);
      pA += n1;
    } while (pA != pF);
    pA -= nn;

    transpose(pA, n2, n1);          // third step

    multiply_matrix(pA, w, 1);      // forth step

    pF = pA+nn;                     // fifth step
    do {
      fft(pA, omega2, n2);
      pA += n2;
    } while (pA != pF);
  } while (pA != pE);
}

void FFT::five_step(const Digit* pE, Digit* pA, const Digit w, const Digit iw) const
// Algorithm:  c.five_step(e, a, w, v)
// Input:      c in FFT, v,w in Digit, a in Digit^(c.n1*c.n2),
//             e in [a_0, a_(c.n1*c.n2)] where e = a_(c.n1*c.n2).
// Output:     a in Digit^(c.n1*c.n2) such that a := F_(c.n1*c.n2)^(-1)(a) ||
{
  const size_t nn = n1*n2;
  do {
    const Digit* pF = pA+nn;        // second step
    do {
      fftinv(pA, omega2, n2);
      pA += n2;
    } while (pA != pF);
    pA -= nn;

    multiply_matrix(pA, w, iw);     // third step

    transpose(pA, n1, n2);          // forth step

    pF = pA+nn;                     // fifth step
    do {
      brevorder(n1, pA, order);
      fftinv(pA, omega, n1);
      pA += n1;
    } while (pA != pF);
    pA -= nn;

    transpose(pA, n2, n1);          // sixth step
    pA += nn;
  } while (pA != pE);
}

void FFT::chinese_remainder()
// Algorithm:  f.chinese_remainder()
// Input:      f in FFT where f.t[0] = f.t[0] mod f.moduli[0],
//                            f.t[1] = f.t[1] mod f.moduli[1],
//                            f.t[2] = f.t[2] mod f.moduli[2].
// Output:     f in FFT such that f.t[2] := (f.moduli[0]*f.moduli[1]*i*f.t[2]
//                                         + f.moduli[0]*f.moduli[2]*j*f.t[1]
//                                         + f.moduli[1]*f.moduli[2]*k*f.t[0])
//                                         mod (f.moduli[0]*f.moduli[1]*f.moduli[2])
//             where i := (f.moduli[0]*f.moduli[1])^(-1) mod f.moduli[2],
//                   j := (f.moduli[0]*f.moduli[2])^(-1) mod f.moduli[1],
//                   k := (f.moduli[1]*f.moduli[2])^(-1) mod f.moduli[0] ||
{
  Digit y;
  Digit m01[2],m02[2],m12[2];
  Digit s[3],x[3],c[3],m[3];
  c[0] = c[1] = c[2] = 0;

  const size_t sT     = t[2].size;
  Digit* pT           = t[2].p+sT;
  const Digit* pA     = t[1].p+sT;
  const Digit* pB     = t[0].p+sT;
  const size_t nbase  = size_t(log2(moduli[0]));
  const Digit base    = (Digit(1) << nbase) - 1;

  digitmul(moduli[0], moduli[1], m01[0], m01[1]);
  digitmul(moduli[0], moduli[2], m02[0], m02[1]);
  digitmul(moduli[1], moduli[2], m12[0], m12[1]);
  digitmod(m01[0], m01[1], moduli[2], y);
  const Digit i = digitinv(y, moduli[2]);
  digitmod(m02[0], m02[1], moduli[1], y);
  const Digit j = digitinv(y, moduli[1]);
  digitmod(m12[0], m12[1], moduli[0], y);
  const Digit k = digitinv(y, moduli[0]);

  digitmul(moduli[2], m01[0], m01[1], m);

  do {
    digitmulmod(i, *--pT, moduli[2], y);

    // s = y * moduli[0] * moduli[1]
    digitmul(y, m01[0], m01[1], s);

    digitmulmod(j, *--pA, moduli[1], y);

    // s += y * moduli[0] * moduli[2]
    digitmul(y, m02[0], m02[1], x);
    bool carry;
    s[2] += x[2];
    if (s[2] < x[2]) {
      s[1] += x[1]+1;
      if (s[1] <= x[1]) {
        s[0] += x[0]+1;
        carry = (s[0] <= x[0]);
      } else {
        s[0] += x[0];
        carry = (s[0] < x[0]);
      }
    } else {
      s[1] += x[1];
      if (s[1] < x[1]) {
        s[0] += x[0]+1;
        carry = (s[0] <= x[0]);
      } else {
        s[0] += x[0];
        carry = (s[0] < x[0]);
      }
    }
    // if (s >= M) s -= M
    if (carry || s[0] > m[0]
        || s[0] == m[0] && (s[1] > m[1] || s[1] == m[1] && s[2] >= m[2])) {
      if (s[2] >= m[2]) {
        s[2] -= m[2];
        if (s[1] >= m[1]) {
          s[1] -= m[1];
          s[0] -= m[0];
        } else {
          s[1] -= m[1];
          s[0] -= m[0]+1;
        }
      } else {
        s[2] -= m[2];
        const Digit c = ~m[1];
        s[1] += c;
        if (s[1] >= c) s[0] -= m[0]+1;
        else s[0] -= m[0];
      }
    }

    digitmulmod(k, *--pB, moduli[0], y);

    // s += y * moduli[1] * moduli[2]
    digitmul(y, m12[0], m12[1], x);
    s[2] += x[2];
    if (s[2] < x[2]) {
      s[1] += x[1]+1;
      if (s[1] <= x[1]) {
        s[0] += x[0]+1;
        carry = (s[0] <= x[0]);
      } else {
        s[0] += x[0];
        carry = (s[0] < x[0]);
      }
    } else {
      s[1] += x[1];
      if (s[1] < x[1]) {
        s[0] += x[0]+1;
        carry = (s[0] <= x[0]);
      } else {
        s[0] += x[0];
        carry = (s[0] < x[0]);
      }
    }
    // if (s >= M) s -= M
    if (carry || s[0] > m[0]
        || s[0] == m[0] && (s[1] > m[1] || s[1] == m[1] && s[2] >= m[2])) {
      if (s[2] >= m[2]) {
        s[2] -= m[2];
        if (s[1] >= m[1]) {
          s[1] -= m[1];
          s[0] -= m[0];
        } else {
          s[1] -= m[1];
          s[0] -= m[0]+1;
        }
      } else {
        s[2] -= m[2];
        const Digit c = ~m[1];
        s[1] += c;
        if (s[1] >= c) s[0] -= m[0]+1;
        else s[0] -= m[0];
      }
    }

    // c += s
    c[2] += s[2];
    if (s[2] > c[2]) {
      c[1] += s[1]+1;
      if (s[1] >= c[1]) c[0] += s[0]+1;
      else c[0] += s[0];
    } else {
      c[1] += s[1];
      if (s[1] > c[1]) c[0] += s[0]+1;
      else c[0] += s[0];
    }

    // *pT = c & base
    *pT = c[2] & base;

    // c >>= nbase
    c[2] >>= nbase;
    c[2] |= c[1] << (BETA-nbase);
    c[1] >>= nbase;
    c[1] |= c[0] << (BETA-nbase);
    c[0] >>= nbase;
  } while (pT != t[2].p);

  const Digit d = c[2];
  if (d) {
    *--t[2].p = d;
    ++(t[2].size);
  } else ++shift;
}

#if defined(_DigitAsm_) && _M_IX86 >= 300 && defined(_MSC_VER)
void FFT::square(const size_t a)
{
  Digit* pT = t[a].p;
  Digit* pE = pT+t[a].size;
  Digit M = m;
  __asm {
    mov edi,pT
    mov ecx,pE
    mov ebx,M
    sub ecx,edi
    mov eax,ecx
    shr ecx,3
    and eax,4
    je  L1
    mov eax,[edi]
    mul eax
    div ebx
    mov [edi],edx
    add edi,4
L1: mov eax,[edi]
    mul eax
    div ebx
    mov [edi],edx
    mov eax,4[edi]
    mul eax
    div ebx
    mov 4[edi],edx
    add edi,8
    dec ecx
    jne L1
  }
}

void FFT::multiply(const size_t a)
{
  Digit* pT = t[a].p;
  Digit* pE = pT+t[a].size;
  Digit* pA = factor->t[0].p;
  Digit M = m;
  __asm {
    mov edi,pT
    mov esi,pA
    mov ecx,pE
    mov ebx,M
    sub ecx,edi
    mov eax,ecx
    shr ecx,3
    and eax,4
    je  L1
    mov edx,[esi]
    mov eax,[edi]
    mul edx
    div ebx
    mov [edi],edx
    add esi,4
    add edi,4
L1: mov edx,[esi]
    mov eax,[edi]
    mul edx
    div ebx
    mov [edi],edx
    mov eax,4[esi]
    mov edx,4[edi]
    mul edx
    div ebx
    mov 4[edi],edx
    add esi,8
    add edi,8
    dec ecx
    jne L1
  }
}
#else
void FFT::square(const size_t a)
// Algorithm:  b.square(a)
// Input:      b in FFT, a in size_t.
// Output:     b in FFT such that (b.t_a)_i := ((b.t_a)_i)^2 (mod b.m)
//             for 0 <= i < b.t_a.size ||
{
  Digit* pT = t[a].p;
  const Digit* pE = pT+t[a].size;
  const Digit M = m;
  do {
    const Digit x = *pT;
    digitmulmod(x, x, M, *pT);
  } while (++pT != pE);
}

void FFT::multiply(const size_t a)
// Algorithm:  b.multiply(a)
// Input:      b in FFT, a in size_t.
// Output:     b in FFT such that (b.t_a)_i := ((b.t_a)_i)*((v.factor->t_0)_i) (mod b.m)
//             for 0 <= i < b.t_a.size ||
{
  CONDITION(factor && factor->t[0].size == t[a].size);

  Digit* pT = t[a].p;
  Digit* pA = factor->t[0].p;
  const Digit* pE = pT+t[a].size;
  const Digit M = m;
  do {
    const Digit x = *pT;
    const Digit y = *pA++;
    digitmulmod(x, y, M, *pT);
  } while (++pT != pE);
}

#endif
void FFT::result(Natural& b) const
// Algorithm:  a.result(b)
// Input:      a in FFT, b in Natural.
// Output:     b in Natural such that b = [a.t_2 / 2^(beta*a.shift)] ||
{
  const size_t nbase = size_t(log2(moduli[0]));
  const size_t sB = (factor)? factor->arg.size+arg.size : 2*arg.size;
  Digit* pF = b.setsize(sB);
  Digit* pB = pF+sB;
  const Digit* pE = t[2].p;
  const Digit* pT = pE+size()-shift;
  size_t i = 0;
  Digit k = 0;
  do {
    Digit j = *--pT;
    k |= j << i;
    const size_t l = BETA-i;
    if (nbase >= l) {
      *--pB = k;
      k = j >> l;
      i -= BETA;
    }
    i += nbase;
  } while (pT != pE);
  if (k) *--pB = k;
  if (pB != pF) {
    b.size -= pB-pF;
    FILL_ZERO(pF, pB);
    b.p = pB;
  }

  NATURALCONDITION(b);
}

void FFT::sqr(Natural& b)
// Algorithm:  a.sqr(b)
// Input:      a in FFT, b in Natural.
// Output:     b in Natural such that b = (a.t_0)^2 ||
{
  t[2] = t[1] = t[0];
  const size_t sT = size();
  if (sT%3 == 0) {
    for (size_t i = 0; i < 3; ++i) {
      setmodulo(moduli[i]);
      Digit w  = pow(primroots[i], (m-1)/(sT/3), m);
      Digit w2 = pow(primroots[i], (m-1)/sT, m);
      init_omega(w);

      // fft:
      Digit* pT = t[i].p;
      innerfft3(pT, sT/3, w2);
      five_step(pT+sT, pT, w);

      square(i);

      // fftinv:
      w = digitinv(w, m);
      w2 = digitinv(w2, m);
      init_omega(w);
      five_step(pT+sT, pT, w, digitinv(sT, m));
      innerfftinv3(pT, sT/3, w2);
    }
  } else {
    for (size_t i = 0; i < 3; ++i) {
      setmodulo(moduli[i]);
      Digit w  = pow(primroots[i], (m-1)/sT, m);
      init_omega(w);

      // fft:
      Digit* pT = t[i].p;
      five_step(pT+sT, pT, w);

      square(i);

      // fftinv:
      w = digitinv(w, m);
      init_omega(w);
      five_step(pT+sT, pT, w, digitinv(sT, m));
    }
  }
  chinese_remainder();
  result(b);
}

void FFT::mul(Natural& b)
// Algorithm:  a.mul(b)
// Input:      a in FFT, b in Natural where not a.factor = 0
//             and a.factor->size = a.size.
// Output:     b in Natural such that b = (a.t_0)*(a.factor->t_0) ||
{
  CONDITION(factor && factor->size() == size());

  t[2] = t[1] = t[0];
  const size_t sT = size();
  if (sT%3 == 0) {
    for (size_t i = 0; i < 3; ++i) {
      setmodulo(moduli[i]);
      Digit w  = pow(primroots[i], (m-1)/(sT/3), m);
      Digit w2 = pow(primroots[i], (m-1)/sT, m);
      init_omega(w);

      // factor->fft:
      if (i >= 1) base_conversion(factor->t[0], factor->arg, factor->shift);
      Digit* pT = factor->t[0].p;
      innerfft3(pT, sT/3, w2);
      five_step(pT+sT, pT, w);

      // fft:
      pT = t[i].p;
      innerfft3(pT, sT/3, w2);
      five_step(pT+sT, pT, w);

      multiply(i);

      // fftinv:
      w = digitinv(w, m);
      w2 = digitinv(w2, m);
      init_omega(w);
      five_step(pT+sT, pT, w, digitinv(sT, m));
      innerfftinv3(pT, sT/3, w2);
    }
  } else {
    for (size_t i = 0; i < 3; ++i) {
      setmodulo(moduli[i]);
      Digit w  = pow(primroots[i], (m-1)/sT, m);
      init_omega(w);

      // factor->fft:
      if (i >= 1) base_conversion(factor->t[0], factor->arg, factor->shift);
      Digit* pT = factor->t[0].p;
      five_step(pT+sT, pT, w);

      // fft:
      pT = t[i].p;
      five_step(pT+sT, pT, w);

      multiply(i);

      // fftinv:
      w = digitinv(w, m);
      init_omega(w);
      five_step(pT+sT, pT, w, digitinv(sT, m));
    }
  }
  chinese_remainder();
  result(b);
}

static void init_order(size_t* a, size_t b)
// Algorithm:  init_order(a, b)
// Input:      a in size_t^b,
//             b in size_t where b > 2 and b = 2^x for x in N.
// Output:     a = [a_0, a_b[ in size_t^b such that
//             a_i = bitreverse(i) for 0 <= i < b ||
{
  CONDITION(b > 2 && (b & (~b+1)) == b);

  const size_t* c = a+b-1;
  b >>= 1;
  size_t j,i = b;
  *a = 0; *++a = i;
  do {
    for (j = b; j <= i; j >>= 1) i -= j;
    *++a = i += j;
  } while (a != c);
}

void FFT::base_conversion(Natural& t, const Natural& a, const size_t shift)
{
  size_t nbase = size_t(log2(moduli[0]));
  const Digit base = (Digit(1) << nbase) - 1;
  
  // base conversion:
  const Digit* pA = a.p;
  size_t sA = a.size;
  const Digit* pE = pA+sA;
  size_t sT = (sA*BETA)/nbase;
  sT += shift;
  Digit* pT = t.setsize(sT)+sT;
  if (shift) {
    Digit* pE = pT-shift;
    FILL_ZERO(pE, pT);
    pT -= shift;
  }
  Digit k = *--pE;
  *--pT = k & base;
  k >>= nbase;
  size_t i = BETA-nbase;
  do {
    Digit j = *--pE;
    k |= j << i;
    *--pT = k & base;
    j >>= BETA-i; k >>= nbase;
    j <<= BETA-nbase; k |= j;
    if (i >= nbase) {
      *--pT = k & base;
      k >>= nbase; i -= nbase;
    }
    i += BETA-nbase;
  } while (pE != pA);
  if (i >= nbase) 
    if (k) *--pT = k;
    else { *(pT-1) = 0; --sT; }
  else if (k) { *--pT = k; ++sT; }
  t.size = sT; t.p = pT;
}

FFT::FFT(const Natural& a, const bool init_static, FFT* b)
// Algorithm:  FFT(a, i, b)
// Input:      a in Natural, i in bool, *b in FFT
//             where 12 < L(a) <= FFT::max_size() and BETA in {32, 64} ||
 : factor(b), arg(a)
{
  NATURALCONDITION(a);
  CONDITION(a.size > 12 && (BETA == 32 || BETA == 64));
  CONDITION(moduli[0] < moduli[1] && moduli[1] < moduli[2]);
  CONDITION(a.size <= FFT::max_size());

  base_conversion(t[0], a, 0);
  size_t sT = t[0].size;
  shift = sT;

  if (init_static) {
    if (b) {
      sT += b->shift;
      const size_t j = size_t(1) << (1+log2(Digit(sT)-1));
      size_t i = (j/4) * 3;
      if (i >= sT) {
        sT = i - b->shift;
        i -= shift;
      } else {
        sT = j - b->shift;
        i = j - shift;
      }
      t[0].lmove(i);
      b->t[0].lmove(sT);
      b->shift = sT;
    } else {
      sT *= 2;
      const size_t j = size_t(1) << (1+log2(Digit(sT)-1));
      const size_t i = (j/4) * 3;
      sT = ((i >= sT)? i : j) - shift;
      t[0].lmove(sT);
    }
    shift = sT-shift;

    const size_t nn = t[0].size;
    n2 = size_t(log2(Digit((nn%3)? nn : (nn/3))));
    n1 = n2/2;
    n2 -= n1;
    n1 = 1 << n1; n2 = 1 << n2;

    omega2 = NOTHROW_NEW Digit[n2];
    if (!omega2) errmsg(2, "(FFT constructor)");
    order = NOTHROW_NEW size_t[n1];
    if (!order) errmsg(2, "(FFT constructor)");

    init_order(order, n1);

    if (n1 != n2) {
      omega = NOTHROW_NEW Digit[n1];
      if (!omega) errmsg(2, "(FFT constructor)");
    } else omega = omega2;


    CONDITION(n1 >= 4 && n2 >= n1);
  }
}

FFT::~FFT()
{
  if (n1 != n2) delete[] omega;
  delete[] order;
  delete[] omega2;
  omega = omega2 = 0;
  order = 0;
  m = 0; n1 = n2 = 0;
}

