/////////////////////////////////
//
// Piologie V 1.3.2
// multi-precision arithmetic
// Digit / NumberBase
//
// (c) 1996-2001 HiPiLib
//               www.hipilib.de
//
// Sebastian Wedeniwski
// 10/21/2001
//

#include "digit.h"
#include "nmbrthry.h"
#ifdef _NEEDS_PIOLOGIE_KEY_
# include "key.h"
#endif

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
  do *a++ = 0; while (a != b)

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



///////////////////// Class NumberBase /////////////////////////////

Digit* NumberBase::k_ary_gcd_inverse = 0;
Digit* NumberBase::k_ary_gcd_linear_id = 0;
bool* NumberBase::k_ary_gcd_linear_id_neg = 0;
Digit* NumberBase::k_ary_gcd_trial_division = 0;
size_t NumberBase::k_ary_gcd_trial_division_size = 0;


class NumberBase_init : private NumberBase {
private:
  static NumberBase_init CheckNumberBase;

public:
  NumberBase_init(char);
  ~NumberBase_init();
};

NumberBase_init NumberBase_init::CheckNumberBase(' ');

NumberBase_init::NumberBase_init(char)
{
  // Testing the bits for low ending
  {
    const Digit d = 1 << 1;
    if (d != 2) errmsg(0, "Bitproblem!");
    
  }
  // Testing size of Digit and SignDigit
  if (sizeof(Digit) != sizeof(SignDigit)) errmsg(0, "Wrong size of SignDigit!");
  
  // precomputation phase for k-ary gcd
  const Digit K_ARY = Digit(1) << ((BETA > 32)? 16 : BETA/2);
  const Digit K_BASE = K_ARY-1;
  
  size_t k = K_ARY >> 1;
  k_ary_gcd_inverse = NOTHROW_NEW Digit[k];
  if (!k_ary_gcd_inverse) errmsg(2, "(precomputation phase for k-ary gcd)");
  size_t i;
  for (i = 0; i < k; ++i) k_ary_gcd_inverse[i] = inverse(Digit(i+i+1), K_ARY);
  k_ary_gcd_linear_id = NOTHROW_NEW Digit[K_ARY];
  k_ary_gcd_linear_id_neg = NOTHROW_NEW bool[K_ARY];
  if (!k_ary_gcd_linear_id || !k_ary_gcd_linear_id_neg) errmsg(2, "(precomputation phase for k-ary gcd)");
  k = sqrt(K_ARY);
  for (i = 1; i < K_ARY; ++i) {
    for (Digit a = 1; a <= k; ++a) {
      const Digit b = (a*i)&K_BASE;
      const Digit c = K_ARY - b;
      if (b <= k || c <= k) { k_ary_gcd_linear_id[i] = a; k_ary_gcd_linear_id_neg[i] = (b <= k); break; }
    }
  }
  
  k_ary_gcd_trial_division = NOTHROW_NEW Digit[++k];
  if (!k_ary_gcd_trial_division) errmsg(2, "(precomputation phase for k-ary gcd)");
  k_ary_gcd_trial_division_size = 0;
  k_ary_gcd_trial_division[k_ary_gcd_trial_division_size] = 1;
  Digit a;
  Primes p;
  p.firstPrime();
  while (p.nextPrime(a) && a <= k) {
    Digit d1,d2;
    digitmul(a, k_ary_gcd_trial_division[k_ary_gcd_trial_division_size], d1, d2);
    if (d1) k_ary_gcd_trial_division[++k_ary_gcd_trial_division_size] = a;
    else k_ary_gcd_trial_division[k_ary_gcd_trial_division_size] = d2;
  }
  Digit* d = NOTHROW_NEW Digit[++k_ary_gcd_trial_division_size];
  if (!d) errmsg(2, "(precomputation phase for k-ary gcd)");
  const Digit* dE = d+k_ary_gcd_trial_division_size;
  COPY(d, k_ary_gcd_trial_division, d, dE);
  k_ary_gcd_trial_division -= k_ary_gcd_trial_division_size;
  delete[] k_ary_gcd_trial_division;
  k_ary_gcd_trial_division = d-k_ary_gcd_trial_division_size;
  
#ifdef _NEEDS_PIOLOGIE_KEY_
  checkPiologieKey();
#endif
}

NumberBase_init::~NumberBase_init()
{
  delete[] k_ary_gcd_inverse;
  delete[] k_ary_gcd_linear_id;
  delete[] k_ary_gcd_linear_id_neg;
  delete[] k_ary_gcd_trial_division;
  k_ary_gcd_inverse = 0;
  k_ary_gcd_linear_id = 0;
  k_ary_gcd_linear_id_neg = 0;
  k_ary_gcd_trial_division = 0;
}

size_t* NumberBase::quad_convergence_sizes(const size_t a, size_t& b) const
// Algorithm:  s := c.quad_convergence_sizes(a, b)
// Input:      c in NumberBase, a,b in size_t where a >= 2.
// Output:     b in size_t, s in size_t^b such that b = log2(a)+2,
//             s[i] = [a/2^i]+epsilon, 0 <= i < b ||
{
  CONDITION(a >= 2);

  size_t* s = NOTHROW_NEW size_t[log2(a)+2];
  if (!s) errmsg(2, "(quad_convergence_sizes)");
  size_t i = 1;
  size_t j = a+2;
  s[0] = j; s[1] = j /= 2;
  do { ++j; s[++i] = j /= 2; } while (j > 1);
  b = i+1;

  CONDITION(b >= 2);

  return s;
}

void default_piologie_error_handler(const int IDerr, const char* msg)
{
  switch(IDerr) {
    case 0: break;
//    case 1: cerr << "Overflow!"; break;
    case 2: cerr << "Out of memory!"; break;
    case 3: cerr << "Result negative (No Natural solution)!"; break;
    case 4: cerr << "Division by zero!"; break;
    case 5: cerr << "Same location error!"; break;
    case 6: cerr << "Discriminant negative (No Real solution)!"; break;

    default: cerr << "Internal error!";
  }
  cerr << msg << endl;
  exit(1);
}

static piologie_error_handler_t
 piologie_error_handler = default_piologie_error_handler;

piologie_error_handler_t set_piologie_error_handler(piologie_error_handler_t a)
{
  piologie_error_handler_t b = piologie_error_handler;
  piologie_error_handler = a;
  return b;
}

void NumberBase::errmsg(const int a, const char* b) const
{
  (*piologie_error_handler)(a, b);
}


#ifndef _DigitAsm_
void NumberBase::digitdiv(const Digit a, const Digit b, const Digit c,
                          Digit& q, Digit& r) const
// Algorithm:  n.digitdiv(a, b, c, q, r)
// Input:      n in NumberBase, a,b,c in Digit where a < c, c >= 2^(BETA-1).
// Output:     q,r in Digit such that q*c + r = a*2^BETA + b ||
{
  CONDITION(a < c && c >= (Digit(1) << (BETA-1)));

  Digit d = c >> BETA/2;
  Digit e = c & GAMMA_LOW;
  Digit x = a/d;
  Digit z = a-x*d;
  Digit y = x*e;
  z = (z << BETA/2) | (b >> BETA/2);
  if (z < y) {
    --x; z += c;
    if (z >= c && z < y) { --x; z += c; }
  }
  q = x << BETA/2;
  z -= y;
  x = z/d;
  z -= x*d;
  y = x*e;
  z = (z << BETA/2) | (b & GAMMA_LOW);
  if (z < y) {
    --x; z += c;
    if (z >= c && z < y) { --x; z += c; }
  }
  q |= x; r = z-y;
}

#endif

Digit pow10(size_t a)
// Algorithm:  b := pow10(a)
// Input:      a in size_t where a <= ALPHA_WIDTH.
// Output:     b in Digit such that b = 10^a ||
{
  CONDITION(a <= ALPHA_WIDTH);

  static const Digit c[10] = {1, 10, 100, Digit(1000), Digit(10000),
                              Digit(100000), Digit(1000000),
                              Digit(10000000), Digit(100000000),
                              Digit(1000000000) };

  if (ALPHA_WIDTH < 10 || a < 10) return c[a];
  Digit b = c[9];
  for (a -= 9; a >= 9; a -= 9) b *= c[9];
  return b*c[a];
}

Digit sqrt(Digit a)
// Algorithm:  c := sqrt(a)
// Input:      a in Digit.
// Output:     c in Digit such that c = [sqrt(a)] ||
{
  Digit b = Digit(1) << (BETA-2);
  Digit c = Digit(1) << (BETA-1) | Digit(1) << (BETA-2);

  do {
    if (a >= b) { a -= b; b |= c; }
    c >>= 2;
    b >>= 1;
    b ^= c;
  } while (c != 3);

  if (a >= b) b |= 2;
  return b >>= 1;
}

void sqrt(Digit a, Digit& x, Digit& y)
// Algorithm:  sqrt(a, x, y)
// Input:      a in Digit.
// Output:     x,y in Digit such that x = [sqrt(a)], y = a - x^2 ||
{
  Digit b = Digit(1) << (BETA-2);
  Digit c = Digit(1) << (BETA-1) | Digit(1) << (BETA-2);

  do {
    if (a >= b) { a -= b; b |= c; }
    c >>= 2;
    b >>= 1;
    b ^= c;
  } while (c != 3);

  if (a >= b) { a -= b; b |= 2; }
  y = a; x = b >>= 1;
}

Digit sqrt(Digit a, Digit b, Digit& x, Digit& y)
// Algorithm:  c := sqrt(a, b, x, y)
// Input:      a,b in Digit.
// Output:     c,x,y in Digit such that c = [sqrt(a*2^BETA+b)],
//             x*2^BETA+y = a*2^BETA+b - c^2 ||
{
  const Digit HIBIT = Digit(1) << (BETA-1);
  Digit c = Digit(1) << (BETA-2);
  Digit d = Digit(1) << (BETA-1) | Digit(1) << (BETA-2);

  do {
    if (a >= c) { a -= c; c |= d; }
    d >>= 2;
    c >>= 1;
    c ^= d;
  } while (d != 3);
  if (a >= c) { a -= c; c |= 2; }
  c >>= 1;

  d = Digit(1) << (BETA-1) | Digit(1) << (BETA-2);
  Digit e = Digit(1) << (BETA-2);
  do {
    if (a > c) {
      a -= c;
      if (b < e) --a;
      b -= e;
      e |= d;
    } else if (a == c && b >= e) { a = 0; b -= e; e |= d; }
    d >>= 2;
    e >>= 1;
    if (c&1) e |= HIBIT;
    c >>= 1;
    e ^= d;
  } while (d);
  x = a; y = b;
  return e;
}

Digit gcd(Digit a, Digit b)
// Algorithm:  c := gcd(a, b)
// Input:      a,b in Digit.
// Output:     c in Digit such that c = gcd(a, b) ||
{
  if (a == 0) return b;
  if (b == 0) return a;

  Digit i,j;
  for (i = 0; (a&1) == 0; ++i) a >>= 1;
  for (j = 0; (b&1) == 0; ++j) b >>= 1;
  while (a != b)
    if (a > b) {
      a -= b;
      do a >>= 1; while ((a&1) == 0);
    } else {
      b -= a;
      do b >>= 1; while ((b&1) == 0);
    }
  return a << ((i > j)? j : i);
}

void NumberBase::gcd(Digit a, Digit b, Digit c, Digit d, Digit& x, Digit& y) const
// Algorithm:  n.gcd(a, b, c, d; x, y)
// Input:      n in NumberBase, a,b,c,d in Digit.
// Output:     x,y in Digit such that x*2^BETA+y = gcd(a*2^BETA+b, c*2^BETA+d) ||
{
  if (a == 0 && b == 0) { x = c; y = d; }
  else if (c == 0 && d == 0) { x = a; y = b; }
  else {
    Digit i = 0;
    Digit j = 0;
    if (b == 0) {
      b = a; a = 0;
      for (i = BETA; (b&1) == 0; ++i) b >>= 1;
    } else if ((b&1) == 0) {
      do { ++i; b >>= 1; } while ((b&1) == 0);
      b |= a << (BETA-i); a >>= i;
    }
    if (d == 0) {
      d = c; c = 0;
      for (j = BETA; (d&1) == 0; ++j) d >>= 1;
    } else if ((d&1) == 0) {
      do { ++j; d >>= 1; } while ((d&1) == 0);
      d |= c << (BETA-j); c >>= j;
    }
    if (j < i) i = j;
    while (b != d && a != c)
      if (a > c) {
        a -= c + (b < d); b -= d;
        j = 0;
        do { ++j; b >>= 1; } while ((b&1) == 0);
        b |= a << (BETA-j); a >>= j;
      } else {
        c -= a + (d < b); d -= b;
        j = 0;
        do { ++j; d >>= 1; } while ((d&1) == 0);
        d |= c << (BETA-j); c >>= j;
      }
    if (b != d || a != c) {
      if (a == c) {
        if (b > d) {
          b -= d;
          digitmod(c, d, b, a);
          a = ::gcd(a, b);
        } else {
          d -= b;
          digitmod(a, b, d, c);
          a = ::gcd(c, d);
        }
      } else {
        if (a > c) {
          a -= c;
          digitmod(c, d, a, b);
          a = ::gcd(a, b);
        } else {
          c -= a;
          digitmod(a, b, c, d);
          a = ::gcd(c, d);
        }
      }
      if (i >= BETA) { x = a << (i-BETA); y = 0; }
      else if (i) { x = a >> (BETA-i); y = a << i; }
      else { x = 0; y = a; }
    } else if (i >= BETA) { x = b << (i-BETA); y = 0; }
    else if (i) { x = a << i | b >> (BETA-i); y = b << i; }
    else { x = a; y = b; }
  }
}

