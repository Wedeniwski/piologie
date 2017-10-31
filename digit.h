/////////////////////////////////
//
// Piologie V 1.3
// multi-precision arithmetic
// Digit / NumberBase
//
// (c) 1996-1999 HiPiLib
//               www.hipilib.de
//
// Sebastian Wedeniwski
// 12/13/1999
//

#ifndef _Include_Digit_H_
#define _Include_Digit_H_

#ifdef _Old_STD_
#include <stdlib.h>

#include <iostream.h>
#include <ctype.h>
#include <assert.h>
#include <limits.h>


# define bool int
# define true 1
# define false 0

# define OSTREAM ostream
# define ISTREAM istream

#else
// New standard

# include <iostream>
# include <cstdlib>
# include <cctype>
# include <cassert>
# include <climits>

# include <algorithm>
# include <utility>

# define OSTREAM std::ostream
# define ISTREAM std::istream

#endif

typedef long SignDigit;

typedef unsigned long Digit;

#define _DigitAsm_          // use assembler if ever it is possible
//#define _Piologie_Debug_    // Debug modus for assert

// error-handling:
typedef void (*piologie_error_handler_t)(const int, const char*);
piologie_error_handler_t set_piologie_error_handler(piologie_error_handler_t);

#ifndef STATIC_VS_INLINE        // some compiler have problems with inline
#define STATIC_VS_INLINE inline
#endif

STATIC_VS_INLINE Digit log2(Digit);
Digit pow10(size_t);
inline void  swap(Digit&, Digit&);
inline Digit sqrt(Digit, Digit);
Digit sqrt(Digit);
void  sqrt(Digit, Digit&, Digit&);
Digit sqrt(Digit, Digit, Digit&, Digit&);
Digit gcd(Digit, Digit);

#ifdef DEFINE_VS_CONST        // some compiler have problems with const
# define BETA       (sizeof(Digit)*CHAR_BIT)
# define GAMMA      (~Digit(0))
# define GAMMA_LOW  (GAMMA >> (BETA/2))
# define GAMMA_HIGH (~GAMMA_LOW)
#else
const size_t BETA        = sizeof(Digit)*CHAR_BIT;

const Digit  GAMMA       = ~Digit(0);                 // 2^BETA - 1

const Digit  GAMMA_LOW   = GAMMA >> (BETA/2);
const Digit  GAMMA_HIGH  = ~GAMMA_LOW;

#endif
const size_t DELTA       = 8;                         // >= 2!

const size_t ALPHA_WIDTH = size_t(0.301029995664*BETA);
const Digit  ALPHA       = pow10(ALPHA_WIDTH);


class NumberBase {
protected:
  NumberBase();                               // default constructor
  ~NumberBase();                              // destructor
  
  size_t* quad_convergence_sizes(const size_t, size_t&) const;
  
  void digitmul(const Digit, const Digit, const Digit, Digit*) const;
  void digitmul(const Digit, const Digit, const Digit, const Digit, Digit*) const;
  void digitsqr(const Digit, const Digit, Digit*) const;
  
public:
  static Digit* k_ary_gcd_inverse;
  static Digit* k_ary_gcd_linear_id;
  static bool* k_ary_gcd_linear_id_neg;
  static Digit* k_ary_gcd_trial_division;
  static size_t k_ary_gcd_trial_division_size;

  void digitmul(const Digit, const Digit, Digit&, Digit&) const;
  void digitdiv(const Digit, const Digit, const Digit, Digit&, Digit&) const;
  void digitmod(Digit a, Digit b, const Digit c, Digit& r) const;
  
  void gcd(Digit, Digit, Digit, Digit, Digit&, Digit&) const;

  void errmsg(const int, const char*) const;
};

/////////////////// Inline-Implementation ///////////////////////

inline NumberBase::NumberBase()
{
}

inline NumberBase::~NumberBase()
{
}

#ifdef _DigitAsm_
#if defined(__386__) && defined(__WATCOMC__)
void _digitmul(const Digit, const Digit, Digit&, Digit&);
#pragma aux _digitmul =               \
   "mul ebx"                          \
   "mov [esi],edx"                    \
   "mov [edi],eax"                    \
   parm [EAX] [EBX] [ESI] [EDI]       \
   modify [ ESI EDI EAX EDX ];

inline void NumberBase::digitmul(const Digit a, const Digit b,
                                 Digit& x, Digit& y) const
{
  _digitmul(a, b, x, y);
}

void _digitdiv(const Digit, const Digit, const Digit, Digit&, Digit&);
#pragma aux _digitdiv =               \
   "div ebx"                          \
   "mov [esi],eax"                    \
   "mov [edi],edx"                    \
   parm [EDX] [EAX] [EBX] [ESI] [EDI] \
   modify [ ESI EDI EAX EDX ];

inline void NumberBase::digitdiv(const Digit a, const Digit b, const Digit c,
                                 Digit& q, Digit& r) const
{
  _digitdiv(a, b, c, q, r);
}

#elif _M_IX86 >= 300 && _MSC_VER
static void __declspec(naked) __fastcall
            _digitmul(Digit* z, const Digit a, const Digit b)
{
  __asm {
    mov eax,[esp+4]
    mul edx
    mov [ecx],edx
    mov 4[ecx],eax
    ret 4
  }
}

inline void NumberBase::digitmul(const Digit a, const Digit b, Digit& x, Digit& y) const
{
  Digit z[2];
  _digitmul(z, a, b);
  x = z[0]; y = z[1];
}

static void __declspec(naked) __fastcall _digitdiv(Digit* q, Digit* r, const Digit a,
                                                   const Digit b, const Digit c)
{
  __asm {
    push edi
    mov edi,edx
    mov edx,[esp+8]
    mov eax,[esp+12]
    div dword ptr [esp+16]
    mov [ecx],eax
    mov [edi],edx
    pop edi
    ret 12
  }
}

inline void NumberBase::digitdiv(const Digit a, const Digit b, const Digit c,
                                 Digit& q, Digit& r) const
{
  _digitdiv(&q, &r, a, b, c);
}

#elif (defined (__i386__) || defined (__i486__)) && defined (__GNUC__)
inline void NumberBase::digitmul(const Digit a, const Digit b,
                                 Digit& x, Digit& y) const
{
  __asm__ ("mull %3"
    : "=a" (y), "=d" (x)
    : "%0" (a), "rm" (b));
}

inline void NumberBase::digitdiv(const Digit a, const Digit b, const Digit c,
                                 Digit& q, Digit& r) const
{
  __asm__ ("divl %4"
    : "=a" (q), "=d" (r)
    : "0"  (b), "1"  (a), "rm" (c));
}

#elif defined (__sparc_v8__) && defined (__GNUC__)
inline void NumberBase::digitmul(const Digit a, const Digit b,
                                 Digit& x, Digit& y) const
{
  __asm__ ("umul %2,%3,%1; rd %%y,%0"
    : "=r" (x), "=r" (y)
    : "r" (a), "r" (b));
}

inline void NumberBase::digitdiv(const Digit a, const Digit b, const Digit c,
                                 Digit& q, Digit& r) const
{
  __asm__ ("mov %2,%%y;nop;nop;nop;udiv %3,%4,%0;umul %0,%4,%1;sub %3,%1,%1"
    : "=&r" (q), "=&r" (r)
    : "r" (a), "r" (b), "r" (c));
}

#else
# undef _DigitAsm_
#endif

#endif // _DigitAsm_

#ifndef _DigitAsm_
inline void NumberBase::digitmul(const Digit a, const Digit b,
                                 Digit& x, Digit& y) const
// Algorithm:  n.digitmul(a, b, x, y)
// Input:      n in NumberBase, a,b in Digit.
// Output:     x,y in Digit such that x*2^BETA + y = a*b ||
{
  const Digit c = a >> BETA/2;
  const Digit d = a & GAMMA_LOW;
  const Digit e = b >> BETA/2;
  const Digit f = b & GAMMA_LOW;

  const Digit z1 = c*f;
  Digit z2 = c*e;
  Digit z3 = d*f;
  Digit z4 = d*e + z1;
  if (z4 < z1) z2 += GAMMA_LOW+1;
  z2 += z4 >> BETA/2;
  z3 += z4 <<= BETA/2;
  if (z3 < z4) ++z2;
  y = z3; x = z2;
}

#endif

inline void NumberBase::digitmod(Digit a, Digit b, const Digit c, Digit& r) const
// Algorithm:  n.digitmod(a, b, c, r)
// Input:      n in NumberBase, a,b,c in Digit where not c = 0.
// Output:     r in Digit such that r = a*2^BETA+b - c*[(a*2^BETA+b)/c] ||
{
#ifdef _DigitAsm_
  if (a < c) {
    Digit d;
    digitdiv(a, b, c, d, r);
  } else {
    Digit d,e = 0;
    digitdiv(e, a, c, d, e);
    digitdiv(e, b, c, d, r);
  }
#else
  if (c > GAMMA/2) {
    Digit d,e = 0;
    digitdiv(e, a, c, d, e);
    digitdiv(e, b, c, d, r);
  } else {
    const Digit n = log2(c)+1;
    const Digit n2 = BETA-n;
    const Digit c2 = c << n2;
    Digit d,z = a >> n;
    digitdiv(z, (a << n2) | (b >> n), c2, d, z);
    digitdiv(z, b << n2, c2, d, z);
    r = z >> n2;
  }
#endif
}

inline void NumberBase::digitmul(const Digit a0, const Digit b0, const Digit b1,
                                 Digit* c) const
// Algorithm:  n.digitmul(a0, b0, b1, d)
// Input:      n in NumberBase, a0,b0,b1 in Digit.
// Output:     d = [d0, d1, d2] in Digit^4
//             such that d0*4^BETA+d1*2^BETA+d2 = a0*(b0*2^BETA+b1) ||
{
  Digit x,y;
  digitmul(a0, b0, c[0], x);
  digitmul(a0, b1, y, c[2]);
  c[1] = y += x;
  if (y < x) ++c[0];
}

inline void NumberBase::digitmul(const Digit a0, const Digit a1, const Digit b0,
                                 const Digit b1, Digit* c) const
// Algorithm:  n.digitmul(a0, a1, b0, b1, d)
// Input:      n in NumberBase, a0,a1,b0,b1 in Digit.
// Output:     d = [d0, d1, d2, d3] in Digit^4 such that
//             d0*8^BETA+d1*4^BETA+d2*2^BETA+d3 = (a0*2^BETA+a1)*(b0*2^BETA+b1) ||
{
  Digit x0,x1,y0,y1,c0,c1,c2;
  digitmul(a0, b0, c0, c1);
  digitmul(a0, b1, x0, x1);
  digitmul(a1, b0, y0, y1);
  digitmul(a1, b1, c2, c[3]);
  x1 += y1;
  if (x1 < y1) {
    x0 += y0+1;
    if (x0 <= y0) ++c0;
  } else {
    x0 += y0;
    if (x0 < y0) ++c0;
  }
  c[2] = c2 += x1;
  if (c2 < x1) {
    c[1] = c1 += x0+1;
    if (c1 <= x0) ++c0;
  } else {
    c[1] = c1 += x0;
    if (c1 < x0) ++c0;
  }
  c[0] = c0;
}

inline void NumberBase::digitsqr(const Digit a0, const Digit a1, Digit* b) const
// Algorithm:  n.digitsqr(a0, a1, b)
// Input:      n in NumberBase, a0,a1 in Digit.
// Output:     b = [b0, b1, b2, b3] in Digit^4
//             such that b0*8^BETA+b1*4^BETA+b2*2^BETA+b3 = (a0*2^BETA+a1)^2 ||
{
  Digit x0,x1,b0,b1,b2;
  digitmul(a0, a0, b0, b1);
  digitmul(a0, a1, x0, x1);
  digitmul(a1, a1, b2, b[3]);
  b0 += x0 >> (BETA-1); x0 <<= 1; x0 |= x1 >> (BETA-1); x1 <<= 1;
  b[2] = b2 += x1;
  if (b2 < x1) {
    b[1] = b1 += x0+1;
    if (b1 <= x0) ++b0;
  } else {
    b[1] = b1 += x0;
    if (b1 < x0) ++b0;
  }
  b[0] = b0;
}

STATIC_VS_INLINE Digit log2(Digit a)
// Algorithm:  b := log2(a)
// Input:      a in Digit.
// Output:     b in Digit
//             such that if a > 0 then b = [log2(a)] else b = 0 ||
{
  Digit b = 0;

#if CHAR_BIT == 8
  static const Digit c[16] = {0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3 };

  if (a > GAMMA_LOW) { b += BETA/2; a >>= BETA/2; }
  if (a >= Digit(1) << BETA/4) { b += BETA/4; a >>= BETA/4; }
  if (a >= Digit(1) << BETA/8) { b += BETA/8; a >>= BETA/8; }
  if (BETA > 32) {
    if (a >= Digit(1) << BETA/16) { b += BETA/16; a >>= BETA/16; }
    if (BETA > 64) {
      while (a >= 16) { b += 4; a >>= 4; }
    }
  }
  return b+c[a];
#else
  while (a >>= 1) ++b;
  return b;
#endif
}

inline void swap(Digit& a, Digit& b)
// Algorithm:  swap(a, b)
//             Let t in Digit.
// Input:      a,b in Digit.
// Output:     a,b in Digit such that t := a, a := b, b := t ||
{
  Digit t = a; a = b; b = t;
}

inline Digit sqrt(Digit a, Digit b)
// Algorithm:  c := sqrt(a, b)
// Input:      a,b in Digit.
// Output:     c in Digit such that c = [sqrt(a*2^BETA+b)] ||
{
  Digit x,y;
  return sqrt(a, b, x, y);
}


#if defined(_Old_STD_) && !defined(__MINMAX_DEFINED) || defined(_MSC_VER)

#ifdef min
# undef min
#endif
#ifdef max
# undef max
#endif

template <class T>
inline const T& min(const T& a, const T& b)
{
  return (a < b)? a : b;
}

template <class T>
inline const T& max(const T& a, const T& b)
{
  return (a < b)? b : a;
}
#else
# define min std::min
# define max std::max
#endif

template <class Arg1, class Arg2, class Tag>
struct binder_arguments {
  const Arg1& x;
  const Arg2& y;

  binder_arguments(const Arg1& a, const Arg2& b) : x(a), y(b) {}
};



#endif
