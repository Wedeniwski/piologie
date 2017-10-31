/////////////////////////////////
//
// Piologie V 1.3
// multi-precision arithmetic
// intensive test
//
// (c) 1996-1999 HiPiLib
//               www.hipilib.de
//
// Sebastian Wedeniwski
// 12/13/1999
//

#include <cstdio>
#include <ctime>
#include <fstream>
#include "rational.h"
#include "nmbrthry.h"

using namespace std;

// BETA >= 16!

size_t ITERATION;
Digit MAX_FACTORIAL;
Digit MAX_FIBONACCI;
Digit MAX_RANDOM;
Digit INNER_ITERATION;



template <class T>
T gcd2(T a, T b)
{
  T c;
  while (b != 0) { c = a%b; a = b; b = c; }
  return abs(a);
}

Natural factorial2(Natural a)
{
  Natural b = 1;
  while (a > 0) { b *= a; --a; }
  return b;
}

Natural fibonacci2(Natural a)
{
  if (a <= 1) return a;
  Natural b = 1;
  Natural t,c = Digit(0);
  do { t = b; b += c; c = t; } while (--a > 1);
  return b;
}

Natural fibonacci3(Natural a)
{
  Natural b = 1;
  Natural c = Digit(0);
  while (a != 0) { c += b; b = c-b; --a; }
  return c;
}

Natural fibonacci4(Natural a)
{
  if (a == 0) return Digit(0);
  Natural b = 1;
  Natural t,c = Digit(0);
  while (--a != 0) {
    t = c+b;
    if (--a == 0) return t;
    c = b+t;
    if (--a == 0) return c;
    b = t+c;
  }
  return b;
}


template <class T, class D>
void identity_generic(T a, D b)
{
  T r,s,t;

  assert((a + b) == (b + a));
  s = a + b; t = b + a;
  assert(s == t);
  r = s;
  s = a; s += b;
  t = b; t += a;
  assert(r == s && s == t);

  assert((a * b) == (b * a));
  s = a * b; t = b * a;
  assert(s == t);
  r = s;
  s = a; s *= b;
  t = b; t *= a;
  assert(r == s && s == t);

  assert(((a + b) - b) == a);
  t = a + b; t = t - b;
  assert(t == a);
  t = a; t += b; t -= b;
  assert(t == a);

  assert(((a + b) - a) == b);
  t = a + b; t = t - a;
  assert(t == b);
  t = a; t += b; t -= a;
  assert(t == b);

  t = b;
  assert((T(b) - t) == 0);    // if T == Natural
  s = T(b) - t;
  assert(s == 0);
  s = b; s -= t;
  assert(s == 0);
  assert((t - b) == 0);
  s = t - b;
  assert(s == 0);
  s = t; s -= b;
  assert(s == 0);

  assert(((a * b) / b) == a);
  assert(((a * b) / a) == b);
  assert(((a * b) / (b * a)) == 1);
}

template <class T, class D>
void identity_Natural_Digit(T a, D b)
{
  T r,s,t;

  identity_generic(a, b);

  assert(((a * b) / b) == a);
  assert(((a * b) % b) == 0);
  t = a * b; div(t, T(b), t, s);
  assert(t == a && s == 0);
  t = a; t *= b; t /= b;
  s = a; s *= b; s %= b;
  assert(t == a && s == 0);

  assert(((a * b) / a) == b);
  assert(((a * b) % a) == 0);
  t = a * b; div(t, a, t, s);
  assert(t == b && s == 0);
  t = a; t *= b; t /= a;
  s = a; s *= b; s %= a;
  assert(t == b && s == 0);

  assert(((a * b) / (b * a)) == 1);
  t = a * b; s = b * a; div(t, s, t, s);
  assert(t == 1 && s == 0);
  t = a; t *= b; s = b; s *= a; r = t; t /= s; s %= r;
  assert(t == 1 && s == 0);

  assert((b * (a / b) + (a % b)) == a);
  div(a, T(b), s, t); s = b * s; s = t + s;
  assert(s == a);
  s = a; s /= b; s *= b; t = a; t %= b; s += t;
  assert(s == a);

  // Logic:
  assert((a | b) == (b | a));
  s = a | b; t = b | a;
  assert(s == t);
  s = a; s |= b; t = b; t |= a;
  assert(s == t);

  assert((a & b) == (b & a));
  s = a & b; t = b & a;
  assert(s == t);
  s = a; s &= b; t = b; t &= a;
  assert(s == t);

  assert((a ^ b) == (b ^ a));
  s = a ^ b; t = b ^ a;
  assert(s == t);
  s = a; s ^= b; t = b; t ^= a;
  assert(s == t);

  assert((a & (a | b)) == a);
  s = a | b; s = a & s;
  assert(s == a);
  s = a; s |= b; s &= a;
  assert(s == a);

  assert((a | (a & b)) == a);
  s = a & b; s = a | s;
  assert(s == a);
  s = a; s &= b; s |= a;
  assert(s == a);
}

template <class T, class D>
void identity_Natural_Digit2(const T a, const D b)
// a >= b if a in Natural.
{
  T s,t;

  s = a-b;
  assert(s == a-b);
  t = a; t -= b;
  assert(s == t);

  assert(((a - b) + b) == a);
  s = a-b; s = s+b;
  assert(s == a);
  t = a; t -= b; t += b;
  assert(s == t);
}

template <class T>
void identity_generic2(const T a, const T b, const T c)
{
  T i,j,r,s,t;

  assert((a + (b + c)) == ((a + b) + c));
  s = b+c; s = a+s; t = a+b; t = t+c;
  assert(s == t);
  r = s;
  s = b; s += c; s += a;
  t = a; t += b; t += c;
  assert(r == s && s == t);

  assert((a * (b * c)) == ((a * b) * c));
  s = b*c; s = s*a; t = a*b; t = c*t;
  assert(s == t);
  r = s;
  s = b; s *= c; s *= a;
  t = a; t *= b; t *= c;
  assert(r == s && s == t);

  assert(((a + b) * (b + c)) == (a*b + b*b + a*c + b*c));
  s = a+b; r = b+c; s = s*r;
  t = a*b; r = b*b; t = t+r;
  r = a*c; t = t+r; r = b*c; t = r+t;
  assert(s == t);
  s = a; s += b; r = b; r += c; s *= r;
  t = a; t *= b; r = b; r *= b; t += r;
  r = a; r *= c; t += r; r = b; r *= c; t += r;
  assert(s == t);

  assert(((a + c) * (c + a)) == (a*a + 2*a*c + c*c));
  s = a+c; r = c+a; s = s*r;
  assert(((a + c) * (c + a)) == s);
  t = a*a; r = a*c; r = r*2; t = r+t;
//  assert(r.even());
  assert(t == (a*a + 2*a*c));
  r = c*c; t = t+r;
  assert(s == t);
  s = a; s += c; r = c; r += a; s *= r;
  t = a; t *= a; r = 2; r *= a; r *= c; t += r;
  r = c; r *= c; t += r;
  assert(s == t);

  assert((a * (b + c)) == ((b * a) + (c * a)));
  s = b+c; s = a*s;
  t = b*a; r = c*a; t = t+r;
  assert(s == t);
  r = b; r += c; s = a; s *= r;
  t = b; t *= a; r = c; r *= a; t += r;
  assert(s == t);

  assert(((a + b) * c) == ((c * a) + (c * b)));
  s = a+b; s = s*c;
  t = c*a; r = c*b; t = t+r;
  assert(s == t);
  s = a; s += b; s *= c;
  t = c; t *= a; r = c; r *= b; t += r;
  assert(s == t);

  assert((a - a) == 0);
  t = a-a;
  assert(t == 0);
  t = a; t -= a;
  assert(t == 0);

  assert((a << 0) == a);
  t = a << 0;
  assert(t == a);
  t = a; t <<= 0;
  assert(t == a);

  assert((a >> 0) == a);
  t = a >> 0;
  assert(t == a);
  t = a; t >>= 0;
  assert(t == a);

  // constructor:
  T t2  = a + b;
  T t3  = a * b;
  T t4  = a * a;
  T t5  = a / b;
  T t7  = a * 103;
  T t8  = a << 103;
  T t9  = a >> 103;
  t = b + a;
  assert(t2 == t);
  t = b * a;
  assert(t3 == t);
  t = a * a;
  assert(t4 == t);
  t = a / b;
  assert(t5 == t);
  t = a * 103;
  assert(t7 == t);
  t = a << 103;
  assert(t8 == t);
  t = a >> 103;
  assert(t9 == t);
}

template <class T>
void identity_Natural(const T& a, const T& b, const T& c)
{
  T i,j,r,s,t;

  identity_generic2(a, b, c);

  s = (a+b) % c;
  t = ((a % c) + (b % c)) % c;
  if (s < 0) s = s+abs(c);
  if (t < 0) t = t+abs(c);
  assert(s < 0 && t > 0 && c+s == t || t < 0 && s > 0 && s == c+t || s == t);
  div(a, c, s, r); div(b, c, s, t); t = r+t;
  div(t, c, t, s); t = a+b; div(t, c, r, t);
  if (s < 0) s = s + abs(c);
  if (t < 0) t = t + abs(c);
  assert(s == t);
  s = a; s += b; s %= c;
  t = a; t %= c; r = b; r %= c; t += r; t %= c;
  if (s < 0) s += abs(c);
  if (t < 0) t += abs(c);
  assert(s == t);

  i = j = 101;
  assert(i == j && i == 101 && j.odd() && (i&15) == 5);
  i = j = i%13;
  assert(i == j && i == 10 && j.even() && (i&3) == 2);
  i = j = 0;
  assert(i == j && i == 0 && j.even() && (i&3) == 0);
  i = j = 1;
  assert(i == j && i == 1 && j.odd() && (i&121) == 1);
  i = j = 111;
  assert(i == j && i == 111 && j.odd() && (i&3) == 3);
  i = j = j/13;
  assert(i == j && i == 8 && j.even() && (i&7) == 0);

  size_t l = size_t(log2(a));
  t = a; s = 1; s <<= l;
  assert((T(1) << l) == s && s <= T(abs(t)));
  for (i = 1, j = 100; i < 100; ++i, j--) {
    r = t; s = t << 1; t <<= 1;
    assert((a << size_t(i.highest())) == s && s == t);
    s = r << 1; r = r << 1;
    assert(r == s && r == t);
    r = a << size_t(i.lowest());
    assert(r == s && r == t);
    r = T(1) << size_t(i.highest()); r *= a;
    assert(r == t);
    swap(r, t);
    ++l; s = 1; s <<= l;
    assert((T(1) << l) == s && s <= T(abs(t)));
    assert(size_t(log2(t)) == l);
    assert(i+j == 101);
  }
  assert(j == 1);

  l = size_t(log2(a));
  t = a; s = 1; s <<= l;
  assert((T(1) << l) == s && s <= T(abs(t)));
  for (i = 1, j = 100; i <= 100; i++, --j) {
    r = t; s = t >> 1; t >>= 1;
    assert((a >> size_t(i.lowest())) == s && s == t);
    s = r >> 1; r = r >> 1;
    assert(r == s);
    r = a >> size_t(i.highest());
    assert(r == s);
    r = T(1) << size_t(i.highest()); s = a/r;
    assert(s == t);
    swap(s, t);
    if (l) {
      --l; s = 1; s <<= l;
      assert((T(1) << l) == s && s <= T(abs(t)));
      assert(size_t(log2(t)) == l);
    }
    assert(i+j == 101);
  }
  assert(j == 0);

  t = s = a;
  j = INNER_ITERATION+1;
  for (i = 1; i < j; ++i) {
    assert(++t == ((s++)+1) && t == s);
    assert(s == (a+i) && s == (a+i%j));
    r = a; r += i.highest();
    assert(s == r);
  }

  for (i = 0; i < 100; i += 3) {
    t = a+(b*(i%101));
    r = a; r += b * (i%101);
    assert(r == t);
    r = b * (i%101); r += a;
    assert(r == t);
    r = b; r *= (i%101); r += a;
    assert(r == t);
  }

  for (i = 0; i < 100; ++i) {
    t = a; t.split(size_t((i%101)), s, t);
    r = a;
    r >>= BETA*size_t((i%101));
    assert(r == s && r == (a >> (BETA*size_t(i%101))));
    s = 1; s <<= BETA*size_t(i%101); --s;
    r = abs(a); r &= s;
    assert(abs(r) == abs(t) && abs(r) == (abs(a)&abs(s)));
    r = T(abs(a)) & s;
    assert(abs(r) == abs(t));
  }

  t = pow(a, 5);
  s = root(t, 5);
  assert(s == a);
  t /= a; s = sqrt(t); s = sqrt(s);
  assert(abs(s) == abs(a));
  t = a*a + 101;
  sqrt(t, r, s);
  assert(r*r+s == t);

  t = gcd(a, b);
  s = gcd(a, c);
  assert(t == gcd2(a, b));
  assert(s == gcd2(a, c));
  assert(T(gcd(b, c)) == gcd2(b, c));
  assert(t*lcm(a, b) == T(abs(T(a*b))));
  assert(lcm(c, a)*s == T(abs(T(c*a))));

  // Logic:
  assert((a & a) == a);
  s = a & a;
  assert(s == a);
  s = a; s &= a;
  assert(s == a);

  assert((a | a) == a);
  s = a | a;
  assert(s == a);
  s = a; s |= a;
  assert(s == a);

  assert((a ^ a) == 0);
  s = a ^ a;
  assert(s == 0);
  s = a; s ^= a;
  assert(s == 0);

  assert((a | (b | c)) == ((a | b) | c));
  s = b | c; s = a | s; t = a | b; t = t | c;
  assert(s == t);
  s = b; s |= c; s |= a;
  t = a; t |= b; t |= c;
  assert(s == t);

  assert((a & (b & c)) == ((a & b) & c));
  s = b & c; s = a & s; t = a & b; t = t & c;
  assert(s == t);
  s = b; s &= c; s &= a;
  t = a; t &= b; t &= c;
  assert(s == t);

  assert((a & (b | c)) == ((a & b) | (a & c)));
  s = b | c; s = a & s;
  t = a & b; r = a & c; t = t | r;
  assert(s == t);
  s = b; s |= c; s &= a;
  r = t = a; t &= b; r &= c; t |= r;
  assert(s == t);

  assert((a | (b & c)) == ((a | b) & (a | c)));
  s = b & c; s = a | s;
  t = a | b; r = a | c; t = t & r;
  assert(s == t);
  s = b; s &= c; s |= a;
  t = a; t |= b; r = a; r |= c; t &= r;
  assert(s == t);

  assert((~(~a) & ~a) == 0);
  s = ~a; s = ~s; t = ~a; s = s & t;
  assert(s == 0);
  s = ~a; s = ~s; s &= ~a;
  assert(s == 0);

  // b can't be a Digit!
  assert((~(a | b) & ~(~a & ~b)) == 0);
  s = a | b; s = ~s;
  r = ~a; t = ~b; r = r & t; r = ~r; s = s & r;
  assert(s == 0);
  s = a; s |= b; s = ~s;
  r = ~a; t = ~b; r &= t; r = ~r; s &= r;
  assert(s == 0);

  // b can't be a Digit!
  assert((~(a & b) & ~(~a | ~b)) == 0);
  s = a & b; s = ~s;
  r = ~a; t = ~b; r = r | t; r = ~r; s = s & r;
  assert(s == 0);
  s = a; s &= b; s = ~s;
  r = ~a; t = ~b; r |= t; r = ~r; s &= r;
  assert(s == 0);

  // b can't be a Digit!
  assert(((a ^ b) & ~((a & ~b) | (~a & b))) == 0);
  s = a ^ b; r = ~b; r = a & r;
  t = ~a; t = t & b; r = r | t; r = ~r; s = s & r;
  assert(s == 0);
  s = a; s ^= b; r = ~b; r &= a;
  t = ~a; t &= b; r |= t; r = ~r; s &= r;
  assert(s == 0);

  s = a;
  r = INNER_ITERATION;
  for (i = 0, j = 3; i < r; i += 3, j += 5) {
    t = s; t.setbit(size_t(j.highest()));
    s |= (T(1) << size_t(j.highest()));
    assert(s == t && s.testbit(size_t(j.highest())) == true);
    s.clearbit(size_t(i.highest()));
    assert(s.testbit(size_t(i.highest())) == false);
  }

  // constructor:
  T t6  = a % b;
  T t10 = a & b;
  T t11 = a | b;
  T t12 = a ^ b;
  T t13 = ~a;
  T t14 = sqrt(T(abs(a)));
  t = a % b;
  assert(t6 == t);
  t = b & a;
  assert(t10 == t);
  t = b | a;
  assert(t11 == t);
  t = b ^ a;
  assert(t12 == t);
  t = ~a;
  assert(t13 == t);
  t = sqrt(T(abs(a)));
  assert(t14 == t);
}

template <class T>
void identity_Natural2(const T& a, const T& b)
// a > b if a,b in Natural.
{
  T r,s,t;

  for (T i = 100; i >= 0; i -= 5) {
    t = a-(b*(i%101));
    r = a; r -= b * (i%101);
    assert(r == t);
    r = b * (i%101); r = a-r;
    assert(r == t);
    r = b; r *= (i%101); s = a; s -= r;
    assert(s == t);
    if (i == 0) break;
  }

  // constructor:
  T t2 = a - b;
  t = a - b;
  assert(t2 == t);
}

template <class T>
void identity_generic3(const T a, const T b, const T c)
// T at least Integer
{
  T r,s,t;

  assert(-(-a) == a);
  s = -a; s = -s;
  assert(s == a);
  s = a; s = -s; s = -s;
  assert(s == a);

  assert((a + (-b)) == (a - b));
  s = -b; s = a+s; t = a-b;
  assert(s == t);
  s = b; s = -s; s += a; t = a; t -= b;
  assert(s == t);
  
  assert((a * (-b)) == -(a * b));
  s = -b; s = a*s; t = a*b; t = -t;
  assert(s == t);
  s = b; s = -s; s *= a;
  t = a; t *= b; t = -t;
  assert(s == t);

  assert((a / (-b)) == -(a / b));

  assert(((a - b) + b) == a);
  s = a-b; s = s+b;
  assert(s == a);
  s = a; s -= b; s += b;
  assert(s == a);

  assert((a / (-b)) == -(a / b));
  r = b; r = -r; s = a; s /= r;
  t = a; t /= b; t = -t;
  assert(s == t);

  // constructor:
  T t2 = a - b;
  T t3 = -a;
  t = a - b;
  assert(t2 == t);
  t = -a;
  assert(t3 == t);
  t = a; t = -t;
  assert(t3 == t);

  assert((a - b) == -(b - a));
  s = a-b; t = b-a; t = -t;
  assert(s == t);
  s = a; s -= b; t = b; t -= a; t = -t;
  assert(s == t);
}

template <class T>
void io_check(const T& a)
{
  ofstream fout("t.tmp");
  fout << print(a);
  fout.close();
  ifstream fin("t.tmp");
  T b;
  assert(b.scan(fin) == true);
  assert(a == b);
  fin.close();

  ofstream fout2("t.tmp");
  fout2 << a;
  fout2.close();
  ifstream fin2("t.tmp");
  ++b;                      // b != a
  fin2 >> b;
  assert(a == b);
}

template <class T, class D>
void check_block_generic(T a, T b, T c, const D d)
// T at least Integer
{
  identity_generic(b, a);
  identity_generic(c, b);
  identity_generic(a, c);
  identity_Natural_Digit2(a, max(a, b));
  identity_Natural_Digit2(b, max(a, b));
  identity_Natural_Digit2(a, max(a, c));
  identity_Natural_Digit2(c, max(a, c));
  identity_Natural_Digit2(b, max(b, c));
  identity_Natural_Digit2(c, max(b, c));
  identity_generic2(b, c, a);
  identity_generic2(c, a, b);
  identity_generic3(a, b, c);
  identity_generic3(c, a, b);
  io_check(a);
  io_check(b);
  io_check(c);
}

template <class T, class D>
void check_generic(T a, T b, T c, const D d)
{
  cout << ';' << flush;
  check_block_generic(T(-a), b, c, d);
  cout << ';' << flush;
  check_block_generic(a, T(-b), c, d);
  cout << ';' << flush;
  check_block_generic(T(-a), T(-b), c, d);
  cout << ';' << flush;
  check_block_generic(a, b, T(-c), d);
  cout << ';' << flush;
  check_block_generic(T(-a), b, T(-c), d);
  cout << ';' << flush;
  check_block_generic(a, T(-b), T(-c), d);
  cout << ';' << flush;
  check_block_generic(T(-a), T(-b), T(-c), d);
}

template <class T, class D>
void check_block(T a, T b, T c, const D d)
{
  identity_Natural_Digit(a, b);
  identity_Natural_Digit(b, c);
  identity_Natural_Digit(c, a);
  identity_Natural_Digit2(max(a, b), a);
  identity_Natural_Digit2(max(a, b), b);
  identity_Natural_Digit2(max(a, c), a);
  identity_Natural_Digit2(max(a, c), c);
  identity_Natural_Digit2(max(b, c), b);
  identity_Natural_Digit2(max(b, c), c);
  identity_Natural2(T((a+b+c)*113), a);
  identity_Natural(a, b, c);
  identity_Natural(c, a, b);
  identity_Natural_Digit(a, d);
  identity_Natural_Digit(b, d);
  identity_Natural_Digit(c, d);
  if (a <= d) { a *= d; a += 241; }
  identity_Natural_Digit2(a, d);
  io_check(a);
  io_check(b);
  io_check(c);
}

template <class T, class D>
void check(T a, T b, T c, const D d)
{
  cout << ':' << flush;
  check_block(T(-a), b, c, d);
  cout << ':' << flush;
  check_block(a, T(-b), c, d);
  cout << ':' << flush;
  check_block(T(-a), T(-b), c, d);
  cout << ':' << flush;
  check_block(a, b, T(-c), d);
  cout << ':' << flush;
  check_block(T(-a), b, T(-c), d);
  cout << ':' << flush;
  check_block(a, T(-b), T(-c), d);
  cout << ':' << flush;
  check_block(T(-a), T(-b), T(-c), d);
}

void checkConversion(const Natural& t)
{
  char* x = new char[t.length()*BETA+1];
  for (Digit j = 2; j <= 36; ++j) {
    Ntoa(t, x, j);
    assert(atoN(x, j) == t);
  }
  delete[] x;
}

int main(int argc, char** argv)
{
  if (argc == 1) {
    ITERATION       = 20;
    MAX_FACTORIAL   = 255;
    MAX_FIBONACCI   = 1000;
    MAX_RANDOM      = 2000;
    INNER_ITERATION = 4000;
  } else if (argc != 6) {
    cout << "arg: ITERATION       (= 10)\n";
    cout << "     MAX_FACTORIAL   (= 255)\n";
    cout << "     MAX_FIBONACCI   (= 1000)\n";
    cout << "     MAX_RANDOM      (= 1000)\n";
    cout << "     INNER_ITERATION (= 2000)\n";
    return 1;
  } else {
    ITERATION       = atoi(argv[1]);
    MAX_FACTORIAL   = atoi(argv[2]);
    MAX_FIBONACCI   = atoi(argv[3]);
    MAX_RANDOM      = atoi(argv[4]);
    INNER_ITERATION = atoi(argv[5]);
    srand((unsigned)time(0));
  }

  Natural a1,a2,a3,a4;
  Natural b1,b2,b3,b4;
  Natural c1,c2,c3,c4;
  size_t i;
  cout << "Check older known bugs" << endl;
  // found by Grant L. Atoyan, September, 2011.
  a1 = atoN("9173994463960286046443283581208347763186259956673124494950355357547691504353939232280074212440502746218496");
  a2 = atoN("1461501637330902918203684832716283019655932542975");
  div(a1, a2, a3, a4);
  assert(a1 == a2*a3+a4);

  cout << "Pass 1";
  for (i = 0; i < ITERATION; ++i) {
    cout << '.' << flush;
    Digit k = rand() % MAX_FIBONACCI + 1;
    a1 = fibonacci(k);
    a2 = fibonacci2(k);
    a3 = fibonacci3(k);
    a4 = fibonacci4(k);
    assert(a1 == a2 && a2 == a3 && a3 == a4);
    checkConversion(a1);
    k = rand() % MAX_FIBONACCI + 1;
    b1 = fibonacci(k);
    b2 = fibonacci2(k);
    b3 = fibonacci3(k);
    b4 = fibonacci4(k);
    assert(b1 == b2 && b2 == b3 && b3 == b4);
    checkConversion(b1);
    k = rand() % MAX_FIBONACCI + 1;
    c1 = fibonacci(k);
    c2 = fibonacci2(k);
    c3 = fibonacci3(k);
    c4 = fibonacci4(k);
    assert(c1 == c2 && c2 == c3 && c3 == c4);
    checkConversion(c1);
    cout << ':' << flush;
    check_block(a1, b1, c1, k);
    check(Integer(a1), Integer(b1), Integer(c1), SignDigit(k));
    a2 = fibonacci(rand() % MAX_FIBONACCI + 1);
    b2 = fibonacci(rand() % MAX_FIBONACCI + 1);
    c2 = fibonacci(rand() % MAX_FIBONACCI + 1);
    check_generic(Integer(a2), Integer(b2), Integer(c2), SignDigit(k));
    check_generic(Rational(a1, a2), Rational(b1, b2), Rational(c1, c2), SignDigit(k));
  }
  cout << "\nPass 2";
  for (i = 0; i < ITERATION; ++i) {
    cout << '.' << flush;
    Digit k = rand() % MAX_FACTORIAL + 1;
    a1 = factorial(k);
    a2 = factorial2(k);
    assert(a1 == a2);
    checkConversion(a1);
    Digit k2 = rand() % MAX_FACTORIAL + 1;
    b1 = factorial(k2);
    b2 = factorial2(k2);
    assert(b1 == b2);
    checkConversion(b1);
    k = rand() % MAX_FACTORIAL + 1;
    c1 = factorial(k);
    c2 = factorial2(k);
    assert(c1 == c2);
    checkConversion(c1);
    c4 = binomial(k, k2);
    Digit k3 = (k >= k2)? k-k2 : k2-k;
    assert(c4 == (c1/(b1*factorial(k3))));
    cout << ':' << flush;
    check_block(a1, b1, c1, k);
    check(Integer(a1), Integer(b1), Integer(c1), SignDigit(k));
    a2 = factorial(rand() % MAX_FACTORIAL + 1);
    b2 = factorial(rand() % MAX_FACTORIAL + 1);
    c2 = factorial(rand() % MAX_FACTORIAL + 1);
    check_generic(Integer(a2), Integer(b2), Integer(c2), SignDigit(k));
    check_generic(Rational(a1, a2), Rational(b1, b2), Rational(c1, c2), SignDigit(k));
  }
  cout << "\nPass 3";
  for (i = 0; i < ITERATION; ++i) {
    cout << '.' << flush;
    Digit k = rand() % MAX_RANDOM + 1;
    a1.rand(size_t(k));
    assert(log2(a1) < k);
    checkConversion(a1);
    k = rand() % MAX_RANDOM + 1;
    b1.rand(size_t(k));
    assert(log2(b1) < k);
    checkConversion(b1);
    k = rand() % MAX_RANDOM + 1;
    c1.rand(size_t(k));
    assert(log2(c1) < k);
    checkConversion(c1);
    cout << ':' << flush;
    check_block(a1, b1, c1, k);
    check(Integer(a1), Integer(b1), Integer(c1), SignDigit(k));
    a2.rand(rand() % MAX_RANDOM + 1);
    b2.rand(rand() % MAX_RANDOM + 1);
    c2.rand(rand() % MAX_RANDOM + 1);
    check_generic(Integer(a2), Integer(b2), Integer(c2), SignDigit(k));
    check_generic(Rational(a1, a2), Rational(b1, b2), Rational(c1, c2), SignDigit(k));
  }
  cout << endl;
  return 0;
}
