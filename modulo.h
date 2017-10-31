///////////////////////////////
//
// modulo arithmetic for Piologie
// multi-precision arithmetic
//
// Sebastian Wedeniwski
// 01/02/1999
//


#ifndef _Include_Modulo_H_
#define _Include_Modulo_H_

template<class T>
T pow(T a, T b, const T& c)
// Algorithm:  d := pow(a, b, c)
// Input:      a,b,c in T where c > 0.
// Output:     d in T such that d == a^b (mod c) ||
{
  a %= c;
  if (b == 1) return a%c;
  else if (b == 0) return 1;
  else if (a == 0) return Digit(0);

  while ((b&1) == 0) { a *= a; a %= c; b >>= 1; }

  T z = a;
  while (--b != 0) {
    while ((b&1) == 0) { a *= a; a %= c; b >>= 1; }
    z *= a; z %= c;
  }
  return z;
}

template<class T>
T inverse(const T& x, const T& m)
// Algorithm:  c := inverse(a, b)
// Input:      a,b in T where b > 0.
// Output:     c in T such that c == a^(-1) (mod b) ||
{
  T a = x;
  T b = m;
  T q = 1;
  T p = 0;
  T s,t;
  do {
    t = a/b;
    a -= s = t*b;
    q += s = t*p;
    if (a == 0) return m - p;
    t = s = b/a;
    b -= s = t*a;
    p += s = t*q;
  } while (b >= 1);
  return q;
}
#include <iostream>
template<class T>
T sqrt(const T& a, const T& b)
// Algorithm:  c := sqrt(a, b)
// Input:      a,b in T where b is prime and (a|b) = 1.
// Output:     c in T such that c^2 == a (mod b) ||
{
  if (a == 0) return a;
  if ((b&3) == 3) {
    T c = b+1;
    c >>= 2;
    return pow(a, c, b);
  } else if ((b&7) == 5) {
    const T c = a << 1;
    const T v = pow<T>(c, (b - 5) >> 3, b);
    T i = c*(v*v);
    --i; i *= a*v;
    return i % b;
  }

  Digit k = 1;
  T d = b-a;
  ++d;
  while (true) {
    const int x = jacobi(d, b);
    if (x == 0 && d != b) return 1;
    else if (x < 0) break;
    d += 2*k+1; ++k;
  }

  T b3 = k;
  d %= b;
  // (k + sqrt(d))^([b/2]+1):
  T c = b >> 1; ++c;
  T t,t2,t3,b2 = 1;
  T x = 1;
  T x2 = 0;

  while (c != 0) {
    while ((c&1) == 0) {
      t = b3*b3; t2 = b2*b2; t += d*t2;
      t2 = b3*b2; t2 <<= 1;
      b3 = t % b; b2 = t2 % b;
      // (b3, b2) := (b3*b3 + d*b2*b2, 2*b2*b3) = (b3, b2)^2
      c >>= 1;
    }
    t = b3+b2; t2 = x+x2; t3 = b2*x2;
    x2 = t*t2; x2 -= t = b3*x; x2 -= t3;
    x2 %= b;
    t += d*t3;
    x = t % b;
    // (x, x2) := (b3*x+d*b2*x2, (b3+b2)*(x+x2)-b3*x-b2*x2) = (b3, b2)*(x, x2)
    --c;
  }
  return x;
}

template<class T>
T chinese(const T& m1, const T& m2, const T& a1, const T& a2)
// Algorithm:  x := chinese(m1, m2, a1, a2)
// Input:      m1,m2,a1,a2 in T where 0 <= a1 < m1, 0 <= a2 < m2.
// Output:     x in T such that x == a1 (mod m1) and x == a2 (mod m2)
//             where 0 <= x < m1*m2 ||
{
  const T d = gcd(m1, m2);
  if (d == 1) {
    const T c = inverse(m1, m2);
    const T m = m1*m2;
    T t;
    if (a2 >= a1) t = a2-a1;
    else { t = m-a1; t += a2; }
    t *= c; t *= m1; t += a1;
    return t %= m;
  }
  T m = m1/d;
  const T b = m2/d;
  const T c = inverse(m, b);
  m *= m2;
  T t = a2%d;
  if (t != a1%d) return T();
  if (a2 >= a1) { t = a2-a1; t /= d; }
  else { t = a1-a2; t /= d; t %= b; t = b-t; }
  t *= c; t *= m1; t += a1;
  return t %= m;
}

template<class Container>
void chinese(const typename Container::value_type& m1,
             const typename Container::value_type& m2,
             const Container& l1, const Container& l2,
             Container& r, typename Container::value_type& m)
// Algorithm:  chinese(m1, m2, l1, l2, r, m)
// Input:      m1,m2 in T, l1,l2 are containers over T
//             where 0 <= a1 < m1, 0 <= a2 < m2.
{
  const typename Container::value_type d = gcd(m1, m2);
  const typename Container::value_type b = m2/d;
  const typename Container::value_type b2 = m1/d;
  const typename Container::value_type c = inverse(b2, b);
  m = b2*m2;
  typename Container::value_type s,t;

  Container* r2 = &r;
  if (r2 == &l1 || r2 == &l2) r2 = new Container();
  r2->clear();
  const typename Container::const_iterator j1 = l1.end();
  for (typename Container::const_iterator i1 = l1.begin(); i1 != j1; ++i1) {
    s = *i1 % d;
    const typename Container::const_iterator j2 = l2.end();
    for (typename Container::const_iterator i2 = l2.begin(); i2 != j2; ++i2)
      if (*i2 % d == s) {       // condition: x1 == x2 (mod gcd(m1, m2))
        if (*i2 >= *i1) { t = *i2 - *i1; t /= d; }
        else { t = *i1 - *i2; t /= d; t %= b; t = b - t; }
        t *= c; t *= m1; t += *i1; t %= m;
        r2->push_back(t);
      }
  }
  if (&r != r2) { r = *r2; delete r2; }
}


#endif
