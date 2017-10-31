/////////////////////////////////
//
// Piologie V 1.3
// multi-precision arithmetic
// Compositeness Test
//
// (c) 1996-2002 HiPiLib
//               www.hipilib.de
//
// Sebastian Wedeniwski
// 02/03/2002
//

#include "nmbrthry.h"
#include "ispprime.h"


bool ispprime(const Natural& n)
{
  if ((n&1) == 0) return (n == 2);

  Natural x,x2;
  sqrt(n, x, x2);
  if (x2 == 0) return false;

  Natural m = n-1;
  Natural l = n >> 1;
  Natural l2 = n+1;
  Natural k = 1;
  x2 = 4;
  Natural x3 = Digit(0);
  while (true) {
    while (true) {
      x2 += k; k += 2;
      if (x2 >= n) x2 -= n;
      const int i = jacobi(x2, n);
      if (i == -1) break;
      if (i == 0) return false;
    }
    if (pow(x2, l, n) != m) return false;
    x2 -= 2;
    if (pow(x2, m, n) != 1) return false;
    x2 += 2;
    x = k >> 1;
    Natural u0 = 0;         // U_0(x, -1) := 0
    Natural u1 = 1;         // U_1(x, -1) := 1
    Digit m2 = log2(l2)-1;
    do {
      // U_{2m-1}(x, -1) = U_m(x, -1)^2 + U_{m-1}(x, -1)^2
      // U_{2m}(x, -1)   = xU_m(x, -1)^2 + 2U_m(x, -1)U_{m-1}(x, -1)
      Natural h0 = u0*u0;
      Natural h1 = u1*u1;
      Natural h3 = u0*u1;
      h3 <<= 1;
      u0 = h0+h1; u0 %= n;
      h1 *= x;
      u1 = h1+h3; u1 %= n;
      if (l2.testbit(m2)) {
        // U_{m+2}(x, -1) = xU_{m+1}(x, -1) + U_m(x, -1)
        h0 = x*u1; h0 += u0;
        u0 = u1;
        u1 = h0%n;
      }
    } while (m2-- > 0);
    if (u1 != 0 || u0 != m) return false;

    if (x3 != 0) break;
    x3 = x;
  }
  k >>= 1;
  return (gcd(x3+k, n) == 1 && gcd(abs(x3, k), n) == 1);
}
