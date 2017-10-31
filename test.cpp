/////////////////////////////////
//
// Piologie V 1.3
// multi-precision arithmetic
// simple test
//
// (c) 1996-1999 HiPiLib
//               www.hipilib.de
//
// Sebastian Wedeniwski
// 12/13/1999
//

#include <time.h>
#include "pi.h"
#include "natural.h"

#ifdef _NEEDS_PIOLOGIE_KEY_
# include "pkey.h"
#endif

#ifndef _Old_STD_
using namespace std;
#endif


int main()
{
  clock_t start,stop;
  double  duration;
  Natural a,b,c,d,e;

  // time(0) is better than clock() on some systems, e.g. PowerPC!

  start = clock();
  a = fibonacci(800000);
  stop = clock();
  duration = double(stop-start)/CLOCKS_PER_SEC;
  cout << "fib1 time [s] = " << duration << ", (" << a%ALPHA << ')' << endl;

  start = clock();
  b = fibonacci(900000);
  stop = clock();
  duration = double(stop-start)/CLOCKS_PER_SEC;
  cout << "fib2 time [s] = " << duration << ", (" << b%ALPHA << ')' << endl;

  start = clock();
  c = sqrt(a);
  stop = clock();
  duration = double(stop-start)/CLOCKS_PER_SEC;
  cout << "sqrt time [s] = " << duration << ", (" << c%ALPHA << ')' << endl;

  start = clock();
  d = a*b;
  stop = clock();
  duration = double(stop-start)/CLOCKS_PER_SEC;
  cout << "mul time [s] = " << duration << ", (" << d%ALPHA << ')' << endl;

  start = clock();
  e = d*d;
  stop = clock();
  duration = double(stop-start)/CLOCKS_PER_SEC;
  cout << "sqr time [s] = " << duration << ", (" << e%ALPHA << ')' << endl;

  start = clock();
  b /= a;
  stop = clock();
  duration = double(stop-start)/CLOCKS_PER_SEC;
  cout << "div time [s] = " << duration << ", (" << b%ALPHA << ')' << endl;

  start = clock();
  Pi pi(20000);
  stop = clock();
  duration = double(stop-start)/CLOCKS_PER_SEC;
  cout << "pi time [s] = " << duration << endl;

  return 0;
}
