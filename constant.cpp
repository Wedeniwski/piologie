///////////////////////////////
//
// Piologie V1.3
//
// (c) Sebastian Wedeniwski
// 08/13/1999
//
// Note: Timings on Pentium II 266 MHz, Windows NT 4.0, MS Visual C++ 6.0
// constant 1 131072
// Computing time [s] = 17.955
// 
// constant 1 1048576
// official: 270 sec
// Computing time [s] = 239.264


#include <time.h>
#include "pi.h"

#ifndef _Old_STD_
# include <fstream>
using namespace std;
#else
# include <fstream.h>
#endif

void calcPi(Digit n)
{
  clock_t start,stop;
  double  duration;

  start = clock();
  Pi pi(n);
  stop = clock();
  duration = double(stop-start)/CLOCKS_PER_SEC;
  cout << "Computing time [s] = " << duration << endl;

  cout << "Writing result to file pi.txt ..." << endl;
  ofstream fout("pi.txt");
  fout << pi;
  stop = clock();
}

void calcExp1(Digit n)
{
  clock_t start,stop;
  double  duration;

  start = clock();
  Exp1 exp(n);
  stop = clock();
  duration = double(stop-start)/CLOCKS_PER_SEC;
  cout << "Computing time [s] = " << duration << endl;

  cout << "Writing result to file exp1.txt ..." << endl;
  ofstream fout("exp1.txt");
  fout << exp;
  stop = clock();
}

void calcZeta3(Digit n)
{
  clock_t start,stop;
  double  duration;

  start = clock();
  Zeta3 zeta(n);
  stop = clock();
  duration = double(stop-start)/CLOCKS_PER_SEC;
  cout << "Computing time [s] = " << duration << endl;

  cout << "Writing result to file zeta3.txt ..." << endl;
  ofstream fout("zeta3.txt");
  fout << zeta;
  stop = clock();
}

void calcGamma(Digit n)
{
  clock_t start,stop;
  double  duration;

  start = clock();
  EulerGamma gamma(n);
  stop = clock();
  duration = double(stop-start)/CLOCKS_PER_SEC;
  cout << "Computing time [s] = " << duration << endl;

  cout << "Writing result to file gamma.txt ..." << endl;
  ofstream fout("gamma.txt");
  fout << gamma;
  stop = clock();
}

void calcLn2(Digit n)
{
  clock_t start,stop;
  double  duration;

  start = clock();
  Ln ln(2, n);
  stop = clock();
  duration = double(stop-start)/CLOCKS_PER_SEC;
  cout << "Computing time [s] = " << duration << endl;

  cout << "Writing result to file ln2.txt ..." << endl;
  ofstream fout("ln2.txt");
  fout << ln;
  stop = clock();
}


int main(int argc, char** argv)
{
  if (argc != 3) {
    cerr << "USAGE: " << argv[0]
         << " <constant> <digits>\n1: pi\n2: exp(1)\n3: zeta(3)\n4: gamma\n5: ln(2)"
         << endl;
    return 1;
  }
  switch (atoi(argv[1])) {
    case 2:   calcExp1(atoi(argv[2])); break;
    case 3:   calcZeta3(atoi(argv[2])); break;
    case 4:   calcGamma(atoi(argv[2])); break;
    case 5:   calcLn2(atoi(argv[2])); break;
    default:  calcPi(atoi(argv[2]));
  }

  return 0;
}
