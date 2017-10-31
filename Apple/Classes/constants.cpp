/*
 *  constants.cpp
 *  Piologie
 *
 *  Created by Sebastian Wedeniwski on 08.08.10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <string>
#include <sstream>
#include <ctime>

#include "pi.h"

extern "C" {
#include "constants.h"
}

using namespace std;


static double tmpCalcDuration = 0.0;
static double tmpConversionDuration = 0.0;
static string strConstantsNumner = "";
static string strFormatted = "";

extern "C"
const char* outputCalcUnformattedNumber() {
  return strConstantsNumner.c_str();
}

extern "C"
const char* outputCalcNumber(int pos, int digitsPerLine, int digitsPerView) {
  char d[10];
  const char* c = strConstantsNumner.c_str();
  for (strFormatted = "        "; *c != '.'; ++c) strFormatted += *c;
  strFormatted += ".\n";
  int length = strlen(++c);
  if (pos > 0 && length > pos) {
    strFormatted = "";
    c += pos;
    length -= pos;
  }
  bool sgn = (pos >= 100000);
  int pos2 = pos%100000;
  for (int i = 1; length > 0 && i <= digitsPerView; i += digitsPerLine) {
    int k = pos2+i;
    if (k >= 100000) {
      k -= 100000;
      sgn = true;
    }
    if (sgn) sprintf(d, "\'%05d:", k);
    else sprintf(d, "%6d:", k);
    strFormatted += d;
    for (int j = 0; length > 0 && j < digitsPerLine; ++j) {
      if (j+10 <= length) {
        strFormatted += ' ';
        strFormatted += c[0];
        strFormatted += c[1];
        strFormatted += c[2];
        strFormatted += c[3];
        strFormatted += c[4];
        strFormatted += c[5];
        strFormatted += c[6];
        strFormatted += c[7];
        strFormatted += c[8];
        strFormatted += c[9];
        j += 9; c += 10; length -= 10;
      } else {
        if (j%10 == 0) strFormatted += ' ';
        strFormatted += *c;
        ++c; --length;
      }
    }
    strFormatted += '\n';
  }
  return strFormatted.c_str();
}

extern "C"
double outputCalcDuration() {
  return tmpCalcDuration;
}

extern "C"
double outputConversionDuration() {
  return tmpConversionDuration;
}
  
extern "C"
void calcPi(unsigned int n) {
  strConstantsNumner = strFormatted = "";
  clock_t start = clock();
  Pi pi(n);
  clock_t stop = clock();
  tmpCalcDuration = double(stop-start)/CLOCKS_PER_SEC;

  string t = "";
  stringstream s(t);
  start = clock();
  s << pi;
  strConstantsNumner = s.str();
  stop = clock();
  tmpConversionDuration = double(stop-start)/CLOCKS_PER_SEC;
}

extern "C"
void calcSqrt2(unsigned int n) {
  strConstantsNumner = strFormatted = "";
  clock_t start = clock();
  Sqrt sqrt(2, n);
  clock_t stop = clock();
  tmpCalcDuration = double(stop-start)/CLOCKS_PER_SEC;
  
  string t = "";
  stringstream s(t);
  start = clock();
  s << sqrt;
  strConstantsNumner = s.str();
  stop = clock();
  tmpConversionDuration = double(stop-start)/CLOCKS_PER_SEC;
}

extern "C"
void calcExp1(unsigned int n) {
  strConstantsNumner = strFormatted = "";
  clock_t start = clock();
  Exp1 exp(n);
  clock_t stop = clock();
  tmpCalcDuration = double(stop-start)/CLOCKS_PER_SEC;
  
  string t = "";
  stringstream s(t);
  start = clock();
  s << exp;
  strConstantsNumner = s.str();
  stop = clock();
  tmpConversionDuration = double(stop-start)/CLOCKS_PER_SEC;
}

extern "C"
void calcZeta3(unsigned int n) {
  strConstantsNumner = strFormatted = "";
  clock_t start = clock();
  Zeta3 zeta(n);
  clock_t stop = clock();
  tmpCalcDuration = double(stop-start)/CLOCKS_PER_SEC;
  
  string t = "";
  stringstream s(t);
  start = clock();
  s << zeta;
  strConstantsNumner = s.str();
  stop = clock();
  tmpConversionDuration = double(stop-start)/CLOCKS_PER_SEC;
}

extern "C"
void calcGamma(unsigned int n) {
  strConstantsNumner = strFormatted = "";
  clock_t start = clock();
  EulerGamma gamma(n);
  clock_t stop = clock();
  tmpCalcDuration = double(stop-start)/CLOCKS_PER_SEC;
  
  string t = "";
  stringstream s(t);
  start = clock();
  s << gamma;
  strConstantsNumner = s.str();
  stop = clock();
  tmpConversionDuration = double(stop-start)/CLOCKS_PER_SEC;
}

extern "C"
void calcLn2(unsigned int n) {
  strConstantsNumner = strFormatted = "";
  clock_t start = clock();
  Ln ln(2, n);
  clock_t stop = clock();
  tmpCalcDuration = double(stop-start)/CLOCKS_PER_SEC;
  
  string t = "";
  stringstream s(t);
  start = clock();
  s << ln;
  strConstantsNumner = s.str();
  stop = clock();
  tmpConversionDuration = double(stop-start)/CLOCKS_PER_SEC;
}
