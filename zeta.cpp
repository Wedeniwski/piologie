/////////////////////////////////
//
// Piologie V 1.3
// multi-precision arithmetic
// Program for Zeta(3) record
//
// (c) 1996-1999 HiPiLib
//               www.hipilib.de
//
// Sebastian Wedeniwski
// 12/13/1999
//

#include "pi.h"
#include "integer.h"

#ifdef _Old_STD_
# include <fstream.h>
#else
# include <fstream>
using namespace std;
#endif

#include <string.h>
#include <time.h>



//#define _Old_Formula_         // Use old formula
//#define _FAST_OUTPUT_         // hex-output




#define LOG(a)  { cout << a << flush; \
                { ofstream lout("../zeta.log", ios::app); \
                  lout << a; } }
#define LOGLN(a)  { cout << a << endl; \
                  { ofstream lout("../zeta.log", ios::app); \
                    lout << a << endl; } }


#ifdef _Old_Formula_
const double ITERATIONS = 3.01;
#else
const double ITERATIONS = 5.03;
#endif


struct Zeta {
  Fixed zeta;

  void series(Digit, const Digit, Integer&, Integer&, Integer&) const;

  Zeta(const size_t, int);
  Zeta(const size_t, const size_t, const size_t, const Digit);

  void fast_output(OSTREAM&);
  void fast_input(ISTREAM&);

  const Fixed& value() const;
};

inline const Fixed& Zeta::value() const
{
  return zeta;
}

inline OSTREAM& operator<<(OSTREAM& out, const Zeta& a)
{
  return out << a.zeta;
}


void Zeta::series(Digit n, const Digit m, Integer& p, Integer& q, Integer& t) const
// n < 2^BETA/72-24
{
  // CONDITION(n < m);
  if (n+1 == m) {

#ifdef _Old_Formula_

    if (n == 0) { p = 1; q = 32; t = 77; }
    else {
      if (n <= GAMMA_LOW/2) { p = n*n; p *= n*n; }
      else { p = n; p = p*p; p = p*p; }
      p *= n; p = -p;
      if (n < GAMMA/414) t = 205*n+250;
      else { t = n; t *= 205; t += 250; }
      t *= n; t += 77; t *= p;
      n = 2*n+1;
      if (n <= GAMMA_LOW/2) { q = n*n; q *= n*n; }
      else { q = n; q = q*q; q = q*q; }
      q *= 32*n;
    }

#else

    if (n == 0) { p = 1; q = 10368; t = 12463; }
    else {
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
    }

#endif

  } else {
    Integer p1,q1,t1;
    Integer p2,q2,t2;

    const Digit l = (n+m)/2;
    series(n, l, p1, q1, t1);
    series(l, m, p2, q2, t2);
    t = t1*q2; t += t2*p1;
    p = p1*p2; q = q1*q2;

    const Digit k = m-n;
    if (k > 8 && k <= 16) {
      p1 = gcd(p, gcd(t, q));
      p /= p1; q /= p1; t /= p1;
    }
  }
}

static void format(bool point)
{
  ifstream fin("zeta3.tmp");
  ofstream fout("zeta3-f.tmp");
  int i = 0;
  char c;
  if (point)
    while (!fin.eof()) {
      if (++i == 2) { fout << '.'; continue; }
      if (!fin.get(c)) break;
      fout.put(c);
      if (i == 78) { fout << "\\\n"; i = 0; break; }
    }
  while (!fin.eof()) {
    if (!fin.get(c)) break;
    fout.put(c);
    if (++i == 78) { fout << "\\\n"; i = 0; }
  }
}

static void reduction(Integer& p, Integer& q, Integer& t)
{
  size_t i = 0;
  if (p != 0)
    while (!abs(p).testbit(i)) ++i;
  size_t j = 0;
  if (q != 0)
    while (!abs(q).testbit(j)) ++j;
  size_t k = 0;
  if (t != 0)
    while (!abs(t).testbit(k)) ++k;
  if (i > j) i = j;
  if (i > k) i = k;
  p >>= i; q >>= i; t >>= i;
}
 
Zeta::Zeta(const size_t n, int idx)
 : zeta(n)
{
  if (idx == 0) {
    ifstream fin("zeta.tmp");
    size_t x;
    fin >> x;
    if (x != zeta.precision()) { LOG("Note (different BETA)!\n"); }
    char c;
    fin.get(c);
    fin >> x;
    if (x != zeta.decimals()) { LOG("Error!\n"); return; }
    fin.get(c);
    if (!zeta.value().scan(fin)) { LOG("Error!\n"); return; }

//    zeta.value() >>= 64;        // BETA: 64 -> 32

#ifdef _Old_Formula_
    zeta.value() >>= 1;
#endif

  } else if (idx == 1) {
    ifstream fin("zeta3.tmp2");
    fast_input(fin);
  } else LOGLN("ERROR!");
}

Zeta::Zeta(const size_t n, const size_t idx, const size_t max, const Digit st)
 : zeta(n)
{
  time_t start,stop;
  const size_t sz = zeta.precision();

  Integer p,q,t;
  if (idx < max) {
    const Digit m = Digit(n/ITERATIONS);
    const Digit k = m/max;
    const Digit i = k*idx;

    start = time(0);
    if (i+2*k > m) {
      LOGLN("process:" << idx << ", from " << ((st)? st : i) << " to " << m);
      series((st)? st : i, m, p, q, t);
    } else {
      LOGLN("process:" << idx << ", from " << ((st)? st : i) << " to " << i+k);
      series((st)? st : i, i+k, p, q, t);
    }
    stop = time(0);
    LOGLN("zeta time [s] = " << difftime(stop, start));

    reduction(p, q, t);

    ofstream foutp("zeta-p.tmp");
    foutp << print(p);
    foutp.close();

    ofstream foutq("zeta-q.tmp");
    foutq << print(q);
    foutq.close();

    ofstream foutt("zeta-t.tmp");
    foutt << print(t);
    foutt.close();
  } else {
    switch (st) {
    case 0: {
      ifstream finq("zeta-q.tmp");
      ifstream fint("zeta-t.tmp");
      q.scan(finq);
      t.scan(fint);

      const bool sign = (q <= 0 && t <= 0 || q >= 0 && t >= 0);

      LOGLN("process division");
      start = time(0);
      size_t i = 0;
      while (!abs(q).testbit(i)) ++i;
      q >>= i; t <<= BETA*sz-i;
      zeta.value() = abs(t) / abs(q);
      stop = time(0);
      LOGLN("zeta time [s] = " << difftime(stop, start));
      ofstream fout("zeta.tmp0");
      fout << zeta.precision() << ',' << zeta.decimals() << ',';
      if (!sign) fout << '-';
      fout << print(zeta.value());
      break;
    }
    case 1: {
      ifstream finq("zeta-q.tmp");
      ifstream finp("zeta-p.tmp");
      q.scan(finq);
      p.scan(finp);

      const bool sign = (q <= 0 && p <= 0 || q >= 0 && p >= 0);

      LOGLN("process division");
      start = time(0);
      size_t i = 0;
      while (!abs(q).testbit(i)) ++i;
      q >>= i; p <<= BETA*sz-i;
      zeta.value() = abs(p) / abs(q);
      stop = time(0);
      LOGLN("zeta time [s] = " << difftime(stop, start));
      ofstream fout("zeta.tmp1");
      fout << zeta.precision() << ',' << zeta.decimals() << ',';
      if (!sign) fout << '-';
      fout << print(zeta.value());
      break;
    }
    default: LOG("ERROR!\n");
    }
  }
}

void Zeta::fast_output(OSTREAM& fout)
{
  const size_t SIZE = zeta.decimals();
  size_t m = zeta.precision();
  Natural c = pow(Natural(10), zeta.decimals());
  zeta.value() *= c;
  zeta.value().rmove(m);
  char* result = new char[SIZE+1];
  Ntoa(zeta.value(), result, 16);
  fout << result;
  delete[] result;
}

void Zeta::fast_input(ISTREAM& fin)
{
  const size_t SIZE = zeta.decimals();
  size_t m = zeta.precision();
  char* a = new char[SIZE+1];
  char* b = a;

  while (!fin.eof()) {
    char c;
    if (!fin.get(c) || c == '\n') break;
    *b++ = c;
  }
  *b = 0;
  zeta.value() = atoN(a, 16);
  delete[] a;
}

static size_t read(const char* path, const char* name,
                   Integer& q, const bool prec, size_t& decimals)
{
  char str[100];
  ifstream fin(strcat(strcpy(str, path), name));
  LOG(name << " read...");

  size_t precision = 0;
  if (prec) {
    char c;
    fin >> precision;
    fin.get(c);
    fin >> decimals;
    fin.get(c);
  }
  if (q.scan(fin)) { LOGLN("...ok!"); }
  else { LOG("...ERROR!\n"); }
  return precision;
}

static void combine(const char* path1, const char* path2, int idx)
{
  time_t start,stop;

  switch (idx) {
  case 1: {
    Integer q,t,p1,t1,z;
    size_t x;
    read(path1, "zeta-p.tmp", p1, false, x);
    read(path1, "zeta-t.tmp", t1, false, x);
    read(path2, "zeta-q.tmp", q, false, x);
    read(path2, "zeta-t.tmp", z, false, x);

    LOGLN("combine1");
    start = time(0);
    t = z*p1;
    LOGLN("combine2");
    t += z = t1*q;
    stop = time(0);
    LOGLN("zeta time [s] = " << difftime(stop, start));

    ofstream fout("zeta-t.tmp");
    fout << print(t);
    break;
  }
  case 2: {
    Integer p,p1;
    size_t x;
    read(path1, "zeta-p.tmp", p1, false, x);
    read(path2, "zeta-p.tmp", p, false, x);

    LOGLN("combine" << idx);
    start = time(0);
    p *= p1;
    stop = time(0);
    LOGLN("zeta time [s] = " << difftime(stop, start));

    ofstream fout("zeta-p.tmp");
    fout << print(p);
    break;
  }
  case 3: {
    Integer q,q1;
    size_t x;
    read(path1, "zeta-q.tmp", q1, false, x);
    read(path2, "zeta-q.tmp", q, false, x);

    start = time(0);
    LOGLN("combine" << idx);
    q *= q1;
    stop = time(0);
    LOGLN("zeta time [s] = " << difftime(stop, start));

    ofstream foutq("zeta-q.tmp");
    foutq << print(q);
    break;
  }
  case 4: {
    Integer pq,tq;
    size_t decimals;
    size_t precision = read(path1, "zeta.tmp1", pq, true, decimals);
    if (precision != read(path2, "", tq, true, decimals)) {
      LOGLN("precision error!");
      return;
    }

    LOGLN("combine" << idx);
    start = time(0);
    pq *= tq;
    pq.rmove(precision);
    stop = time(0);
    LOGLN("zeta time [s] = " << difftime(stop, start));

    ofstream foutq("zeta.tmp2");
    foutq << precision << ',' << decimals << ',' << print(pq);
    break;
  }
  case 5: {
    Integer s1,s2;
    size_t decimals;
    size_t precision = read(path1, "zeta.tmp0", s1, true, decimals);
    if (precision != read(path2, "zeta.tmp2", s2, true, decimals)) {
      LOGLN("precision error!");
      return;
    }

    LOGLN("combine" << idx);
    start = time(0);
    s1 += s2;
    stop = time(0);
    LOGLN("zeta time [s] = " << difftime(stop, start));

    ofstream foutq("zeta.tmp");
    foutq << precision << ',' << decimals << ',' << print(s1);
    break;
  }
  default: LOG("ERROR!\n");
  }
}

static void pre_output(ostream& out, const Natural& a, size_t n)
{
  cout << "level " << n << endl;
  if (n <= 1000) {
    out.width(n);
    out.fill('0');
    out << a;
  } else {
    const size_t m = n/2;
    Natural q;
    cout << "pow" << endl;
    Natural r = pow(Natural(10), m);
    cout << "div" << endl;
    div(a, r, q, r);
    pre_output(out, q, m + (n&1));
    pre_output(out, r, m);
  }
}


int main(int argc, char** argv)
{
  size_t i,j,n;
  Digit  m = 0;

  for (i = 0; i < argc; ++i) LOG(argv[i] << ' ');

  if (argc == 4) {
    n = atoi(argv[1]);
    if (n >= 1 && n <= 5) {
      combine(argv[2], argv[3], n);
      return 0;
    }
    i = atoi(argv[2]);
    j = atoi(argv[3]);
  } else if (argc == 5) {
    n = atoi(argv[1]);
    i = atoi(argv[2]);
    j = atoi(argv[3]);
    m = atoi(argv[4]);
  } else if (argc == 3) {
    i = atoi(argv[2]);
    Zeta z(atoi(argv[1]), i);
    time_t start,stop;
    if (i == 0) {
      ofstream fout("zeta3.tmp2");
      start = time(0);
      z.fast_output(fout);
      stop = time(0);
    } else if (i == 1) {
      ofstream fout("zeta3.tmp");
      start = time(0);
      pre_output(fout, z.value().value(), z.value().decimals());
      stop = time(0);
      fout.close();
      format(true);
    } else LOGLN("ERROR!\n");
    LOGLN("zeta time [s] = " << difftime(stop, start));
    return 0;
  } else if (argc == 2) {
    Zeta z(atoi(argv[1]), 0);
    time_t start,stop;
    ofstream fout("zeta3.tmp");
    start = time(0);
    fout << z << endl;
    stop = time(0);
    LOGLN("zeta time [s] = " << difftime(stop, start));
    fout.close();
    format(false);
    return 0;
  } else {
    cout << "useage:  zeta <decimals> <step> <parts> [start]\n";
    cout << "         zeta <decimals> <step> <parts> <idx>\n";
    cout << "         zeta <decimals> <idx>\n";
    cout << "         zeta <decimals>\n";
    cout << "      or zeta <idx> <path1> <path2>\n";
    return 1;
  }

  Zeta z(n, i, j, m);

  return 0;
}
