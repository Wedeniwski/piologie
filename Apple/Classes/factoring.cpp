/*
 *  factoring.cpp
 *  Piologie
 *
 *  Created by Sebastian Wedeniwski on 10.08.10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <sstream>
#include <ctime>
#include <list>

#include "nmbrthry.h"
#include "ispprime.h"
#include "natural.h"

extern "C" {
#include "factoring.h"
}

using namespace std;

static double tmpFactoringDuration = 0.0;
static string strFactoring = "";
static bool isFactoringCompleteFlag = true;


extern "C"
double outputFactoringDuration() {
  return tmpFactoringDuration;
}


void insert(const Natural& a, list<Natural>& p)
{
  list<Natural>::iterator i = p.end();
  if (i == p.begin()) {
    p.push_back(a);
  } else {
    do {
      if (*--i <= a) {
        p.insert(++i, a);
        return;
      }
    } while (i != p.begin());
    p.insert(i, a);
  }
}

void insert(const list<Natural>& a, list<Natural>& p)
{
  for (list<Natural>::const_iterator i = a.begin(); i != a.end(); ++i) {
    insert(*i, p);
  }
}

Natural pollard(const Natural& n)
{
  Natural a = 3;
  Natural x = 17;
  Natural y = x;
  Natural q = 1;
  for (Digit i = 1; i < 10000; ++i) {
    x = abs(x*x, a); x %= n;
    y = abs(y*y, a); y %= n; y = abs(y*y, a); y %= n;
    q *= abs(x, y);
    if ((i&31) == 0) {
      Natural g = gcd(q, n);
      if (g > 1) return g;
    }
  }
  return 1;
}

#ifndef _Old_STD_
# define NOTHROW_NEW  new(nothrow)
#else
# define NOTHROW_NEW  new
#endif


Natural cfrac(const Natural& a, const Digit m, const Digit c)
// Algorithm:  q := cfrac(a, b, c)
// Input:      a in Natural, b,c in Digit where b >= 50, c >= 1.
// Output:     b in Natural such that q | a ||
{
  struct prim_t {
    Natural* x;
    char*    c;
  };
  Natural s,s2,s3,r,v,y;
  Natural t = 1;
  const Digit pm = 3;
  Digit* p = NOTHROW_NEW Digit[m];
  if (!p) a.errmsg(2, "(cfrac)");
  p[0] = p[1] = 0;
  
  prim_t* d = NOTHROW_NEW prim_t[m+1];
  if (!d) a.errmsg(2, "(cfrac)");
  Digit i;
  for (i = 0; i <= m; ++i) {
    d[i].x = NOTHROW_NEW Natural[3];
    if (!d[i].x) a.errmsg(2, "(cfrac)");
    d[i].c = NOTHROW_NEW char[m];
    if (!d[i].c) a.errmsg(2, "(cfrac)");
    memset(d[i].c, 0, m);
  }
  Primes primes;
  for (Natural b = c*a; t == 1; b += a) {
    sqrt(b, r, v);
    Natural w = r;
    Natural z = 1;
    Natural u = 1;
    r <<= 1;
    Natural x = r;
    Digit prim = primes.firstPrime();
    Digit n = 2;
    while (primes.nextPrime(prim) && n < m) {
      Digit q = b % prim;
      Digit q2 = 1;
      for (i = prim/2; i > 1; i /= 2) {
        if (i&1) { q2 *= q; q2 %= prim; }
        q *= q; q %= prim;
      }
      q2 *= q; q2 %= prim;
      if (q2 <= 1) p[n++] = prim;
    }
    Natural trial = 1;
    for (i = 2; i < pm; ++i);
    while (i < n) trial *= p[i++];
    
    Digit nd = 0;
    Digit iteration = 0;
    while (true) {
      ++iteration;
      t = v;
      if (t == 1) break;
      while ((t&3) == 0) t >>= 2;
      d[nd].c[1] = char(t.even());
      if (d[nd].c[1]) t >>= 1;
      for (i = 2; i < pm && t > 1; ++i) {
        const Digit k = p[i]*p[i];
        const Digit l = k*p[i];
        Digit m;
        while (true) {
          div(t, k*k, s, m);
          if (m) break;
          swap(t, s);
        }
        if (m%l == 0) {
          d[nd].c[i] = 1;
          t /= l;
        } else if (m%k == 0) {
          d[nd].c[i] = 0;
          t /= k;
        } else if (m%p[i] == 0) {
          d[nd].c[i] = 1;
          t /= p[i];
        } else d[nd].c[i] = 0;
      }
      if (t > 1) {
        div(trial, t, s, s2);
        if (s2 != 0) {
          s2 *= s2; s2 %= t;
          if (s2 != 0) { s2 *= s2; s2 %= t; }
        }
      }
      if (t == 1 || s2 == 0) {
        if (t > 1) {
          s = t;
          Digit j,k;
          do {
            k = p[i]*p[i];
            div(s, k, t, j);
            if (j == 0) {
              swap(s, t);
              --i;
            } else if (j%p[i] == 0) {
              d[nd].c[i] = 1;
              s /= p[i];
            }
            ++i;
          } while (s > k);
          if (s != 1) {
            while (s > p[i]) ++i;
            d[nd].c[i] = 1;
          }
        }
        d[nd].c[0] = char(iteration&1);
        d[nd].x[0] = w; d[nd].x[1] = v; d[nd].x[2] = 1;
        Digit j;
        for (j = 0; j < n && d[nd].c[j] == 0; ++j);
        Digit k = 0;
        for (i = 0; i < nd && j < n; ++i) {
          while (d[i].c[k] == 0) ++k;
          if (k == j) {
            d[nd].c[k] = 0;
            for (Digit l = ++k; l < n; ++l) d[nd].c[l] ^= d[i].c[l];
            d[nd].x[0] *= d[i].x[0];
            d[nd].x[0] %= a;
            t = gcd(d[nd].x[1], d[i].x[1]);
            d[nd].x[1] *= d[i].x[1];
            d[nd].x[1] /= t*t;
            d[nd].x[2] *= t;
            d[nd].x[2] *= d[i].x[2];
            d[nd].x[2] %= a;
            while (j < n && d[nd].c[j] == 0) ++j;
          } else if (k > j) {
            prim_t tP = d[nd];
            Digit l;
            for (l = nd; l > i; --l) d[l] = d[l-1];
            d[l] = tP;
            ++i;
            while (true) {
              if (d[l].c[k]) {
                d[l].c[k] = 0;
                for (Digit l2 = k+1; l2 < n; ++l2) d[l].c[l2] ^= d[i].c[l2];
                d[l].x[0] *= d[i].x[0];
                d[l].x[0] %= a;
                t = gcd(d[l].x[1], d[i].x[1]);
                d[l].x[1] *= d[i].x[1];
                d[l].x[1] /= t*t;
                d[l].x[2] *= t;
                d[l].x[2] *= d[i].x[2];
                d[l].x[2] %= a;
              }
              if (++i > nd) break;
              while (d[i].c[++k] == 0);
            }
            for (i = 0; i < l; ++i)
              if (d[i].c[j]) {
                d[i].c[j] = 0;
                for (Digit l2 = j+1; l2 < n; ++l2) d[i].c[l2] ^= d[l].c[l2];
                d[i].x[0] *= d[l].x[0];
                d[i].x[0] %= a;
                t = gcd(d[l].x[1], d[i].x[1]);
                d[i].x[1] *= d[l].x[1];
                d[i].x[1] /= t*t;
                d[i].x[2] *= t;
                d[i].x[2] *= d[l].x[2];
                d[i].x[2] %= a;
              }
            i = ++nd;
          } else ++k;
        }
        if (j == n) {   // Solution?
          Natural y = sqrt(d[nd].x[1]) * d[nd].x[2];
          y %= a;
          t = a-y;
          if (d[nd].x[0] != y && d[nd].x[0] != t) {
            y += d[nd].x[0];
            t = gcd(y, a);
            v = 2;
            break;
          } else memset(d[nd].c, 0, n);
        } else if (i == nd) {
          for (i = 0; i < nd; ++i)
            if (d[i].c[j]) {
              d[i].c[j] = 0;
              for (Digit l2 = j+1; l2 < n; ++l2) d[i].c[l2] ^= d[nd].c[l2];
              d[i].x[0] *= d[nd].x[0];
              d[i].x[0] %= a;
              t = gcd(d[nd].x[1], d[i].x[1]);
              d[i].x[1] *= d[nd].x[1];
              d[i].x[1] /= t*t;
              d[i].x[2] *= t;
              d[i].x[2] *= d[nd].x[2];
              d[i].x[2] %= a;
            }
          ++nd;
        }
      }
      div(x, v, s, t);                  // Kettenbruchentwicklung
      y = r-t; z += w*s; z %= a;
      if (x >= y) {
        t = x-y; t *= s; u += t;
      } else {
        t = y-x; t *= s; u -= t;
      }
      swap(v, u); swap(x, y); swap(w, z);
    }
    for (i = 0; i < nd; ++i) memset(d[i].c, 0, n);
    //cout << "Iterationen=" << iteration << endl;
    //cout << "Matrix size=" << nd << endl;
  }
  for (i = 0; i <= m; ++i) {
    delete[] d[i].x;
    delete[] d[i].c;
  }
  delete[] d;
  delete[] p;
  return t;
}

bool fact2(Natural a, list<Natural>& p)
// a contains larger prime numbers
{
  Natural t = pollard(a);
  if (t == 1) t = cfrac(a, 400, 257);
  if (t > 1) {
    a /= t;
    bool complete = false;
    if (ispprime(t)) {
      insert(t, p);
      complete = true;
    } else {
      list<Natural> p2;
      complete = fact2(t, p2);
      insert(p2, p);
    }
    if (ispprime(a)) {
      insert(a, p);
    } else {
      list<Natural> p2;
      if (!fact2(a, p2) && complete) complete = false;
      insert(p2, p);
    }
    return complete;
  }
  return false;
}

bool fact(Natural a, list<Natural>& p)
{
  p.erase(p.begin(), p.end());
  if (a == 0) return true;
  while (a.even()) {                      // less factors
    p.push_back(2);
    a >>= 1;
  }
  if (a == 1) return true;
  Natural t;
  Digit i;
  Primes primes;
  Digit prim = primes.firstPrime();            // didn't need 2
  while (primes.nextPrime(prim)) {
    while (true) {
      div(a, prim, t, i);
      if (i) break;
      p.push_back(prim);
      if (t == 1) return true;
      swap(a, t);
    }
    if (t < prim) { p.push_back(a); return true; }
  }
  if (ispprime(a)) { p.push_back(a); return true; }   // greater
  return fact2(a, p);
}

void factLong(Natural a, list<Natural>& p)
// bruteforst for long prime factors
{
  Primes primes;
  Natural s = primes.lastPrime();
  const Digit n[8] = { 4, 2, 4, 2, 4, 6, 2, 6 };
  Digit i = 0;
  switch (s % 30) {
    case 1:  ++i;
    case 29: ++i;
    case 23: ++i;
    case 19: ++i;
    case 17: ++i;
    case 13: ++i;
    case 11: ++i;
  }
  Natural q,r;
  Natural t = root(a, 3);
  for (s += n[i]; s.length() == 1 && s <= t; s += n[i&7]) {
    if (a%s.highest() == 0) {
      p.push_back(s);
      a /= s.highest();
      if (isprime(a)) { p.push_back(a); return; }
      t = root(a, 3);
    }
    ++i;
  }
  while (s <= t) {
    div(a, s, q, r);
    if (r == 0) {
      p.push_back(s);
      if (isprime(q)) { p.push_back(q); return; }
      swap(a, q);
      t = root(a, 3);
    }
    ++i;
    s += n[i&7];
  }
  
  Natural w;                                        // large factors
  sqrt(a, s, w);
  if (w == 0) { p.push_back(s); p.push_back(s); return; }
  s = root(a, 6);
  q = a << 2;
  t *= a; t <<= 2;
  
  Natural x,y,z,d = 4;
  Natural e = 2;
  r = s >> 2;
  const char c[32] = { 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 };
  Natural k = q;
  for (; k <= t; k += q) {
    if (e == 0) {
      d += 4;
      if (d > s) break;
      r = d >> 1;
      div(s, d, r, x);
    } else --e;
    sqrt(k, x, y);
    z = x << 1; ++x; ++z; z -= y;
    while (true) {
      if (c[z.lowest() & 31]) {
        sqrt(z, y, w);
        if (w == 0) {
          y = gcd(x+y, a);
          a /= y;
          if (y < a) { p.push_back(y); p.push_back(a); }
          else { p.push_back(a); p.push_back(y); }
          return;
        }
      }
      if (r == 0) break;
      --r; z = x << 1; ++x; ++z;
    }
  } 
  while (k <= t) {
    sqrt(k, x, z);
    y = x << 1; ++x; ++y; y -= z;
    if (c[y.lowest() & 31]) {
      sqrt(y, y, w);
      if (w == 0) {
        y = gcd(x+y, a);
        a /= y;
        if (y < a) { p.push_back(y); p.push_back(a); }
        else { p.push_back(a); p.push_back(y); }
        return;
      }
    } 
    k += q;
  }
  p.push_back(a);
}

extern "C"
const char* factoring(const char* a)
{
  list<Natural> p;
  Natural n = atoN(a, 10);

  clock_t start = clock();
  bool complete = fact(n, p);
  clock_t stop = clock();
  isFactoringCompleteFlag = complete;
  tmpFactoringDuration = double(stop-start)/CLOCKS_PER_SEC;
  
  strFactoring = "";
  string t = "";
  stringstream s(t);
  s << "<html><head><title>Factorization</title></head><body>";
  for (list<Natural>::iterator j = p.begin(); j != p.end();) {
    list<Natural>::iterator k = j;
    Digit i = 1;
    if (!complete) n /= *k;
    while (++j != p.end() && *j == *k) {
      ++i;
      if (!complete) n /= *k;
    }
    if (i == 1) s << *k << ' ';
    else s << *k << "<sup>" << i << "</sup> ";
    if (j != p.end()) s << "<sup><b>.</b></sup> "; //\n ";
  }
  if (!complete && n > 1) {
    if (p.begin() != p.end()) s << "<sup><b>.</b></sup> ";
    string strLengthOfC;
    stringstream sLengthOfC(strLengthOfC);
    sLengthOfC << n;
    const int l = strlen(sLengthOfC.str().c_str());
    s << "C<sub>" << l << "</sub><p>where C<sub>" << l << "</sub> is a composite number with " << l << " digits.</p>";
  }
  s << "</body></html>";
  strFactoring = s.str();
  return strFactoring.c_str();
}

extern "C"
int isFactoringComplete()
{
  return isFactoringCompleteFlag;
}

extern "C"
const char* longFactoring(const char* a)
{
  list<Natural> p;
  Natural n = atoN(a, 10);
  
  clock_t start = clock();
  if (!fact(n, p)) factLong(n, p);
  clock_t stop = clock();
  tmpFactoringDuration = double(stop-start)/CLOCKS_PER_SEC;
  
  strFactoring = "";
  string t = "";
  stringstream s(t);
  s << "<html><head><title>Factorization</title><body>";
  for (list<Natural>::iterator j = p.begin(); j != p.end();) {
    list<Natural>::iterator k = j;
    Digit i = 1;
    while (++j != p.end() && *j == *k) ++i;
    if (i == 1) s << *k << " ";
    else s << *k << "<sup>" << i << "</sup> ";
    if (j != p.end()) s << "<sup><b>.</b></sup> "; //\n ";
  }
  s << "</body></html>";
  strFactoring = s.str().c_str();
  return strFactoring.c_str();
}

