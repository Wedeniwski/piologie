/////////////////////////////////
//
// Piologie V 1.3.2
// multi-precision arithmetic
// Calculation of Constants
//
// (c) 1996-2001 HiPiLib
//               www.hipilib.de
//
// Sebastian Wedeniwski
// 10/21/2001
//

#ifndef _Include_Pi_H_
#define _Include_Pi_H_

#include "integer.h"
#include "nmbrthry.h"

class Fixed : private Natural {
private:
  size_t prec,dec;

public:
  Fixed(const size_t);
  Fixed(const Fixed&);

  void      output(OSTREAM&);
  void      set_decimals(const size_t);
  size_t    precision() const;
  size_t    decimals() const;
  Natural&  value() const;
};


inline OSTREAM& operator<<(OSTREAM&, Fixed);

class Pi {
private:
  Fixed pi;

  void stoermer(const size_t, Natural&);
  void schoenhage(const size_t, Natural&);
  void chudnovsky(Digit, const Digit, Integer&, Integer&, Integer&) const;
  void chudnovsky2(const Digit, const Digit, Integer&, Integer&) const;
  void sqrt_series(Digit, const Digit, Integer&, Integer&, Integer&) const;
  void sqrt_series2(const Digit, const Digit, Integer&, Integer&) const;

public:
  Pi(const size_t);

  inline friend OSTREAM& operator<<(OSTREAM&, const Pi&);
};

class Sqrt {
private:
  Digit d;
  Fixed sqrt;
  Natural u,v;

  void series(Digit, const Digit, Natural&, Natural&, Integer&) const;
  void series2(const Digit, const Digit, Natural&, Integer&) const;

public:
  Sqrt(const Digit, const size_t);

  inline friend OSTREAM& operator<<(OSTREAM&, const Sqrt&);
};

class Zeta3 {
private:
  Fixed zeta;

  void linear(const size_t, Natural&) const;
  void series(Digit, const Digit, Integer&, Integer&, Integer&) const;
  void series2(const Digit, const Digit, Integer&, Integer&) const;

public:
  Zeta3(const size_t);

  inline friend OSTREAM& operator<<(OSTREAM&, const Zeta3&);
};

class Exp1 : protected Natural {
private:
  Fixed exp;

  void linear(const size_t, Natural&) const;
  void series(Digit, const Digit, Natural&, Natural&) const;

public:
  Exp1(const size_t);

  inline friend OSTREAM& operator<<(OSTREAM&, const Exp1&);
};

class Ln : protected Natural {
private:
  Fixed ln;
  Digit u,v,x;

  Natural atanh_inv_linear(const Digit, const size_t) const;
  void    atanh_inv_series(Digit, const Digit, Natural&, Natural&, Natural&) const;
  Natural atanh_inv_series(const Digit);
  Natural ln2(const size_t);

  void    linear(const size_t, Natural&) const;
  void    series(Digit, const Digit, Natural&, Integer&, Digit&, Integer&) const;

public:
  Ln(const Digit, const size_t);

  inline const Natural& value() const;

  inline friend OSTREAM& operator<<(OSTREAM&, const Ln&);
};

class EulerGamma : protected Natural {
private:
  Fixed euler;

  void linear(const size_t, Natural&);
  void series(Digit, const Digit,
              Natural&, Natural&, Natural&, Natural&, Natural&, Natural&) const;
  void series2(Digit, const Digit, Natural&, Natural&, Natural&, Natural&) const;

public:
  EulerGamma(const size_t);

  inline friend OSTREAM& operator<<(OSTREAM&, const EulerGamma&);
};

/////////////////// Inline-Implementation ///////////////////////

inline Natural& Fixed::value() const
{
  return (Natural&)*this;
}

inline Fixed::Fixed(const size_t n)
 : prec(NumberOfDecimals(n+ALPHA_WIDTH)), dec(n)
{
  RestoreSize();
}

inline Fixed::Fixed(const Fixed& a)
 : Natural(a.value()), prec(a.prec), dec(a.dec)
{
}

inline void Fixed::set_decimals(const size_t n)
{
  prec = NumberOfDecimals(n+ALPHA_WIDTH);
  dec = n;
  RestoreSize();
}

inline size_t Fixed::precision() const
{
  return prec;
}

inline size_t Fixed::decimals() const
{
  return dec;
}

inline OSTREAM& operator<<(OSTREAM& out, Fixed a)
// Algorithm:  o := operator<<(o, a)
// Input:      o in ostream, a in Fixed.
// Output:     o in ostream ||
//
// Note:       puts Fixed a on output stream.
{
  a.output(out);
  return out;
}

inline OSTREAM& operator<<(OSTREAM& out, const Pi& a)
// Algorithm:  o := o << a
// Input:      o in ostream, a in Pi.
// Output:     o in ostream ||
//
// Note:       puts Pi a on output stream.
{
  return out << a.pi;
}

inline OSTREAM& operator<<(OSTREAM& out, const Sqrt& a)
// Algorithm:  o := o << a
// Input:      o in ostream, a in Sqrt.
// Output:     o in ostream ||
//
// Note:       puts Sqrt a on output stream.
{
  return out << a.sqrt;
}

inline OSTREAM& operator<<(OSTREAM& out, const Zeta3& a)
// Algorithm:  o := o << a
// Input:      o in ostream, a in Zeta3.
// Output:     o in ostream ||
//
// Note:       puts Zeta3 a on output stream.
{
  return out << a.zeta;
}

inline OSTREAM& operator<<(OSTREAM& out, const Exp1& a)
// Algorithm:  o := o << a
// Input:      o in ostream, a in Exp1.
// Output:     o in ostream ||
//
// Note:       puts Exp1 a on output stream.
{
  return out << a.exp;
}

inline const Natural& Ln::value() const
{
  return ln.value();
}
inline OSTREAM& operator<<(OSTREAM& out, const Ln& a)
// Algorithm:  o := o << a
// Input:      o in ostream, a in Ln.
// Output:     o in ostream ||
//
// Note:       puts Ln a on output stream.
{
  return out << a.ln;
}

inline OSTREAM& operator<<(OSTREAM& out, const EulerGamma& a)
// Algorithm:  o := o << a
// Input:      o in ostream, a in EulerGamma.
// Output:     o in ostream ||
//
// Note:       puts EulerGamma a on output stream.
{
  return out << a.euler;
}


#endif
