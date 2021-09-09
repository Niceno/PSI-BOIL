#ifndef RATIONAL_H
#define RATIONAL_H

///////////////////////
//                   //
// Rational numbers  //
//                   //
///////////////////////
class Rational {
  public:
    Rational() : n(0), d(1) {}
    Rational(const int N, const int D) : n(N), d(D) {}

    real evaluate() const { return real(n)/real(d); }

    int & num() { return n; }
    int & den() { return d; }

  private:
    int n,d;
};

inline int operator * (Rational lhs, const int & rhs) {
  return ( lhs.num() * rhs )/lhs.den();
}

#endif
