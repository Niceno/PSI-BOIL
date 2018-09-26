#ifndef SCALAR_ACC_H
#define SCALAR_ACC_H

/* forward declaration needed for operators */
class Matrix;
class Scalar;

/*-----------+
|  alfa * x  |
+-----------*/
struct alfa_x {
  const real & alfa; 
  const Scalar & x;
  alfa_x(const real & a_, const Scalar & x_) : alfa(a_), x(x_) {};
};

inline alfa_x operator * 
 (const real & alfa, const Scalar & x) 
   {return alfa_x(alfa, x);}

/*-----------+
|  x * alfa  |
+-----------*/
struct x_alfa {
  const Scalar & x;
  const real & alfa; 
  x_alfa(const Scalar & x_, const real & a_) : x(x_), alfa(a_) {};
};

inline x_alfa operator * 
 (const Scalar & x, const real & alfa) 
   {return x_alfa(x, alfa);}

/*-----------+
|  x / alfa  |
+-----------*/
struct x_d_alfa {
  const Scalar & x;
  const real & alfa; 
  x_d_alfa(const Scalar & x_, const real & a_) : x(x_), alfa(a_) {};
};

inline x_d_alfa operator / 
 (const Scalar & x, const real & alfa) 
   {return x_d_alfa(x, alfa);}

/*---------------+
|  y + alfa * x  |
+---------------*/
struct y_p_alfa_x {
  const Scalar & y;
  const alfa_x & c; // c stands for compound
  y_p_alfa_x(const Scalar & y_, const alfa_x & c_) : y(y_), c(c_) {}; 
};

inline y_p_alfa_x operator +
 (const Scalar & y, const alfa_x & ax)
   {return y_p_alfa_x(y, ax);}

/*---------------+
|  y - alfa * x  |
+---------------*/
struct y_m_alfa_x {
  const Scalar & y;
  const alfa_x & c; // c stands for compound
  y_m_alfa_x(const Scalar & y_, const alfa_x & c_) : y(y_), c(c_) {}; 
};

inline y_m_alfa_x operator -
 (const Scalar & y, const alfa_x & ax)
   {return y_m_alfa_x(y, ax);}

/*--------+
|  A * x  |
+--------*/
struct A_x {
  const Matrix & A; 
  const Scalar & x;
  A_x(const Matrix & A_, const Scalar & x_) : A(A_), x(x_) {};
};

inline A_x operator * 
 (const Matrix & A, const Scalar & x) 
   {return A_x(A, x);}

/*------------+
|  y + A * x  |
+------------*/
struct y_p_A_x {
  const Scalar & y;
  const A_x    & c; // c stands for compound
  y_p_A_x(const Scalar & y_, const A_x & c_) : y(y_), c(c_) {}; 
};

inline y_p_A_x operator +
 (const Scalar & y, const A_x & ax)
   {return y_p_A_x(y, ax);}

/*------------+
|  y - A * x  |
+------------*/
struct y_m_A_x {
  const Scalar & y;
  const A_x    & c; // c stands for compound
  y_m_A_x(const Scalar & y_, const A_x & c_) : y(y_), c(c_) {}; 
};

inline y_m_A_x operator -
 (const Scalar & y, const A_x & ax)
   {return y_m_A_x(y, ax);}

/*--------+
|  x * y  | -> this is NOT a dot vector product
+--------*/
struct x_y {
  const Scalar & x; 
  const Scalar & y;
  x_y(const Scalar & x_, const Scalar & y_) : x(x_), y(y_) {};
};

inline x_y operator * 
 (const Scalar & x, const Scalar & y) 
   {return x_y(x, y);}

/*--------+
|  x / y  | 
+--------*/
struct x_d_y {
  const Scalar & x; 
  const Scalar & y;
  x_d_y(const Scalar & x_, const Scalar & y_) : x(x_), y(y_) {};
};

inline x_d_y operator / 
 (const Scalar & x, const Scalar & y) 
   {return x_d_y(x, y);}

/*--------+
|  x + y  | 
+--------*/
struct x_p_y {
  const Scalar & x; 
  const Scalar & y;
  x_p_y(const Scalar & x_, const Scalar & y_) : x(x_), y(y_) {};
};

inline x_p_y operator + 
 (const Scalar & x, const Scalar & y) 
   {return x_p_y(x, y);}

/*--------+
|  x - y  | 
+--------*/
struct x_m_y {
  const Scalar & x; 
  const Scalar & y;
  x_m_y(const Scalar & x_, const Scalar & y_) : x(x_), y(y_) {};
};

inline x_m_y operator - 
 (const Scalar & x, const Scalar & y) 
   {return x_m_y(x, y);}

#endif

/*-----------------------------------------------------------------------------+
 '$Id: scalar_acc.h,v 1.10 2011/05/29 09:44:20 niceno Exp $'/
+-----------------------------------------------------------------------------*/
