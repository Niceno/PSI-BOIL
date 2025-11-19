#include "scalar.h"
#include "../../Matrix/matrix.h"

/******************************************************************************/
const Scalar & Scalar::operator = (const char * c) {

  Formula F;

  for_aijk(i,j,k) {
    std::stringstream x, y, z, f;
    x << "x=" << xc(i); F.evaluate(x);
    y << "y=" << yc(j); F.evaluate(y);
    z << "z=" << zc(k); F.evaluate(z);
    f << c;

    val[i][j][k] = F.evaluate(f);
  }
  return *this;
}

/******************************************************************************/
const Scalar & Scalar::operator = (const y_m_A_x & c) {

  c.c.x.exchange();

  /* central coefficient */
  for_ijk(i,j,k)
    val[i][j][k] = c.y[i][j][k] - c.c.A.c[i][j][k] * c.c.x[i][j][k];

  /* w - e coefficient */
  for_ijk(i,j,k)
    val[i][j][k] += c.c.A.w[i][j][k] * c.c.x[i-1][j][k]
                  + c.c.A.e[i][j][k] * c.c.x[i+1][j][k];

  /* s - n coefficient */
  for_ijk(i,j,k)
    val[i][j][k] += c.c.A.s[i][j][k] * c.c.x[i][j-1][k]
                  + c.c.A.n[i][j][k] * c.c.x[i][j+1][k];

  /* b - t coefficient */
  for_ijk(i,j,k)
    val[i][j][k] += c.c.A.b[i][j][k] * c.c.x[i][j][k-1]
                  + c.c.A.t[i][j][k] * c.c.x[i][j][k+1];

  return *this;
}

/******************************************************************************/
const Scalar & Scalar::operator = (const A_x & c) {

  c.x.exchange();

  /* central coefficient */
  for_ijk(i,j,k) val[i][j][k] = c.A.c[i][j][k] * c.x[i][j][k];

  /* w - e coefficient */
  for_ijk(i,j,k) 
    val[i][j][k] -= (  c.A.w[i][j][k] * c.x[i-1][j][k] 
                     + c.A.e[i][j][k] * c.x[i+1][j][k]);

  /* s - n coefficient */
  for_ijk(i,j,k) 
    val[i][j][k] -= (  c.A.s[i][j][k] * c.x[i][j-1][k]
                     + c.A.n[i][j][k] * c.x[i][j+1][k]);

  /* b - t coefficient */
  for_ijk(i,j,k) 
    val[i][j][k] -= (  c.A.b[i][j][k] * c.x[i][j][k-1]
                     + c.A.t[i][j][k] * c.x[i][j][k+1]);

  return *this;
}
