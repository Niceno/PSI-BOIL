#include "vector.h"                

/******************************************************************************/
void Vector::copy(const Vector & other, 
                  const Dir & dir, const real & y_off, const real & z_off) {

  int cnt = 0;

  if(dir == Dir::imin() && dom->coord(Comp::i())==0) /* check on coord needed? */
    for_m(m)
      for_mjk(m,j,k)
        for_vmjk(other,m,J,K) 
          if( approx( yc(m,j), other.yc(m,J)+y_off ) &&
              approx( zc(m,k), other.zc(m,K)+z_off ) ) {
            vec[m][si(m)][j][k] = other[m][other.ei(m)][J][K];
            cnt++;
          }

  if(cnt) APR(cnt);

  exchange();
}

/******************************************************************************/
void Vector::copy(const Vector & other, 
                  const real & x_off, const Dir & dir, const real & z_off) {

  int cnt = 0;

  if(dir == Dir::jmin() && dom->coord(Comp::j())==0) /* check on coord needed? */
    for_m(m)
      for_mik(m,i,k)
        for_vmik(other,m,I,K) 
          if( approx( xc(m,i), other.xc(m,I)+x_off ) &&
              approx( zc(m,k), other.zc(m,K)+z_off ) ) {
            vec[m][i][sj(m)][k] = other[m][I][other.ej(m)][K];
            cnt++;
          }

  if(cnt) APR(cnt);

  exchange();
}

/******************************************************************************/
void Vector::copy(const Vector & other, 
                  const real & x_off, const real & y_off, const Dir & dir) {

  int cnt = 0;

  if(dir == Dir::kmin() && dom->coord(Comp::k())==0) /* check on coord needed? */
    for_m(m)
      for_mij(m,i,j)
        for_vmij(other,m,I,J) 
          if( approx( xc(m,i), other.xc(m,I)+x_off ) &&
              approx( yc(m,j), other.yc(m,J)+y_off ) ) {
            vec[m][i][j][sk(m)] = other[m][I][J][other.ek(m)];
            cnt++;
          }

  if(cnt) APR(cnt);

  exchange();
}

/*-----------------------------------------------------------------------------+
 '$Id: vector_copy.cpp,v 1.10 2011/07/07 07:00:54 niceno Exp $'/
+-----------------------------------------------------------------------------*/
