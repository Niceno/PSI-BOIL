#include "momentum.h"

/******************************************************************************/
void Momentum::vanishing_derivative_outlet(Vector & veloc) {

  for_m(m)
    for( int b=0; b<veloc.bc(m).count(); b++ ) {

      if( veloc.bc(m).type(b) == BndType::outlet() ) {

        Dir d = veloc.bc(m).direction(b);

        if( d == Dir::imin() ) 
          for_vjk(veloc.bc(m).at(b),j,k) 
            for(int ii=1; ii<=boil::BW; ++ii)
              veloc[m][si(m)-ii][j][k] = veloc[m][si(m)][j][k]
                                       * dSx(m,si(m),j,k)/dSx(m,si(m)-ii,j,k);
        if( d == Dir::imax() ) 
          for_vjk(veloc.bc(m).at(b),j,k)
            for(int ii=1; ii<=boil::BW; ++ii)
              veloc[m][ei(m)+ii][j][k] = veloc[m][ei(m)][j][k]
                                       * dSx(m,ei(m),j,k)/dSx(m,ei(m)+ii,j,k);

        if( d == Dir::jmin() ) 
          for_vik(veloc.bc(m).at(b),i,k)
            for(int jj=1; jj<=boil::BW; ++jj)
              veloc[m][i][sj(m)-jj][k] = veloc[m][i][sj(m)][k]
                                       * dSy(m,i,sj(m),k)/dSy(m,i,sj(m)-jj,k);
        if( d == Dir::jmax() ) 
          for_vik(veloc.bc(m).at(b),i,k)
            for(int jj=1; jj<=boil::BW; ++jj)
              veloc[m][i][ej(m)+jj][k] = veloc[m][i][ej(m)][k]
                                       * dSy(m,i,ej(m),k)/dSy(m,i,ej(m)+jj,k);

        if( d == Dir::kmin() ) 
          for_vij(veloc.bc(m).at(b),i,j)
            for(int kk=1; kk<=boil::BW; ++kk)
              veloc[m][i][j][sk(m)-kk] = veloc[m][i][j][sk(m)]
                                       * dSz(m,i,j,sk(m))/dSz(m,i,j,sk(m)-kk);
        if( d == Dir::kmax() ) 
          for_vij(veloc.bc(m).at(b),i,j)
            for(int kk=1; kk<=boil::BW; ++kk)
              veloc[m][i][j][ek(m)+kk] = veloc[m][i][j][ek(m)]
                                       * dSz(m,i,j,ek(m))/dSz(m,i,j,ek(m)+kk);
      }
    } /* m & b */

  return;
}
