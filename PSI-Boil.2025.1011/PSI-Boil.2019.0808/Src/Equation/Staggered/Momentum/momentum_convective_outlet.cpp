#include "momentum.h"

/******************************************************************************/
void Momentum::convective_outlet(Vector & veloc, const real fubo) {
  for_m(m)
    for( int b=0; b<veloc.bc(m).count(); b++ ) {

      if( veloc.bc(m).type(b) == BndType::outlet() ) {

        Dir d = veloc.bc(m).direction(b);

        if( d == Dir::imin() ) 
          for_vjk(veloc.bc(m).at(b),j,k) {
            veloc[m][si(m)-1][j][k] += fubo
                                 * (veloc[m][si(m)][j][k] - veloc[m][si(m)-1][j][k])
                                  / veloc.dxw(m,si(m))
                                  * time->dt();
            /* value in further buffers is only extended */
            for(int ii=2; ii<=boil::BW; ++ii)
              veloc[m][si(m)-ii][j][k] = veloc[m][si(m)-1][j][k]
                                       * dSx(m,si(m)-1,j,k)/dSx(m,si(m)-ii,j,k);
          }
        if( d == Dir::imax() ) 
          for_vjk(veloc.bc(m).at(b),j,k) {
             veloc[m][ei(m)+1][j][k] -= fubo 
                                  * (veloc[m][ei(m)+1][j][k]-veloc[m][ei(m)][j][k]) 
                                  / veloc.dxe(m,ei(m)) 
                                  * time->dt();
            /* value in further buffers is only extended */
            for(int ii=2; ii<=boil::BW; ++ii)
              veloc[m][ei(m)+ii][j][k] = veloc[m][ei(m)+1][j][k]
                                       * dSx(m,ei(m)+1,j,k)/dSx(m,ei(m)+ii,j,k);
          }

        if( d == Dir::jmin() ) 
          for_vik(veloc.bc(m).at(b),i,k) {
             veloc[m][i][sj(m)-1][k] += fubo
                                  * (veloc[m][i][sj(m)][k] - veloc[m][i][sj(m)-1][k])
                                  / veloc.dys(m,sj(m))
                                  * time->dt();
            /* value in further buffers is only extended */
              for(int jj=2; jj<=boil::BW; ++jj)
            veloc[m][i][sj(m)-jj][k] = veloc[m][i][sj(m)-1][k]
                                     * dSy(m,i,sj(m)-1,k)/dSy(m,i,sj(m)-jj,k);
          }
        if( d == Dir::jmax() ) 
          for_vik(veloc.bc(m).at(b),i,k) {
             veloc[m][i][ej(m)+1][k] -= fubo
                                  * (veloc[m][i][ej(m)+1][k] - veloc[m][i][ej(m)][k])
                                  / veloc.dyn(m,ej(m))
                                  * time->dt();
            /* value in further buffers is only extended */
            for(int jj=2; jj<=boil::BW; ++jj)
              veloc[m][i][ej(m)+jj][k] = veloc[m][i][ej(m)+1][k]
                                       * dSy(m,i,ej(m)+1,k)/dSy(m,i,ej(m)+jj,k);
          }

        if( d == Dir::kmin() ) 
          for_vij(veloc.bc(m).at(b),i,j) {
             veloc[m][i][j][sk(m)-1] += fubo 
                                  * (veloc[m][i][j][sk(m)] - veloc[m][i][j][sk(m)-1])
                                  / veloc.dzb(m,sk(m))
                                  * time->dt();
            /* value in further buffers is only extended */
            for(int kk=2; kk<=boil::BW; ++kk)
              veloc[m][i][j][sk(m)-kk] = veloc[m][i][j][sk(m)-1]
                                       * dSz(m,i,j,sk(m)-1)/dSz(m,i,j,sk(m)-kk);
          }
        if( d == Dir::kmax() ) 
          for_vij(veloc.bc(m).at(b),i,j) {
             veloc[m][i][j][ek(m)+1] -= fubo
                                  * (veloc[m][i][j][ek(m)+1] - veloc[m][i][j][ek(m)])
                                  / veloc.dzt(m,ek(m))
                                  * time->dt();
            /* value in further buffers is only extended */
            for(int kk=2; kk<=boil::BW; ++kk)
              veloc[m][i][j][ek(m)+kk] = veloc[m][i][j][ek(m)+1]
                                       * dSz(m,i,j,ek(m)+1)/dSz(m,i,j,ek(m)+kk);
          }
      }
    } /* m & b */

  return;
}
