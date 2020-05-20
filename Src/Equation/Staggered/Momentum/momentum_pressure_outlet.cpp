#include "momentum.h"

/******************************************************************************/
void Momentum::pressure_outlet(const Scalar & frc) {
 
  Comp m;
  real rho;

  /*----------------------------------------------------+
  |  project the velocities on the computational cells  |
  +- - - - - - - - - - - - - - - - - - - - - - - - - - -+
  |     not the most efficient -> too many divisions    |
  +- - - - - - - - - - - - - - - - - - - - - - - - - - -+
  |      do not try to introduce some integral form     |
  |       here gradients are constant in each cell      |
  |          and depend on boundary values only         |
  +----------------------------------------------------*/

  for_m(m) {
    for(int b=0; b<u.bc(m).count(); b++) {
      if( u.bc(m).type(b) == BndType::outlet() ) {

        Dir d = u.bc(m).direction(b);

        if( m == Comp::u() && d == Dir::imin() ) {
          for_vjk(u.bc(m).at(b),j,k) {
            rho = fluid()->rho(m,si(m),j,k);
            real corr = (frc[si(m)-2][j][k]-frc[si(m)-1][j][k])
                        / rho * time->dt()
                        / (u.dxc(m,si(m)-1));
            for(int ii=1; ii<=boil::BW; ii++)
              u[m][si(m)-ii][j][k] += corr;
          }
        }
        if( m == Comp::u() && d == Dir::imax() ) {
          for_vjk(u.bc(m).at(b),j,k) {
            rho = fluid()->rho(m,ei(m),j,k);
            real corr = (frc[ei(m)][j][k]-frc[ei(m)+1][j][k])
                        / rho * time->dt()
                        / (u.dxc(m,ei(m)+1));
            for(int ii=1; ii<=boil::BW; ii++)
              u[m][ei(m)+ii][j][k] += corr;
          }
        }

        if( m == Comp::v() && d == Dir::jmin() ) {
          for_vik(u.bc(m).at(b),i,k) {
            rho = fluid()->rho(m,i,sj(m),k);
            real corr = (frc[i][sj(m)-2][k]-frc[i][sj(m)-1][k])
                        / rho * time->dt()
                        / (u.dyc(m,sj(m)-1)); 
            for(int jj=1; jj<=boil::BW; jj++)
              u[m][i][sj(m)-jj][k] += corr;
          }
        }
        if( m == Comp::v() && d == Dir::jmax() ) {
          for_vik(u.bc(m).at(b),i,k) {
            rho = fluid()->rho(m,i,ej(m),k);
            real corr = (frc[i][ej(m)][k]-frc[i][ej(m)+1][k])
                        / rho * time->dt()
                        / (u.dyc(m,ej(m)+1)); 
            for(int jj=1; jj<=boil::BW; jj++)
              u[m][i][ej(m)+jj][k] += corr;
          }
        }

        if( m == Comp::w() && d == Dir::kmin() ) {
          for_vij(u.bc(m).at(b),i,j) {
            rho = fluid()->rho(m,i,j,sk(m));
            real corr = (frc[i][j][sk(m)-2]-frc[i][j][sk(m)-1])
                        / rho * time->dt()
                        / (u.dzc(m,sk(m)-1));
            for(int kk=1; kk<=boil::BW; kk++)
              u[m][i][j][sk(m)-kk] += corr;
          }
        }
        if( m == Comp::w() && d == Dir::kmax() ) {
          for_vij(u.bc(m).at(b),i,j) {
            rho = fluid()->rho(m,i,j,ek(m));
            real corr = (frc[i][j][ek(m)]-frc[i][j][ek(m)+1])
                        / rho * time->dt()
                        / (u.dzc(m,ek(m)+1));
            for(int kk=1; kk<=boil::BW; kk++)
              u[m][i][j][ek(m)+kk] += corr;
          }
        }
      } /* outlet */
    } /* bc */
  } /* m */

  u.exchange_all();

  return;
}
