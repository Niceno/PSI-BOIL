#include "momentum.h"

/******************************************************************************/
void Momentum::extrapolate_outlet_velocity(const real fubo, const real ratio) {
  /*-------------------------------------------------+ 
  |  set range in which convective outlet will work  |
  +-------------------------------------------------*/
  Range<real> valid(0.97, 1.03);

  /*-------------------------+
  |  use convective outflow  |
  +-------------------------*/
  if( valid.contains(ratio) ) {
    boil::oout << "using convective outflow " << boil::endl;
    for_m(m)
      for( int b=0; b<u.bc(m).count(); b++ ) {

        if( u.bc(m).type(b) == BndType::outlet() ) {

          Dir d = u.bc(m).direction(b);

          if( d == Dir::imin() ) 
            for_vjk(u.bc(m).at(b),j,k) {
              u[m][si(m)-1][j][k] += fubo
                                   * (u[m][si(m)][j][k] - u[m][si(m)-1][j][k])
                                    / u.dxw(m,si(m))
                                    * time->dt();
              /* value in further buffers is only extended */
              for(int ii=2; ii<=boil::BW; ++ii)
                u[m][si(m)-ii][j][k] = u[m][si(m)-1][j][k];
            }
          if( d == Dir::imax() ) 
            for_vjk(u.bc(m).at(b),j,k) {
               u[m][ei(m)+1][j][k] -= fubo 
                                    * (u[m][ei(m)+1][j][k]-u[m][ei(m)][j][k]) 
                                    / u.dxe(m,ei(m)) 
                                    * time->dt();
              /* value in further buffers is only extended */
              for(int ii=2; ii<=boil::BW; ++ii)
                u[m][ei(m)+ii][j][k] = u[m][ei(m)+1][j][k];
            }

          if( d == Dir::jmin() ) 
            for_vik(u.bc(m).at(b),i,k) {
               u[m][i][sj(m)-1][k] += fubo
                                    * (u[m][i][sj(m)][k] - u[m][i][sj(m)-1][k])
                                    / u.dys(m,sj(m))
                                    * time->dt();
              /* value in further buffers is only extended */
              for(int jj=2; jj<=boil::BW; ++jj)
                u[m][i][sj(m)-jj][k] = u[m][i][sj(m)-1][k];
            }
          if( d == Dir::jmax() ) 
            for_vik(u.bc(m).at(b),i,k) {
               u[m][i][ej(m)+1][k] -= fubo
                                    * (u[m][i][ej(m)+1][k] - u[m][i][ej(m)][k])
                                    / u.dyn(m,ej(m))
                                    * time->dt();
              /* value in further buffers is only extended */
              for(int jj=2; jj<=boil::BW; ++jj)
                u[m][i][ej(m)+jj][k] = u[m][i][ej(m)+1][k];
            }

          if( d == Dir::kmin() ) 
            for_vij(u.bc(m).at(b),i,j) {
               u[m][i][j][sk(m)-1] += fubo 
                                    * (u[m][i][j][sk(m)] - u[m][i][j][sk(m)-1])
                                    / u.dzb(m,sk(m))
                                    * time->dt();
              /* value in further buffers is only extended */
              for(int kk=2; kk<=boil::BW; ++kk)
                u[m][i][j][sk(m)-kk] = u[m][i][j][sk(m)-1];
            }
          if( d == Dir::kmax() ) 
            for_vij(u.bc(m).at(b),i,j) {
               u[m][i][j][ek(m)+1] -= fubo
                                    * (u[m][i][j][ek(m)+1] - u[m][i][j][ek(m)])
                                    / u.dzt(m,ek(m))
                                    * time->dt();
              /* value in further buffers is only extended */
              for(int kk=2; kk<=boil::BW; ++kk)
                u[m][i][j][ek(m)+kk] = u[m][i][j][ek(m)+1];
            }
        }
      } /* m & b */

  /*-----------------------------------+
  |  use vanishing derivative outflow  |
  +-----------------------------------*/
  } else {
    boil::oout << "using vanishing derivative outflow " << boil::endl;
    for_m(m)
      for( int b=0; b<u.bc(m).count(); b++ ) {
  
        if( u.bc(m).type(b) == BndType::outlet() ) {
  
          Dir d = u.bc(m).direction(b);
  
          if( d == Dir::imin() ) 
            for_vjk(u.bc(m).at(b),j,k) 
              for(int ii=1; ii<=boil::BW; ++ii)
                u[m][si(m)-ii][j][k] = u[m][si(m)][j][k];
          if( d == Dir::imax() ) 
            for_vjk(u.bc(m).at(b),j,k)
              for(int ii=1; ii<=boil::BW; ++ii)
                u[m][ei(m)+ii][j][k] = u[m][ei(m)][j][k];
  
          if( d == Dir::jmin() ) 
            for_vik(u.bc(m).at(b),i,k)
              for(int jj=1; jj<=boil::BW; ++jj)
                u[m][i][sj(m)-jj][k] = u[m][i][sj(m)][k];
          if( d == Dir::jmax() ) 
            for_vik(u.bc(m).at(b),i,k)
              for(int jj=1; jj<=boil::BW; ++jj)
                u[m][i][ej(m)+jj][k] = u[m][i][ej(m)][k];
  
          if( d == Dir::kmin() ) 
            for_vij(u.bc(m).at(b),i,j)
              for(int kk=1; kk<=boil::BW; ++kk)
                u[m][i][j][sk(m)-kk] = u[m][i][j][sk(m)];
          if( d == Dir::kmax() ) 
            for_vij(u.bc(m).at(b),i,j)
              for(int kk=1; kk<=boil::BW; ++kk)
                u[m][i][j][ek(m)+kk] = u[m][i][j][ek(m)];
        }
      } /* m & b */
  }
}
