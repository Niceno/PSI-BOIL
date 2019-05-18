#include "momentum.h"

/******************************************************************************/
void Momentum::outlet() {
/*----------------------------------+ 
|  set bulk velocity at the outlet  |
+----------------------------------*/

  /*---------------------------------+
  |  if outlet does not exist, exit  |
  +---------------------------------*/
  if( 
      !u.bc( Comp::u() ).exists( BndType::outlet() ) &&
      !u.bc( Comp::v() ).exists( BndType::outlet() ) && 
      !u.bc( Comp::w() ).exists( BndType::outlet() ) 
    ) return;

  real ao=0.0, aox=0.0, aoy=0.0, aoz=0.0, ubo=0.0;

  /*---------------------------------------+
  |  get volume flux the inlet and outlet  |
  +---------------------------------------*/
  /* v_phase_change() should be called from main.cpp, if phase change occurs */
  const real volf_in  = volf_bct( BndType::inlet() ) + v_phase_change;
  const real volf_out = volf_bct( BndType::outlet(), &aox, &aoy, &aoz ); 

  /* compute bulk velocity */
  ubo = volf_out / (aox + aoy + aoz);

  /* outlet bulk velocity is negative if it leaves the domain.  
     let's make it positive from here on */
  ubo = fabs(ubo);

  real ratio;
  if(volf_out==0.0 && volf_in==0.0){
    ratio=1.0;
  } else if(volf_out==0.0 && volf_in!=0.0){
    ratio=1.0e+12;
  } else {
    ratio= -volf_in / volf_out;
  }

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
              u[m][si(m)-1][j][k] += ubo
                                   * (u[m][si(m)][j][k] - u[m][si(m)-1][j][k])
                                    / u.dxw(m,si(m))
                                    * time->dt();
              /* value in further buffers is only extended */
              for(int ii=2; ii<=boil::BW; ++ii)
                u[m][si(m)-ii][j][k] = u[m][si(m)-1][j][k];
            }
          if( d == Dir::imax() ) 
            for_vjk(u.bc(m).at(b),j,k) {
               u[m][ei(m)+1][j][k] -= ubo 
                                    * (u[m][ei(m)+1][j][k]-u[m][ei(m)][j][k]) 
                                    / u.dxe(m,ei(m)) 
                                    * time->dt();
              /* value in further buffers is only extended */
              for(int ii=2; ii<=boil::BW; ++ii)
                u[m][ei(m)+ii][j][k] = u[m][ei(m)+1][j][k];
            }

          if( d == Dir::jmin() ) 
            for_vik(u.bc(m).at(b),i,k) {
               u[m][i][sj(m)-1][k] += ubo
                                    * (u[m][i][sj(m)][k] - u[m][i][sj(m)-1][k])
                                    / u.dys(m,sj(m))
                                    * time->dt();
              /* value in further buffers is only extended */
              for(int jj=2; jj<=boil::BW; ++jj)
                u[m][i][sj(m)-jj][k] = u[m][i][sj(m)-1][k];
            }
          if( d == Dir::jmax() ) 
            for_vik(u.bc(m).at(b),i,k) {
               u[m][i][ej(m)+1][k] -= ubo
                                    * (u[m][i][ej(m)+1][k] - u[m][i][ej(m)][k])
                                    / u.dyn(m,ej(m))
                                    * time->dt();
              /* value in further buffers is only extended */
              for(int jj=2; jj<=boil::BW; ++jj)
                u[m][i][ej(m)+jj][k] = u[m][i][ej(m)+1][k];
            }

          if( d == Dir::kmin() ) 
            for_vij(u.bc(m).at(b),i,j) {
               u[m][i][j][sk(m)-1] += ubo 
                                    * (u[m][i][j][sk(m)] - u[m][i][j][sk(m)-1])
                                    / u.dzb(m,sk(m))
                                    * time->dt();
              /* value in further buffers is only extended */
              for(int kk=2; kk<=boil::BW; ++kk)
                u[m][i][j][sk(m)-kk] = u[m][i][j][sk(m)-1];
            }
          if( d == Dir::kmax() ) 
            for_vij(u.bc(m).at(b),i,j) {
               u[m][i][j][ek(m)+1] -= ubo
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
