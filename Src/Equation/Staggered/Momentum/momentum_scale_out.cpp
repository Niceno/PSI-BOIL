#include "momentum.h"

/******************************************************************************/
void Momentum::scale_out() {
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

  real aox=0.0, aoy=0.0, aoz=0.0;

  /*---------------------------------------+
  |  get volume flux the inlet and outlet  |
  +---------------------------------------*/
  /* v_phase_change() should be called from main.cpp, if phase change occurs */
  const real volf_in  = volf_bct( BndType::inlet() )
                      + volf_bct( BndType::insert() )
                      + v_phase_change;
  const real volf_out = volf_bct( BndType::outlet(), &aox, &aoy, &aoz ); 

  /* compute bulk velocity */
  const real ub = volf_in/ (aox + aoy + aoz);

  real ratio;
  if(volf_out==0.0 && volf_in==0.0){
    ratio=1.0;
  } else if(volf_out==0.0 && volf_in!=0.0){
    ratio=1.0e+12;
  } else {
    ratio= -volf_in / volf_out;
  }
/*
  boil::oout << "inlet / outlet volume flux ratio: " 
             << ratio
             << boil::endl;
*/
  Range<real> valid(0.5, 2.0);

  /*------------------------+
  |  nothing at the outlet  |
  +------------------------*/
  if( !valid.contains(ratio) ) {
    /* bad, code duplication, this has been coppied from insert_bc */
    for_m(m){
    for( int b=0; b<u.bc(m).count(); b++ ) {

      if( u.bc(m).type(b) == BndType::outlet() ) {

        Dir d = u.bc(m).direction(b);

        if( m == Comp::u() && d == Dir::imin() )
          for_vjk(u.bc(m).at(b),j,k) u[m][si(m)-1][j][k] = -ub;
        if( m == Comp::u() && d == Dir::imax() )
          for_vjk(u.bc(m).at(b),j,k) u[m][ei(m)+1][j][k] = +ub;

        if( m == Comp::v() && d == Dir::jmin() )
          for_vik(u.bc(m).at(b),i,k) u[m][i][sj(m)-1][k] = -ub;
        if( m == Comp::v() && d == Dir::jmax() )
          for_vik(u.bc(m).at(b),i,k) u[m][i][ej(m)+1][k] = +ub;

        if( m == Comp::w() && d == Dir::kmin() )
          for_vij(u.bc(m).at(b),i,j) u[m][i][j][sk(m)-1] = -ub;
        if( m == Comp::w() && d == Dir::kmax() )
          for_vij(u.bc(m).at(b),i,j) u[m][i][j][ek(m)+1] = +ub;
      }
    } /* b */
    }
  }

  /*--------------------------+
  |  something at the outlet  |
  +--------------------------*/
  else {

    for_m(m) {              
      /* again code duplication, this has been coppied from insert_bc */
      for( int b=0; b<u.bc(m).count(); b++ ) {

        if( u.bc(m).type(b) == BndType::outlet() ) {

          Dir d = u.bc(m).direction(b);

          if( d == Dir::imin() ) 
            for_vjk(u.bc(m).at(b),j,k) u[m][si(m)-1][j][k] *= ratio;    
          if( d == Dir::imax() ) 
            for_vjk(u.bc(m).at(b),j,k) u[m][ei(m)+1][j][k] *= ratio;    

          if( d == Dir::jmin() ) 
            for_vik(u.bc(m).at(b),i,k) u[m][i][sj(m)-1][k] *= ratio;    
          if( d == Dir::jmax() ) 
            for_vik(u.bc(m).at(b),i,k) u[m][i][ej(m)+1][k] *= ratio;    

          if( d == Dir::kmin() ) 
            for_vij(u.bc(m).at(b),i,j) u[m][i][j][sk(m)-1] *= ratio;    
          if( d == Dir::kmax() ) 
            for_vij(u.bc(m).at(b),i,j) u[m][i][j][ek(m)+1] *= ratio;    
        }
      } /* b */
    } /* m */
  }
}

/*-----------------------------------------------------------------------------+
 '$Id: momentum_scale_out.cpp,v 1.18 2015/01/05 17:30:26 sato Exp $'/
+-----------------------------------------------------------------------------*/
