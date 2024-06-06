#include "momentum.h"

/******************************************************************************/
void Momentum::scale_outlet_velocity(const real ubo, const real ratio) {

  Range<real> valid(0.5, 2.0);

  /*---------------------+
  |  ratio out of range  |
  +---------------------*/
  if( !valid.contains(ratio) ) {
    /* bad, code duplication, this has been copied from insert_bc */
    //boil::oout<<"scale_outlet_velocity:ratio_out_of_range "<<ratio<<"\n";
    for_m(m){
    for( int b=0; b<u.bc(m).count(); b++ ) {

      if( u.bc(m).type(b) == BndType::outlet() ) {

        Dir d = u.bc(m).direction(b);

        if( m == Comp::u() && d == Dir::imin() )
          for_vjk(u.bc(m).at(b),j,k)
            for(int ii=1; ii<=boil::BW; ii++) {
              if (approx(dom->ibody().fSw(si(m)-1,j,k),0.0) ) {
                u[m][si(m)-ii][j][k] = 0.0;
              }else {
                u[m][si(m)-ii][j][k] = -ubo;
              }
	    }
        if( m == Comp::u() && d == Dir::imax() )
          for_vjk(u.bc(m).at(b),j,k)
            for(int ii=1; ii<=boil::BW; ii++) {
              if (approx(dom->ibody().fSe(ei(m)+1,j,k),0.0) ) {
                u[m][ei(m)+ii][j][k] = 0.0;
              } else {
                u[m][ei(m)+ii][j][k] = +ubo;
              }
            }
        if( m == Comp::v() && d == Dir::jmin() )
          for_vik(u.bc(m).at(b),i,k) 
            for(int jj=1; jj<=boil::BW; jj++) {
              if (approx(dom->ibody().fSs(i,sj(m)-1,k),0.0) ) {
                u[m][i][sj(m)-jj][k] = 0.0;
              } else {
                u[m][i][sj(m)-jj][k] = -ubo;
              }
            }
        if( m == Comp::v() && d == Dir::jmax() )
          for_vik(u.bc(m).at(b),i,k)
            for(int jj=1; jj<=boil::BW; jj++) {
              if (approx(dom->ibody().fSn(i,ej(m)+1,k),0.0) ) {
                u[m][i][ej(m)+jj][k] = 0.0;
              } else {
                u[m][i][ej(m)+jj][k] = +ubo;
              }
            }
        if( m == Comp::w() && d == Dir::kmin() )
          for_vij(u.bc(m).at(b),i,j) 
            for(int kk=1; kk<=boil::BW; kk++) {
              if (approx(dom->ibody().fSb(i,j,sk(m)-1),0.0) ) {
                u[m][i][j][sk(m)-kk] = 0.0;
	      } else {
                u[m][i][j][sk(m)-kk] = -ubo;
              }
            }
        if( m == Comp::w() && d == Dir::kmax() )
          for_vij(u.bc(m).at(b),i,j) 
            for(int kk=1; kk<=boil::BW; kk++) {
              if (approx(dom->ibody().fSt(i,j,ek(m)+1),0.0) ) {
                u[m][i][j][ek(m)+kk] = 0.0;
              } else {
                u[m][i][j][ek(m)+kk] = +ubo;
              }
            }
      }
    } /* b */
    }
  }

  /*-----------------+
  |  ratio in range  |
  +-----------------*/
  else {

    for_m(m) {              
      /* again code duplication, this has been copied from insert_bc */
      for( int b=0; b<u.bc(m).count(); b++ ) {

        if( u.bc(m).type(b) == BndType::outlet() ) {

          Dir d = u.bc(m).direction(b);

          if( d == Dir::imin() ) 
            for_vjk(u.bc(m).at(b),j,k)
              for(int ii=1; ii<=boil::BW; ii++)
                u[m][si(m)-ii][j][k] *= ratio;    
          if( d == Dir::imax() ) 
            for_vjk(u.bc(m).at(b),j,k) 
              for(int ii=1; ii<=boil::BW; ii++)
                u[m][ei(m)+ii][j][k] *= ratio;    

          if( d == Dir::jmin() ) 
            for_vik(u.bc(m).at(b),i,k) 
              for(int jj=1; jj<=boil::BW; jj++)
                u[m][i][sj(m)-jj][k] *= ratio;    
          if( d == Dir::jmax() ) 
            for_vik(u.bc(m).at(b),i,k)
              for(int jj=1; jj<=boil::BW; jj++)
                u[m][i][ej(m)+jj][k] *= ratio;    

          if( d == Dir::kmin() ) 
            for_vij(u.bc(m).at(b),i,j)
              for(int kk=1; kk<=boil::BW; kk++)
                u[m][i][j][sk(m)-kk] *= ratio;    
          if( d == Dir::kmax() ) 
            for_vij(u.bc(m).at(b),i,j)
              for(int kk=1; kk<=boil::BW; kk++)
                u[m][i][j][ek(m)+kk] *= ratio;    
        }
      } /* b */
    } /* m */
  }

  return;
}
