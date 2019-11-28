#include "momentum.h"

/******************************************************************************/
void Momentum::scale_outlet_velocity(const real ubo, const real ratio) {

  Range<real> valid(0.5, 2.0);

  /*---------------------+
  |  ratio out of range  |
  +---------------------*/
  if( !valid.contains(ratio) ) {
    /* bad, code duplication, this has been copied from insert_bc */
    for_m(m){
    for( int b=0; b<u.bc(m).count(); b++ ) {

      if( u.bc(m).type(b) == BndType::outlet() ) {

        Dir d = u.bc(m).direction(b);

        if( m == Comp::u() && d == Dir::imin() )
          for_vjk(u.bc(m).at(b),j,k)
            for(int ii=1; ii<=boil::BW; ii++)
              u[m][si(m)-ii][j][k] = -ubo;
        if( m == Comp::u() && d == Dir::imax() )
          for_vjk(u.bc(m).at(b),j,k)
            for(int ii=1; ii<=boil::BW; ii++)
              u[m][ei(m)+ii][j][k] = +ubo;

        if( m == Comp::v() && d == Dir::jmin() )
          for_vik(u.bc(m).at(b),i,k) 
            for(int jj=1; jj<=boil::BW; jj++)
              u[m][i][sj(m)-jj][k] = -ubo;
        if( m == Comp::v() && d == Dir::jmax() )
          for_vik(u.bc(m).at(b),i,k) 
            for(int jj=1; jj<=boil::BW; jj++)
              u[m][i][ej(m)+jj][k] = +ubo;

        if( m == Comp::w() && d == Dir::kmin() )
          for_vij(u.bc(m).at(b),i,j) 
            for(int kk=1; kk<=boil::BW; kk++)
              u[m][i][j][sk(m)-kk] = -ubo;
        if( m == Comp::w() && d == Dir::kmax() )
          for_vij(u.bc(m).at(b),i,j) 
            for(int kk=1; kk<=boil::BW; kk++)
              u[m][i][j][ek(m)+kk] = +ubo;
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
