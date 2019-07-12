#include "cipcsl2.h"

/******************************************************************************/
void CIPCSL2::insert_bc_norm() {
/***************************************************************************//**
*  \brief Set boundary value for normal vector of interface
*******************************************************************************/

  Formula form;

  int i,j,k;

  for( int b=0; b<phi.bc().count(); b++ ) {

    if(phi.bc().type_decomp(b)) continue;

    if( phi.bc().type(b) == BndType::wall()
      ||phi.bc().type(b) == BndType::dirichlet()
      ||phi.bc().type(b) == BndType::inlet()
      ||phi.bc().type(b) == BndType::neumann()
      ||phi.bc().type(b) == BndType::outlet()
      ||phi.bc().type(b) == BndType::insert() ) {

      int iof=0, jof=0, kof=0;
      Dir d = phi.bc().direction(b);

      if(d != Dir::undefined()) {

        if(d == Dir::imin()) iof++; if(d == Dir::imax()) iof--;
        if(d == Dir::jmin()) jof++; if(d == Dir::jmax()) jof--;
        if(d == Dir::kmin()) kof++; if(d == Dir::kmax()) kof--;

        for_vijk( phi.bc().at(b), i,j,k ){
          nx[i][j][k]=nx[i+iof][j+jof][k+kof];
          ny[i][j][k]=ny[i+iof][j+jof][k+kof];
          nz[i][j][k]=nz[i+iof][j+jof][k+kof];
        }
      }
    } else if( phi.bc().type(b) == BndType::symmetry() ) {

      Dir d = phi.bc().direction(b);

      if(d != Dir::undefined()) {

        int iof=0, jof=0, kof=0;
        real xsign=1.0, ysign=1.0, zsign=1.0;

        if (d == Dir::imin()) {
          iof++; xsign=-1.0; 
        }
        if (d == Dir::imax()) {
          iof--; xsign=-1.0;
        }
        if (d == Dir::jmin()) {
          jof++; ysign=-1.0;
        }
        if (d == Dir::jmax()) {
          jof--; ysign=-1.0;
        }
        if (d == Dir::kmin()) {
          kof++; zsign=-1.0;
        }
        if (d == Dir::kmax()) {
          kof--; zsign=-1.0;
        }

        for_vijk( phi.bc().at(b), i,j,k ){
          if (i<=si()-2) continue;
          if (i>=ei()+2) continue;
          if (j<=sj()-2) continue;
          if (j>=ej()+2) continue;
          if (k<=sk()-2) continue;
          if (k>=ek()+2) continue;
          nx[i][j][k] = xsign*nx[i+iof][j+jof][k+kof];
          ny[i][j][k] = ysign*ny[i+iof][j+jof][k+kof];
          nz[i][j][k] = zsign*nz[i+iof][j+jof][k+kof];
        }

      }
    }
  }
}
