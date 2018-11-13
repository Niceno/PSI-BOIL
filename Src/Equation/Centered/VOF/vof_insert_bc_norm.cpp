#include "vof.h"

/******************************************************************************/
void VOF::insert_bc_norm() {
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
      ||phi.bc().type(b) == BndType::outlet()
      ||phi.bc().type(b) == BndType::insert() ) {

      int iof=0, jof=0, kof=0;
      Dir d      = phi.bc().direction(b);

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
    } else if( phi.bc().type(b) == BndType::symmetry()
      ||phi.bc().type(b) == BndType::neumann() ) {

      int iof=0, jof=0, kof=0;
      Dir d      = phi.bc().direction(b);

      if(d != Dir::undefined()) {

        if(d == Dir::imin()) iof++; if(d == Dir::imax()) iof--;
        if(d == Dir::jmin()) jof++; if(d == Dir::jmax()) jof--;
        if(d == Dir::kmin()) kof++; if(d == Dir::kmax()) kof--;

        for_vijk( phi.bc().at(b), i,j,k ){
          nx[i][j][k]=nx[i+iof][j+jof][k+kof];
          ny[i][j][k]=ny[i+iof][j+jof][k+kof];
          nz[i][j][k]=nz[i+iof][j+jof][k+kof];
        }

        if(d == Dir::imin()||d == Dir::imax()){
          for_vijk( phi.bc().at(b), i,j,k )
            //nx[i][j][k]=0.0;
            nx[i][j][k]=-nx[i+iof][j+jof][k+kof];
        }
        if(d == Dir::jmin()||d == Dir::jmax()){
          for_vijk( phi.bc().at(b), i,j,k )
            //ny[i][j][k]=0.0;
            ny[i][j][k]=-ny[i+iof][j+jof][k+kof];
        }
        if(d == Dir::kmin()||d == Dir::kmax()){
          for_vijk( phi.bc().at(b), i,j,k )
            //nz[i][j][k]=0.0;
            nz[i][j][k]=-nz[i+iof][j+jof][k+kof];
        }
      }
    }
  }
}

