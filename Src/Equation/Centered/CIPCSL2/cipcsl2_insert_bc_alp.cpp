#include "cipcsl2.h"

/******************************************************************************/
void CIPCSL2::insert_bc_alp(const Scalar & val) {
/***************************************************************************//**
*  \brief Set boundary value for alp.
*         dirichlet boundary condition is neglected.
*         scalar_exchange should take account of periodic condition.
*           1st: insert_bc_alp(phi);
*           2nd: phi.exchange();
*******************************************************************************/

  int i,j,k;

  for( int b=0; b<val.bc().count(); b++ ) {

    if( val.bc().type_decomp(b) ) continue;

    if( val.bc().type(b) == BndType::neumann()
      ||val.bc().type(b) == BndType::symmetry()
      ||val.bc().type(b) == BndType::wall()
      ||val.bc().type(b) == BndType::dirichlet()
      ||val.bc().type(b) == BndType::inlet()
      ||val.bc().type(b) == BndType::outlet()
      ||val.bc().type(b) == BndType::insert() ) {

      int iof=0, jof=0, kof=0;

      Dir d      = val.bc().direction(b);

      /*------------+
      |  direction  |
      +------------*/
      if(d != Dir::undefined()) {

        if(d == Dir::imin()) iof++; if(d == Dir::imax()) iof--;
        if(d == Dir::jmin()) jof++; if(d == Dir::jmax()) jof--;
        if(d == Dir::kmin()) kof++; if(d == Dir::kmax()) kof--;

        for_vijk( val.bc().at(b), i,j,k )
          val[i][j][k]=val[i+iof][j+jof][k+kof];
      }
    }
  }
}
