#include "levelset.h"

/******************************************************************************/
void LevelSet::insert_bc_kappa2(Scalar & val) {
/***************************************************************************//**
*  \brief Set boundary value for curvature except wall boundary
*******************************************************************************/

  int i,j,k;

  for( int b=0; b<val.bc().count(); b++ ) {

    if( val.bc().type_decomp(b) ) continue;

    if( val.bc().type(b) == BndType::neumann()
      ||val.bc().type(b) == BndType::symmetry()
      ||val.bc().type(b) == BndType::dirichlet()
      ||val.bc().type(b) == BndType::inlet()
      ||val.bc().type(b) == BndType::insert()
      ||val.bc().type(b) == BndType::outlet() ) {

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
/*-----------------------------------------------------------------------------+
 '$Id: levelset_insert_bc_kappa2.cpp,v 1.3 2015/01/05 17:16:38 sato Exp $'/
+-----------------------------------------------------------------------------*/
