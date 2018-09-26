#include "colorcip.h"

/******************************************************************************/
void ColorCIP::insert_bc_tanh(const Scalar & val, const Scalar & sca) {
/***************************************************************************//**
*  \brief Boundary condition for a color function of tanh profile.
*         Value next to boundary is calculated from distance function.
*         val: color function of tanh profile
*         sca: distance function
*******************************************************************************/

  int i,j,k;

  for( int b=0; b<val.bc().count(); b++ ) {

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
          val[i+iof][j+jof][k+kof] = tanh(clrn[i+iof][j+jof][k+kof]
                                                     /(sqrt(2.0)*ww));
      }
    }
  }
}
/*-----------------------------------------------------------------------------+
 '$Id: colorcip_insert_bc_tanh.cpp,v 1.2 2015/01/05 17:11:57 sato Exp $'/
+-----------------------------------------------------------------------------*/
