#include "cipcsl2.h"

/******************************************************************************/
void CIPCSL2::insert_bc_gradphi(const Scalar & val) {
/***************************************************************************//**
*  \brief Wall boundary condition for grad(val).
*         val is supporsed to be color function 
*******************************************************************************/

  int i,j,k;

  for( int b=0; b<val.bc().count(); b++ ) {

    if(val.bc().type_decomp(b)) continue;

    if(  val.bc().type(b) == BndType::symmetry()
      || val.bc().type(b) == BndType::neumann() ) {

      int iof=0, jof=0, kof=0;

      Dir d      = val.bc().direction(b);

      /*------------+
      |  direction  |
      +------------*/
      if(d != Dir::undefined()) {
        real ni,nj,nk,magn;
        if(d == Dir::imin() || d == Dir::imax()){
          int ii=si();
          if(d == Dir::imax()) ii=ei()+1;
          int jst=val.bc().at(b).sj();
          int jed=val.bc().at(b).ej()+1;
          int kst=val.bc().at(b).sk();
          int ked=val.bc().at(b).ek()+1;
          for(int j = jst; j<=jed; j++){
          for(int k = kst; k<=ked; k++){
            nx[ii][j][k] = 0.0;
          } }
        }
        if(d == Dir::jmin() || d == Dir::jmax()){
          int jj=sj();
          if(d == Dir::jmax()) jj=ej()+1;
          int ist=val.bc().at(b).si();
          int ied=val.bc().at(b).ei()+1;
          int kst=val.bc().at(b).sk();
          int ked=val.bc().at(b).ek()+1;
          for(int i = ist; i<=ied; i++){
          for(int k = kst; k<=ked; k++){
            ny[i][jj][k] = 0.0;
          } }
        }
        if(d == Dir::kmin() ||d == Dir::kmax()){
          int kk=sk();
          if(d == Dir::kmax()) kk=ek()+1;
          int ist=val.bc().at(b).si();
          int ied=val.bc().at(b).ei()+1;
          int jst=val.bc().at(b).sj();
          int jed=val.bc().at(b).ej()+1;
          for(int i = ist; i<=ied; i++){
          for(int j = jst; j<=jed; j++){
            nz[i][j][kk] = 0.0;
          } }
        }
      }
    }
  }
}

/*-----------------------------------------------------------------------------+
 '$Id: cipcsl2_insert_bc_gradphi.cpp,v 1.3 2015/02/16 08:44:45 sato Exp $'/
+-----------------------------------------------------------------------------*/
