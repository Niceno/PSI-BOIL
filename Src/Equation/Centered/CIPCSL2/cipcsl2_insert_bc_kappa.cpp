#include "cipcsl2.h"

/******************************************************************************/
void CIPCSL2::insert_bc_kappa(Scalar & val) {
/***************************************************************************//**
*  \brief Set boundary value for curvature except wall boundary
*******************************************************************************/

  for_aijk(i,j,k) {
    if(abs(iflag[i][j][k])>=7)  // crude code: iflag should be computed before
      val[i][j][k]=0.0;
  }

  int i,j,k;

  for( int b=0; b<val.bc().count(); b++ ) {

    if( val.bc().type_decomp(b) ) continue;

    if( val.bc().type(b) == BndType::neumann()
      ||val.bc().type(b) == BndType::symmetry()
      //||val.bc().type(b) == BndType::wall()
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

        // set loop range  (necessary for corner)
        int ist,ied,jst,jed,kst,ked;
        if(val.bc().at(b).si()==val.bc().at(b).ei()){
          ist=ied=val.bc().at(b).si();
        } else {
          if(val.bc().at(b).si()==si()){ist=si()-1;}
          else{ist=val.bc().at(b).si();}
          if(val.bc().at(b).ei()==ei()){ied=ei()+1;}
          else{ied=val.bc().at(b).ei();}
        }
        if(val.bc().at(b).sj()==val.bc().at(b).ej()){
          jst=jed=val.bc().at(b).sj();
        } else {
          if(val.bc().at(b).sj()==sj()){jst=sj()-1;}
          else{jst=val.bc().at(b).sj();}
          if(val.bc().at(b).ej()==ej()){jed=ej()+1;}
          else{jed=val.bc().at(b).ej();}
        }
        if(val.bc().at(b).sk()==val.bc().at(b).ek()){
          kst=ked=val.bc().at(b).sk();
        } else {
          if(val.bc().at(b).sk()==sk()){kst=sk()-1;}
          else{kst=val.bc().at(b).sk();}
          if(val.bc().at(b).ek()==ek()){ked=ek()+1;}
          else{ked=val.bc().at(b).ek();}
        }

        //for_vijk( val.bc().at(b), i,j,k )
        for(i=ist; i<=ied; i++)
          for(j=jst; j<=jed; j++)
            for(k=kst; k<=ked; k++) {
          val[i][j][k]=val[i+iof][j+jof][k+kof];
        }
      }
    }
  }
}
/*-----------------------------------------------------------------------------+
 '$Id: cipcsl2_insert_bc_kappa.cpp,v 1.4 2015/05/06 07:52:09 sato Exp $'/
+-----------------------------------------------------------------------------*/
