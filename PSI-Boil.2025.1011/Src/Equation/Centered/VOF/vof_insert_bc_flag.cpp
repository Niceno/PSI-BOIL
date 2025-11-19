#include "vof.h"

/******************************************************************************/
void VOF::insert_bc_flag(ScalarInt & val, const bool adj) {
/***************************************************************************//**
*  \brief Set boundary value for flag
*  if adj == true, flag in wall adjacent cells = -1001
*******************************************************************************/

  int iproc=boil::cart.iam();

  for( int b=0; b<val.bc().count(); b++ ) {

    if(val.bc().type_decomp(b)) continue;

    if( val.bc().type(b) == BndType::neumann()
      ||val.bc().type(b) == BndType::symmetry() 
      ||val.bc().type(b) == BndType::dirichlet()
      ||val.bc().type(b) == BndType::inlet()
      ||val.bc().type(b) == BndType::outlet()
      ||val.bc().type(b) == BndType::pseudo()
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

        int ist=val.bc().at(b).si()-1;
        int ied=val.bc().at(b).ei()+1;
        int jst=val.bc().at(b).sj()-1;
        int jed=val.bc().at(b).ej()+1;
        int kst=val.bc().at(b).sk()-1;
        int ked=val.bc().at(b).ek()+1;
        if(d == Dir::imin() || d == Dir::imax()){
          ist=ied=val.bc().at(b).si();
        }
        if(d == Dir::jmin() || d == Dir::jmax()){
          jst=jed=val.bc().at(b).sj();
        }
        if(d == Dir::kmin() || d == Dir::kmax()){
          kst=ked=val.bc().at(b).sk();
        }
        if(ist==-1 || jst==-1 || kst ==-1) continue;

        //for_vijk( val.bc().at(b), i,j,k ){
        for(int i=ist; i<=ied; i++){
        for(int j=jst; j<=jed; j++){
        for(int k=kst; k<=ked; k++){
          val[i][j][k]=val[i+iof][j+jof][k+kof];
          //do not use extrapolation. it causes trouble in bdcurv!
        } } }

      }
    } else if( val.bc().type(b) == BndType::wall()) {

      int iof=0, jof=0, kof=0;

      Dir d      = val.bc().direction(b);

      /*------------+
      |  direction  |
      +------------*/
      if(d != Dir::undefined()) {

        if(d == Dir::imin()) iof++; if(d == Dir::imax()) iof--;
        if(d == Dir::jmin()) jof++; if(d == Dir::jmax()) jof--;
        if(d == Dir::kmin()) kof++; if(d == Dir::kmax()) kof--;

        int ist=val.bc().at(b).si()-1;
        int ied=val.bc().at(b).ei()+1;
        int jst=val.bc().at(b).sj()-1;
        int jed=val.bc().at(b).ej()+1;
        int kst=val.bc().at(b).sk()-1;
        int ked=val.bc().at(b).ek()+1;
        if(d == Dir::imin() || d == Dir::imax()){
          ist=ied=val.bc().at(b).si();
        }
        if(d == Dir::jmin() || d == Dir::jmax()){
          jst=jed=val.bc().at(b).sj();
        }
        if(d == Dir::kmin() || d == Dir::kmax()){
          kst=ked=val.bc().at(b).sk();
        }

        //for_vijk( val.bc().at(b), i,j,k ){
        if(ist!=-1 && jst!=-1 && kst!=-1) {
          for(int i=ist; i<=ied; i++){
            for(int j=jst; j<=jed; j++){
              for(int k=kst; k<=ked; k++){
                if (adj) {
                  val[i][j][k]=-1001;
                  val[i+iof][j+jof][k+kof]=-1001;
                } else {
                  val[i][j][k]=-1001;
                  //val[i][j][k]=val[i+iof][j+jof][k+kof];
                }
              }
            }
          }
        }
      }
    }
  }
}
