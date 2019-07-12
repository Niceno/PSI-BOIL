#include "cipcsl2.h"

/******************************************************************************/
void CIPCSL2::insert_bc_dist(Scalar & val) {
/***************************************************************************//**
*  \brief Set boundary value for a distance function, etc.
*         Dirichlet boundary condition is neglected.
*         scalar_exchange(_all) should take account of periodic condition.
*           1st: insert_bc_dist(phi);
*           2nd: phi.exchange_all();
*******************************************************************************/

  int i,j,k;

  for( int b=0; b<val.bc().count(); b++ ) {

    if(val.bc().type_decomp(b)) continue;

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
        //if(ist==si()-1 || jst==sj()-1 || kst ==sk()-1) continue;

        for(int i=ist; i<=ied; i++){
        for(int j=jst; j<=jed; j++){
        for(int k=kst; k<=ked; k++){
          val[i][j][k]=val[i+iof][j+jof][k+kof];
          //do not use extrapolation. it causes trouble in bdcurv!
        } } }

      }
    }
  }
#ifdef IB
  // extrapolate clr  //crude 2014.01.16, Ver1.1.10.23R2
  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);

    /* set direction */
    // (ux,uy,uz) points liquid to solid 
    // crude code!!!
    real ux=dom->ibody().nwx(i,j,k);
    real uy=dom->ibody().nwy(i,j,k);
    real uz=dom->ibody().nwz(i,j,k);
    Dir d = Dir::undefined();
    if (fabs(uz)>0.707) {
      d = Dir::kmin();
      if(uz>0.707){
        d = Dir::kmax();
      }
    } else if (fabs(ux)>0.707) {
      d = Dir::imin();
      if(ux>0.707){
        d = Dir::imax();
      }
    } else if (fabs(uy)>0.707) {
      d = Dir::jmin();
      if(uy>0.707){
        d = Dir::jmax();
      }
    } else {
      std::cout<<"phasechange_micro: Underdevelopment!!!\n";
      exit(0);
    }

    int iof=0, jof=0, kof=0;
    if(d == Dir::imin()) iof--; if(d == Dir::imax()) iof++;
    if(d == Dir::jmin()) jof--; if(d == Dir::jmax()) jof++;
    if(d == Dir::kmin()) kof--; if(d == Dir::kmax()) kof++;

    /*----------------------------------------------------------+
    |  Note:  (i    , j    , k    ) is in computational domain  |
    |         (i+iof, j+jof, k+kof) is on wall                  |
    +----------------------------------------------------------*/
    if(dom->ibody().off(i+iof,j+jof,k+kof)){
      val[i+iof][j+jof][k+kof] = val[i][j][k];
    }
  }
#endif
}
