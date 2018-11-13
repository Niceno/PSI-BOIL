#include "vof.h"

/******************************************************************************/
void VOF::calc_sxyz() {

  for(int i=1; i<=ei()+1; i++){
  for(int j=1; j<=ej()+1; j++){
  for(int k=1; k<=ek()+1; k++){
    real dy=phi.dyc(j);
    real dz=phi.dzc(k);
    /* initial condition for face density */
    sxyz[Comp::i()][i][j][k]=(phi[i][j][k]+phi[i-1][j][k])/2.0*dz*dy;
  }}}

  for(int i=1; i<=ei()+1; i++){
  for(int j=1; j<=ej()+1; j++){
  for(int k=1; k<=ek()+1; k++){
    real dx=phi.dxc(i);
    real dz=phi.dzc(k);
    /* initial condition for face density */
    sxyz[Comp::j()][i][j][k]=(phi[i][j][k]+phi[i][j-1][k])/2.0*dz*dx;
  }}}

  for(int i=1; i<=ei()+1; i++){
  for(int j=1; j<=ej()+1; j++){
  for(int k=1; k<=ek()+1; k++){
    real dx=phi.dxc(i);
    real dy=phi.dyc(j);
    /* initial condition for face density */
    sxyz[Comp::k()][i][j][k]=(phi[i][j][k]+phi[i][j][k-1])/2.0*dx*dy;
  }}}

  for_m(m)
    bdsxyz(sxyz,m,phi);

#ifdef IB
  if(phi.domain()->ibody().ncall() == 0);
  else {
    for(int cc=0; cc<dom->ibody().nccells(); cc++){
      int i,j,k;
      dom->ibody().ijk(cc,&i,&j,&k);

      /* set direction */    // crude code!!!
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
        std::cout<<"cipcsl2_ib_insert_bc: Underdevelopment!!!\n";
        exit(0);
      }

      real area;
      Comp m;
      if(d == Dir::imin() || d == Dir::imax()) {
        m = Comp::u();
        area = phi.dyc(j)*phi.dzc(k);
      } else if(d == Dir::jmin() || d == Dir::jmax()) {
        m = Comp::v();
        area = phi.dxc(i)*phi.dzc(k);
      } else {
        m = Comp::w();
        area = phi.dxc(i)*phi.dyc(j);
      }

      int iof=0, jof=0, kof=0;
      if(d == Dir::imin()) iof--; if(d == Dir::imax()) iof++;
      if(d == Dir::jmin()) jof--; if(d == Dir::jmax()) jof++;
      if(d == Dir::kmin()) kof--; if(d == Dir::kmax()) kof++;

      /*----------------------------------------------------------+
      |  Note:  (i    , j    , k    ) is in fluid domain  |
      |         (i+iof, j+jof, k+kof) is in wall                  |
      +----------------------------------------------------------*/
      if(dom->ibody().off(i+iof,j+jof,k+kof)){
        iof=max(0,iof);
        jof=max(0,jof);
        kof=max(0,kof);
        sxyz[m][i+iof][j+jof][k+kof] = phi[i][j][k]*area;
      }
    }
  }
#endif

  sxyz.exchange_all();

  return;
}
