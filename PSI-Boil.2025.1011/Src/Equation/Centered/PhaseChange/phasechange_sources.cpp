#include "phasechange.h"

/******************************************************************************/
void PhaseChange::sources() {
/***************************************************************************//**
*  \brief calculate source terms
*******************************************************************************/

  for_vijk(tpr,i,j,k){
    real mdotc=phi[i][j][k];
    real volc =dV(i,j,k);
    real dt = time->dt();
    fext[i][j][k] = volc/dt*(1.0/rhov-1.0/rhol)*mdotc;
    clrs[i][j][k]= -(1.0/rhol)*mdotc;
    //tprs[i][j][k] = -volc*mdotc*latent;
  }
  fext.exchange();
  clrs.exchange();

#ifdef IB
  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);

    if(phi[i][j][k]==0.0) continue;

    /* set direction */  // crude code!!!
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
    real volb = dV(i+iof, j+jof, k+kof);
    real mdotc=phi[i][j][k];
    tprs[i+iof][j+jof][k+kof]  -= volb*mdotc*latent;
  }
  tprs.exchange();
#endif

  // sum positive-mdot and negative-mdot
  // mdot: positive for boiling, negative for condensation
  // Note that boiling due to micro-region model is not included.
  smdot_pos = 0.0;
  smdot_neg = 0.0;
  for_vijk(phi,i,j,k){
#ifdef IB
    if(dom->ibody().off(i,j,k))continue;
#endif
    real volc =dV(i,j,k);
    if(phi[i][j][k]>0.0) smdot_pos += phi[i][j][k]*volc;
    if(phi[i][j][k]<0.0) smdot_neg += phi[i][j][k]*volc;
  }
  boil::cart.sum_real(&smdot_pos);
  boil::cart.sum_real(&smdot_neg);

  return;
}
