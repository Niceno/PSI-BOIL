#include "phasechange.h"

/******************************************************************************/
void PhaseChange::sources_tprs() {
/***************************************************************************//**
*  \brief calculate source terms
*   tprs[i][j][k] = -volc*mdotc*latent;
*******************************************************************************/

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
    |  Note:  (i    , j    , k    ) is in fluid domain          |
    |         (i+iof, j+jof, k+kof) is in solid domain          |
    +----------------------------------------------------------*/
    /* solid */
    //real volb = dV(i+iof, j+jof, k+kof);
    real volc = dV(i, j, k);
    real mdotc=phi[i][j][k];
    tprs[i+iof][j+jof][k+kof]  -= volc*mdotc*latent;

    /* liquid */ //correction 
    tprs[i][j][k] -= cpv *(tpr[i][j][k]-tsat)*(1.0/rhov-1.0/rhol)*mdotc*volc;
    
  }
  tprs.exchange();
#endif

  return;
}
