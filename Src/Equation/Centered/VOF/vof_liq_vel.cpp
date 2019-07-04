#include "vof.h"

/******************************************************************************/
/* calculate liquid velocity */
/******************************************************************************/
void VOF::cal_liq_vel(Vector * umass, Vector * uliq) {

  for_ijk(i,j,k) {

    /* calculate normal vector  
     * n points to the liquid
     * m is the real space normal vector  */
    real mmx = -mx[i][j][k];
    real mmy = -my[i][j][k];
    real mmz = -mz[i][j][k];
 
    /* cell centre velocity */
    real uxc = 0.5 * ((*umass)[Comp::u()][i][j][k] + (*umass)[Comp::u()][i+1][j][k]);
    real uyc = 0.5 * ((*umass)[Comp::v()][i][j][k] + (*umass)[Comp::v()][i][j+1][k]);
    real uzc = 0.5 * ((*umass)[Comp::w()][i][j][k] + (*umass)[Comp::w()][i][j][k+1]);

    /* normal velocity amplitude */
    real normvel = uxc*mmx+uyc*mmy+uzc*mmz;

    /* tangential velocity vector */
    uxc -= normvel*mmx;
    uyc -= normvel*mmy;
    uzc -= normvel*mmz;

    utx[i][j][k] = uxc;
    uty[i][j][k] = uyc;
    utz[i][j][k] = uzc;

    unliq[i][j][k] = normvel;

    //boil::oout<<i<<" "<<j<<" "<<k<<" "<<mmx<<" "<<mmy<<" "<<mmz<<boil::endl;
  }

  unliq.exchange();

  /* extrapolation flag */
  for_aijk(i,j,k)
    stmp[i][j][k]=0.;
  /* step one: pure liquid cells and interfacial cells */
  for_ijk(i,j,k) {
    if(dom->ibody().off(i,j,k))  continue; /* fix value in solid */
    if(phi[i][j][k]>phisurf&&adens[i][j][k]<boil::pico) {
      stmp[i][j][k]=1.; /* liquid cell without interface */
    } else if(adens[i][j][k]>boil::pico) {
      stmp[i][j][k]=2.; /* interfacial cell */
    }
  }
  stmp.exchange();

  /* step two: neigbours of above */
  for_ijk(i,j,k) {
    if(dom->ibody().off(i,j,k))  continue; /* fix value in solid */
    if(stmp[i][j][k]==0.) {
      if(stmp[i+1][j][k]>0.) {
        stmp[i][j][k]=-2.;
        continue;
      } 
      if(stmp[i-1][j][k]>0.) {
        stmp[i][j][k]=-2.;
        continue;
      } 
      if(stmp[i][j+1][k]>0.) {
        stmp[i][j][k]=-2.;
        continue;
      } 
      if(stmp[i][j-1][k]>0.) {
        stmp[i][j][k]=-2.;
        continue;
      } 
      if(stmp[i][j][k+1]>0.) {
        stmp[i][j][k]=-2.;
        continue;
      } 
      if(stmp[i][j][k-1]>0.) {
        stmp[i][j][k]=-2.;
        continue;
      } 
    }
  }
  stmp.exchange();

#if 0
  if(time->current_step()==4000) {
    for_vi(stmp,i) {
      boil::aout<<i<<" "<<stmp[i][1][1]<<" "<<unliq[i][1][1]<<boil::endl;
    }  
    exit(0);
  }
#endif

  /* stmp = -2, 2 -> extrapolate */
  ext_vel(unliq,stmp,-1);

#if 0
  for_vi(stmp,i) {
    boil::oout<<i<<" "<<stmp[i][1][1]<<" "<<unliq[i][1][1]<<boil::endl;
  }  
  exit(0);
#endif

  /* calculate cell-centred velocity */
  for_ijk(i,j,k) {
    utx[i][j][k] += -mx[i][j][k]*unliq[i][j][k];
    uty[i][j][k] += -my[i][j][k]*unliq[i][j][k];
    utz[i][j][k] += -mz[i][j][k]*unliq[i][j][k];
  }
  utx.exchange();
  uty.exchange();
  utz.exchange();
 
  /* calculate face-centred velocity */
  int ofx(0),ofy(0),ofz(0);
  for_m(m) {
    if       (m==Comp::i()) {
      ofx = -1;
      ofy = 0;
      ofz = 0;
    } else if(m==Comp::j()) {
      ofx = 0;
      ofy = -1;
      ofz = 0;
    } else {
      ofx = 0;
      ofy = 0;
      ofz = -1;
    }

    for_avmijk((*uliq),m,i,j,k) {
        (*uliq)[m][i][j][k] = (*umass)[m][i][j][k];
    }
    for_vmijk((*uliq),m,i,j,k) {
      int ii = i+ofx;
      int jj = j+ofy;
      int kk = k+ofz;

      real flagp = stmp[i ][j ][k ];
      real flagm = stmp[ii][jj][kk];     
 
      if(fabs(flagp)>1.0&&fabs(flagm)>1.0) {
        (*uliq)[m][i][j][k] = -0.5*( ofx*(utx[i][j][k]+utx[ii][jj][kk])
                                   +ofy*(uty[i][j][k]+uty[ii][jj][kk]) 
                                   +ofz*(utz[i][j][k]+utz[ii][jj][kk]) ); 
      } 
      //boil::oout<<"vof liq vel: "<<m<<" "<<i<< " "<<j<<" "<<k<<" "<<(*u)[m][i][j][k]<<" "<<uliq[m][i][j][k]<<" "<<vel_value(m,i,j,k)<<boil::endl;
    }
  }
  uliq->exchange();
#if 0
  if(time->current_step()==4000) {
    boil::plot->plot(uliq,stmp,phi,nx,ny,nz, "un-flag-phi-nx-ny-nz", time->current_step());
    exit(0);
  }
#endif

  return;
}
