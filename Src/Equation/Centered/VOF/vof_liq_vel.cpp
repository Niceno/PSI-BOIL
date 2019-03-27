#include "vof.h"

/******************************************************************************/
/* calculate liquid velocity */
/******************************************************************************/
void VOF::cal_liq_vel() {
  
  for_ijk(i,j,k) {

    /* calculate normal vector  
     * n points to the liquid
     * m is the real space normal vector  */
    real mmx = -mx[i][j][k];
    real mmy = -my[i][j][k];
    real mmz = -mz[i][j][k];
 
    //boil::oout<<i<<" "<<j<<" "<<k<<" "<<mmx<<" "<<mmy<<" "<<mmz<<boil::endl;

    /* cell centre velocity */
    real uxc = 0.5 * ((*u)[Comp::u()][i][j][k] + (*u)[Comp::u()][i+1][j][k]);
    real uyc = 0.5 * ((*u)[Comp::v()][i][j][k] + (*u)[Comp::v()][i][j+1][k]);
    real uzc = 0.5 * ((*u)[Comp::w()][i][j][k] + (*u)[Comp::w()][i][j][k+1]);

    /* normal velocity amplitude */
    real normvel = uxc*mmx+uyc*mmy+uzc*mmz;

    /* tangential velocity vector */
    uxc -= normvel*mmx;
    uyc -= normvel*mmy;
    uzc -= normvel*mmz;

    utx[i][j][k] = uxc;
    uty[i][j][k] = uyc;
    utz[i][j][k] = uzc;

    /* tangential velocity amplitude */
    real tangvel = sqrt(uxc*uxc+uyc*uyc+uzc*uzc);
 
    unliq[i][j][k] = normvel;
    utliq[i][j][k] = tangvel;
  }

  unliq.exchange();
  utliq.exchange();
#if 0
  utx.exchange();
  uty.exchange();
  utz.exchange();

  /* calculate divergence of tangential velocity */
  for_ijk(i,j,k) {
    divutliq[i][j][k] = 0.5*( (utx[i+1][j][k]-utx[i-1][j][k])/(dxw(i)+dxe(i))
                             +(uty[i][j+1][k]-uty[i][j-1][k])/(dys(j)+dyn(j))
                             +(utz[i][j][k+1]-utz[i][j][k-1])/(dzb(k)+dzt(k))
                            );
  }
  /* it is necessary to give fext neumann/periodic bcs */
  /* however, this might still possibly be incorrect at boundaries */
  divutliq.bnd_update();
  divutliq.exchange();
#endif

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

#if 0
  /* calculate divergence of velocity */
  for_ijk(i,j,k) {
    real uxp = utx[i+1][j][k] + nx[i+1][j][k]*unliq[i+1][j][k];
    real uxm = utx[i-1][j][k] + nx[i-1][j][k]*unliq[i-1][j][k];
    real uyp = uty[i][j+1][k] + ny[i][j+1][k]*unliq[i][j+1][k];
    real uym = uty[i][j-1][k] + ny[i][j-1][k]*unliq[i][j-1][k];
    real uzp = utz[i][j+1][k] + nz[i][j][k+1]*unliq[i][j][k+1];
    real uzm = utz[i][j-1][k] + nz[i][j][k-1]*unliq[i][j][k-1];

    divutliq[i][j][k] = 0.5*( (uxp-uxm)/(dxw(i)+dxe(i))
                             +(uyp-uym)/(dys(j)+dyn(j))
                             +(uzp-uzm)/(dzb(k)+dzt(k))
                            );
  }
  /* it is necessary to give fext neumann/periodic bcs */
  /* however, this might still possibly be incorrect at boundaries */
  divutliq.bnd_update();
  divutliq.exchange();

  boil::plot->plot(unliq,stmp,divutliq, "un-flag-div", time->current_step());
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

    for_avmijk(uliq,m,i,j,k) {
        uliq[m][i][j][k] = (*u)[m][i][j][k];
    }
    for_vmijk(uliq,m,i,j,k) {
      int ii = i+ofx;
      int jj = j+ofy;
      int kk = k+ofz;

      real flagp = stmp[i ][j ][k ];
      real flagm = stmp[ii][jj][kk];     
 
      if(fabs(flagp)>1.0&&fabs(flagm)>1.0) {
        uliq[m][i][j][k] = -0.5*( ofx*(utx[i][j][k]+utx[ii][jj][kk])
                                 +ofy*(uty[i][j][k]+uty[ii][jj][kk]) 
                                 +ofz*(utz[i][j][k]+utz[ii][jj][kk]) ); 
      } 
      //boil::oout<<"vof liq vel: "<<m<<" "<<i<< " "<<j<<" "<<k<<" "<<(*u)[m][i][j][k]<<" "<<uliq[m][i][j][k]<<" "<<vel_value(m,i,j,k)<<boil::endl;
    }
  }
  uliq.exchange();
#if 0
  if(time->current_step()==4000) {
    boil::plot->plot(uliq,stmp,phi,nx,ny,nz, "un-flag-phi-nx-ny-nz", time->current_step());
    exit(0);
  }
#endif

  return;
}
