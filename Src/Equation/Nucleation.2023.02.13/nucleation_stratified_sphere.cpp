#include "nucleation.h"

/***************************************************************************//**
*  approximate volume fraction in a cell (part of sphere)
*******************************************************************************/
real Nucleation::stratified_sphere(const int i, const int j, const int k,
                         const real xcent, const real ycent, const real zcent) {

  real adist = sqrt( pow(vf->xc(i)-xcent,2.0)
                   + pow(vf->yc(j)-ycent,2.0)
                   + pow(vf->zc(k)-zcent,2.0))
                   - rseed;

  real cseed = 1.0;
  if(adist<-eps) {
     cseed=0.0;
  } else if(adist<eps) {
#if 0
     cseed= 0.5 + adist/(2.0*eps) 
                + 1.0/(2.0*boil::pi)*sin(boil::pi*adist/eps);
#else
     real mm=8;
     real x0=vf->xn(i);
     real y0=vf->yn(j);
     real z0=vf->zn(k);
     real ddx=vf->dxc(i)/real(mm);
     real ddy=vf->dyc(j)/real(mm);
     real ddz=vf->dzc(k)/real(mm);
     int itmp=0;
     for (int ii=0; ii<mm; ii++){
       for (int jj=0; jj<mm; jj++){
         for (int kk=0; kk<mm; kk++){
           real xxc=x0+0.5*ddx+real(ii)*ddx;
           real yyc=y0+0.5*ddy+real(jj)*ddy;
           real zzc=z0+0.5*ddz+real(kk)*ddz;
           real dist=sqrt(pow(xxc-xcent,2.0)
                         +pow(yyc-ycent,2.0)+pow(zzc-zcent,2.0));
           if (dist>rseed){
             itmp=itmp+1;
           }
         }
       }
     }
     cseed=real(itmp)/real(mm*mm*mm);
#endif
   }

   return cseed;
}
