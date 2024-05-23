#include "nucleation.h"

/***************************************************************************//**
*  approximate volume fraction in a cell (part of sphere)
*******************************************************************************/
real Nucleation::stratified_sphere(const int i, const int j, const int k,
                         const real xcent, const real ycent, const real zcent) {

  //real adist = sqrt( pow(vf->xc(i)-xcent,2.0)
  //                 + pow(vf->yc(j)-ycent,2.0)
  //                 + pow(vf->zc(k)-zcent,2.0))
  //                 - rseed;
  real dist_min=boil::yotta;
  real dist_max=0.0;
  for(int ii=0; ii<=1; ii++)
    for(int jj=0; jj<=1; jj++)
      for(int kk=0; kk<=1; kk++){
        // distance between vertex of cell (i,j,k) and (xcent,ycent,zcent)
        real dist = sqrt ( pow(vf->xn(i+ii)-xcent,2.0)
                         + pow(vf->yn(j+jj)-ycent,2.0)
                         + pow(vf->zn(k+kk)-zcent,2.0));
        if (dist_min>dist) dist_min=dist;
        if (dist_max<dist) dist_max=dist;
      }

  real cseed = 1.0;
  if(dist_max<rseed) {
    // whole cell is inside bubble
    cseed=0.0;
    //std::cout<<"whole inside bubble "<<cseed<<"\n";
  } else if(dist_min>rseed) {
    // whole cell outside of bubble
    cseed=1.0;
    //std::cout<<"whole outside bubble "<<cseed<<"\n";
  } else {
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
    //std::cout<<"stratified_sphere:itmp= "<<itmp<<"\n";
  }

  //if (i==18&&j==18) {
  //  std::cout<<"stratified_sphere:k "<<k<<" "<<dist_max<<" "<<dist_min<<" "<<rseed<<" "<<cseed<<"\n";
  //}

   return cseed;
}
