#include "custom.h"

namespace boil {
  real setup_sphere(Scalar & c, const real radius, 
                    const real xcent, const real ycent, const real zcent,
                    const real mm) {  

    int neles(0);
    for_vijk(c,i,j,k) {
      real wsb_x = c.xc(i) - c.dxc(i)*0.5;
      real wst_x = c.xc(i) - c.dxc(i)*0.5;
      real wnb_x = c.xc(i) - c.dxc(i)*0.5;
      real wnt_x = c.xc(i) - c.dxc(i)*0.5;
      real esb_x = c.xc(i) + c.dxc(i)*0.5;
      real est_x = c.xc(i) + c.dxc(i)*0.5;
      real enb_x = c.xc(i) + c.dxc(i)*0.5;
      real ent_x = c.xc(i) + c.dxc(i)*0.5;

      real wsb_y = c.yc(j) - c.dyc(j)*0.5;
      real wst_y = c.yc(j) - c.dyc(j)*0.5;
      real wnb_y = c.yc(j) + c.dyc(j)*0.5;
      real wnt_y = c.yc(j) + c.dyc(j)*0.5;
      real esb_y = c.yc(j) - c.dyc(j)*0.5;
      real est_y = c.yc(j) - c.dyc(j)*0.5;
      real enb_y = c.yc(j) + c.dyc(j)*0.5;
      real ent_y = c.yc(j) + c.dyc(j)*0.5;

      real wsb_z = c.zc(k) - c.dzc(k)*0.5;
      real wst_z = c.zc(k) + c.dzc(k)*0.5;
      real wnb_z = c.zc(k) - c.dzc(k)*0.5;
      real wnt_z = c.zc(k) + c.dzc(k)*0.5;
      real esb_z = c.zc(k) - c.dzc(k)*0.5;
      real est_z = c.zc(k) + c.dzc(k)*0.5;
      real enb_z = c.zc(k) - c.dzc(k)*0.5;
      real ent_z = c.zc(k) + c.dzc(k)*0.5;

      real wsb_dist = sqrt(pow(wsb_x-xcent,2.0)+pow(wsb_y-ycent,2.0)+pow(wsb_z-zcent,2.0));
      real wst_dist = sqrt(pow(wst_x-xcent,2.0)+pow(wst_y-ycent,2.0)+pow(wst_z-zcent,2.0));
      real wnb_dist = sqrt(pow(wnb_x-xcent,2.0)+pow(wnb_y-ycent,2.0)+pow(wnb_z-zcent,2.0));
      real wnt_dist = sqrt(pow(wnt_x-xcent,2.0)+pow(wnt_y-ycent,2.0)+pow(wnt_z-zcent,2.0));
      real esb_dist = sqrt(pow(esb_x-xcent,2.0)+pow(esb_y-ycent,2.0)+pow(esb_z-zcent,2.0));
      real est_dist = sqrt(pow(est_x-xcent,2.0)+pow(est_y-ycent,2.0)+pow(est_z-zcent,2.0));
      real enb_dist = sqrt(pow(enb_x-xcent,2.0)+pow(enb_y-ycent,2.0)+pow(enb_z-zcent,2.0));
      real ent_dist = sqrt(pow(ent_x-xcent,2.0)+pow(ent_y-ycent,2.0)+pow(ent_z-zcent,2.0));

      if(wsb_dist<radius&&wst_dist<radius&&wnb_dist<radius&&wnt_dist<radius&&
         esb_dist<radius&&est_dist<radius&&enb_dist<radius&&ent_dist<radius) {
         c[i][j][k] = 1.0;
      } else if(wsb_dist<=radius||wst_dist<=radius||wnb_dist<=radius||wnt_dist<=radius||
                esb_dist<=radius||est_dist<=radius||enb_dist<=radius||ent_dist<=radius) {
         neles++;
         real x0=c.xn(i);
         real y0=c.yn(j);
         real z0=c.zn(k);
         real ddx=c.dxc(i)/real(mm);
         real ddy=c.dyc(j)/real(mm);
         real ddz=c.dzc(k)/real(mm);
         int itmp=0;
         for (int ii=0; ii<mm; ii++){
           for (int jj=0; jj<mm; jj++){
             for (int kk=0; kk<mm; kk++){
               real xxc=x0+0.5*ddx+real(ii)*ddx;
               real yyc=y0+0.5*ddy+real(jj)*ddy;
               real zzc=z0+0.5*ddz+real(kk)*ddz;
               real dist=sqrt(pow(xxc-xcent,2.0)
                             +pow(yyc-ycent,2.0)+pow(zzc-zcent,2.0));
               if (dist<radius){
                 itmp=itmp+1;
               }
             }
           }
         }
         c[i][j][k]=real(itmp)/real(mm*mm*mm);
      }
    }

    c.exchange_all();
    c.bnd_update();

    return neles;
  }
}

