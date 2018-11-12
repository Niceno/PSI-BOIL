#include "colorcip.h"
#include <cmath>

/******************************************************************************/void ColorCIP::advance() {
/*-----------------------------------------------+
|  computes convective equation with CIP method  |
+-----------------------------------------------*/

#ifdef USE_TAN
  for_aijk(i,j,k)
     clr[i][j][k]=tan((clr[i][j][k]-0.5)*tanfac*pi);
#endif

  insert_bc(clr);
  clr.exchange_all();

  real phim, phip;
  real umf, upf, vmf, vpf, wmf, wpf;

  const real dt= time->dt();
  real dtdx, dtdy, dtdz, dtdx2, dtdy2, dtdz2;
  real dx, dy, dz;
  real dx2, dx3, dy2, dy3, dz2, dz3;
  real xx, yy, zz, x1, y1, z1;
  real a01, a02, a03, a04, a05, a06, a07, a08, a09, a10
     , a11, a12, a13, a14, a15, a16;
  real b1, b2, b3;
  real cx, cy, cz;
  int im1, jm1, km1;
  Comp mc;

  for_ijk(i,j,k){

    dx=clr.dxc(i);
    dy=clr.dyc(j);
    dz=clr.dzc(k);

    dtdx=dt/dx;
    dtdy=dt/dy;
    dtdz=dt/dz;

    dx2 =dx*dx;
    dx3 =dx2*dx;
    dy2 =dy*dy;
    dy3 =dy2*dy;
    dz2 =dz*dz;
    dz3 =dz2*dz;

    /* Velocity at cell center */
    mc=Comp::u();
    cx=0.5*((*u)[mc][i][j][k]+(*u)[mc][i+1][j][k]);
    mc=Comp::v();
    cy=0.5*((*u)[mc][i][j][k]+(*u)[mc][i][j+1][k]);
    mc=Comp::w();
    cz=0.5*((*u)[mc][i][j][k]+(*u)[mc][i][j][k+1]);

    xx=-cx*dt;
    yy=-cy*dt;
    zz=-cz*dt;

    x1=-copysign(1.0,cx);
    y1=-copysign(1.0,cy);
    z1=-copysign(1.0,cz);

    im1=i+int(x1);
    jm1=j+int(y1);
    km1=k+int(z1);

    a01 = ( ( gpx[i][j][k]+gpx[im1][j][k] )*dx*x1
             +2.0*( clr[i][j][k]-clr[im1][j][k] ) )/(dx3*x1);
    a11 = ( 3.0*( clr[im1][j][k]-clr[i][j][k] )
           - ( gpx[im1][j][k]+2.*gpx[i][j][k] )*dx*x1 )/dx2;
    a02 = ( ( gpy[i][j][k]+gpy[i][jm1][k] )*dy*y1
           +2.0*( clr[i][j][k]-clr[i][jm1][k] ) )/(dy3*y1);
    a12 = ( 3.0*( clr[i][jm1][k]-clr[i][j][k] )
           - ( gpy[i][jm1][k]+2.*gpy[i][j][k] )*dy*y1 )/dy2;
    a03 = ( ( gpz[i][j][k]+gpz[i][j][km1] )*dz*z1
           +2.0*( clr[i][j][k]-clr[i][j][km1] ) )/(dz3*z1);
    a13 = ( 3.0*( clr[i][j][km1]-clr[i][j][k] )
           - ( gpz[i][j][km1]+2.*gpz[i][j][k] )*dz*z1 )/dz2;
    b1  = clr[i][j][k]-clr[im1][j][k]-clr[i][jm1][k]+clr[im1][jm1][k];
    b2  = clr[i][j][k]-clr[i][jm1][k]-clr[i][j][km1]+clr[i][jm1][km1];
    b3  = clr[i][j][k]-clr[im1][j][k]-clr[i][j][km1]+clr[im1][j][km1];
    a05 = ( b1-(-gpy[i][j][k]+gpy[im1][j][k])*dy*y1)/(dx*x1*dy2);
    a04 = ( b1-(-gpx[i][j][k]+gpx[i][jm1][k])*dx*x1)/(dx2*dy*y1);
    a14 = (-b1+(-gpx[i][j][k]+gpx[i][jm1][k])*dx*x1
              +(-gpy[i][j][k]+gpy[im1][j][k])*dy*y1)/(dx*x1*dy*y1);
    a09 = ( b2-(-gpz[i][j][k]+gpz[i][jm1][k])*dz*z1)/(dy*y1*dz2);
    a08 = ( b2-(-gpy[i][j][k]+gpy[i][j][km1])*dy*y1)/(dy2*dz*z1);
    a15 = (-b2+(-gpy[i][j][k]+gpy[i][j][km1])*dy*y1
              +(-gpz[i][j][k]+gpz[i][jm1][k])*dz*z1)/(dy*y1*dz*z1);
    a06 = ( b3-(-gpz[i][j][k]+gpz[im1][j][k])*dz*z1)/(dx*x1*dz2);
    a07 = ( b3-(-gpx[i][j][k]+gpx[i][j][km1])*dx*x1)/(dx2*dz*z1);
    a16 = (-b3+(-gpz[i][j][k]+gpz[im1][j][k])*dz*z1
              +(-gpx[i][j][k]+gpx[i][j][km1])*dx*x1)/(dx*x1*dz*z1);
    a10 = ( -clr[i][j][k] + (clr[im1][j][k]+clr[i][jm1][k]+clr[i][j][km1])
            - (clr[im1][jm1][k]+clr[i][jm1][km1]+clr[im1][j][km1])
            + clr(im1,jm1,km1) ) / (dx*x1*dy*y1*dz*z1);

    clrn[i][j][k] = ((a01*xx+a04*yy+a07*zz+a11)*xx+a14*yy+gpx[i][j][k])*xx
               +((a02*yy+a05*xx+a08*zz+a12)*yy+a15*zz+gpy[i][j][k])*yy
               +((a03*zz+a06*xx+a09*yy+a13)*zz+a16*xx+gpz[i][j][k])*zz
               +a10*xx*yy*zz+clr[i][j][k];

#if 0
    if((i==1||i==20)&&j==1&&k==11){
      std::cout<<"a01="<<a01<<",xx="<<xx<<"\n";
      std::cout<<"cx="<<cx<<"cy="<<cy<<"cz="<<cz<<",dt="<<dt<<"\n";
      std::cout<<"clr[i][j][k]"<<clr[i][j][k]<<"\n";
      std::cout<<"clrn[i][j][k]"<<clrn[i][j][k]<<"\n";
    }
#endif

    gpxn[i][j][k]= (3.*a01*xx+2.*(a04*yy+a07*zz+a11))*xx
               +(a05*yy+a10*zz+a14)*yy+(a06*zz+a16)*zz+gpx[i][j][k];
    gpyn[i][j][k]= (3.*a02*yy+2.*(a05*xx+a08*zz+a12))*yy
               +(a09*zz+a10*xx+a15)*zz+(a04*xx+a14)*xx+gpy[i][j][k];
    gpzn[i][j][k]= (3.*a03*zz+2.*(a06*xx+a09*yy+a13))*zz
              +(a07*xx+a10*yy+a16)*xx+(a08*yy+a15)*yy+gpz[i][j][k];
  }

  insert_bc(clrn);
  clrn.exchange_all();

  gpxn.bnd_grad_update(Comp::i());
  gpyn.bnd_grad_update(Comp::j());
  gpzn.bnd_grad_update(Comp::k());

  gpxn.exchange();
  gpyn.exchange();
  gpzn.exchange();

  real um1,up1,vm1,vp1,wm1,wp1;
  for_ijk(i,j,k){

    dx=clr.dxc(i);
    dy=clr.dyc(j);
    dz=clr.dzc(k);

    dtdx=dt/dx;
    dtdy=dt/dy;
    dtdz=dt/dz;

    dtdx2=dt/(dxw(i)+dxe(i));
    dtdy2=dt/(dys(j)+dyn(j));
    dtdz2=dt/(dzb(k)+dzt(k));

    clr[i][j][k]=clrn[i][j][k];

     mc=Comp::u();
     um1=(*u)[mc][i][j][k];
     up1=(*u)[mc][i+1][j][k];
     mc=Comp::v();
     vm1=0.5*((*u)[mc][i-1][j][k]+(*u)[mc][i-1][j+1][k]);
     vp1=0.5*((*u)[mc][i+1][j][k]+(*u)[mc][i+1][j+1][k]);
     mc=Comp::w();
     wm1=0.5*((*u)[mc][i-1][j][k]+(*u)[mc][i-1][j][k+1]);
     wp1=0.5*((*u)[mc][i+1][j][k]+(*u)[mc][i+1][j][k+1]);
     gpx[i][j][k] = gpxn[i][j][k]
          - (gpxn[i][j][k]*(up1-um1)*dtdx
            +gpyn[i][j][k]*(vp1-vm1)*dtdx2
            +gpzn[i][j][k]*(wp1-wm1)*dtdx2);

     mc=Comp::u();
     um1=0.5*((*u)[mc][i][j-1][k]+(*u)[mc][i+1][j-1][k]);
     up1=0.5*((*u)[mc][i][j+1][k]+(*u)[mc][i+1][j+1][k]);
     mc=Comp::v();
     vm1=(*u)[mc][i][j][k];
     vp1=(*u)[mc][i][j+1][k];
     mc=Comp::w();
     wm1=0.5*((*u)[mc][i][j-1][k]+(*u)[mc][i][j-1][k+1]);
     wp1=0.5*((*u)[mc][i][j+1][k]+(*u)[mc][i][j+1][k+1]);
     gpy[i][j][k] = gpyn[i][j][k]
          - (gpxn[i][j][k]*(up1-um1)*dtdy2
            +gpyn[i][j][k]*(vp1-vm1)*dtdy
            +gpzn[i][j][k]*(wp1-wm1)*dtdy2);

     mc=Comp::u();
     um1=0.5*((*u)[mc][i][j][k-1]+(*u)[mc][i+1][j][k-1]);
     up1=0.5*((*u)[mc][i][j][k+1]+(*u)[mc][i+1][j][k+1]);
     mc=Comp::v();
     vm1=0.5*((*u)[mc][i][j][k-1]+(*u)[mc][i][j+1][k-1]);
     vp1=0.5*((*u)[mc][i][j][k+1]+(*u)[mc][i][j+1][k+1]);
     mc=Comp::w();
     wm1=(*u)[mc][i][j][k];
     wp1=(*u)[mc][i][j][k+1];
     gpz[i][j][k] = gpzn[i][j][k]
          - (gpxn[i][j][k]*(up1-um1)*dtdz2
            +gpyn[i][j][k]*(vp1-vm1)*dtdz2
            +gpzn[i][j][k]*(wp1-wm1)*dtdz);
  }

#ifdef USE_TAN
  for_aijk(i,j,k)
     clr[i][j][k]=atan(clr[i][j][k])/(tanfac*pi)+0.5;
#endif

  insert_bc(clr);
  clr.exchange_all();

  gpx.bnd_grad_update(Comp::i());
  gpy.bnd_grad_update(Comp::j());
  gpz.bnd_grad_update(Comp::k());

  gpx.exchange();
  gpy.exchange();
  gpz.exchange();

  /********************************************************+
  | Final treatment                                        |
  +********************************************************/
  for_aijk(i,j,k)
    fext[i][j][k]=clr[i][j][k];

#if 1
  smooth(clr, phi, 0);
#else
  smooth_ls(clr,clrn,8);
  for_aijk(i,j,k)
    phi[i][j][k]=0.5+0.5*tanh(clrn[i][j][k]/(sqrt(2.0)*ww));
#endif

  return;
}
