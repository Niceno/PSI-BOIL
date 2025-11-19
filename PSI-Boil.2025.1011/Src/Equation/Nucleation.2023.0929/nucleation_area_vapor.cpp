#include "nucleation.h"

real area_triangle( real x1, real y1, real z1, real x2, real y2, real z2
                  , real x3, real y3, real z3);
/******************************************************************************/
real Nucleation::area_vapor(const int i
                          , const int j
                          , const int k
                          , const Dir d) {
/***************************************************************************//**
*  \brief calculate iso-surface length on wall
*******************************************************************************/
  real area=0.0;
  real phisurf=0.5;

  int i1,j1,k1; //node-ID 1
  int i2,j2,k2; //node-ID 2
  int i3,j3,k3; //node-ID 3
  int i4,j4,k4; //node-ID 4
  real clrc,clr1,clr2,clr3,clr4;

  clrc=(*clr)[i][j][k];

  i1=i; j1=j; k1=k;

  if (d == Dir::imin() || d == Dir::imax()) {

    i2=i; j2=j+1; k2=k;
    i3=i; j3=j  ; k3=k+1;
    i4=i; j4=j+1; k4=k+1;
    clr1 = 0.25*((*clr)[i1][j1-1][k1-1]+(*clr)[i1  ][j1  ][k1-1]
                +(*clr)[i1][j1-1][k1  ]+(*clr)[i1  ][j1  ][k1  ]);
    clr2 = 0.25*((*clr)[i2][j2-1][k2-1]+(*clr)[i2  ][j2  ][k2-1]
                +(*clr)[i2][j2-1][k2  ]+(*clr)[i2  ][j2  ][k2  ]);
    clr3 = 0.25*((*clr)[i3][j3-1][k3-1]+(*clr)[i3  ][j3  ][k3-1]
                +(*clr)[i3][j3-1][k3  ]+(*clr)[i3  ][j3  ][k3  ]);
    clr4 = 0.25*((*clr)[i4][j4-1][k4-1]+(*clr)[i4  ][j4  ][k4-1]
                +(*clr)[i4][j4-1][k4  ]+(*clr)[i4  ][j4  ][k4  ]);
    area=clr->dSx(i,j,k);

  } else if (d == Dir::jmin() || d == Dir::jmax()) {

    i2=i+1; j2=j  ; k2=k;
    i3=i;   j3=j  ; k3=k+1;
    i4=i+1; j4=j  ; k4=k+1;
    clr1 = 0.25*((*clr)[i1-1][j1][k1-1]+(*clr)[i1  ][j1][k1-1]
                +(*clr)[i1-1][j1][k1  ]+(*clr)[i1  ][j1][k1  ]);
    clr2 = 0.25*((*clr)[i2-1][j2][k2-1]+(*clr)[i2  ][j2][k2-1]
                +(*clr)[i2-1][j2][k2  ]+(*clr)[i2  ][j2][k2  ]);
    clr3 = 0.25*((*clr)[i3-1][j3][k3-1]+(*clr)[i3  ][j3][k3-1]
                +(*clr)[i3-1][j3][k3  ]+(*clr)[i3  ][j3][k3  ]);
    clr4 = 0.25*((*clr)[i4-1][j4][k4-1]+(*clr)[i4  ][j4][k4-1]
                +(*clr)[i4-1][j4][k4  ]+(*clr)[i4  ][j4][k4  ]);
    area=clr->dSy(i,j,k);

  } else if (d == Dir::kmin() || d == Dir::kmax()) {

    i2=i+1; j2=j  ; k2=k;
    i3=i  ; j3=j+1; k3=k;
    i4=i+1; j4=j+1; k4=k;
    clr1 = 0.25*((*clr)[i1-1][j1-1][k1]+(*clr)[i1  ][j1-1][k1]
                +(*clr)[i1-1][j1  ][k1]+(*clr)[i1  ][j1  ][k1]);
    clr2 = 0.25*((*clr)[i2-1][j2-1][k2]+(*clr)[i2  ][j2-1][k2]
                +(*clr)[i2-1][j2  ][k2]+(*clr)[i2  ][j2  ][k2]);
    clr3 = 0.25*((*clr)[i3-1][j3-1][k3]+(*clr)[i3  ][j3-1][k3]
                +(*clr)[i3-1][j3  ][k3]+(*clr)[i3  ][j3  ][k3]);
    clr4 = 0.25*((*clr)[i4-1][j4-1][k4]+(*clr)[i4  ][j4-1][k4]
                +(*clr)[i4-1][j4  ][k4]+(*clr)[i4  ][j4  ][k4]);
    area=clr->dSz(i,j,k);

  }

  real eps_clr=1.0e-8;
  if(fabs(clr1-phisurf)<eps_clr) {
    clr1=phisurf+copysign(1.0,clr1-phisurf)*eps_clr;
  }
  if(fabs(clr2-phisurf)<eps_clr) {
    clr2=phisurf+copysign(1.0,clr2-phisurf)*eps_clr;
  }
  if(fabs(clr3-phisurf)<eps_clr) {
    clr3=phisurf+copysign(1.0,clr3-phisurf)*eps_clr;
  }
  if(fabs(clr4-phisurf)<eps_clr) {
    clr4=phisurf+copysign(1.0,clr4-phisurf)*eps_clr;
  }

  int iflag=0;
  if (clr1 >= phisurf) iflag++;
  if (clr2 >= phisurf) iflag++;
  if (clr3 >= phisurf) iflag++;
  if (clr4 >= phisurf) iflag++;

  if (iflag==0) return area;
  if (iflag==4) return 0.0;


  int icross=0;  // number of cross point
  real x[4],y[4],z[4];  // cordinates of cross points
  int icross_v=0;  // number of cross point
  real x_v[6],y_v[6],z_v[6];  // cordinates of cross points
  int icross_l=0;  // number of cross point
  real x_l[6],y_l[6],z_l[6];  // cordinates of cross points

  //  1--2
  //  |  |
  //  3--4
  
  if (clr1<phisurf) {
    x_v[icross_v]=clr->xn(i1);
    y_v[icross_v]=clr->yn(j1);
    z_v[icross_v]=clr->zn(k1);
    icross_v++;
  } else {
    x_l[icross_l]=clr->xn(i1);
    y_l[icross_l]=clr->yn(j1);
    z_l[icross_l]=clr->zn(k1);
    icross_l++;
  }

  if ((clr1-phisurf)*(clr2-phisurf)<0.0) {
    real x1=clr->xn(i1);
    real y1=clr->yn(j1);
    real z1=clr->zn(k1);
    real x2=clr->xn(i2);
    real y2=clr->yn(j2);
    real z2=clr->zn(k2);
    real coef1=fabs((clr2-phisurf)/(clr2-clr1));
    real coef2=1.0-coef1;
    x[icross] = coef1*x1 + coef2*x2;
    y[icross] = coef1*y1 + coef2*y2;
    z[icross] = coef1*z1 + coef2*z2;
    icross++;
    x_v[icross_v] = coef1*x1 + coef2*x2;
    y_v[icross_v] = coef1*y1 + coef2*y2;
    z_v[icross_v] = coef1*z1 + coef2*z2;
    icross_v++;
    x_l[icross_l] = coef1*x1 + coef2*x2;
    y_l[icross_l] = coef1*y1 + coef2*y2;
    z_l[icross_l] = coef1*z1 + coef2*z2;
    icross_l++;
  }

  if (clr2<phisurf) {
    x_v[icross_v]=clr->xn(i2);
    y_v[icross_v]=clr->yn(j2);
    z_v[icross_v]=clr->zn(k2);
    icross_v++;
  } else {
    x_l[icross_l]=clr->xn(i2);
    y_l[icross_l]=clr->yn(j2);
    z_l[icross_l]=clr->zn(k2);
    icross_l++;
  }

  if ((clr2-phisurf)*(clr4-phisurf)<0.0) {
    real x2=clr->xn(i2);
    real y2=clr->yn(j2);
    real z2=clr->zn(k2);
    real x4=clr->xn(i4);
    real y4=clr->yn(j4);
    real z4=clr->zn(k4);
    real coef1=fabs((clr4-phisurf)/(clr4-clr2));
    real coef2=1.0-coef1;
    x[icross] = coef1*x2 + coef2*x4;
    y[icross] = coef1*y2 + coef2*y4;
    z[icross] = coef1*z2 + coef2*z4;
    icross++;
    x_v[icross_v] = coef1*x2 + coef2*x4;
    y_v[icross_v] = coef1*y2 + coef2*y4;
    z_v[icross_v] = coef1*z2 + coef2*z4;
    icross_v++;
    x_l[icross_l] = coef1*x2 + coef2*x4;
    y_l[icross_l] = coef1*y2 + coef2*y4;
    z_l[icross_l] = coef1*z2 + coef2*z4;
    icross_l++;
  }

  if (clr4<phisurf) {
    x_v[icross_v]=clr->xn(i4);
    y_v[icross_v]=clr->yn(j4);
    z_v[icross_v]=clr->zn(k4);
    icross_v++;
  } else {
    x_l[icross_l]=clr->xn(i4);
    y_l[icross_l]=clr->yn(j4);
    z_l[icross_l]=clr->zn(k4);
    icross_l++;
  }

  if ((clr3-phisurf)*(clr4-phisurf)<0.0) {
    real x3=clr->xn(i3);
    real y3=clr->yn(j3);
    real z3=clr->zn(k3);
    real x4=clr->xn(i4);
    real y4=clr->yn(j4);
    real z4=clr->zn(k4);
    real coef1=fabs((clr4-phisurf)/(clr4-clr3));
    real coef2=1.0-coef1;
    x[icross] = coef1*x3 + coef2*x4;
    y[icross] = coef1*y3 + coef2*y4;
    z[icross] = coef1*z3 + coef2*z4;
    icross++;
    x_v[icross_v] = coef1*x3 + coef2*x4;
    y_v[icross_v] = coef1*y3 + coef2*y4;
    z_v[icross_v] = coef1*z3 + coef2*z4;
    icross_v++;
    x_l[icross_l] = coef1*x3 + coef2*x4;
    y_l[icross_l] = coef1*y3 + coef2*y4;
    z_l[icross_l] = coef1*z3 + coef2*z4;
    icross_l++;
  }

  if (clr3<phisurf) {
    x_v[icross_v]=clr->xn(i3);
    y_v[icross_v]=clr->yn(j3);
    z_v[icross_v]=clr->zn(k3);
    icross_v++;
  } else {
    x_l[icross_l]=clr->xn(i3);
    y_l[icross_l]=clr->yn(j3);
    z_l[icross_l]=clr->zn(k3);
    icross_l++;
  }

  if ((clr1-phisurf)*(clr3-phisurf)<0.0) {
    real x1=clr->xn(i1);
    real y1=clr->yn(j1);
    real z1=clr->zn(k1);
    real x3=clr->xn(i3);
    real y3=clr->yn(j3);
    real z3=clr->zn(k3);
    real coef1=fabs((clr3-phisurf)/(clr3-clr1));
    real coef2=1.0-coef1;
    x[icross] = coef1*x1 + coef2*x3;
    y[icross] = coef1*y1 + coef2*y3;
    z[icross] = coef1*z1 + coef2*z3;
    icross++;
    x_v[icross_v] = coef1*x1 + coef2*x3;
    y_v[icross_v] = coef1*y1 + coef2*y3;
    z_v[icross_v] = coef1*z1 + coef2*z3;
    icross_v++;
    x_l[icross_l] = coef1*x1 + coef2*x3;
    y_l[icross_l] = coef1*y1 + coef2*y3;
    z_l[icross_l] = coef1*z1 + coef2*z3;
    icross_l++;
  }


  real areav=0.0;
  real areal=0.0;
  if (icross_v==3) {
    areav = area_triangle( x_v[0], y_v[0], z_v[0]
                         , x_v[1], y_v[1], z_v[1]
                         , x_v[2], y_v[2], z_v[2]);
    return areav;
  } else if (icross_v==4) {
    areav = area_triangle( x_v[0], y_v[0], z_v[0]
                         , x_v[1], y_v[1], z_v[1]
                         , x_v[2], y_v[2], z_v[2]);
    areav+= area_triangle( x_v[0], y_v[0], z_v[0]
                         , x_v[2], y_v[2], z_v[2]
                         , x_v[3], y_v[3], z_v[3]);
    return areav;
  } else if (icross_v==5) {
    if (icross_l!=3) {
      std::cout<<"phasechange_area_vapor: Error!!!\n ";
    }
    areal = area_triangle( x_l[0], y_l[0], z_l[0]
                         , x_l[1], y_l[1], z_l[1]
                         , x_l[2], y_l[2], z_l[2]);
    return area - areal;
  } else if (icross_v==6) {
    if (clrc<phisurf) {
      //cout<<"center is vapor\n";
      areav = area_triangle( x_v[0], y_v[0], z_v[0]
                           , x_v[1], y_v[1], z_v[1]
                           , x_v[2], y_v[2], z_v[2]);
      areav+= area_triangle( x_v[2], y_v[2], z_v[2]
                           , x_v[3], y_v[3], z_v[3]
                           , x_v[4], y_v[4], z_v[4]);
      areav+= area_triangle( x_v[4], y_v[4], z_v[4]
                           , x_v[5], y_v[5], z_v[5]
                           , x_v[0], y_v[0], z_v[0]);
      areav+= area_triangle( x_v[0], y_v[0], z_v[0]
                           , x_v[2], y_v[2], z_v[2]
                           , x_v[4], y_v[4], z_v[4]);
      //cout<<"icross_v=6 "<<areav<<"\n";
      return areav;
    } else {
      //cout<<"center is liquid\n";
      areal = area_triangle( x_l[0], y_l[0], z_l[0]
                           , x_l[1], y_l[1], z_l[1]
                           , x_l[2], y_l[2], z_l[2]);
      areal+= area_triangle( x_l[2], y_l[2], z_l[2]
                           , x_l[3], y_l[3], z_l[3]
                           , x_l[4], y_l[4], z_l[4]);
      areal+= area_triangle( x_l[4], y_l[4], z_l[4]
                           , x_l[5], y_l[5], z_l[5]
                           , x_l[0], y_l[0], z_l[0]);
      areal+= area_triangle( x_l[0], y_l[0], z_l[0]
                           , x_l[2], y_l[2], z_l[2]
                           , x_l[4], y_l[4], z_l[4]);
      //cout<<"icross_l=6 "<<areal<<"\n";
      return area - areal;
    }
  } else {
    std::cout<<"nucleation_area_vapor:Error!!!\n";
    exit(0);
  }
  return 0;
}

/******************************************************************************/
real area_triangle( real x1, real y1, real z1, real x2, real y2, real z2
                  , real x3, real y3, real z3){
  real x21 = x2 - x1;
  real y21 = y2 - y1;
  real z21 = z2 - z1;
  real x31 = x3 - x1;
  real y31 = y3 - y1;
  real z31 = z3 - z1;
  real outx = y21*z31 - z21*y31;
  real outy = z21*x31 - x21*z31;
  real outz = x21*y31 - y21*x31;
  real area=0.5*sqrt(outx*outx+outy*outy+outz*outz);

  return area;
}
