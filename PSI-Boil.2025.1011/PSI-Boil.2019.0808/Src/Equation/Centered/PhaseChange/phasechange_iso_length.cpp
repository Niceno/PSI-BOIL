#include "phasechange.h"

/******************************************************************************/
real PhaseChange::iso_length(const int i
                           , const int j
                           , const int k
                           , const Dir d) {
/***************************************************************************//**
*  \brief calculate iso-surface length on wall
*******************************************************************************/
  real alen=0.0;
  int i1,j1,k1; //node-ID 1
  int i2,j2,k2; //node-ID 2
  int i3,j3,k3; //node-ID 3
  int i4,j4,k4; //node-ID 4
  real clr1,clr2,clr3,clr4;

  i1=i; j1=j; k1=k;

  if (d == Dir::imin() || d == Dir::imax()) {

    i2=i; j2=j+1; k2=k;
    i3=i; j3=j  ; k3=k+1;
    i4=i; j4=j+1; k4=k+1;
    clr1 = 0.25*(clr[i1][j1-1][k1-1]+clr[i1  ][j1  ][k1-1]
                +clr[i1][j1-1][k1  ]+clr[i1  ][j1  ][k1  ]);
    clr2 = 0.25*(clr[i2][j2-1][k2-1]+clr[i2  ][j2  ][k2-1]
                +clr[i2][j2-1][k2  ]+clr[i2  ][j2  ][k2  ]);
    clr3 = 0.25*(clr[i3][j3-1][k3-1]+clr[i3  ][j3  ][k3-1]
                +clr[i3][j3-1][k3  ]+clr[i3  ][j3  ][k3  ]);
    clr4 = 0.25*(clr[i4][j4-1][k4-1]+clr[i4  ][j4  ][k4-1]
                +clr[i4][j4-1][k4  ]+clr[i4  ][j4  ][k4  ]);

  } else if (d == Dir::jmin() || d == Dir::jmax()) {

    i2=i+1; j2=j  ; k2=k;
    i3=i;   j3=j  ; k3=k+1;
    i4=i+1; j4=j  ; k4=k+1;
    clr1 = 0.25*(clr[i1-1][j1][k1-1]+clr[i1  ][j1][k1-1]
                +clr[i1-1][j1][k1  ]+clr[i1  ][j1][k1  ]);
    clr2 = 0.25*(clr[i2-1][j2][k2-1]+clr[i2  ][j2][k2-1]
                +clr[i2-1][j2][k2  ]+clr[i2  ][j2][k2  ]);
    clr3 = 0.25*(clr[i3-1][j3][k3-1]+clr[i3  ][j3][k3-1]
                +clr[i3-1][j3][k3  ]+clr[i3  ][j3][k3  ]);
    clr4 = 0.25*(clr[i4-1][j4][k4-1]+clr[i4  ][j4][k4-1]
                +clr[i4-1][j4][k4  ]+clr[i4  ][j4][k4  ]);

  } else if (d == Dir::kmin() || d == Dir::kmax()) {

    i2=i+1; j2=j  ; k2=k;
    i3=i  ; j3=j+1; k3=k;
    i4=i+1; j4=j+1; k4=k;
    clr1 = 0.25*(clr[i1-1][j1-1][k1]+clr[i1  ][j1-1][k1]
                +clr[i1-1][j1  ][k1]+clr[i1  ][j1  ][k1]);
    clr2 = 0.25*(clr[i2-1][j2-1][k2]+clr[i2  ][j2-1][k2]
                +clr[i2-1][j2  ][k2]+clr[i2  ][j2  ][k2]);
    clr3 = 0.25*(clr[i3-1][j3-1][k3]+clr[i3  ][j3-1][k3]
                +clr[i3-1][j3  ][k3]+clr[i3  ][j3  ][k3]);
    clr4 = 0.25*(clr[i4-1][j4-1][k4]+clr[i4  ][j4-1][k4]
                +clr[i4-1][j4  ][k4]+clr[i4  ][j4  ][k4]);
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

  if (iflag==0 || iflag==4) return alen;

  int icross=0;  // number of cross point
  real x[4],y[4],z[4];  // cordinates of cross points

  //  1--2
  //  |  |
  //  3--4

  if ((clr1-phisurf)*(clr2-phisurf)<0.0) {
    real x1=clr.xn(i1);
    real y1=clr.yn(j1);
    real z1=clr.zn(k1);
    real x2=clr.xn(i2);
    real y2=clr.yn(j2);
    real z2=clr.zn(k2);
    real coef1=fabs((clr2-phisurf)/(clr2-clr1));
    real coef2=1.0-coef1;
    x[icross] = coef1*x1 + coef2*x2;
    y[icross] = coef1*y1 + coef2*y2;
    z[icross] = coef1*z1 + coef2*z2;
    icross++;
  }

  if ((clr2-phisurf)*(clr4-phisurf)<0.0) {
    real x2=clr.xn(i2);
    real y2=clr.yn(j2);
    real z2=clr.zn(k2);
    real x4=clr.xn(i4);
    real y4=clr.yn(j4);
    real z4=clr.zn(k4);
    real coef1=fabs((clr4-phisurf)/(clr4-clr2));
    real coef2=1.0-coef1;
    x[icross] = coef1*x2 + coef2*x4;
    y[icross] = coef1*y2 + coef2*y4;
    z[icross] = coef1*z2 + coef2*z4;
    icross++;
  }

  if ((clr3-phisurf)*(clr4-phisurf)<0.0) {
    real x3=clr.xn(i3);
    real y3=clr.yn(j3);
    real z3=clr.zn(k3);
    real x4=clr.xn(i4);
    real y4=clr.yn(j4);
    real z4=clr.zn(k4);
    real coef1=fabs((clr4-phisurf)/(clr4-clr3));
    real coef2=1.0-coef1;
    x[icross] = coef1*x3 + coef2*x4;
    y[icross] = coef1*y3 + coef2*y4;
    z[icross] = coef1*z3 + coef2*z4;
    icross++;
  }

  if ((clr1-phisurf)*(clr3-phisurf)<0.0) {
    real x1=clr.xn(i1);
    real y1=clr.yn(j1);
    real z1=clr.zn(k1);
    real x3=clr.xn(i3);
    real y3=clr.yn(j3);
    real z3=clr.zn(k3);
    real coef1=fabs((clr3-phisurf)/(clr3-clr1));
    real coef2=1.0-coef1;
    x[icross] = coef1*x1 + coef2*x3;
    y[icross] = coef1*y1 + coef2*y3;
    z[icross] = coef1*z1 + coef2*z3;
    icross++;
  }

  if (icross==2) {
    alen = sqrt( pow(x[1]-x[0],2.0)
                +pow(y[1]-y[0],2.0)
                +pow(z[1]-z[0],2.0));
  } else if (icross==4) {
    if (clr1>=phisurf ) {
      alen = sqrt( pow(x[1]-x[0],2.0)
                  +pow(y[1]-y[0],2.0)
                  +pow(z[1]-z[0],2.0));
      alen += sqrt( pow(x[3]-x[2],2.0)
                   +pow(y[3]-y[2],2.0)
                   +pow(z[3]-z[2],2.0));
    } else {
      alen = sqrt( pow(x[3]-x[0],2.0)
                  +pow(y[3]-y[0],2.0)
                  +pow(z[3]-z[0],2.0));
      alen += sqrt( pow(x[2]-x[1],2.0)
                   +pow(y[2]-y[1],2.0)
                   +pow(z[2]-z[1],2.0));
    }
  } else {
    std::cout<<"phasechange_iso_length: Error! icross= "<<icross<<"\n";
    std::cout<<"i,j,k,d= "<<i<<" "<<j<<" "<<k<<" "<<d<<"\n";
    std::cout<<"clr= "<<clr1<<" "<<clr2<<" "<<clr3<<" "<<clr4<<"\n";
    std::cout<<"iflag= "<<iflag<<"\n";
    std::cout<<(clr1-phisurf)*(clr2-phisurf)<<" "
	     <<(clr2-phisurf)*(clr4-phisurf)<<" "
	     <<(clr3-phisurf)*(clr4-phisurf)<<" "
	     <<(clr1-phisurf)*(clr3-phisurf)<<"\n";
    exit(0);
  }
  
  return alen;
}
