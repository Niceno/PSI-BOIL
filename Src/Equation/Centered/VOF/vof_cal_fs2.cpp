#include "vof.h"
#include <string>
#include <sstream>
#include <iostream>
//#define OUTPUT
using namespace std;

real frontPosition(real xyz1, real xyz2, real phi1, real phi2);

/******************************************************************************/
void VOF::cal_fs2() {
/***************************************************************************//**
*  \brief Calculate free-surface position between cell centers
*     if there is no interface in the cell, yotta (=1e+24) is stored.
*     plane: vm1*x + vm2*y + vm3*z = alpha
*     output: fsx, fsy, fsz
*******************************************************************************/

  // tolerance is necessary because the linear-planes are not closed
  //real tol = 0.0e-3;  

  // initialize
  fsx=boil::yotta;
  fsy=boil::yotta;
  fsz=boil::yotta;

  for_ijk(i,j,k) {
    real cc = phi[i][j][k];
    real cw = phi[i-1][j][k];
    real ce = phi[i+1][j][k];
    real cs = phi[i][j-1][k];
    real cn = phi[i][j+1][k];
    real cb = phi[i][j][k-1];
    real ct = phi[i][j][k+1];
    if ((cw-0.5)*(cc-0.5)<0.0) {
      fsx[i][j][k]=frontPosition(phi.xc(i-1),phi.xc(i),cw,cc);
    }
    if ((cc-0.5)*(ce-0.5)<0.0) {
      fsx[i][j][k]=frontPosition(phi.xc(i),phi.xc(i+1),cc,ce);
    }
    if ((cs-0.5)*(cc-0.5)<0.0) {
      fsy[i][j][k]=frontPosition(phi.yc(j-1),phi.yc(j),cs,cc);
    }
    if ((cc-0.5)*(cn-0.5)<0.0) {
      fsy[i][j][k]=frontPosition(phi.yc(j),phi.yc(j+1),cc,cn);
    }
    if ((cb-0.5)*(cc-0.5)<0.0) {
      fsz[i][j][k]=frontPosition(phi.zc(k-1),phi.zc(k),cb,cc);
    }
    if ((cc-0.5)*(ct-0.5)<0.0) {
      fsz[i][j][k]=frontPosition(phi.zc(k),phi.zc(k+1),cc,ct);
    }
  }


  fsx.exchange_all();
  fsy.exchange_all();
  fsz.exchange_all();

  /* SPECIAL TREATMENT IS NECESSARY FOR WALL BOUNDARY */

#ifdef OUTPUT
  real xpos=0.5, ypos=0.5, zpos=0.5;
  //if (time->current_step()==350) {
    std::string fname_x = name_file("fsx", ".dat", time->current_step(), boil::cart.iam());
    std::string fname_y = name_file("fsy", ".dat", time->current_step(), boil::cart.iam());
    std::string fname_z = name_file("fsz", ".dat", time->current_step(), boil::cart.iam());
    std::ofstream foutx,fouty,foutz;
    foutx.open(fname_x.c_str());
    fouty.open(fname_y.c_str());
    foutz.open(fname_z.c_str());
    foutx<<"VARIABLES=\"X\" \"Y\" \"Z\"\n";
    fouty<<"VARIABLES=\"X\" \"Y\" \"Z\"\n";
    foutz<<"VARIABLES=\"X\" \"Y\" \"Z\"\n";
  
    for_ijk(i,j,k) {
      if (fsx[i][j][k]<boil::zetta) {
        foutx<<fsx[i][j][k]<<" "
             <<phi.yn(j)+xpos*phi.dyc(j)<<" "
             <<phi.zn(k)+zpos*phi.dzc(k)<<"\n";
      }
      if (fsy[i][j][k]<boil::zetta) {
        fouty<<phi.xn(i)+xpos*phi.dxc(i)<<" "
             <<fsy[i][j][k]<<" "
             <<phi.zn(k)+zpos*phi.dzc(k)<<"\n";
      }
      if (fsz[i][j][k]<boil::zetta) {
        foutz<<phi.xn(i)+xpos*phi.dxc(i)<<" "
             <<phi.yn(j)+ypos*phi.dyc(j)<<" "
             <<fsz[i][j][k]<<"\n";
      }
    }
    foutx.close();
    fouty.close();
    foutz.close();
  //}
    exit(0);
#endif
#if 0
  int ii=26, jj=23, kk=22;
  std::cout<<"cal_fs: "<<fsx[ii][jj][kk]<<" "<<fsy[ii][jj][kk]<<" "
           <<fsx[ii][jj][kk]<<"\n";
#endif
  return;
}

real frontPosition(real xyz1, real xyz2, real phi1, real phi2){
   real xyzfront;
   if (phi1 != phi2) {
      xyzfront=xyz1+(0.5-phi1)*(xyz2-xyz1)/(phi2-phi1);
   } else {
      xyzfront=0.5*(xyz1+xyz2);
   }
   return xyzfront;
}

/*-----------------------------------------------------------------------------+
