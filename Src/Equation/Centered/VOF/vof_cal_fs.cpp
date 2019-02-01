#include "vof.h"
#include <string>
#include <sstream>
#include <iostream>
//#define OUTPUT
using namespace std;

/******************************************************************************/
void VOF::cal_fs() {
/***************************************************************************//**
*  \brief Calculate free-surface position between cell centers
*     if there is no interface in the cell, yotta (=1e+24) is stored.
*     plane: vm1*x + vm2*y + vm3*z = alpha
*     output: fsx, fsy, fsz
*******************************************************************************/

  // tolerance is necessary because the linear-planes are not closed
  //real tolf = 0.0e-3;  
  //real tolf = 2.0e-2;  
  real tolf = 5.0e-2;  

  // initialize
  fsx=boil::yotta;
  fsy=boil::yotta;
  fsz=boil::yotta;
  iflagx=0;
  iflagy=0;
  iflagz=0;

  // 0 < xpos, ypos, zpos < 1
  real xpos = 0.5;
  real ypos = 0.5;
  real zpos = 0.5;

  for_ijk(i,j,k) {
    real c = phi[i][j][k];
#if 0
      if (i==26&&j==1&&k==56) {
        cout<<"cal_fs: "<<c<<"\n";
      }
#endif
    if (abs(iflag[i][j][k])>=1 ) {
    //if (abs(iflag[i][j][k])>=0 ) {
      continue;
    }

#if 0
    if (i==26&&j==23&&k==22) {
      cout<<"cal_fs: "<<c<<"\n";
    }
#endif

    if ( c < boil::pico) {
      if (phi[i-1][j][k]>1.0-boil::pico) {
        fsx[i][j][k]=phi.xn(i);
      } 
      if (phi[i+1][j][k]>1.0-boil::pico) {
        fsx[i][j][k]=phi.xn(i+1);
      }
      if (phi[i][j-1][k]>1.0-boil::pico) {
        fsy[i][j][k] = phi.yn(j);
      } 
      if (phi[i][j+1][k]>1.0-boil::pico) {
        fsy[i][j][k] = phi.yn(j+1);
      } 
      if (phi[i][j][k-1]>1.0-boil::pico) {
        fsz[i][j][k] = phi.zn(k);
      } 
      if (phi[i][j][k+1]>1.0-boil::pico) {
        fsz[i][j][k] = phi.zn(k+1);
      } 
    } else if (1.0-boil::pico<c) {
      if (phi[i-1][j][k]<boil::pico) {
        fsx[i][j][k]=phi.xn(i);
      }
      if (phi[i+1][j][k]<boil::pico) {
        fsx[i][j][k]=phi.xn(i+1);
      }
      if (phi[i][j-1][k]<boil::pico) {
        fsy[i][j][k] = phi.yn(j);
      }
      if (phi[i][j+1][k]<boil::pico) {
        fsy[i][j][k] = phi.yn(j+1);
      }
      if (phi[i][j][k-1]<boil::pico) {
        fsz[i][j][k] = phi.zn(k);
      }
      if (phi[i][j][k+1]<boil::pico) {
        fsz[i][j][k] = phi.zn(k+1);
      }
    } else {

#if 0
    if (i==26&&j==1&&k==50) {
      cout<<"cal_fs:B "<<c<<"\n";
    }
#endif
      // calculate vn1, vn2, vn3: normal vector at face center
      real vn1 = -nx[i][j][k];
      real vn2 = -ny[i][j][k];
      real vn3 = -nz[i][j][k];

      real vm1 = fabs(vn1);
      real vm2 = fabs(vn2);
      //real vm3 = fabs(vn3)+boil::pico;
      real vm3 = fabs(vn3);
      //real qa = 1.0/(vm1+vm2+vm3+boil::pico);
      real qa = 1.0/(vm1+vm2+vm3);
      vm1 *= qa;
      vm2 *= qa;
      vm3 *= qa;
      real alpha = calc_alpha(c, vm1, vm2, vm3);

#if 1
      /* special case: xpos=0.5 & ypos=0.5 & zpos=0.5 */
      // cal xuni
      real xuni = (alpha-vm2*0.5-vm3*0.5)/(vm1+boil::pico);
      //real xuni = (alpha-vm2*0.5-vm3*0.5)/vm1;
      if (vn1<0) xuni = 1.0-xuni;

      // cal yuni
      real yuni = (alpha-vm1*0.5-vm3*0.5)/(vm2+boil::pico);
      //real yuni = (alpha-vm1*0.5-vm3*0.5)/vm2;
      if (vn2<0) yuni = 1.0-yuni;

      // cal zuni
      real zuni = (alpha-vm1*0.5-vm2*0.5)/(vm3+boil::pico);
      //real zuni = (alpha-vm1*0.5-vm2*0.5)/vm3;
      if (vn3<0) zuni = 1.0-zuni;
#else
      /* general case: 0 <= xpos <= zpos */
      real xtmp,ytmp,ztmp;
      // cal xuni
      ytmp=ypos;
      ztmp=zpos;
      if (vn2<0) ytmp=1.0-ytmp;
      if (vn3<0) ztmp=1.0-ztmp;
      real xuni = (alpha-vm2*ytmp-vm3*ztmp)/(vm1+boil::pico);
      if (vn1<0) xuni = 1.0-xuni;

      // cal yuni
      xtmp=xpos;
      ztmp=zpos;
      if (vn1<0) xtmp=1.0-xtmp;
      if (vn3<0) ztmp=1.0-ztmp;
      real yuni = (alpha-vm1*xtmp-vm3*ztmp)/(vm2+boil::pico);
      if (vn2<0) yuni = 1.0-yuni;

      // cal zuni
      xtmp=xpos;
      ytmp=ypos;
      if (vn1<0) xtmp=1.0-xtmp;
      if (vn2<0) ytmp=1.0-ytmp;
      real zuni = (alpha-vm1*xtmp-vm2*ytmp)/vm3;
      if (vn3<0) zuni = 1.0-zuni;
#endif
#if 0
      if (i==26&&j==22&&k==22) {
        cout<<"cal_fs:C "<<c<<" "<<alpha<<" "
                 <<xuni<<" "<<yuni<<" "<<zuni<<"\n";
        cout<<"cal_fs:vn "<<vn1<<" "<<vn2<<" "<<vn3<<"\n";
        cout<<"cal_fs:vm "<<vm1<<" "<<vm2<<" "<<vm3<<"\n";
      }
#endif
      /* store data */
      if (0.0-tolf <= xuni && xuni <= 1.0+tolf) {
        //xuni = min(1.0,max(0.0,xuni));
        fsx[i][j][k] = phi.xn(i) + phi.dxc(i) * xuni;
        if (0.0<=xuni && xuni<=1.0) {
          iflagx[i][j][k]=1;         // correctly calculated
        } else if (xuni<0.0) {
          iflagx[i][j][k]=-2;        // torelance in negative
        } else {
          iflagx[i][j][k]= 2;        // torelance in positive 
        }
      }

      if (0.0-tolf <= yuni && yuni <= 1.0+tolf) {
        //yuni = min(1.0,max(0.0,yuni));
        fsy[i][j][k] = phi.yn(j) + phi.dyc(j) * yuni;
        if (0.0<=yuni && yuni<=1.0) {
          iflagy[i][j][k]=1;         // correctly calculated
        } else if (yuni<0.0) {
          iflagy[i][j][k]=-2;        // torelance in negative
        } else {
          iflagy[i][j][k]= 2;        // torelance in positive
        }
      }

      if (0.0-tolf <= zuni && zuni <= 1.0+tolf) {
        //zuni = min(1.0,max(0.0,zuni));
        fsz[i][j][k] = phi.zn(k) + phi.dzc(k) * zuni;
        if (0.0<=zuni && zuni<=1.0) {
          iflagz[i][j][k]=1;         // correctly calculated
        } else if (zuni<0.0) {
          iflagz[i][j][k]=-2;        // torelance in negative
        } else {
          iflagz[i][j][k]= 2;        // torelance in positive
        }
      }
#if 0
      if (i==26&&j==22&&k==22) {
        cout<<"cal_fs:calend "<<fsx[i][j][k]<<" "<<fsy[i][j][k]<<" "
                 <<fsz[i][j][k]<<"\n";
      }
#endif
    }
  }

  fsx.exchange_all();
  fsy.exchange_all();
  fsz.exchange_all();

  /* treatment for trelance */
  for_ijk(i,j,k) {
    if (iflagx[i][j][k]==-2) {
      if (iflagx[i-1][j][k]==1) {
        fsx[i][j][k]=fsx[i-1][j][k];  // i-1 is correctly calculated
      } else if (iflagx[i-1][j][k]==2) {
        fsx[i][j][k]=phi.xn(i);
      } else if (iflagx[i-1][j][k]==0) {
        // use calculated value
      } else {
        cout<<"fs_cal: Error! Trelance treatment west\n";
        cout<<"at"<<boil::cart.iam()<<" "<<i<<" "<<j<<" "<<k<<"\n";
        cout<<"iflag= "<<iflagx[i][j][k]<<" "<<iflagx[i-1][j][k]<<"\n";
      }
    }
    if (iflagx[i][j][k]==2) {
      if (iflagx[i+1][j][k]==1) {
        fsx[i][j][k]=fsx[i+1][j][k];  // i-1 is correctly calculated
      } else if (iflagx[i+1][j][k]==-2) {
        fsx[i][j][k]=phi.xn(i+1);
      } else if (iflagx[i+1][j][k]==0) {
        // use calculated value
      } else {
        cout<<"fs_cal: Error! Trelance treatment east\n";
        cout<<"at "<<boil::cart.iam()<<" "<<i<<" "<<j<<" "<<k<<"\n";
      }
    }
    if (iflagy[i][j][k]==-2) {
      if (iflagy[i][j-1][k]==1) {
        fsy[i][j][k]=fsy[i][j-1][k];  // i-1 is correctly calculated
      } else if (iflagy[i][j-1][k]==2) {
        fsy[i][j][k]=phi.yn(j);
      } else if (iflagy[i][j-1][k]==0) {
        // use calculated value
      } else {
        cout<<"fs_cal: Error! Trelance treatment south\n";
        cout<<"at "<<boil::cart.iam()<<" "<<i<<" "<<j<<" "<<k<<"\n";
      }
    }
    if (iflagy[i][j][k]==2) {
      if (iflagy[i][j+1][k]==1) {
        fsy[i][j][k]=fsy[i][j+1][k];  // i-1 is correctly calculated
      } else if (iflagy[i][j+1][k]==-2) {
        fsy[i][j][k]=phi.yn(j+1);
      } else if (iflagy[i][j+1][k]==0) {
        // use calculated value
      } else {
        cout<<"fs_cal: Error! Trelance treatment north\n";
        cout<<"at "<<boil::cart.iam()<<" "<<i<<" "<<j<<" "<<k<<"\n";
      }
    }
    if (iflagz[i][j][k]==-2) {
      if (iflagz[i][j][k-1]==1) {
        fsz[i][j][k]=fsz[i][j][k-1];  // i-1 is correctly calculated
      } else if (iflagz[i][j][k-1]==2) {
        fsz[i][j][k]=phi.zn(k);
      } else if (iflagz[i][j][k-1]==0) {
        // use calculated value
      } else {
        cout<<"fs_cal: Error! Trelance treatment bottom\n";
        cout<<"at "<<boil::cart.iam()<<" "<<i<<" "<<j<<" "<<k<<"\n";
      }
    }
    if (iflagz[i][j][k]==2) {
      if (iflagz[i][j][k+1]==1) {
        fsz[i][j][k]=fsz[i][j][k+1];  // i-1 is correctly calculated
      } else if (iflagz[i][j][k+1]==-2) {
        fsz[i][j][k]=phi.zn(k+1);
      } else if (iflagz[i][j][k+1]==0) {
        // use calculated value
      } else {
        cout<<"fs_cal: Error! Trelance treatment top\n";
        cout<<"at "<<boil::cart.iam()<<" "<<i<<" "<<j<<" "<<k<<"\n";
      }
    }
  }

#ifdef OUTPUT
  if (time->current_step()==350) {
    string fname_x, fname_y, fname_z;
    fname_x = name_file("fsx",".dat",time->current_step(),boil::cart.iam());
    fname_y = name_file("fsy",".dat",time->current_step(),boil::cart.iam());
    fname_z = name_file("fsz",".dat",time->current_step(),boil::cart.iam());
    ofstream foutx,fouty,foutz;
    foutx.open(fname_x.c_str());
    fouty.open(fname_y.c_str());
    foutz.open(fname_z.c_str());
  
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
  }
#endif
#if 0
  int ii=26, jj=23, kk=22;
  cout<<"cal_fs: "<<fsx[ii][jj][kk]<<" "<<fsy[ii][jj][kk]<<" "
           <<fsx[ii][jj][kk]<<"\n";
#endif
  return;
}
