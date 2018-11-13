#include "vof.h"

/******************************************************************************/
void VOF::ext_fs() {
/***************************************************************************//**
*  \brief Extrapolate free-surface position
*     if there is no interface in the cell, 1.0e+300 is stored.
*     output: fsx, fsy, fsz
*******************************************************************************/
  for (int loop=1; loop<=n_ext_fs; loop++) {
    stmp=fsx;
    for_ijk(i,j,k) {
#if 0
      if(i==27&&j==1&&k==56){
        std::cout<<fsx[i][j][k]<<"\n";
      }
#endif
      if (fsx[i][j][k]>boil::zetta) {  // no fsx data
        int ii=0;
        if (fsx[i-1][j][k]<boil::zetta) {
          stmp[i][j][k]=fsx[i-1][j][k];
          ii++;
        }
        if (fsx[i+1][j][k]<boil::zetta) {
          stmp[i][j][k]=fsx[i+1][j][k];
          ii++;
        }
        //if (ii==2) {
        //  stmp[i][j][k]=boil::yotta;
        //}
#if 0
        if(i==27&&j==1&&k==56){
          std::cout<<fsx[i][j][k]<<" "<<fsx[i+1][j][k]<<"\n";
        }
#endif
      }
    }
    fsx=stmp;
    fsx.exchange_all();
  }

  for (int loop=1; loop<=n_ext_fs; loop++) {
    stmp=fsy;
    for_ijk(i,j,k) {
      if (fsy[i][j][k]>boil::zetta) {  // no fsy data
        int jj=0;
        if (fsy[i][j-1][k]<boil::zetta) {
          stmp[i][j][k]=fsy[i][j-1][k];
          jj++;
        }
        if (fsy[i][j+1][k]<boil::zetta) {
          stmp[i][j][k]=fsy[i][j+1][k];
          jj++;
        }
        if (jj==2) {
          stmp[i][j][k]=boil::yotta;
        }
      }
    }
    fsy=stmp;
    fsy.exchange_all();
  }

  for (int loop=1; loop<=n_ext_fs; loop++) {
    stmp=fsz;
    for_ijk(i,j,k) {
      if (fsz[i][j][k]>boil::zetta) {  // no fsz data
        int kk=0;
        if (fsz[i][j][k-1]<boil::zetta) {
          stmp[i][j][k]=fsz[i][j][k-1];
          kk++;
        }
        if (fsz[i][j][k+1]<boil::zetta) {
          stmp[i][j][k]=fsz[i][j][k+1];
          kk++;
        }
        if (kk==2) {
          stmp[i][j][k]=boil::yotta;
        }
      }
    }
    fsz=stmp;
    fsz.exchange_all();
  }


#if 0
  boil::plot->plot(phi,fsx,fsy,fsz, "phi-fsx-fsy-fsz", time->current_step());
  exit(0);
#endif

  return;
}

