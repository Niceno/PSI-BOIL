#include "vof.h"
using namespace std;

/******************************************************************************/
void VOF::norm_young(const Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate normal vector at interface
*  Young's method in E.Aulisa,JCP,225(2007),2301-2319, Section 2.2.2.1
*  Young's method is used by J.Lopez et al., An improved height function...
*  Comput. Methods Appl. Mech. Engrg. 198 (2009) 2555-2564
*         Results: nx, ny, nz : normal vector in the VOF space
*******************************************************************************/
  nx=0.0;
  ny=0.0;
  nz=0.0;

  for_ijk(i,j,k){

    int jflag, sum_jflag=0;
    real nxave=0.0, nyave=0.0, nzave=0.0;

    for (int ii=0; ii<=1; ii++) {
    for (int jj=0; jj<=1; jj++) {
    for (int kk=0; kk<=1; kk++) {
      int iii=i+ii;
      int jjj=j+jj;
      int kkk=k+kk;
      real nxtmp,nytmp,nztmp; // normal at corner
      nxtmp = 0.25 * (sca[iii  ][jjj][kkk  ]+sca[iii  ][jjj-1][kkk  ]
                    + sca[iii  ][jjj][kkk-1]+sca[iii  ][jjj-1][kkk-1])
             -0.25 * (sca[iii-1][jjj][kkk  ]+sca[iii-1][jjj-1][kkk  ]
                    + sca[iii-1][jjj][kkk-1]+sca[iii-1][jjj-1][kkk-1]);
      nxtmp /= (sca.xc(iii)-sca.xc(iii-1));
      nytmp = 0.25 * (sca[iii][jjj  ][kkk  ]+sca[iii-1][jjj  ][kkk  ]
                    + sca[iii][jjj  ][kkk-1]+sca[iii-1][jjj  ][kkk-1])
             -0.25 * (sca[iii][jjj-1][kkk  ]+sca[iii-1][jjj-1][kkk  ]
                    + sca[iii][jjj-1][kkk-1]+sca[iii-1][jjj-1][kkk-1]);
      nytmp /= (sca.yc(jjj)-sca.yc(jjj-1));
      nztmp = 0.25 * (sca[iii][jjj  ][kkk  ]+sca[iii-1][jjj  ][kkk  ]
                    + sca[iii][jjj-1][kkk  ]+sca[iii-1][jjj-1][kkk  ])
             -0.25 * (sca[iii][jjj  ][kkk-1]+sca[iii-1][jjj  ][kkk-1]
                    + sca[iii][jjj-1][kkk-1]+sca[iii-1][jjj-1][kkk-1]);
      nztmp /= (sca.zc(kkk)-sca.zc(kkk-1));

      //normalize(nxtmp,nytmp,nztmp);  // comment out 2019.07.03

      /* check boundary condition */
      jflag=1;
      if (iii  ==si() && iminp==false) jflag=0;
      if (iii-1==ei() && imaxp==false) jflag=0;
      if (jjj  ==sj() && jminp==false) jflag=0;
      if (jjj-1==ej() && jmaxp==false) jflag=0;
      if (kkk  ==sk() && kminp==false) jflag=0;
      if (kkk-1==ek() && kmaxp==false) jflag=0;

      sum_jflag += jflag;

      nxave += real(jflag) * nxtmp;
      nyave += real(jflag) * nytmp;
      nzave += real(jflag) * nztmp;

    }}}

    /* average norm at coner */
    nxave /= sum_jflag;
    nyave /= sum_jflag;
    nzave /= sum_jflag;

#if 0
    if (sum_jflag==0) {
      boil::aout<<"sum_jflag=0: "<<i<<" "<<j<<" "<<k<<"\n";
    }
#endif

    normalize(nxave,nyave,nzave);

    nx[i][j][k]=nxave;
    ny[i][j][k]=nyave;
    nz[i][j][k]=nzave;
  }

  nx.bnd_update();
  ny.bnd_update();
  nz.bnd_update();
  nx.exchange_all();
  ny.exchange_all();
  nz.exchange_all();

  //boil::plot->plot(sca,nx,ny,nz, "clr-nx-ny-nz", time->current_step());
  //exit(0);

  return;
}
