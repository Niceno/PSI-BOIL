#include "vof.h"
using namespace std;

/******************************************************************************/
void VOF::norm_young(const Scalar & sca) {
/***************************************************************************//**
*  \brief Calculate normal vector at interface
*  Young's method in E.Aulisa,JCP,225(2007),2301-2319, Section 2.2.2.1
*  Young's method is used by J.Lopez et al., An improved height function...
*  Comput. Methods Appl. Mech. Engrg. 198 (2009) 2555-2564
*         Results: mx, my, mz : true normal vector
*                  nx, ny, nz : normal vector in the VOF space
*******************************************************************************/
  mx=0.0;
  my=0.0;
  mz=0.0;
#if 0
  nx=0.0;
  ny=0.0;
  nz=0.0;
#endif

  for_ijk(i,j,k){

    int jflag, sum_jflag=0;
    real mxave=0.0, myave=0.0, mzave=0.0;

    for (int ii=0; ii<=1; ii++) {
    for (int jj=0; jj<=1; jj++) {
    for (int kk=0; kk<=1; kk++) {
      int iii=i+ii;
      int jjj=j+jj;
      int kkk=k+kk;
      real mxtmp,mytmp,mztmp; // normal at corner
      mxtmp = 0.25 * (sca[iii  ][jjj][kkk  ]+sca[iii  ][jjj-1][kkk  ]
                    + sca[iii  ][jjj][kkk-1]+sca[iii  ][jjj-1][kkk-1])
             -0.25 * (sca[iii-1][jjj][kkk  ]+sca[iii-1][jjj-1][kkk  ]
                    + sca[iii-1][jjj][kkk-1]+sca[iii-1][jjj-1][kkk-1]);
      mxtmp /= (sca.xc(iii)-sca.xc(iii-1));
      mytmp = 0.25 * (sca[iii][jjj  ][kkk  ]+sca[iii-1][jjj  ][kkk  ]
                    + sca[iii][jjj  ][kkk-1]+sca[iii-1][jjj  ][kkk-1])
             -0.25 * (sca[iii][jjj-1][kkk  ]+sca[iii-1][jjj-1][kkk  ]
                    + sca[iii][jjj-1][kkk-1]+sca[iii-1][jjj-1][kkk-1]);
      mytmp /= (sca.yc(jjj)-sca.yc(jjj-1));
      mztmp = 0.25 * (sca[iii][jjj  ][kkk  ]+sca[iii-1][jjj  ][kkk  ]
                    + sca[iii][jjj-1][kkk  ]+sca[iii-1][jjj-1][kkk  ])
             -0.25 * (sca[iii][jjj  ][kkk-1]+sca[iii-1][jjj  ][kkk-1]
                    + sca[iii][jjj-1][kkk-1]+sca[iii-1][jjj-1][kkk-1]);
      mztmp /= (sca.zc(kkk)-sca.zc(kkk-1));

      normalize(mxtmp,mytmp,mztmp);

      /* check boundary condition */
      jflag=1;
#if 0
      if (iii  ==si() && iminp==false) jflag=0;
      if (iii-1==ei() && imaxp==false) jflag=0;
      if (jjj  ==sj() && jminp==false) jflag=0;
      if (jjj-1==ej() && jmaxp==false) jflag=0;
      if (kkk  ==sk() && kminp==false) jflag=0;
      if (kkk-1==ek() && kmaxp==false) jflag=0;
#elif 1
      /* planes */
      if (iii  ==si() && iminw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp,-1,0,0);
      if (iii-1==ei() && imaxw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp, 1,0,0);
      if (jjj  ==sj() && jminw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp,0,-1,0);
      if (jjj-1==ej() && jmaxw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp,0, 1,0);
      if (kkk  ==sk() && kminw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp,0,0,-1);
      if (kkk-1==ek() && kmaxw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp,0,0, 1);
      
      /* lines */

      /* imin and ... */
      if (iii  ==si() && iminw==true &&
          jjj  ==sj() && jminw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp,-0.5*sqrt(2.0),-0.5*sqrt(2.0),0);
      if (iii  ==si() && iminw==true &&
          jjj-1==ej() && jmaxw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp,-0.5*sqrt(2.0), 0.5*sqrt(2.0),0);
      if (iii  ==si() && iminw==true &&
          kkk  ==sk() && kminw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp,-0.5*sqrt(2.0),0,-0.5*sqrt(2.0));
      if (iii  ==si() && iminw==true &&
          kkk-1==ek() && kmaxw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp,-0.5*sqrt(2.0),0, 0.5*sqrt(2.0));
      /* imax and ... */
      if (iii-1==ei() && imaxw==true &&
          jjj  ==sj() && jminw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp, 0.5*sqrt(2.0),-0.5*sqrt(2.0),0);
      if (iii-1==ei() && imaxw==true &&
          jjj-1==ej() && jmaxw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp, 0.5*sqrt(2.0), 0.5*sqrt(2.0),0);
      if (iii-1==ei() && imaxw==true &&
          kkk  ==sk() && kminw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp, 0.5*sqrt(2.0),0,-0.5*sqrt(2.0));
      if (iii-1==ei() && imaxw==true &&
          kkk-1==ek() && kmaxw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp, 0.5*sqrt(2.0),0, 0.5*sqrt(2.0));
      /* jmin and k.. */
      if (jjj  ==sj() && jminw==true &&
          kkk  ==sk() && kminw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp,0,-0.5*sqrt(2.0),-0.5*sqrt(2.0));
      if (jjj  ==sj() && jminw==true &&
          kkk-1==ek() && kmaxw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp,0,-0.5*sqrt(2.0), 0.5*sqrt(2.0));
      /* jmax and k.. */
      if (jjj-1==ej() && jmaxw==true &&
          kkk  ==sk() && kminw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp,0, 0.5*sqrt(2.0),-0.5*sqrt(2.0));
      if (jjj-1==ej() && jmaxw==true &&
          kkk-1==ek() && kmaxw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp,0, 0.5*sqrt(2.0), 0.5*sqrt(2.0));

      /* corners */

      /* imin and ... */
      if (iii  ==si() && iminw==true &&
          jjj  ==sj() && jminw==true &&
          kkk  ==sk() && kminw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp,-sqrt(1./3.),-sqrt(1./3.),-sqrt(1./3.));
      if (iii  ==si() && iminw==true &&
          jjj-1==ej() && jmaxw==true &&
          kkk  ==sk() && kminw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp,-sqrt(1./3.), sqrt(1./3.),-sqrt(1./3.));
      if (iii  ==si() && iminw==true &&
          jjj  ==sj() && jminw==true &&
          kkk-1==ek() && kmaxw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp,-sqrt(1./3.),-sqrt(1./3.), sqrt(1./3.));
      if (iii  ==si() && iminw==true &&
          jjj-1==ej() && jmaxw==true &&
          kkk-1==ek() && kmaxw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp,-sqrt(1./3.), sqrt(1./3.), sqrt(1./3.));
      /* imax and ... */
      if (iii-1==ei() && imaxw==true &&
          jjj  ==sj() && jminw==true &&
          kkk  ==sk() && kminw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp, sqrt(1./3.),-sqrt(1./3.),-sqrt(1./3.));
      if (iii-1==ei() && imaxw==true &&
          jjj-1==ej() && jmaxw==true &&
          kkk  ==sk() && kminw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp, sqrt(1./3.), sqrt(1./3.),-sqrt(1./3.));
      if (iii-1==ei() && imaxw==true &&
          jjj  ==sj() && jminw==true &&
          kkk-1==ek() && kmaxw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp, sqrt(1./3.),-sqrt(1./3.), sqrt(1./3.));
      if (iii-1==ei() && imaxw==true &&
          jjj-1==ej() && jmaxw==true &&
          kkk-1==ek() && kmaxw==true)
        wall_adhesion_norm(mxtmp,mytmp,mztmp, sqrt(1./3.), sqrt(1./3.), sqrt(1./3.));
#endif

      sum_jflag += jflag;

      mxave += real(jflag) * mxtmp;
      myave += real(jflag) * mytmp;
      mzave += real(jflag) * mztmp;

    }}}

    /* average norm at coner */
    mxave /= sum_jflag;
    myave /= sum_jflag;
    mzave /= sum_jflag;

#if 0
    if (sum_jflag==0) {
      boil::aout<<"sum_jflag=0: "<<i<<" "<<j<<" "<<k<<"\n";
    }
#endif

    normalize(mxave,myave,mzave);

    mx[i][j][k]=mxave;
    my[i][j][k]=myave;
    mz[i][j][k]=mzave;

  }

  mx.exchange_all();
  my.exchange_all();
  mz.exchange_all();

#if 0
  /* nx, ny, nz is calculated */
  standardized_norm_vect();
#endif

#if 0
  boil::plot->plot(sca,mx,my,mz, "clr-mx-my-mz", time->current_step());
  boil::plot->plot(sca,nx,ny,nz, "clr-nx-ny-nz", time->current_step());
  exit(0);
#endif

  return;
}

/* ancillary function */
void VOF::wall_adhesion_norm(real & nx, real & ny, real & nz,
                             const real nwx, const real nwy, const real nwz) {
  /* n is a normal vector tangent to wall obtained by mirroring clr to wall
   * Brackbill, 1991; p. 341 (section D) */

  /* n_wall = nwl * cos(cangle) + ntan * sin(cangle) */
  real nnx = nwx*cos(cangle) + nx*sin(cangle);
  real nny = nwy*cos(cangle) + ny*sin(cangle);
  real nnz = nwz*cos(cangle) + nz*sin(cangle);

  normalize(nnx,nny,nnz);

  nx = nnx;
  ny = nny;
  nz = nnz;

  return;
}