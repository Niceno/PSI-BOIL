#include "cipcsl2.h"
using namespace std;

real trilinear_interpolate(const Scalar & sca, real xx, real yy, real zz, 
                           int ii, int jj, int kk);

/******************************************************************************/
void CIPCSL2::curv_interface() {
/***************************************************************************//**
*  \brief Calculate at interface cells using interpolation
*     WARNING: iflag must be set before (e.g. in distfunc)
*     output: kappa
*******************************************************************************/
  /* normal vector based on distance function*/
  gradphic(dist);

  for_ijk(i,j,k) {
    if (iflag[i][j][k]==0) {

      /*-----------------------------------------+
      |  nearest interface position (xx, yy, zz) |
      +-----------------------------------------*/
      real dd = dist[i][j][k];
      real xx = phi.xc(i) -  dd * nx[i][j][k];
      real yy = phi.yc(j) -  dd * ny[i][j][k];
      real zz = phi.zc(k) -  dd * nz[i][j][k];

      /*-------------------+
      |  find (ii, jj, kk) |
      +-------------------*/
      /* xc(ii)< xx < xc(ii+1) */
      int ii = -1001;
      for (int itmp=-1; itmp<=0; itmp++) {
         if ((phi.xc(i+itmp)-xx)*(phi.xc(i+itmp+1)-xx)<=0.0) {
           ii = i+itmp;
         }
      }
      /* irregular */
      if (ii<-1000) {
        if( xx<phi.xc(i-1) ) {
          ii = i-1;
          xx = phi.xc(i-1);
        } else if ( xx> phi.xc(i+1)) {
          ii = i;
          xx = phi.xc(i+1);
        }
      }

      /* yc(jj)< yy < yc(jj+1) */
      int jj = -1001;
      for (int jtmp=-1; jtmp<=0; jtmp++) {
         if ((phi.yc(j+jtmp)-yy)*(phi.yc(j+jtmp+1)-yy)<=0.0) {
           jj = j+jtmp;
         }
      }
      /* irregular */
      if (jj<-1000) {
        if( yy<phi.yc(j-1) ) {
          jj = j-1;
          yy = phi.yc(j-1);
        } else if ( yy> phi.yc(j+1)) {
          jj = j;
          yy = phi.yc(j+1);
        }
      }

      /* zc(kk)< zz < zc(kk+1) */
      int kk = -1001;
      for (int ktmp=-1; ktmp<=0; ktmp++) {
         if ((phi.zc(k+ktmp)-zz)*(phi.zc(k+ktmp+1)-zz)<=0.0) {
           kk = k+ktmp;
         }
      }
      /* irregular */
      if (kk<-1000) {
        if( zz<phi.zc(k-1) ) {
          kk = k-1;
          zz = phi.zc(k-1);
        } else if ( zz> phi.zc(k+1)) {
          kk = k;
          zz = phi.zc(k+1);
        }
      }

      /* Error message */
      if (ii+jj+kk<-500) {
        boil::aout<<"# cipcsl2_curv_interface: Error!!!!\n";
        boil::aout<<"# cannot find (ii,jj,kk).\n";
        boil::aout<<"proc: "<<boil::cart.iam()<<" (i,j,k)= ( "<<i<<" , "
                  <<j<<" , "<<k<<" ).\n";
        boil::aout<<"dd= "<<dd<<" (xx,yy,zz)= ( "<<xx<<" , "
                  <<yy<<" , "<<zz<<" ).\n";
        exit(0);
      }

      /*-----------------------------------+
      |  interpolate kappa at (xx, yy, zz) |
      +-----------------------------------*/
      fn[i][j][k] = trilinear_interpolate(kappa, xx, yy, zz, ii, jj, kk);

    }
  }

  /* update */
  for_ijk(i,j,k) {
    if (iflag[i][j][k]==0)
    kappa[i][j][k] = fn[i][j][k];
  }

  return;
}

/*----------------------------------------------------------------------------*/
real trilinear_interpolate(const Scalar & sca, real xx, real yy, real zz, 
                           int ii, int jj, int kk) {

  /* check sign */
  int sa = copysign( 1.0, sca[ii  ][jj  ][kk  ]);
  int sb = copysign( 1.0, sca[ii+1][jj  ][kk  ]);
  int sc = copysign( 1.0, sca[ii  ][jj+1][kk  ]);
  int sd = copysign( 1.0, sca[ii+1][jj+1][kk  ]);
  int se = copysign( 1.0, sca[ii  ][jj  ][kk+1]);
  int sf = copysign( 1.0, sca[ii+1][jj  ][kk+1]);
  int sg = copysign( 1.0, sca[ii  ][jj+1][kk+1]);
  int sh = copysign( 1.0, sca[ii+1][jj+1][kk+1]);

  bool harmonic = false;
  int s_sum = sa + sb + sc + sd + se + sf + sg + sh;
  if ( abs(s_sum)==8 ) {
    harmonic = true;
  }

  /* coordinates and volume */
  real x0 = sca.xc(ii);
  real y0 = sca.yc(jj);
  real z0 = sca.zc(kk);
  real x1 = sca.xc(ii+1);
  real y1 = sca.yc(jj+1);
  real z1 = sca.zc(kk+1);
  real vol = (x1-x0)*(y1-y0)*(z1-z0);

  real va, vb, vc, vd, ve, vf, vg, vh;

  if (harmonic) {
    /* harmonic mean */
    va = vb = vc = vd = ve = vf = vg = vh = 1.0e+300;
    if (sca[ii  ][jj  ][kk  ]!=0.0) va = 1.0/sca[ii  ][jj  ][kk  ];
    if (sca[ii+1][jj  ][kk  ]!=0.0) vb = 1.0/sca[ii+1][jj  ][kk  ]; 
    if (sca[ii  ][jj+1][kk  ]!=0.0) vc = 1.0/sca[ii  ][jj+1][kk  ];
    if (sca[ii+1][jj+1][kk  ]!=0.0) vd = 1.0/sca[ii+1][jj+1][kk  ];
    if (sca[ii  ][jj  ][kk+1]!=0.0) ve = 1.0/sca[ii  ][jj  ][kk+1];
    if (sca[ii+1][jj  ][kk+1]!=0.0) vf = 1.0/sca[ii+1][jj  ][kk+1];
    if (sca[ii  ][jj+1][kk+1]!=0.0) vg = 1.0/sca[ii  ][jj+1][kk+1];
    if (sca[ii+1][jj+1][kk+1]!=0.0) vh = 1.0/sca[ii+1][jj+1][kk+1];
  } else {
    /* arithmetic mean */
    va = sca[ii  ][jj  ][kk  ];
    vb = sca[ii+1][jj  ][kk  ];
    vc = sca[ii  ][jj+1][kk  ];
    vd = sca[ii+1][jj+1][kk  ];
    ve = sca[ii  ][jj  ][kk+1];
    vf = sca[ii+1][jj  ][kk+1];
    vg = sca[ii  ][jj+1][kk+1];
    vh = sca[ii+1][jj+1][kk+1];
  }

  if (vol==0.0) {
    cout<<"curv_interface: ERROR!!! vol = 0\n";
    cout<<"proc= "<<boil::cart.iam()<<" (i,j,k)= "<<ii<<","<<jj<<","<<kk<<"\n";
    exit(0);
  }

  real ca = (x1-xx)*(y1-yy)*(z1-zz)/vol;
  real cb = (xx-x0)*(y1-yy)*(z1-zz)/vol;
  real cc = (x1-xx)*(yy-y0)*(z1-zz)/vol;
  real cd = (xx-x0)*(yy-y0)*(z1-zz)/vol;

  real ce = (x1-xx)*(y1-yy)*(zz-z0)/vol;
  real cf = (xx-x0)*(y1-yy)*(zz-z0)/vol;
  real cg = (x1-xx)*(yy-y0)*(zz-z0)/vol;
  real ch = (xx-x0)*(yy-y0)*(zz-z0)/vol;

  real vv = ca*va + cb*vb + cc*vc + cd*vd + ce*ve + cf*vf + cg*vg + ch*vh;

  if (harmonic) {
    if (vv==0.0) {
      cout<<"curv_interface: ERROR!!! vv = 0\n";
      cout<<"proc= "<<boil::cart.iam()<<"(i,j,k)= "<<ii<<","<<jj<<","<<kk<<"\n";
      exit(0);
    }
    return 1.0/vv;
  } else {
    return vv;
  }
}
