#include "cipcsl2.h"
using namespace std;

real linear_interpolate(const real c1, const real c2, const real c3, 
                        const real c4); 

/******************************************************************************/
void CIPCSL2::bdcurv_interface() {
/***************************************************************************//**
*  \brief Interpolate curvature at interface in wall adjacent cells.
*     normal vector is not used in this function.
*     output: kappa
*******************************************************************************/
  /* store kappa in fn */
  for_aijk(i,j,k)
    fn[i][j][k]=kappa[i][j][k];

  /* wall boundary */
  for( int b=0; b<phi.bc().count(); b++ ) {
    if(phi.bc().type_decomp(b)) continue;
    if( phi.bc().type(b) == BndType::wall() ) {
      Dir d      = phi.bc().direction(b);

      /*------------+
      |  direction  |
      +------------*/
      if(d != Dir::undefined()) {
        int iof=0, jof=0, kof=0;
        if(d == Dir::imin()) iof++; if(d == Dir::imax()) iof--;
        if(d == Dir::jmin()) jof++; if(d == Dir::jmax()) jof--;
        if(d == Dir::kmin()) kof++; if(d == Dir::kmax()) kof--;

        for_vijk( phi.bc().at(b), i,j,k ){
          int ii=i+iof;
          int jj=j+jof;
          int kk=k+kof;
          //cout<<d<<" "<<ii<<" "<<jj<<" "<<kk<<"\n";
          if (dom->ibody().on(i,j,k)) {
            if (wflag[ii][jj][kk]==0) {

              real kw, ke, ks, kn, kb, kt;
              real dw, de, ds, dn, db, dt;
              kw = ke = ks = kn = kb = kt = 0.0;       // curvature
              dw = de = ds = dn = db = dt = 1.0e+300;  // distance
              real clrw = clr[ii-1][jj][kk];
              real clre = clr[ii+1][jj][kk];
              real clrs = clr[ii][jj-1][kk];
              real clrn = clr[ii][jj+1][kk];
              real clrb = clr[ii][jj][kk-1];
              real clrt = clr[ii][jj][kk+1];
              real clrc = clr[ii][jj][kk];

              if (wflag[ii-1][jj][kk]==0) {
                if ((clrw-phisurf)*(clrc-phisurf)<0.0) {
                  kw = linear_interpolate(clrw, clrc,
                         kappa[ii-1][jj][kk], kappa[ii][jj][kk]);
                  dw = (phisurf-clrc)/(clrw-clrc)*dxw(ii);
                }
              }
              if (wflag[ii+1][jj][kk]==0) {
                if ((clre-phisurf)*(clrc-phisurf)<0.0) {
                  ke = linear_interpolate(clre, clrc,
                         kappa[ii+1][jj][kk], kappa[ii][jj][kk]);
                  de = (phisurf-clrc)/(clre-clrc)*dxe(ii);
                }
              }
              if (wflag[ii][jj-1][kk]==0) {
                if ((clrs-phisurf)*(clrc-phisurf)<0.0) {
                  ks = linear_interpolate(clrs, clrc,
                         kappa[ii][jj-1][kk], kappa[ii][jj][kk]);
                  ds = (phisurf-clrc)/(clrs-clrc)*dys(jj);
                }
              }
              if (wflag[ii][jj+1][kk]==0) {
                if ((clrn-phisurf)*(clrc-phisurf)<0.0) {
                  kn = linear_interpolate(clrn, clrc,
                         kappa[ii][jj+1][kk], kappa[ii][jj][kk]);
                  dn = (phisurf-clrc)/(clrn-clrc)*dyn(jj);
                }
              }
              if (wflag[ii][jj][kk-1]==0) {
                if ((clrb-phisurf)*(clrc-phisurf)<0.0) {
                  kb = linear_interpolate(clrb, clrc,
                         kappa[ii][jj][kk-1], kappa[ii][jj][kk]);
                  db = (phisurf-clrc)/(clrb-clrc)*dzb(kk);
                }
              }
              if (wflag[ii][jj][kk+1]==0) {
                if ((clrt-phisurf)*(clrc-phisurf)<0.0) {
                  kt = linear_interpolate(clrt, clrc,
                         kappa[ii][jj][kk+1], kappa[ii][jj][kk]);
                  dt = (phisurf-clrc)/(clrt-clrc)*dzt(kk);
                }
              }
              real ww, we, ws, wn, wb, wt;
              ww = we = ws = wn = wb = wt = 1.0e+300;
              if (dw!=0.0) ww = 1.0/dw;
              if (de!=0.0) we = 1.0/de;
              if (ds!=0.0) ws = 1.0/ds;
              if (dn!=0.0) wn = 1.0/dn;
              if (db!=0.0) wb = 1.0/db;
              if (dt!=0.0) wt = 1.0/dt;
              real wsum = ww + we + ws + wn + wb + wt;
              fn[ii][jj][kk] = (ww*kw + we*ke + ws*ks + wn*kn + wb*kb + wt*kt)
                             / (wsum+boil::atto);
            }
          }
        }
      }
    }
  }

#if 0
  boil::plot->plot(clr,kappa,fn, "clr-kappa-fn", time->current_step());
  boil::plot->plot(clr,wflag, "clr-wflag", time->current_step());
  exit(0);
#endif

#ifdef IB
  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    dom->ibody().ijk(cc,&i,&j,&k);

    if (wflag[i][j][k]==0) {

      real kw, ke, ks, kn, kb, kt;
      real dw, de, ds, dn, db, dt;
      kw = ke = ks = kn = kb = kt = 0.0;       // curvature
      dw = de = ds = dn = db = dt = 1.0e+300;  // distance
      real clrw = clr[i-1][j][k];
      real clre = clr[i+1][j][k];
      real clrs = clr[i][j-1][k];
      real clrn = clr[i][j+1][k];
      real clrb = clr[i][j][k-1];
      real clrt = clr[i][j][k+1];
      real clrc = clr[i][j][k];

      if (wflag[i-1][j][k]==0) {
        if ((clrw-phisurf)*(clrc-phisurf)<0.0) {
          kw = linear_interpolate(clrw, clrc,
                 kappa[i-1][j][k], kappa[i][j][k]);
          dw = (phisurf-clrc)/(clrw-clrc)*dxw(i);
        }
      }
      if (wflag[i+1][j][k]==0) {
        if ((clre-phisurf)*(clrc-phisurf)<0.0) {
          ke = linear_interpolate(clre, clrc,
                 kappa[i+1][j][k], kappa[i][j][k]);
          de = (phisurf-clrc)/(clre-clrc)*dxe(i);
        }
      }
      if (wflag[i][j-1][k]==0) {
        if ((clrs-phisurf)*(clrc-phisurf)<0.0) {
          ks = linear_interpolate(clrs, clrc,
                 kappa[i][j-1][k], kappa[i][j][k]);
          ds = (phisurf-clrc)/(clrs-clrc)*dys(j);
        }
      }
      if (wflag[i][j+1][k]==0) {
        if ((clrn-phisurf)*(clrc-phisurf)<0.0) {
          kn = linear_interpolate(clrn, clrc,
                 kappa[i][j+1][k], kappa[i][j][k]);
          dn = (phisurf-clrc)/(clrn-clrc)*dyn(j);
        }
      }
      if (wflag[i][j][k-1]==0) {
        if ((clrb-phisurf)*(clrc-phisurf)<0.0) {
          kb = linear_interpolate(clrb, clrc,
                 kappa[i][j][k-1], kappa[i][j][k]);
          db = (phisurf-clrc)/(clrb-clrc)*dzb(k);
        }
      }
      if (wflag[i][j][k+1]==0) {
        if ((clrt-phisurf)*(clrc-phisurf)<0.0) {
          kt = linear_interpolate(clrt, clrc,
                 kappa[i][j][k+1], kappa[i][j][k]);
          dt = (phisurf-clrc)/(clrt-clrc)*dzt(k);
        }
      }

      real ww, we, ws, wn, wb, wt;
      ww = we = ws = wn = wb = wt = 1.0e+300;
      if (dw!=0.0) ww = 1.0/dw;
      if (de!=0.0) we = 1.0/de;
      if (ds!=0.0) ws = 1.0/ds;
      if (dn!=0.0) wn = 1.0/dn;
      if (db!=0.0) wb = 1.0/db;
      if (dt!=0.0) wt = 1.0/dt;
      real wsum = ww + we + ws + wn + wb + wt;
      fn[i][j][k] = (ww*kw + we*ke + ws*ks + wn*kn + wb*kb + wt*kt)
                     / (wsum+boil::atto);
    }
  }
#endif

  /* update */
  for_ijk(i,j,k) {
    if (wflag[i][j][k]==0)
      kappa[i][j][k] = fn[i][j][k];
  }
  insert_bc_kappa(kappa);  // boundary condition
  kappa.exchange_all();

#if 0
  boil::plot->plot(clr,kappa,fn, "clr-kappa-fn", time->current_step());
  boil::plot->plot(clr,wflag, "clr-wflag", time->current_step());
  exit(0);
#endif

  return;
}

/*----------------------------------------------------------------------------*/
real linear_interpolate(const real c1, const real c2, const real k1, 
                        const real k2) {
/***************************************************************************//**
*  \brief Interpolate curvature at interface in wall adjacent cells.
*    |                 |                 |
*    +--------x--------+--------x--------+
*    |        1   |    |        2        |
*                 |    cell face
*                 interface
*******************************************************************************/
  int sign1 = copysign(1.0,k1);
  int sign2 = copysign(1.0,k2);
  int s_sum = sign1 + sign2;

  bool harmonic = false;
  if (abs(s_sum)==2) {
    harmonic = true;
  }

  if (c2==c1) return k1;

  real weight1 = (c2-0.5)/(c2-c1);
  real weight2 = (0.5-c1)/(c2-c1);
  
  if (harmonic) {
    return k1 * k2 / (weight1*k2 + weight2*k1);
  } else {
    return k1 * weight1 + k2 * weight2;
  }
}

/*-----------------------------------------------------------------------------+
 '$Id: cipcsl2_bdcurv_interface.cpp,v 1.2 2018/07/11 14:13:19 guo_w Exp $'/
+-----------------------------------------------------------------------------*/

