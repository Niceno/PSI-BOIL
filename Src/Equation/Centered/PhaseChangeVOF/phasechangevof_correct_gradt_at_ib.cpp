#include "phasechangevof.h"
using namespace std;

static real grad3(real ww, real dm, real dp, real tm, real tc, real tp, real epsl);
/******************************************************************************/
void PhaseChangeVOF::correct_gradt_at_ib(const Scalar * diff_eddy) {
/***************************************************************************//**
*  \brief calculate gradient of temperature near immersed bodies
*******************************************************************************/
  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    /* cell[i][j][k] is wall adjacent cell in fluid domain */
    dom->ibody().ijk(cc,&i,&j,&k);
    real t_c  = tpr[i][j][k];
    real clrc = clr[i][j][k];

    /* x direction */
    Comp m = Comp::u();
    real dx_w = phi.dxw(i);
    real dx_e = phi.dxe(i);
    real t_w  = tpr[i-1][j][k];
    real t_e  = tpr[i+1][j][k];

    /* west is in wall & there is no interface in west */
    //if(dom->ibody().off(m,i,j,k)&&!Interface(-1,m,i,j,k)) {
    if(dom->ibody().off(i-1,j,k)&&!Interface(-1,m,i,j,k)) {
      dx_w *= dom->ibody().fdxw(cc);
#if 0
      real len_s = phi.dxw(i) - 0.5*phi.dxc(i);
      real lam_s = solid()->lambda(i-1,j,k);
      real tmp_s = tpr[i-1][j][k];
      real len_f = 0.5*phi.dxc(i);
      real lam_f;
      if(clrc<clrsurf) {
        lam_f=lambdav;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
      } else {
        lam_f=lambdal;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      }
      real tmp_f = tpr[i][j][k];
      t_w = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);
#else
      t_w = bndtpr[m][i][j][k];
      assert(boil::realistic(t_w));
#endif
      real dtdx = (t_e-t_w)/(dx_w+dx_e);

      /* interface in east */
      if(Interface(+1,m,i,j,k)) {
        real ww = 1.0/dx_e;
        dx_e = distance_x(i,j,k,+1,t_e);
        ww *= dx_e;
	dtdx = grad3(ww, dx_w, dx_e, t_w, t_c, t_e, epsl);
      }

      /* update grad */
      if(clrc>=clrsurf) {
        txl[i][j][k]=dtdx;
      } else {
        txv[i][j][k]=dtdx;
      }
    }

    /* east is in wall & there is no interface in east */
    //if(dom->ibody().off(m,i+1,j,k)&&!Interface(+1,m,i,j,k)) {
    if(dom->ibody().off(i+1,j,k)&&!Interface(+1,m,i,j,k)) {
      dx_e *= dom->ibody().fdxe(cc);
#if 0
      real len_s = phi.dxe(i) - 0.5*phi.dxc(i);
      real lam_s = solid()->lambda(i+1,j,k);
      real tmp_s = tpr[i+1][j][k];
      real len_f = 0.5*phi.dxc(i);
      real lam_f;
      if(clrc<clrsurf) {
        lam_f=lambdav;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
      } else {
        lam_f=lambdal;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      }
      real tmp_f = tpr[i][j][k];
      t_e = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);
#else
      t_e = bndtpr[m][i+1][j][k];
      assert(boil::realistic(t_e));
#endif
      real dtdx = (t_e-t_w)/(dx_w+dx_e);

      /* interface in west */
      if(Interface(-1,m,i,j,k)) {
        real ww = 1.0/dx_w;
        dx_w = distance_x(i,j,k,-1,t_w);
        ww *= dx_w;
	dtdx = grad3(ww, dx_w, dx_e, t_w, t_c, t_e, epsl);
      }

      /* update grad */
      if(clrc>=clrsurf) {
        txl[i][j][k]=dtdx;
      } else {
        txv[i][j][k]=dtdx;
      }
    }

    /* y direction */
    m = Comp::v();
    real dy_s = phi.dys(j);
    real dy_n = phi.dyn(j);
    real t_s  = tpr[i][j-1][k];
    real t_n  = tpr[i][j+1][k];

    /* south is in wall & there is no interface in south */
    //if(dom->ibody().off(m,i,j,k)&&!Interface(-1,m,i,j,k)) {
    if(dom->ibody().off(i,j-1,k)&&!Interface(-1,m,i,j,k)) {
      dy_s *= dom->ibody().fdys(cc);
#if 0
      real len_s = phi.dys(j) - 0.5*phi.dyc(j);
      real lam_s = solid()->lambda(i,j-1,k);
      real tmp_s = tpr[i][j-1][k];
      real len_f = 0.5*phi.dyc(j);
      real lam_f;
      if(clrc<clrsurf) {
        lam_f=lambdav;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
      } else {
        lam_f=lambdal;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      }
      real tmp_f = tpr[i][j][k];
      t_s = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);
#else
      t_s = bndtpr[m][i][j][k];
      assert(boil::realistic(t_s));
#endif
      real dtdy = (t_n-t_s)/(dy_s+dy_n);

      /* interface in north */
      if(Interface(+1,m,i,j,k)) {
        real ww = 1.0/dy_n;
        dy_n = distance_y(i,j,k,+1,t_n);
        ww *= dy_n;
	dtdy = grad3(ww, dy_s, dy_n, t_s, t_c, t_n, epsl);
      }

      /* update grad */
      if(clrc>=clrsurf) {
        tyl[i][j][k]=dtdy;
      } else {
        tyv[i][j][k]=dtdy;
      }
    }

    /* north is in wall & there is no interface in north */
    //if(dom->ibody().off(m,i,j+1,k)&&!Interface(+1,m,i,j,k)) {
    if(dom->ibody().off(i,j+1,k)&&!Interface(+1,m,i,j,k)) {
      dy_n *= dom->ibody().fdyn(cc);
#if 0
      real len_s = phi.dyn(j) - 0.5*phi.dyc(j);
      real lam_s = solid()->lambda(i,j+1,k);
      real tmp_s = tpr[i][j+1][k];
      real len_f = 0.5*phi.dyc(j);
      real lam_f;
      if(clrc<clrsurf) {
        lam_f=lambdav;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
      } else {
        lam_f=lambdal;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      }
      real tmp_f = tpr[i][j][k];
      t_n = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);
#else
      t_n = bndtpr[m][i][j+1][k];
      assert(boil::realistic(t_n));
#endif
      real dtdy = (t_n-t_s)/(dy_s+dy_n);

      /* interface in south */
      if(Interface(-1,m,i,j,k)) {
        real ww = 1.0/dy_s;
        dy_s = distance_y(i,j,k,-1,t_s);
        ww *= dy_s;
	dtdy = grad3(ww, dy_s, dy_n, t_s, t_c, t_n, epsl);
      }

      /* update grad */
      if(clrc>=clrsurf) {
        tyl[i][j][k]=dtdy;
      } else {
        tyv[i][j][k]=dtdy;
      }
    }

    /* z direction */
    m = Comp::w();
    real dz_b = phi.dzb(k);
    real dz_t = phi.dzt(k);
    real t_b  = tpr[i][j][k-1];
    real t_t  = tpr[i][j][k+1];

    /* bottom is in wall & there is no interface in bottom */
    //if(dom->ibody().off(m,i,j,k)&&!Interface(-1,m,i,j,k)) {
    if(dom->ibody().off(i,j,k-1)&&!Interface(-1,m,i,j,k)) {
      dz_b *= dom->ibody().fdzb(cc);
#if 0
      real len_s = phi.dzb(k) - 0.5*phi.dzc(k);
      real lam_s = solid()->lambda(i,j,k-1);
      real tmp_s = tpr[i][j][k-1];
      real len_f = 0.5*phi.dzc(k);
      real lam_f;
      if(clrc<clrsurf) {
        lam_f=lambdav;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
      } else {
        lam_f=lambdal;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      }
      real tmp_f = tpr[i][j][k];
      t_b = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);
#else
      t_b = bndtpr[m][i][j][k];
      assert(boil::realistic(t_b));
#endif
      real dtdz = (t_t-t_b)/(dz_b+dz_t);

      /* interface in top */
      if(Interface(+1,m,i,j,k)) {
        real ww = 1.0/dz_t;
        dz_t = distance_z(i,j,k,+1,t_t);
        ww *= dz_t;
	dtdz = grad3(ww, dz_b, dz_t, t_b, t_c, t_t, epsl);
      }

      /* update grad */
      if(clrc>=clrsurf) {
        tzl[i][j][k]=dtdz;
      } else {
        tzv[i][j][k]=dtdz;
      }
    }

    /* top is in wall & there is no interface in top */
    //if(dom->ibody().off(m,i,j,k+1)&&!Interface(+1,m,i,j,k)) {
    if(dom->ibody().off(i,j,k+1)&&!Interface(+1,m,i,j,k)) {
      dz_t *= dom->ibody().fdzt(cc);
#if 0
      real len_s = phi.dzt(k) - 0.5*phi.dzc(k);
      real lam_s = solid()->lambda(i,j,k+1);
      real tmp_s = tpr[i][j][k+1];
      real len_f = 0.5*phi.dzc(k);
      real lam_f;
      if(clrc<clrsurf) {
        lam_f=lambdav;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
      } else {
        lam_f=lambdal;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      }
      real tmp_f = tpr[i][j][k];
      t_t = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);
#else
      t_t = bndtpr[m][i][j][k+1];
      assert(boil::realistic(t_t));
#endif
      real dtdz = (t_t-t_b)/(dz_b+dz_t);

      /* interface in bottom */
      if(Interface(-1,m,i,j,k)) {
        real ww = 1.0/dz_b;
        dz_b = distance_z(i,j,k,-1,t_b);
        ww *= dz_b;
	dtdz = grad3(ww, dz_b, dz_t, t_b, t_c, t_t, epsl);
      }

      /* update grad */
      if(clrc>=clrsurf) {
        tzl[i][j][k]=dtdz;
      } else {
        tzv[i][j][k]=dtdz;
      }
    }
  }

#if 0
  if(dom->ibody().nccells() > 0) {
    for_ijk(i,j,k) {
      if(dom->ibody().off(i,j,k)) {
        txl[i][j][k]=0.0;
        tyl[i][j][k]=0.0;
        tzl[i][j][k]=0.0;
        txv[i][j][k]=0.0;
        tyv[i][j][k]=0.0;
        tzv[i][j][k]=0.0;
      }
    }
  }
#endif

  return;
}

/******************************************************************************/
real grad3(real ww, real dm, real dp, real tm, real tc, real tp, real epsl){
   real dtdx;
/***************************************************************************//**
*  \brief calculate gradient from three point
*               dm           dp
*         *-------------*------------*
*        tm            tc           tp 
*******************************************************************************/
  if (ww>epsl) {
    real a = dm;
    real b = dp;
    dtdx = b*b*(tc - tm) +a*a*(tp - tc);
    dtdx /= (a*b*(a+b));
   } else {
    dtdx = (tp - tm)/(dm + dp);
   }
  return dtdx;
}

