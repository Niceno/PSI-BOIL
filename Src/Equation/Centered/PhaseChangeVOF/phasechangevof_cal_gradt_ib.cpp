#include "phasechangevof.h"

/******************************************************************************/
void PhaseChangeVOF::cal_gradt_ib(const Scalar * diff_eddy) {
/***************************************************************************//**
*  \brief calculate fluid temperature gradient in immersed body 
*******************************************************************************/

  /* normal cell, 2nd order */
  for_vijk(tpr,i,j,k){

    if(dom->ibody().on(i,j,k)) continue;

    real dtdx,dtdy,dtdz;
    tpr.grad(i,j,k,&dtdx,&dtdy,&dtdz);

    real factl = solid()->lambda(i,j,k)/lambdal;
    real factv = solid()->lambda(i,j,k)/lambdav;

    txl[i][j][k]=dtdx*factl;
    tyl[i][j][k]=dtdy*factl;
    tzl[i][j][k]=dtdz*factl;
    txv[i][j][k]=dtdx*factv;
    tyv[i][j][k]=dtdy*factv;
    tzv[i][j][k]=dtdz*factv;

  }

  /* near interfaces */
  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    /* cell[i][j][k] is wall adjacent cell in fluid domain */
    dom->ibody().ijk(cc,&i,&j,&k);
    real clrc = clr[i][j][k];

    /* x direction */
    Comp m = Comp::u();

    /* west is in wall & there is interface in west */
    if(dom->ibody().off(i-1,j,k)&&Interface(-1,m,i,j,k)) {
      real len_s = 0.5*phi.dxc(i-1);
      real lam_s = solid()->lambda(i-1,j,k);
      real tmp_s = tpr[i-1][j][k];

      real tmp_f;
      real len_f = 0.5*phi.dxc(i)-distance_x(i,j,k,-1,tmp_f);
      real lam_f;
      /* note the inversion */
      if(clrc<clrsurf) {
        lam_f=lambdal;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      } else {
        lam_f=lambdav;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
      }

      real tmp_w = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);

      real factl = lam_s/lambdal;
      real factv = lam_s/lambdav;
      txl[i-1][j][k] = factl * (tmp_w-tmp_s)/len_s;
      txv[i-1][j][k] = factv * (tmp_w-tmp_s)/len_s;
    }

    /* east is in wall & there is interface in east */
    if(dom->ibody().off(i+1,j,k)&&Interface(+1,m,i,j,k)) {
      real len_s = 0.5*phi.dxc(i+1);
      real lam_s = solid()->lambda(i+1,j,k);
      real tmp_s = tpr[i+1][j][k];

      real tmp_f;
      real len_f = 0.5*phi.dxc(i)-distance_x(i,j,k,+1,tmp_f);
      real lam_f;
      /* note the inversion */
      if(clrc<clrsurf) {
        lam_f=lambdal;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      } else {
        lam_f=lambdav;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
      } 

      real tmp_w = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);

      real factl = lam_s/lambdal;
      real factv = lam_s/lambdav;
      txl[i+1][j][k] = factl * (tmp_s-tmp_w)/len_s;
      txv[i+1][j][k] = factv * (tmp_s-tmp_w)/len_s;
    }
   
    /* y direction */
    m = Comp::v();
 
    /* south is in wall & there is interface in south */
    if(dom->ibody().off(i,j-1,k)&&Interface(-1,m,i,j,k)) {
      real len_s = 0.5*phi.dyc(j-1);
      real lam_s = solid()->lambda(i,j-1,k);
      real tmp_s = tpr[i][j-1][k];

      real tmp_f;
      real len_f = 0.5*phi.dyc(j)-distance_y(i,j,k,-1,tmp_f);
      real lam_f;
      /* note the inversion */
      if(clrc<clrsurf) {
        lam_f=lambdal;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      } else {
        lam_f=lambdav;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
      }

      real tmp_w = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);

      real factl = lam_s/lambdal;
      real factv = lam_s/lambdav;
      txl[i][j-1][k] = factl * (tmp_w-tmp_s)/len_s;
      txv[i][j-1][k] = factv * (tmp_w-tmp_s)/len_s;
    }

    /* north is in wall & there is interface in north */
    if(dom->ibody().off(i,j+1,k)&&Interface(+1,m,i,j,k)) {
      real len_s = 0.5*phi.dyc(j+1);
      real lam_s = solid()->lambda(i,j+1,k);
      real tmp_s = tpr[i][j+1][k];

      real tmp_f;
      real len_f = 0.5*phi.dyc(j)-distance_y(i,j,k,+1,tmp_f);
      real lam_f;
      /* note the inversion */
      if(clrc<clrsurf) {
        lam_f=lambdal;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      } else {
        lam_f=lambdav;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
      }

      real tmp_w = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);

      real factl = lam_s/lambdal;
      real factv = lam_s/lambdav;
      txl[i][j+1][k] = factl * (tmp_s-tmp_w)/len_s;
      txv[i][j+1][k] = factv * (tmp_s-tmp_w)/len_s;
    }

    /* z direction */
    m = Comp::w();

    /* bottom is in wall & there is interface in bottom */
    if(dom->ibody().off(i,j,k-1)&&Interface(-1,m,i,j,k)) {
      real len_s = 0.5*phi.dzc(k-1);
      real lam_s = solid()->lambda(i,j,k-1);
      real tmp_s = tpr[i][j][k-1];

      real tmp_f;
      real len_f = 0.5*phi.dzc(k)-distance_z(i,j,k,-1,tmp_f);
      real lam_f;
      /* note the inversion */
      if(clrc<clrsurf) {
        lam_f=lambdal;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      } else {
        lam_f=lambdav;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
      }

      real tmp_w = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);

      real factl = lam_s/lambdal;
      real factv = lam_s/lambdav;
      txl[i][j][k-1] = factl * (tmp_w-tmp_s)/len_s;
      txv[i][j][k-1] = factv * (tmp_w-tmp_s)/len_s;
    }

    /* top is in wall & there is interface in top */
    if(dom->ibody().off(i,j,k+1)&&Interface(+1,m,i,j,k)) {
      real len_s = 0.5*phi.dzc(k+1);
      real lam_s = solid()->lambda(i,j,k+1);
      real tmp_s = tpr[i][j][k+1];

      real tmp_f;
      real len_f = 0.5*phi.dzc(k)-distance_z(i,j,k,+1,tmp_f);
      real lam_f;
      /* note the inversion */
      if(clrc<clrsurf) {
        lam_f=lambdal;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      } else {
        lam_f=lambdav;
        if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
      }

      real tmp_w = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);

      real factl = lam_s/lambdal;
      real factv = lam_s/lambdav;
      txl[i][j][k+1] = factl * (tmp_w-tmp_s)/len_s;
      txv[i][j][k+1] = factv * (tmp_w-tmp_s)/len_s;
    }
  } /* ib cells */
    
  return;
}
