#include "phasechangevof.h"

/******************************************************************************/
void PhaseChangeVOF::subgrid_gradt() {
/***************************************************************************//**
*  \brief correct gradient of temperature near walls
*******************************************************************************/

  //boil::plot->plot(clr,txl,tyl,tzl, "ib+II+clr-tx-ty-tz", time->current_step());

  /* wall */
  for(int b = 0; b < clr.bc().count(); b++) {

    if(clr.bc().type_decomp(b))
      continue;

    if(clr.bc().type(b) == BndType::wall()) {

      Dir d = clr.bc().direction(b);
      if(d != Dir::undefined()) {
        Comp mcomp;
        int dir(0), ofx(0), ofy(0), ofz(0);
        if (d == Dir::imin()) {
          mcomp = Comp::i();
          dir = -1;
          ofx = +1;
        } else if (d == Dir::imax()) {
          mcomp = Comp::i();
          dir = +1;
          ofx = -1;
        } else if (d == Dir::jmin()) {
          mcomp = Comp::j();
          dir = -1;
          ofy = +1;
        } else if (d == Dir::jmax()) {
          mcomp = Comp::j();
          dir = +1;
          ofy = -1;
        } else if (d == Dir::kmin()) {
          mcomp = Comp::k();
          dir = -1;
          ofz = +1;
        } else if (d == Dir::kmax()) {
          mcomp = Comp::k();
          dir = +1;
          ofz = -1;
        } else {
          continue;
        }

        /* in the following code, subgrid cell temperature gradients are calculated */
        for_vijk( clr.bc().at(b), i,j,k ) {
          int ii = i+ofx;
          int jj = j+ofy;
          int kk = k+ofz;

          /* is there an interface between cell centre and wall? */
          if( (tempflag[ii][jj][kk]==3) || (tempflag[ii][jj][kk]==-3) ) continue;

          /* is there an interface between cell centre and wall? */
          if(Interface(dir,mcomp,ii,jj,kk)) {
            /* wall temperature */
            real tw = tpr[i][j][k];
            /* interface temperature and distance wall-interface */
            real ti;
            real dist;

            /* the temperature gradient for inverse phase is set
               (the other two components are extrapolated) */
            if       (mcomp==Comp::i()) {
              dist = distance_x(ii,jj,kk,dir,ti);
              dist = clr.dxc(ii)/2.0 - dist;
              if(clr[ii][jj][kk]>=clrsurf) {
                //txv[ii][jj][kk] = (tw-ti)/dist * real(dir);
                tnv[ii][jj][kk] = (tw-ti)/dist * (+1);
              } else {
                //txl[ii][jj][kk] = (tw-ti)/dist * real(dir);
                tnl[ii][jj][kk] = (tw-ti)/dist * (-1);
              }
            } else if(mcomp==Comp::j()) {
              dist = distance_y(ii,jj,kk,dir,ti);
              dist = clr.dyc(jj)/2.0 - dist;
              if(clr[ii][jj][kk]>=clrsurf) {
                //tyv[ii][jj][kk] = (tw-ti)/dist * real(dir);
                tnv[ii][jj][kk] = (tw-ti)/dist * (+1);
              } else {
                //tyl[ii][jj][kk] = (tw-ti)/dist * real(dir);
                tnl[ii][jj][kk] = (tw-ti)/dist * (-1);
              }
            } else {
              dist = distance_z(ii,jj,kk,dir,ti);
              dist = clr.dzc(kk)/2.0 - dist;
              if(clr[ii][jj][kk]>=clrsurf) {
                //tzv[ii][jj][kk] = (tw-ti)/dist * real(dir);
                tnv[ii][jj][kk] = (tw-ti)/dist * (+1);
              } else {
                //tzl[ii][jj][kk] = (tw-ti)/dist * real(dir);
                tnl[ii][jj][kk] = (tw-ti)/dist * (-1);
              }
            }
           
          } /* is interface? */
        } /* loop over cells */

      } /* dir not undefined */
    } /* is wall? */
  } /* loop over bcs */
        

#ifdef IB
  /* immersed body */
  for(int cc=0; cc<dom->ibody().nccells(); cc++) {
    int i,j,k;
    /* cell[i][j][k] is wall adjacent cell in fluid domain */
    dom->ibody().ijk(cc,&i,&j,&k);

    if( (tempflag[i][j][k]==3) || (tempflag[i][j][k]==-3) ) continue;

    real clrc = clr[i][j][k];

    /* x direction */
    Comp m = Comp::u();

    /* west is in wall & there is interface in west */
    if(dom->ibody().off(i-1,j,k)&&Interface(-1,m,i,j,k)) {
      real tmp_s = tpr[i-1][j][k];
      real len_s = phi.dxw(i) - 0.5*phi.dxc(i);
      real lam_s = solid()->lambda(i-1,j,k);

      real tmp_f;
      real len_f = distance_x(i,j,k,-1,tmp_f);
      len_f = 0.5*phi.dxc(i) - len_f;
      real lam_f;
      if(clrc>=clrsurf) { /* note the inversion! */
        lam_f=lambdav;
        /* diff eddy omitted in subgrid layers */
        //if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
      } else {
        lam_f=lambdal;
        //if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      }

      real tmp_w = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);

      if(clrc>=clrsurf) {
        //txv[i][j][k] = (tmp_f-tmp_w)/len_f;
        tnv[i][j][k] = (tmp_f-tmp_w)/len_f * (+1);
      } else {
        //txl[i][j][k] = (tmp_f-tmp_w)/len_f;
        tnl[i][j][k] = (tmp_f-tmp_w)/len_f * (-1);
      }
    }

    /* east is in wall & there is interface in east */
    if(dom->ibody().off(i+1,j,k)&&Interface(+1,m,i,j,k)) {
      real tmp_s = tpr[i+1][j][k];
      real len_s = phi.dxe(i) - 0.5*phi.dxc(i);
      real lam_s = solid()->lambda(i+1,j,k);

      real tmp_f;
      real len_f = distance_x(i,j,k,+1,tmp_f);
      len_f = 0.5*phi.dxc(i) - len_f;
      real lam_f;
      if(clrc>=clrsurf) { /* note the inversion! */
        lam_f=lambdav;
        /* diff eddy omitted in subgrid layers */
        //if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
      } else {
        lam_f=lambdal;
        //if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      }

      real tmp_w = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);

      if(clrc>=clrsurf) {
        //txv[i][j][k] = (tmp_w-tmp_f)/len_f;
        tnv[i][j][k] = (tmp_f-tmp_w)/len_f * (+1);
      } else {
        //txl[i][j][k] = (tmp_w-tmp_f)/len_f;
        tnl[i][j][k] = (tmp_f-tmp_w)/len_f * (-1);
      }
    }

    /* y direction */
    m = Comp::v();

    /* south is in wall & there is interface in south */
    if(dom->ibody().off(i,j-1,k)&&Interface(-1,m,i,j,k)) {
      real tmp_s = tpr[i][j-1][k];
      real len_s = phi.dys(j) - 0.5*phi.dyc(j);
      real lam_s = solid()->lambda(i,j-1,k);

      real tmp_f;
      real len_f = distance_y(i,j,k,-1,tmp_f);
      len_f = 0.5*phi.dyc(j) - len_f;
      real lam_f;
      if(clrc>=clrsurf) { /* note the inversion! */
        lam_f=lambdav;
        /* diff eddy omitted in subgrid layers */
        //if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
      } else {
        lam_f=lambdal;
        //if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      }

      real tmp_w = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);

      if(clrc>=clrsurf) {
        //tyv[i][j][k] = (tmp_f-tmp_w)/len_f;
        tnv[i][j][k] = (tmp_f-tmp_w)/len_f * (+1);
      } else {
        //tyl[i][j][k] = (tmp_f-tmp_w)/len_f;
        tnl[i][j][k] = (tmp_f-tmp_w)/len_f * (-1);
      }
    }

    /* north is in wall & there is interface in north */
    if(dom->ibody().off(i,j+1,k)&&Interface(+1,m,i,j,k)) {
      real tmp_s = tpr[i][j+1][k];
      real len_s = phi.dyn(j) - 0.5*phi.dyc(j);
      real lam_s = solid()->lambda(i,j+1,k);

      real tmp_f;
      real len_f = distance_y(i,j,k,+1,tmp_f);
      len_f = 0.5*phi.dyc(j) - len_f;
      real lam_f;
      if(clrc>=clrsurf) { /* note the inversion! */
        lam_f=lambdav;
        /* diff eddy omitted in subgrid layers */
        //if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
      } else {
        lam_f=lambdal;
        //if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      }

      real tmp_w = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);

      if(clrc>=clrsurf) {
        //tyv[i][j][k] = (tmp_w-tmp_f)/len_f;
        tnv[i][j][k] = (tmp_f-tmp_w)/len_f * (+1);
      } else {
        //tyl[i][j][k] = (tmp_w-tmp_f)/len_f;
        tnl[i][j][k] = (tmp_f-tmp_w)/len_f * (-1);
      }
    }

    /* z direction */
    m = Comp::w();

    /* bottom is in wall & there is interface in bottom */
    if(dom->ibody().off(i,j,k-1)&&Interface(-1,m,i,j,k)) {
      real tmp_s = tpr[i][j][k-1];
      real len_s = phi.dzb(k) - 0.5*phi.dzc(k);
      real lam_s = solid()->lambda(i,j,k-1);

      real tmp_f;
      real len_f = distance_z(i,j,k,-1,tmp_f);
      len_f = 0.5*phi.dzc(k) - len_f;
      real lam_f;
      if(clrc>=clrsurf) { /* note the inversion! */
        lam_f=lambdav;
        /* diff eddy omitted in subgrid layers */
        //if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
      } else {
        lam_f=lambdal;
        //if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      }

      real tmp_w = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);

      if(clrc>=clrsurf) {
        //tzv[i][j][k] = (tmp_f-tmp_w)/len_f;
        tnv[i][j][k] = (tmp_f-tmp_w)/len_f * (+1);
      } else {
        //tzl[i][j][k] = (tmp_f-tmp_w)/len_f;
        tnl[i][j][k] = (tmp_f-tmp_w)/len_f * (-1);
      }
    }

    /* top is in wall & there is interface in top */
    if(dom->ibody().off(i,j,k+1)&&Interface(+1,m,i,j,k)) {
      real tmp_s = tpr[i][j][k+1];
      real len_s = phi.dzt(k) - 0.5*phi.dzc(k);
      real lam_s = solid()->lambda(i,j,k+1);

      real tmp_f;
      real len_f = distance_z(i,j,k,+1,tmp_f);
      len_f = 0.5*phi.dzc(k) - len_f;
      real lam_f;
      if(clrc>=clrsurf) { /* note the inversion! */
        lam_f=lambdav;
        /* diff eddy omitted in subgrid layers */
        //if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpv/rhov/turbP;
      } else {
        lam_f=lambdal;
        //if(diff_eddy) lam_f += (*diff_eddy)[i][j][k]*cpl/rhol/turbP;
      }

      real tmp_w = temperature_node(len_s, lam_s, tmp_s, len_f, lam_f, tmp_f);

      if(clrc>=clrsurf) {
        //tzv[i][j][k] = (tmp_w-tmp_f)/len_f;
        tnv[i][j][k] = (tmp_f-tmp_w)/len_f * (+1);
      } else {
        //tzl[i][j][k] = (tmp_w-tmp_f)/len_f;
        tnl[i][j][k] = (tmp_f-tmp_w)/len_f * (-1);
      }
    }

  } /* loop over ib cells */
#endif

#if 0 /* a simple 1D version is being used */
  /* normal component of gradient is calculated */
  for_vijk(clr,i,j,k) {
    if        (tempflag[i][j][k]== 1) { /* subgrid cell with liquid majority */
      tnv[i][j][k] = txv[i][j][k]*nx[i][j][k]
                   + tyv[i][j][k]*ny[i][j][k]
                   + tzv[i][j][k]*nz[i][j][k];
    } else if (tempflag[i][j][k]==-1) { /* subgrid cell with vapour majority */
      tnl[i][j][k] = txl[i][j][k]*nx[i][j][k]
                   + tyl[i][j][k]*ny[i][j][k]
                   + tzl[i][j][k]*nz[i][j][k];
    }
  }
#endif

  //boil::plot->plot(clr,txl,tyl,tzl, "ib+III+clr-tx-ty-tz", time->current_step());

  return;
}
