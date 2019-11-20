#include "phasechangevof.h"

/******************************************************************************/
void PhaseChangeVOF::subgrid_setflag() {
/***************************************************************************//**
*  \brief set flag for subgrid cells 
         vapor cell = -3
         liquid cell = 3
         subgrid cell with centre being vapor = -1
         subgrid cell with centre being liquid = 1
*******************************************************************************/
  for_vijk(clr,i,j,k) {
    if(clr[i][j][k]>=clrsurf){
      tempflag[i][j][k]=3;
    } else {
      tempflag[i][j][k]=-3;
    }
  }

#ifdef IB
  for_ijk(i,j,k) {
    if(dom->ibody().off(i,j,k))
      tempflag[i][j][k]=-1001;
  }
#endif

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

        /* in the following code, subgrid cells are flagged */
        for_vijk( clr.bc().at(b), i,j,k ) {
          int ii = i+ofx;
          int jj = j+ofy;
          int kk = k+ofz;

          /* is there an interface between cell centre and wall? */
          if(Interface(dir,mcomp,ii,jj,kk)) {
            if(clr[ii][jj][kk]>=clrsurf) {
              tempflag[ii][jj][kk] = +1;
            } else {
              tempflag[ii][jj][kk] = -1;
            }
          }
        }

      } /* dir not undefined */
    } /* is wall? */
  } /* loop over bcs */
        

#ifdef IB
  /* immersed body */
  for(int cc=0; cc<dom->ibody().nccells(); cc++){
    int i,j,k;
    /* cell[i][j][k] is wall adjacent cell in fluid domain */
    dom->ibody().ijk(cc,&i,&j,&k);

    /* west is in wall & there is interface in west etc */
    if(  (dom->ibody().off(i-1,j,k)&&Interface(-1,Comp::i(),i,j,k))
       ||(dom->ibody().off(i+1,j,k)&&Interface(+1,Comp::i(),i,j,k))
       ||(dom->ibody().off(i,j-1,k)&&Interface(-1,Comp::j(),i,j,k))
       ||(dom->ibody().off(i,j+1,k)&&Interface(+1,Comp::j(),i,j,k))
       ||(dom->ibody().off(i,j,k-1)&&Interface(-1,Comp::k(),i,j,k))
       ||(dom->ibody().off(i,j,k+1)&&Interface(+1,Comp::k(),i,j,k))) {

      if(clr[i][j][k]>=clrsurf) {
        tempflag[i][j][k] = +1;
      } else {
        tempflag[i][j][k] = -1;
      }
    }
  }
#endif

  tempflag.exchange_all();

  //boil::plot->plot(clr,tempflag, "ib+setflag+clr-dflag", time->current_step());
  //boil::plot->plot(clr,nx,ny,nz, "clr-nx-ny-nz", time->current_step());
  //exit(0);

  return;
}

