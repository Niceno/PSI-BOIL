#include "vof.h"

/******************************************************************************/
/* calculate volume fraction in staggered cells */
/******************************************************************************/
void VOF::cal_bndclr(const Scalar & scp) {
  int ofx(0),ofy(0),ofz(0);
  for_m(m) {
    if       (m==Comp::i()) {
      ofx = -1;
      ofy = 0;
      ofz = 0;
    } else if(m==Comp::j()) {
      ofx = 0;
      ofy = -1;
      ofz = 0;
    } else {
      ofx = 0;
      ofy = 0;
      ofz = -1;
    }

    /* wvmijk should take care of boundaries, e.g. walls */
    for_wvmijk((*bndclr),m,i,j,k) {
      int ii = i+ofx;
      int jj = j+ofy;
      int kk = k+ofz;

      real scpp = scp[i ][j ][k ];
      real scpm = scp[ii][jj][kk];     
     
      bool liqp = scpp-1.>-boil::pico;
      bool gasp = scpp<boil::pico;
      bool liqm = scpm-1.>-boil::pico;
      bool gasm = scpm<boil::pico;

      real stagscpp, stagscpm;

      /* degenerate cases */
      if( (liqp|gasp) & (liqm|gasm) ) {
        (*bndclr)[m][i][j][k] = 0.5 * (scpp + scpm);
        continue;
      }

      if(liqp||gasp) {
        stagscpp = scpp*0.5*dV(i,j,k); /* this only works for cart. grid...*/
      } else {
        /* we are looking in the negative direction */
        real g = -0.5;
        real c = scpp;

        real vn1 = -nx[i][j][k];
        real vn2 = -ny[i][j][k];
        real vn3 = -nz[i][j][k];

        real absg = fabs(g);
        real vm1 = fabs(vn1);
        real vm2 = fabs(vn2);
        real vm3 = fabs(vn3)+boil::pico;
        real qa = 1.0/(vm1+vm2+vm3);

        vm1 *= qa;
        vm2 *= qa;
        vm3 *= qa;
        real alpha = qa*nalpha[i][j][k];

        real * alig_vm;
        real * alig_vn;
        if       (m==Comp::i()) {
	  alig_vm = &vm1;
          alig_vn = &vn1;
        } else if(m==Comp::j()) {
          alig_vm = &vm2;
          alig_vn = &vn2;
        } else {
          alig_vm = &vm3;
          alig_vn = &vn3;
        }

        real ra = (*alig_vm) * (1.0-absg);
        qa = 1.0/(1.0-ra);
        if((*alig_vn)*g>0.0) alpha -= ra; 
        (*alig_vm) *= absg;

        /* here we use absg to avoid negative sign */
        stagscpp = calc_v(alpha*qa, vm1*qa, vm2*qa, vm3*qa)*absg*dV(i,j,k);
      }

      if(liqm||gasm) {
        stagscpm = scpm*0.5*dV(ii,jj,kk);
      } else {
        real g = 0.5;
        real c = scpm;

        real vn1 = -nx[ii][jj][kk];
        real vn2 = -ny[ii][jj][kk];
        real vn3 = -nz[ii][jj][kk];

        real vm1 = fabs(vn1);
        real vm2 = fabs(vn2);
        real vm3 = fabs(vn3)+boil::pico;
        real qa = 1.0/(vm1+vm2+vm3);

        vm1 *= qa;
        vm2 *= qa;
        vm3 *= qa;
        //real alpha = calc_alpha(c, vm1, vm2, vm3);
        real alpha = qa*nalpha[ii][jj][kk];

        real * alig_vm;
        real * alig_vn;
        if       (m==Comp::i()) {
          alig_vm = &vm1;
          alig_vn = &vn1;
        } else if(m==Comp::j()) {
          alig_vm = &vm2;
          alig_vn = &vn2;
        } else {
          alig_vm = &vm3;
          alig_vn = &vn3;
        }

        real ra = (*alig_vm) * (1.0-g);
        qa = 1.0/(1.0-ra);
        if((*alig_vn)>0.0) alpha -= ra;
        (*alig_vm) *= g;

        stagscpm = calc_v(alpha*qa, vm1*qa, vm2*qa, vm3*qa)*g*dV(ii,jj,kk);
      }

#if 0
      if(j==i-1&&m==Comp::i())
      boil::oout<<"VOF::bdclr "<<i<<" "<<stagscpp<<" "<<stagscpm<<" "<<(*bndclr).dV(m,i,j,k)<<boil::endl;
#endif

    
      (*bndclr)[m][i][j][k] = (stagscpp+stagscpm)/(*bndclr).dV(m,i,j,k);
    } /* for_wvmijk */
  } /* for_m */
}
