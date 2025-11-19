#include "phasechange4.h"
#include "../../../Parallel/Out/out.h"

/******************************************************************************/
void PhaseChange4::mdot() {
/***************************************************************************//**
*         Calculate mdot, using redistribution of cut-mdot.
*         output : phi [kg/m3s]
*         temporary: A
*******************************************************************************/

#if 0
  boil::oout<<"PCVOF::mdot"<<boil::endl;
#endif

  /* step one: calculated mdot and cut-mdot (=A) */
  for_ijk(i,j,k) {
    A.c[i][j][k] = 0.0;
    A.w[i][j][k] = 0.0;
    A.e[i][j][k] = 0.0;
    A.s[i][j][k] = 0.0;
    A.n[i][j][k] = 0.0;
    A.b[i][j][k] = 0.0;
    A.t[i][j][k] = 0.0;
    if(cht.interface(i,j,k)) {
      if(dom->ibody().on(i,j,k)) {
        real mflxc=M[i][j][k];

        /* iso-surface area */
        phi[i][j][k] = mflxc * cht.topo->get_adens()[i][j][k];
        phi[i][j][k] = mdot_cut(phi[i][j][k],vf[i][j][k],
                                cht.rhol(i,j,k),A.c[i][j][k]);

        /* calculate directional cut-mdot. normvector components squared
           are used as weighs as they naturally sum to one. redistribution 
           direction is decided on phi and normvector sign considerations */
        if(A.c[i][j][k] != 0.0) {
          real wx = nx[i][j][k]*nx[i][j][k];
          real wy = ny[i][j][k]*ny[i][j][k];
          real wz = nz[i][j][k]*nz[i][j][k];

          real wc(1.0);
          if(A.c[i][j][k]<0.0) {
            wc = -1.0;
          }

          A.w[i][j][k] = ((wc*nx[i][j][k])<0.)*wx*A.c[i][j][k];
          A.e[i][j][k] = ((wc*nx[i][j][k])>0.)*wx*A.c[i][j][k];
          A.s[i][j][k] = ((wc*ny[i][j][k])<0.)*wy*A.c[i][j][k];
          A.n[i][j][k] = ((wc*ny[i][j][k])>0.)*wy*A.c[i][j][k];
          A.b[i][j][k] = ((wc*nz[i][j][k])<0.)*wz*A.c[i][j][k];
          A.t[i][j][k] = ((wc*nz[i][j][k])>0.)*wz*A.c[i][j][k];
        }
      } else {
        phi[i][j][k] = 0.0;
      }
    } else {
      phi[i][j][k] = 0.0;
    }
  }
#if 1
  A.w.exchange();
  A.e.exchange();
  A.s.exchange();
  A.n.exchange();
  A.b.exchange();
  A.t.exchange();

  /* step two: redistribute cut-mdot */
  for_ijk(i,j,k) {
    if(dom->ibody().on(i,j,k)) {
      phi[i][j][k] += A.w[i+1][j][k];
      phi[i][j][k] += A.e[i-1][j][k];
      phi[i][j][k] += A.s[i][j+1][k];
      phi[i][j][k] += A.n[i][j-1][k];
      phi[i][j][k] += A.b[i][j][k+1];
      phi[i][j][k] += A.t[i][j][k-1];
      /* cut-again */
      real dummy;
      phi[i][j][k] = mdot_cut(phi[i][j][k],vf[i][j][k],cht.rhol(i,j,k),dummy);
    }
  }
#endif

  phi.bnd_update();
  phi.exchange_all();

  return;
}

