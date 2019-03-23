#include "vof.h"

void VOF::sharpen() {

  for_ijk(i,j,k) {
    real clr = phi[i][j][k];
    real adns = adens[i][j][k];
    if(clr<phisurf&&clr>boil::pico&&adns<boil::micro) {
      real nnx = nx[i][j][k];
      real nny = ny[i][j][k];
      real nnz = nz[i][j][k];

      real snnx = (nnx>0.0)-(nnx<0.0);
      real snny = (nny>0.0)-(nny<0.0);
      real snnz = (nnz>0.0)-(nnz<0.0);

      int ofx(0),ofy(0),ofz(0);
      real adensx(-1.0),adensy(-1.0),adensz(-1.0);
      if       (snnx>0.0) {
        ofx++;
        adensx = adens[i+ofx][j][k];
      } else if(snnx<0.0) {
        ofx--;
        adensx = adens[i+ofx][j][k];
      }
      if       (snny>0.0) {
        ofy++;
        adensy = adens[i][j+ofy][k];
      } else if(snny<0.0) {
        ofy--;
        adensy = adens[i][j+ofy][k];
      }
      if       (snnz>0.0) {
        ofz++;
        adensz = adens[i][j][k+ofz];
      } else if(snnz<0.0) {
        ofz--;
        adensz = adens[i][j][k+ofz];
      }

      real denom(1.0);
   
      if(adensx<=0.0)
        denom -= nnx*nnx;
      if(adensy<=0.0)
        denom -= nny*nny;
      if(adensz<=0.0)
        denom -= nnz*nnz;
 
      if(denom>boil::pico) {
        if(adensx>0.0) {
          phi[i+ofx][j][k] += clr*nnx*nnx/denom * dV(i,j,k)/dV(i+ofx,j,k);
          phi[i    ][j][k] -= clr*nnx*nnx/denom * dV(i,j,k)/dV(i+ofx,j,k);
        }
        if(adensy>0.0) {
          phi[i][j+ofy][k] += clr*nny*nny/denom * dV(i,j,k)/dV(i,j+ofy,k);
          phi[i][j    ][k] -= clr*nny*nny/denom * dV(i,j,k)/dV(i,j+ofy,k);
        }
        if(adensz>0.0) {
          phi[i][j][k+ofz] += clr*nnz*nnz/denom * dV(i,j,k)/dV(i,j,k+ofz);
          phi[i][j][k    ] -= clr*nnz*nnz/denom * dV(i,j,k)/dV(i,j,k+ofz);
        }

        boil::aout<<i<<" "<<j<<" "<<k<<" "<<clr<<" "<<nnx<<" "<<nny<<" "<<nnz<<" "<<adensx<<" "<<adensy<<" "<<adensz<<boil::endl;
      } /* denominator nonzero */
    } /* no adens and nonzero clr */
  } /* for ijk */

  phi.exchange();

  return;
}
