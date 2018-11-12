#include "colorfunction.h"

#define CONSERVATIVE true

/******************************************************************************/
void ColorFunction::sharpen() {
/*-----------------------------------------------------------------------------+
|  the weakness of this subroutine is in the computation of derivatives close  |
|  to boundaries. calls to bnd_update simply copy the values from interior to  |
|  the boundary and derivatives (and fluxes) are therefore underestimated.     |
|                                                                              |
|  note: bnd_update() and exchange_all() are called each time phi changes.     |
+-----------------------------------------------------------------------------*/

  real maxd = 0.0;

  /*----------------------------+
  |  compute normals only once  |
  +----------------------------*/
  phi.bnd_update(); 
  phi.exchange_all();
  for_ijk(i,j,k) {

    real ni = nx[i][j][k] = (phi[i+1][j][k]-phi[i-1][j][k])/(dxw(i)+dxe(i));
    real nj = ny[i][j][k] = (phi[i][j+1][k]-phi[i][j-1][k])/(dys(j)+dyn(j));
    real nk = nz[i][j][k] = (phi[i][j][k+1]-phi[i][j][k-1])/(dzb(k)+dzt(k));

    real magn = sqrt(ni*ni + nj*nj + nk*nk) + boil::micro;

    nx[i][j][k] /= magn;
    ny[i][j][k] /= magn;
    nz[i][j][k] /= magn;
  }
  nx.bnd_update(); 
  ny.bnd_update(); 
  nz.bnd_update(); 
  nx.exchange_all();
  ny.exchange_all();
  nz.exchange_all();

  //debug:boil::plot->plot(nx,ny,nz, "normals-a", it);

  /*------------+
  |             |
  |  time-loop  | 
  |             |
  +------------*/
  for(int it=0; it<10; it++) {

    phi_old = phi;

    /*------------+
    |  diffusion  | 
    +------------*/
    fold = 0.0;

    /* compute diffusive fluxes */
    for_ijk(i,j,k) {
      fold[i][j][k] += (A.w[i][j][k] * (phi[i-1][j][k] - phi[i][j][k]) 
                      + A.e[i][j][k] * (phi[i+1][j][k] - phi[i][j][k])); 
      fold[i][j][k] += (A.s[i][j][k] * (phi[i][j-1][k] - phi[i][j][k]) 
                      + A.n[i][j][k] * (phi[i][j+1][k] - phi[i][j][k])); 
      fold[i][j][k] += (A.b[i][j][k] * (phi[i][j][k-1] - phi[i][j][k])
                      + A.t[i][j][k] * (phi[i][j][k+1] - phi[i][j][k]));
    }

    #if !CONSERVATIVE
    /* advance */
    for_ijk(i,j,k)
      phi[i][j][k] = phi[i][j][k] + ldt/dV(i,j,k) * fold[i][j][k];
    phi.bnd_update(); 
    phi.exchange_all();

    //debug:boil::plot->plot(phi, "diff-phi", it);

    /*-------------+
    |  convection  |
    +-------------*/
    fold = 0.0;
    #endif

    /* compute fluxes */
    real tm, tc, tp;
    for_ijk(i,j,k) {
      tm = phi[i-1][j][k] * (1.0-phi[i-1][j][k]) * nx[i-1][j][k];
      tc = phi[i]  [j][k] * (1.0-phi[i]  [j][k]) * nx[i]  [j][k];
      tp = phi[i+1][j][k] * (1.0-phi[i+1][j][k]) * nx[i+1][j][k];
      fold[i][j][k] += 0.5 * (tm+tc) * dSx(i,j,k);
      fold[i][j][k] -= 0.5 * (tp+tc) * dSx(i,j,k);

      tm = phi[i][j-1][k] * (1.0-phi[i][j-1][k]) * ny[i][j-1][k];
      tc = phi[i][j]  [k] * (1.0-phi[i][j]  [k]) * ny[i][j]  [k];
      tp = phi[i][j+1][k] * (1.0-phi[i][j+1][k]) * ny[i][j+1][k];
      fold[i][j][k] += 0.5 * (tm+tc) * dSy(i,j,k);
      fold[i][j][k] -= 0.5 * (tp+tc) * dSy(i,j,k);

      tm = phi[i][j][k-1] * (1.0-phi[i][j][k-1]) * nz[i][j][k-1];
      tc = phi[i][j][k]   * (1.0-phi[i][j][k]  ) * nz[i][j][k];
      tp = phi[i][j][k+1] * (1.0-phi[i][j][k+1]) * nz[i][j][k+1];
      fold[i][j][k] += 0.5 * (tm+tc) * dSz(i,j,k);
      fold[i][j][k] -= 0.5 * (tp+tc) * dSz(i,j,k);
    }

    /* advance */
    for_ijk(i,j,k)
      phi[i][j][k] = phi[i][j][k] + ldt/dV(i,j,k) * fold[i][j][k];
    phi.bnd_update(); 
    phi.exchange_all();

    //debug:boil::plot->plot(phi, "conv-phi", it);

    /*--------------------+  
    |  check convergence  |
    +--------------------*/
    maxd = 0.0;
    for_ijk(i,j,k) 
      maxd = std::max( maxd, fabs( phi_old[i][j][k] - phi[i][j][k] ) );
    boil::cart.max_real(&maxd);
    OPR(maxd);

    if( maxd < 0.01 ) {
      for_ijk(i,j,k) {
        if( phi[i][j][k] < 0.0 ) phi[i][j][k] = 0.0;
        if( phi[i][j][k] > 1.0 ) phi[i][j][k] = 1.0;
      }
      phi.bnd_update(); 
      phi.exchange_all();
      return;
    }
  }
}
