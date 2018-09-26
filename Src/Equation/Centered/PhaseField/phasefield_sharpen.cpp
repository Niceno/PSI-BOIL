#include "phasefield.h"

#define NEW true

/******************************************************************************/
void PhaseField::sharpen() {
/*-----------------------------------------------------------------------------+
|                                                                              |
+-----------------------------------------------------------------------------*/

  real maxd = 0.0;

  const real ldt = 1e-5;

  /* width = 3.0 * sqrt(2.0) * w 
     w = width / 3.0 / sqrt(2.0) */
 
  const real width = phi.dxc(1) * 6.0;    /* width is 6 cells, let's say */
  const real w = width / 3.0 / sqrt(2.0); /* 0.01 */
  const real b = 1.0;

  const real max_value = 1.0;
  #if NEW==true
    const real min_value = 0.0;
  #else
    const real min_value = -1.0;
  #endif


  OPR( width );
  OPR( w );

  /*------------+
  |             |
  |  time-loop  | 
  |             |
  +------------*/
  for(int it=0; it<10; it++) {

    phi_old = phi;

    /*-----------------+
    |  diffusion-like  | 
    +-----------------*/
    fold = 0.0;

    /* compute diffusive fluxes using finite differences */
    for_ijk(i,j,k) {
      /* i */
      fold[i][j][k] += b * ( (phi[i+1][j][k] - phi[i][j][k]) / phi.dxe(i) -
                             (phi[i][j][k] - phi[i-1][j][k]) / phi.dxw(i) ) 
                         / phi.dxc(i);
      /* j */
      fold[i][j][k] += b * ( (phi[i][j+1][k] - phi[i][j][k]) / phi.dyn(j) -
                             (phi[i][j][k] - phi[i][j-1][k]) / phi.dys(j) ) 
                         / phi.dyc(j);
      /* k */
      fold[i][j][k] += b * ( (phi[i][j][k+1] - phi[i][j][k]) / phi.dzt(k) -
                             (phi[i][j][k] - phi[i][j][k-1]) / phi.dzb(k) ) 
                         / phi.dzc(k);
    }

    /* advance */
    for_ijk(i,j,k)
      phi[i][j][k] = phi[i][j][k] + ldt * fold[i][j][k];
    phi.bnd_update(); 
    phi.exchange_all();

    //debug:boil::plot->plot(phi, "diff-phi", it);

    /*-----------+
    |  source 1  |
    +-----------*/
    fold = 0.0;

    for_ijk(i,j,k) {
      #if NEW==true
        fold[i][j][k] -= 2.0 * b / w / w 
                       * (phi[i][j][k] - phi[i][j][k]*phi[i][j][k]) 
                       * (1.0 - 2.0*phi[i][j][k]);
      #else
        fold[i][j][k] += b / w / w * phi[i][j][k] 
                       * (1.0 - phi[i][j][k]*phi[i][j][k]); 
      #endif
    }

    /* advance */
    for_ijk(i,j,k)
      phi[i][j][k] = phi[i][j][k] + ldt * fold[i][j][k];
    phi.bnd_update(); 
    phi.exchange_all();

    //debug:boil::plot->plot(phi, "src-1-phi", it);

    /*-----------+
    |  source 2  |
    +-----------*/
    fold = 0.0;

    real phix, phiy, phiz;
    real abs_nabla_phi_w, abs_nabla_phi_e, 
         abs_nabla_phi_s, abs_nabla_phi_n, 
         abs_nabla_phi_b, abs_nabla_phi_t;

    /* compute \nabla \cdot (\nabla \phi / |\nabla \phi|) */
    for_ijk(i,j,k) {

      /* w */
      phix = (phi[i][j][k] - phi[i-1][j][k]) / phi.dxw(i);
      phiy = 0.5 * (   phi[i][j+1][k] + phi[i-1][j+1][k]
                     - phi[i][j-1][k] - phi[i-1][j-1][k] ) / (phi.dys(j)+phi.dyn(j));
      phiz = 0.5 * (   phi[i][j][k+1] + phi[i-1][j][k+1]
                     - phi[i][j][k-1] - phi[i-1][j][k-1] ) / (phi.dzb(k)+phi.dzt(k));
      abs_nabla_phi_w = sqrt( phix*phix + phiy*phiy + phiz*phiz );

      /* e */
      phix = (phi[i+1][j][k] - phi[i][j][k]) / phi.dxe(i);
      phiy = 0.5 * (   phi[i+1][j+1][k] + phi[i][j+1][k]
                     - phi[i+1][j-1][k] - phi[i][j-1][k] ) / (phi.dys(j)+phi.dyn(j));
      phiz = 0.5 * (   phi[i+1][j][k+1] + phi[i][j][k+1]
                     - phi[i+1][j][k-1] - phi[i][j][k-1] ) / (phi.dzb(k)+phi.dzt(k));
      abs_nabla_phi_e = sqrt( phix*phix + phiy*phiy + phiz*phiz );
   
      /* s */
      phiy = (phi[i][j][k] - phi[i][j-1][k]) / phi.dys(j);
      phix = 0.5 * (   phi[i+1][j][k] + phi[i+1][j-1][k] 
                     - phi[i-1][j][k] - phi[i-1][j-1][k] ) / (phi.dxw(i)+phi.dxe(i));
      phiz = 0.5 * (   phi[i][j][k+1] + phi[i][j-1][k+1] 
                     - phi[i][j][k-1] - phi[i][j-1][k-1] ) / (phi.dzb(k)+phi.dzt(k));
      abs_nabla_phi_s = sqrt( phix*phix + phiy*phiy + phiz*phiz );

      /* n */
      phiy = (phi[i][j+1][k] - phi[i][j][k]) / phi.dyn(j);
      phix = 0.5 * (   phi[i+1][j+1][k] + phi[i+1][j][k] 
                     - phi[i-1][j+1][k] - phi[i-1][j][k] ) / (phi.dxw(i)+phi.dxe(i));
      phiz = 0.5 * (   phi[i][j+1][k+1] + phi[i][j][k+1] 
                     - phi[i][j+1][k-1] - phi[i][j][k-1] ) / (phi.dzb(k)+phi.dzt(k));
      abs_nabla_phi_n = sqrt( phix*phix + phiy*phiy + phiz*phiz );

      /* b */
      phiz = (phi[i][j][k] - phi[i][j][k-1]) / phi.dzb(k);
      phix = 0.5 * (   phi[i+1][j][k] + phi[i+1][j][k-1] 
                     - phi[i-1][j][k] - phi[i-1][j][k-1] ) / (phi.dxw(i)+phi.dxe(i));
      phiy = 0.5 * (   phi[i][j+1][k] + phi[i][j+1][k-1]
                     - phi[i][j-1][k] - phi[i][j-1][k-1] ) / (phi.dys(j)+phi.dyn(j));
      abs_nabla_phi_b = sqrt( phix*phix + phiy*phiy + phiz*phiz );

      /* t */
      phiz = (phi[i][j][k+1] - phi[i][j][k]) / phi.dzt(k);
      phix = 0.5 * (   phi[i+1][j][k+1] + phi[i+1][j][k] 
                     - phi[i-1][j][k+1] - phi[i-1][j][k] ) / (phi.dxw(i)+phi.dxe(i));
      phiy = 0.5 * (   phi[i][j+1][k+1] + phi[i][j+1][k]
                     - phi[i][j-1][k+1] - phi[i][j-1][k] ) / (phi.dys(j)+phi.dyn(j));
      abs_nabla_phi_t = sqrt( phix*phix + phiy*phiy + phiz*phiz );

      real term = 0.0;
      /* i */
      term += (   (phi[i+1][j][k] - phi[i][j][k]) 
                 / phi.dxe(i) / (abs_nabla_phi_e+boil::pico) 
                - (phi[i][j][k] - phi[i-1][j][k]) 
                 / phi.dxw(i) / (abs_nabla_phi_w+boil::pico) ) 
           / phi.dxc(i);
      /* j */
      term += (   (phi[i][j+1][k] - phi[i][j][k]) 
                / phi.dyn(j) / (abs_nabla_phi_n+boil::pico)
                - (phi[i][j][k] - phi[i][j-1][k]) 
                / phi.dys(j) / (abs_nabla_phi_s+boil::pico) ) 
           / phi.dyc(j);
      /* j */
      term += (   (phi[i][j][k+1] - phi[i][j][k]) 
                / phi.dzt(k) / (abs_nabla_phi_t+boil::pico)
                - (phi[i][j][k] - phi[i][j][k-1]) 
                / phi.dzb(k) / (abs_nabla_phi_b+boil::pico) ) 
           / phi.dzc(k);

      /* compute |\nabla \phi| (stored in abs_nabla_phi) */
      phix = (phi[i+1][j][k] - phi[i-1][j][k]) / (phi.dxw(i)+phi.dxe(i));
      phiy = (phi[i][j+1][k] - phi[i][j-1][k]) / (phi.dys(j)+phi.dyn(j));
      phiz = (phi[i][j][k+1] - phi[i][j][k-1]) / (phi.dzb(k)+phi.dzt(k));
      real absnablaphi = sqrt( phix*phix + phiy*phiy + phiz*phiz );

      fold[i][j][k] = - b * absnablaphi * term; /* keep the minus ;-) */
    }

    /* advance */
    for_ijk(i,j,k)
      phi[i][j][k] = phi[i][j][k] + ldt * fold[i][j][k];
    phi.bnd_update(); 
    phi.exchange_all();

    /*-------+
    |  trim  |
    for_ijk(i,j,k) {
      if( phi[i][j][k] < min_value+boil::micro ) phi[i][j][k] =  min_value;
      if( phi[i][j][k] > max_value-boil::micro ) phi[i][j][k] =  max_value;
    }
    phi.bnd_update(); 
    phi.exchange_all();
    +-------*/

    /*--------------------+  
    |  check convergence  |
    +--------------------*/
    maxd = 0.0;
    for_ijk(i,j,k) 
      maxd = std::max( maxd, fabs( phi_old[i][j][k] - phi[i][j][k] ) );
    boil::cart.max_real(&maxd);

    OPR(maxd);
    if( maxd < boil::micro ) return;
  }
}

/*-----------------------------------------------------------------------------+
 '$Id: phasefield_sharpen.cpp,v 1.7 2014/10/15 13:39:38 niceno Exp $'/
+-----------------------------------------------------------------------------*/
