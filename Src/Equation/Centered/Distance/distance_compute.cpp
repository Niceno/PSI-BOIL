#include "distance.h"

/******************************************************************************/
void Distance::compute() {

  boil::timer.start("distance");

  /*----------------------------+
  |  create discretized system  |
  +----------------------------*/
  v_fluid.lambda( 1.0 );
  /* create only diffusion matrix ... */
  create_system_diffusive(flu->lambda());
  create_system_bnd();   

  /*--------------------+
  |  create the source  |
  +--------------------*/
  for_ijk(i,j,k) 
    fnew[i][j][k] = dom->dV(i,j,k);

  if( dom->ibody().nccells() > 0 )
    for_ijk(i,j,k) {
      real fV = dom->ibody().fV(i,j,k);  /* fraction in fluid */ 
      fnew[i][j][k] *= fV; 
      if(dom->ibody().off(i,j,k)) fnew[i][j][k] = 0.0;
    } 

  /*
  for_ijk(i,j,k) {
    if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j,k-1) ) 
      fnew[i][j][k] += dom->dV(i,j,k-1) * dom->ibody().fV(i,j,k-1);
    
    if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j,k+1) ) 
      fnew[i][j][k] += dom->dV(i,j,k+1) * dom->ibody().fV(i,j,k+1);
    
    if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j-1,k) ) 
      fnew[i][j][k] += dom->dV(i,j-1,k) * dom->ibody().fV(i,j-1,k);
    
    if( dom->ibody().on(i,j,k) && dom->ibody().off(i,j+1,k) ) 
      fnew[i][j][k] += dom->dV(i,j+1,k) * dom->ibody().fV(i,j+1,k);
    
    if( dom->ibody().on(i,j,k) && dom->ibody().off(i-1,j,k) ) 
      fnew[i][j][k] += dom->dV(i-1,j,k) * dom->ibody().fV(i-1,j,k);
    
    if( dom->ibody().on(i,j,k) && dom->ibody().off(i+1,j,k) ) 
      fnew[i][j][k] += dom->dV(i+1,j,k) * dom->ibody().fV(i+1,j,k);
  }
  */

  /*--------+
  |  solve  |
  +--------*/
  solver->solve(A, phi, fnew, min_iter,
                MaxIter(1000), "distance", 
                ResRat(boil::micro), ResTol(boil::pico));

  phi.bnd_update();

  /*-------------------+
  |  compute distance  |
  +-------------------*/
  Scalar dist( *(phi.domain()) );
  for_ijk(i,j,k) {
    real phi_x = (phi[i+1][j][k]-phi[i-1][j][k])/(phi.dxe(i)+phi.dxw(i));
    real phi_y = (phi[i][j+1][k]-phi[i][j-1][k])/(phi.dys(j)+phi.dyn(j));
    real phi_z = (phi[i][j][k+1]-phi[i][j][k-1])/(phi.dzb(k)+phi.dzt(k));

    real nabla_phi_dot = phi_x * phi_x + phi_y * phi_y + phi_z * phi_z;
    real nabla_phi_abs = sqrt( nabla_phi_dot );

    dist[i][j][k] = sqrt( nabla_phi_dot + 2.0 * phi[i][j][k] ) - nabla_phi_abs;
  }
  phi = dist;

  phi.exchange();

  boil::timer.stop("distance");
}
