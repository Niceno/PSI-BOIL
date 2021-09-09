#include "vof.h"

/***************************************************************************//**
*  \brief Different advection schemes. Naive simply advances in the directions
*  independently, without constraints on over-fluxing. Reconstructed tries to 
*  limit this by reconstructing the geometry after each step. This requires di-
*  rection alternation to avoid systematic bias. Bounded introduces a corrective
*  term to prevent over-fluxing under the condition CFL<0.5. See the reference
*  provided in the definition of the AdvectionMethod class.
*******************************************************************************/

/* !!!!! cold is updated inside those advance_. functions !!!!! */

/******************************************************************************/
void VOF::advect_naive(Scalar & scp) {
/******************************************************************************/

  /* advance in x-direction */
  if(bflag_struct.ifull)
    advance_x(scp,cold);

  /* advance in y-direction */
  if(bflag_struct.jfull)
    advance_y(scp,cold);
    
  /* advance in z-direction */
  if(bflag_struct.kfull)
    advance_z(scp,cold);

  /* convert volume to vf */
  update_phi(cold,scp);

  /* reconstruct geometry */
  reconstruct_geometry(scp);

  return;
}

/******************************************************************************/
void VOF::advect_reconstructed(Scalar & scp) {
/******************************************************************************/

  /* alternate directions to minimize bias */
  /* xyz */
  if       (time->current_step()%(2*bflag_struct.dim)==label_adv[0]) {
    if(bflag_struct.ifull) 
      { advance_x(scp,cold); update_phi(cold,scp); reconstruct_geometry(scp); }
    if(bflag_struct.jfull) 
      { advance_y(scp,cold); update_phi(cold,scp); reconstruct_geometry(scp); }
    if(bflag_struct.kfull) 
      { advance_z(scp,cold); update_phi(cold,scp); reconstruct_geometry(scp); }

  /* yzx */
  } else if(time->current_step()%(2*bflag_struct.dim)==label_adv[1]) {
    if(bflag_struct.jfull) 
      { advance_y(scp,cold); update_phi(cold,scp); reconstruct_geometry(scp); }
    if(bflag_struct.kfull) 
      { advance_z(scp,cold); update_phi(cold,scp); reconstruct_geometry(scp); }
    if(bflag_struct.ifull) 
      { advance_x(scp,cold); update_phi(cold,scp); reconstruct_geometry(scp); }

  /* zxy */
  } else if(time->current_step()%(2*bflag_struct.dim)==label_adv[2]) {
    if(bflag_struct.kfull) 
      { advance_z(scp,cold); update_phi(cold,scp); reconstruct_geometry(scp); }
    if(bflag_struct.ifull) 
      { advance_x(scp,cold); update_phi(cold,scp); reconstruct_geometry(scp); }
    if(bflag_struct.jfull) 
      { advance_y(scp,cold); update_phi(cold,scp); reconstruct_geometry(scp); }

  /* xzy */
  } else if(time->current_step()%(2*bflag_struct.dim)==label_adv[3]) {
    if(bflag_struct.ifull) 
      { advance_x(scp,cold); update_phi(cold,scp); reconstruct_geometry(scp); }
    if(bflag_struct.kfull) 
      { advance_z(scp,cold); update_phi(cold,scp); reconstruct_geometry(scp); }
    if(bflag_struct.jfull) 
      { advance_y(scp,cold); update_phi(cold,scp); reconstruct_geometry(scp); }

  /* yxz */
  } else if(time->current_step()%(2*bflag_struct.dim)==label_adv[4]) {
    if(bflag_struct.jfull) 
      { advance_y(scp,cold); update_phi(cold,scp); reconstruct_geometry(scp); }
    if(bflag_struct.ifull) 
      { advance_x(scp,cold); update_phi(cold,scp); reconstruct_geometry(scp); }
    if(bflag_struct.kfull) 
      { advance_z(scp,cold); update_phi(cold,scp); reconstruct_geometry(scp); }

  /* zyx */
  } else if(time->current_step()%(2*bflag_struct.dim)==label_adv[5]) {
    if(bflag_struct.kfull) 
      { advance_z(scp,cold); update_phi(cold,scp); reconstruct_geometry(scp); }
    if(bflag_struct.jfull) 
      { advance_y(scp,cold); update_phi(cold,scp); reconstruct_geometry(scp); }
    if(bflag_struct.ifull) 
      { advance_x(scp,cold); update_phi(cold,scp); reconstruct_geometry(scp); }

  /* error */
  } else {
    OMS(Catastrophic error. Exiting);
  }

  return;
}

/******************************************************************************/
void VOF::advect_bounded(Scalar & scp) {
/******************************************************************************/

  /* store marker function */
  for_ijk(i,j,k) {
    tempflag2[i][j][k] = color()[i][j][k] > 0.5;
  }

  /* divergence correction to preserve volume if div.u != 0 */
  real dt = time->dt();
  for_ijk(i,j,k) {
    if(tempflag2[i][j][k])
      cold[i][j][k] -= u->divergence(i,j,k)*dt*dV(i,j,k);
  }
  cold.exchange();

  /* alternate directions to minimize bias */
  /* xyz */
  if       (time->current_step()%(2*bflag_struct.dim)==label_adv[0]) {
    //OPR(label_adv[0]);
    if(bflag_struct.ifull) 
    { advance_x(scp,cold); 
      divergence_skew(Comp::u(),tempflag2,cold);
      update_phi(cold,scp); 
      reconstruct_geometry(scp);
    }
    if(bflag_struct.jfull) 
    { advance_y(scp,cold); 
      divergence_skew(Comp::v(),tempflag2,cold);
      update_phi(cold,scp); 
      reconstruct_geometry(scp);
    }
    if(bflag_struct.kfull) 
    { advance_z(scp,cold);
      divergence_skew(Comp::w(),tempflag2,cold);
      update_phi(cold,scp);
      reconstruct_geometry(scp);
    }

  /* yzx */
  } else if(time->current_step()%(2*bflag_struct.dim)==label_adv[1]) {
    //OPR(label_adv[1]);
    if(bflag_struct.jfull) 
    { advance_y(scp,cold);
      divergence_skew(Comp::v(),tempflag2,cold);
      update_phi(cold,scp);
      reconstruct_geometry(scp);
    }
    if(bflag_struct.kfull) 
    { advance_z(scp,cold);
      divergence_skew(Comp::w(),tempflag2,cold);
      update_phi(cold,scp);
      reconstruct_geometry(scp);
    }
    if(bflag_struct.ifull) 
    { advance_x(scp,cold);
      divergence_skew(Comp::u(),tempflag2,cold);
      update_phi(cold,scp);
      reconstruct_geometry(scp);
    }

  /* zxy */
  } else if(time->current_step()%(2*bflag_struct.dim)==label_adv[2]) {
    //OPR(label_adv[2]);
    if(bflag_struct.kfull) 
    { advance_z(scp,cold);
      divergence_skew(Comp::w(),tempflag2,cold);
      update_phi(cold,scp);
      reconstruct_geometry(scp);
    }
    if(bflag_struct.ifull) 
    { advance_x(scp,cold);
      divergence_skew(Comp::u(),tempflag2,cold);
      update_phi(cold,scp);
      reconstruct_geometry(scp);
    }
    if(bflag_struct.jfull) 
    { advance_y(scp,cold);
      divergence_skew(Comp::v(),tempflag2,cold);
      update_phi(cold,scp);
      reconstruct_geometry(scp);
    }

  /* xzy */
  } else if(time->current_step()%(2*bflag_struct.dim)==label_adv[3]) {
    //OPR(label_adv[3]);
    if(bflag_struct.ifull) 
    { advance_x(scp,cold);
      divergence_skew(Comp::u(),tempflag2,cold);
      update_phi(cold,scp);
      reconstruct_geometry(scp);
    }
    if(bflag_struct.kfull) 
    { advance_z(scp,cold);
      divergence_skew(Comp::w(),tempflag2,cold);
      update_phi(cold,scp);
      reconstruct_geometry(scp);
    }
    if(bflag_struct.jfull) 
    { advance_y(scp,cold);
      divergence_skew(Comp::v(),tempflag2,cold);
      update_phi(cold,scp);
      reconstruct_geometry(scp);
    }

  /* yxz */
  } else if(time->current_step()%(2*bflag_struct.dim)==label_adv[4]) {
    //OPR(label_adv[4]);
    if(bflag_struct.jfull) 
    { advance_y(scp,cold);
      divergence_skew(Comp::v(),tempflag2,cold);
      update_phi(cold,scp);
      reconstruct_geometry(scp);
    }
    if(bflag_struct.ifull) 
    { advance_x(scp,cold);
      divergence_skew(Comp::u(),tempflag2,cold);
      update_phi(cold,scp);
      reconstruct_geometry(scp);
    }
    if(bflag_struct.kfull) 
    { advance_z(scp,cold);
      divergence_skew(Comp::w(),tempflag2,cold);
      update_phi(cold,scp);
      reconstruct_geometry(scp);
    }

  /* zyx */
  } else if(time->current_step()%(2*bflag_struct.dim)==label_adv[5]) {
    //OPR(label_adv[5]);
    if(bflag_struct.kfull) 
    { advance_z(scp,cold);
      divergence_skew(Comp::w(),tempflag2,cold);
      update_phi(cold,scp);
      reconstruct_geometry(scp);
    }
    if(bflag_struct.jfull) 
    { advance_y(scp,cold);
      divergence_skew(Comp::v(),tempflag2,cold);
      update_phi(cold,scp);
      reconstruct_geometry(scp);
    }
    if(bflag_struct.ifull) 
    { advance_x(scp,cold);
      divergence_skew(Comp::u(),tempflag2,cold);
      update_phi(cold,scp);
      reconstruct_geometry(scp);
    }

  /* error */
  } else {
    OMS(Catastrophic error. Exiting);
  }

  return;
}

/******************************************************************************/
void VOF::divergence_skew(const Comp & m, const ScalarInt & marker,
                          Scalar & cellvol) {
/******************************************************************************/
  real dt = time->dt();

  for_ijk(i,j,k) {
    if(marker[i][j][k])
      cellvol[i][j][k] += u->divergence(m,i,j,k)*dt*dV(i,j,k);
  }
 
  return;
}
