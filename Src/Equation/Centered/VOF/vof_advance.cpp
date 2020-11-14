#include "vof.h"

/******************************************************************************/
void VOF::advance(const bool anci) {
/***************************************************************************//**
*  \brief Advance VOF in time.
*******************************************************************************/
  advance(phi,anci);
}

void VOF::advance(const Scalar & scp, const bool anci) {
 
  boil::timer.start("vof advance");
  
  /* reset phi */
  phi = scp;

  /*---------------------+
  |  phase change shift  |
  +---------------------*/
  advance_phase_change(phi);

  /*-----------------------+
  |  reconstruct geometry  |
  +-----------------------*/
  reconstruct_geometry(phi);

  /*--------------------+
  |  geometric advance  |
  +--------------------*/
  advance_geometric(phi);

  boil::timer.stop("vof advance");

  /*---------------------------+
  |  ancillary ITM parameters  |
  +---------------------------*/
  if(anci)
    ancillary(color(),false);

  return;
}

/******************************************************************************/
void VOF::advance_with_extrapolation(const bool anci, const ResTol & restol,
                                     const Vector & umixed, const Scalar & fext, 
                                     const Matter * fluid_1, Vector * uvw_1, 
                                     const Matter * fluid_2, Vector * uvw_2) {
/***************************************************************************//**
*  \brief Advance VOF, extrapolating velocities in the process.
*******************************************************************************/
  advance_with_extrapolation(phi,anci,restol,umixed,fext,
                             fluid_1,uvw_1,fluid_2,uvw_2);
}

void VOF::advance_with_extrapolation(const Scalar & scp, 
                                     const bool anci, const ResTol & restol,
                                     const Vector & umixed, const Scalar & fext, 
                                     const Matter * fluid_1, Vector * uvw_1,
                                     const Matter * fluid_2, Vector * uvw_2) {

  boil::timer.start("vof advance");

  /* reset phi */
  phi = scp;

  /*---------------------+
  |  phase change shift  |
  +---------------------*/
  advance_phase_change(phi);

  /*-----------------------+
  |  reconstruct geometry  |
  +-----------------------*/
  reconstruct_geometry(phi);

  boil::timer.stop("vof advance");

  /*-------------------------+
  |  extrapolate velocities  |
  +-------------------------*/
  if(fluid_1) {
    extrapolate_velocity(color(),fext,fluid_1,umixed,(*uvw_1),
                         restol,Sign::pos(),true); /* true for int. flagging */
  }
  if(fluid_2) {
    bool iflagging = false;
    if(!fluid_1) iflagging = true;

    extrapolate_velocity(color(),fext,fluid_2,umixed,(*uvw_2),
                         restol,Sign::neg(),iflagging);
  }

  boil::timer.start("vof advance");

  /*--------------------+
  |  geometric advance  |
  +--------------------*/
  advance_geometric(phi);

  boil::timer.stop("vof advance");

  /*---------------------------+
  |  ancillary ITM parameters  |
  +---------------------------*/
  if(anci)
    ancillary(color(),false);

  return;
}

/******************************************************************************/
void VOF::advance_with_extrapolation(const bool anci, const ResTol & restol,
                                     const Vector & umixed,
                                     const Matter * fluid_1, Vector * uvw_1, 
                                     const Matter * fluid_2, Vector * uvw_2) {
/***************************************************************************//**
*  \brief Advance VOF, extrapolating velocities in the process.
*******************************************************************************/
  advance_with_extrapolation(phi,anci,restol,umixed,
                             fluid_1,uvw_1,fluid_2,uvw_2);
}

void VOF::advance_with_extrapolation(const Scalar & scp, 
                                     const bool anci, const ResTol & restol,
                                     const Vector & umixed,
                                     const Matter * fluid_1, Vector * uvw_1,
                                     const Matter * fluid_2, Vector * uvw_2) {

  boil::timer.start("vof advance");

  /* reset phi */
  phi = scp;

  /*---------------------+
  |  phase change shift  |
  +---------------------*/
  advance_phase_change(phi);

  /*-----------------------+
  |  reconstruct geometry  |
  +-----------------------*/
  reconstruct_geometry(phi);

  boil::timer.stop("vof advance");

  /*-------------------------+
  |  extrapolate velocities  |
  +-------------------------*/
  if(fluid_1) {
    extrapolate_velocity(color(),fluid_1,umixed,(*uvw_1),
                         restol,Sign::pos(),true); /* true for int. flagging */
  }
  if(fluid_2) {
    bool iflagging = false;
    if(!fluid_1) iflagging = true;

    extrapolate_velocity(color(),fluid_2,umixed,(*uvw_2),
                         restol,Sign::neg(),iflagging);
  }

  boil::timer.start("vof advance");

  /*--------------------+
  |  geometric advance  |
  +--------------------*/
  advance_geometric(phi);

  boil::timer.stop("vof advance");

  /*---------------------------+
  |  ancillary ITM parameters  |
  +---------------------------*/
  if(anci)
    ancillary(color(),false);

  return;
}

/******************************************************************************/
void VOF::advance_phase_change(Scalar & scp) {
  /*------------------------------+
  |  source term for phase change |
  +------------------------------*/
  for_ijk(i,j,k){
    scp[i][j][k]=scp[i][j][k]+time->dt()*fext[i][j][k];
  }
  scp.bnd_update();
  scp.exchange_all();

  return;
}

void VOF::advance_geometric(Scalar & scp) {

  /* liquid volume */
  for_aijk(i,j,k){
    cold[i][j][k] = scp[i][j][k] * dV(i,j,k);
  }

  for_m(m)
    vflow(m) = 0.;

  /* different methods */
  if(advect_method==AdvectionMethod::NaiveSplit()) {
    advect_naive(scp);
  } else if(advect_method==AdvectionMethod::ReconstructedSplit()) {
    advect_reconstructed(scp);
  } else if(advect_method==AdvectionMethod::BoundedSplit()) {
    advect_bounded(scp);
  }

  return;
}

/******************************************************************************/
void VOF::update_phi(const Scalar & cellvol, Scalar & scp) {

  /* update phi */
  if(limit_color) {
    for_ijk(i,j,k){
      real phi_tmp = cellvol[i][j][k] / dV(i,j,k);
      scp[i][j][k] = std::min(1.0,std::max(0.0,phi_tmp));
    }
  } else {
    int ierr=0;
    for_ijk(i,j,k){
      scp[i][j][k] = cellvol[i][j][k] / dV(i,j,k);
      if (scp[i][j][k]<-1.0||scp[i][j][k]>2.0) {
        boil::aout<<"Error!!! Too small or too large phi in vof_advance. phi= "
                  <<scp[i][j][k]<<"\n";
        boil::aout<<"proc= "<<boil::cart.iam()
                  <<" (i,j,k)= "<<i<<" "<<j<<" "<<k<<"\n";
        boil::aout<<"(x,y,z)= "<<phi.xc(i)<<" "<<phi.yc(j)<<" "<<phi.zc(k)<<"\n";
        ierr=1;
      }
    }
    boil::cart.sum_int(&ierr);
    if(ierr>=1) {
      boil::plot->plot(vflow,scp,color(),cellvol,"flux-phi-clr-vol", time->current_step());
      exit(0);
    }
  }

  scp.bnd_update();
  scp.exchange_all();
  
  return;
}
