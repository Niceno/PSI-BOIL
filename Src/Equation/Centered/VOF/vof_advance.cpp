#include "vof.h"

/******************************************************************************/
void VOF::advance() {
 
  boil::timer.start("vof advance");

  /* calculate liquid velocity, uses stmp and stmp2 */
  cal_liq_vel();

/* this is the startpoint of using 
 * stmp (phi*V),stmp2 (residual fext),stmp3 (flux multiplier due to vol exp) */

  /*------------------------------+
  |  source term for phase change |
  +------------------------------*/
  for_aijk(i,j,k){   //must be aijk for insert boundary
   #if 1
    /* fext cut */
    real ftot = fext[i][j][k];
    real fval = fext_cut(i,j,k,ftot);

    real phival = phi[i][j][k]+time->dt()*fval;
    /* maybe we can re-apply fext after advance -> stmp2 stores the value */
    stmp2[i][j][k] = ftot-fval;

    /* liquid content */
    stmp[i][j][k] = phival * dV(i,j,k);

    /* intermediate phi value */
    /* fext is related to mdot as follows:
     * fext = -m'''/rhol  */
    real denscoef = 1.0-rhol/rhov;
    real denom = std::max(1.0 + denscoef*fval*time->dt(),boil::pico);
    phival /= denom;
    phi  [i][j][k] = std::max(0.0,std::min(1.0,phival));
    stmp3[i][j][k] = denom;

  //if(fabs(ftot)>boil::atto&&j==2)
  //  boil::aout<<"VOF::advance_I: "<<i<<" "<<phi[i][j][k]<<" "<<ftot<<" "<<fval<<boil::endl;
   #else
    phi[i][j][k]=phi[i][j][k]+time->dt()*fext[i][j][k];
   #endif
  }
  //phi.bnd_update();
  //update_at_walls();

  phi.exchange_all();
  stmp2.exchange();
  stmp3.exchange();

  //boil::aout<<"VOF::advance: "<<stmp[13][2][2]<<" "<<stmp[14][2][2]<<" "<<stmp[15][2][2]<<" | "<<(*u)[Comp::u()][14][2][2]<<" "<<(*u)[Comp::u()][15][2][2]<<" | "<<uliq[Comp::u()][14][2][2]<<" "<<uliq[Comp::u()][15][2][2]<<"\n";

#if 0
  adens.exchange();
  /* calculate phi in staggered cells */
  if(bndclr)
    cal_bndclr();
#endif
  // advance in x-direction
  advance_x();

  //boil::aout<<"VOF::advance: "<<stmp[13][2][2]<<" "<<stmp[14][2][2]<<" "<<stmp[15][2][2]<<" | "<<(*u)[Comp::u()][14][2][2]<<" "<<(*u)[Comp::u()][15][2][2]<<" | "<<uliq[Comp::u()][14][2][2]<<" "<<uliq[Comp::u()][15][2][2]<<"\n";
#if 1
  // advance in y-direction
  advance_y();
  
  // advance in z-direction
  advance_z();
#endif
  
  /* superpose fluxes */
  superpose();

  // update phi
  for_ijk(i,j,k){
    real phi_tmp = stmp[i][j][k] / dV(i,j,k);
#if 0
    // limit C
    phi[i][j][k] = std::min(1.0,std::max(0.0,phi_tmp));
    if(phi_tmp>1.0+boil::pico || phi_tmp< -boil::pico){
      std::cout.setf(std::ios_base::scientific);
      std::cout<<"limit phi "<<phi_tmp<<" i "<<i<<" j "<<j<<" k "<<k<<"\n";
      std::cout.unsetf(std::ios_base::floatfield);
    }
#else
    // unlimit C
    phi[i][j][k] = phi_tmp;

  #if 1
    /* fext cut II */
    real ftot = stmp2[i][j][k];
    real fval = fext_cut(i,j,k,ftot);
    phi[i][j][k]=phi[i][j][k]+time->dt()*fval;
  //if(fabs(ftot)>boil::atto&&j==2)
  //  boil::aout<<"VOF::advance_II: "<<i<<" "<<phi[i][j][k]<<" "<<ftot<<" "<<fval<<boil::endl;
  #endif
#endif
  }

/* this is the endpoint of using stmp, stmp2, stmp3 */

#if 0  
  vf_limiter();
#endif

  phi.bnd_update();
  phi.exchange_all();

  boil::timer.stop("vof advance");
#if 1
  boil::timer.start("vof ancillary");

  #if 0
  sharpen();
  #endif

  /*-------------------------------+
  |  normal vector at cell center  |
  +-------------------------------*/
  norm_cc(phi);

  /* calculate alpha in cells */
  extract_alpha();

  /* calculate free surface position */
  cal_fs3();

  /* prerequisite for marching cubes */
  update_at_walls();

  /* calculate the real-space normal vector */
  true_norm_vect(); 

  /* calculate area */
  cal_adens();
  //cal_adens_geom(adens);
  //set_adens(adensgeom);

  /* calculate phi in staggered cells */
  if(bndclr)
    cal_bndclr();

  boil::timer.stop("vof ancillary");

  /* calculate curvature */
  //curv_HF();
#endif
  return;
}

