#include "momentum.h"

/******************************************************************************/
void Momentum::outlet() {
/*----------------------------------+ 
|  set bulk velocity at the outlet  |
+----------------------------------*/

  /*---------------------------------+
  |  if outlet does not exist, exit  |
  +---------------------------------*/
  if( 
      !u.bc( Comp::u() ).exists( BndType::outlet() ) &&
      !u.bc( Comp::v() ).exists( BndType::outlet() ) && 
      !u.bc( Comp::w() ).exists( BndType::outlet() ) 
    ) return;

  real aox=0.0, aoy=0.0, aoz=0.0, ubo=0.0, fubo=0.0;

  /*---------------------------------------+
  |  get volume flux the inlet and outlet  |
  +---------------------------------------*/
  /* v_phase_change() should be called from main.cpp, if phase change occurs */
  const real volf_in  = volf_bct( BndType::inlet() )
                      + volf_bct( BndType::insert() )
                      + v_phase_change;
  real volf_out = volf_bct( BndType::outlet(), &aox, &aoy, &aoz ); 

  /* compute bulk velocity */
  ubo = volf_out / (aox + aoy + aoz);

  /* outlet bulk velocity is negative if it leaves the domain.  
     let's make it positive from here on */
  fubo = fabs(ubo);

  real ratio;
  if(volf_out==0.0 && volf_in==0.0){
    ratio=1.0;
  } else if(volf_out==0.0 && volf_in!=0.0){
    ratio=1.0e+12;
  } else {
    ratio= -volf_in / volf_out;
  }

  /*-------------------------------+
  |  step 1: extrapolate velocity  |
  +-------------------------------*/
  extrapolate_outlet_velocity(fubo,ratio);

  /*------------------------------------------------------+
  |  step 2: scale velocity to satisfy mass conservation  |
  +------------------------------------------------------*/

  /* due to step 1, volf_out has changed */
  volf_out = volf_bct( BndType::outlet(), &aox, &aoy, &aoz ); 
  ubo = volf_in / (aox + aoy + aoz);
  if(volf_out==0.0 && volf_in==0.0){
    ratio=1.0;
  } else if(volf_out==0.0 && volf_in!=0.0){
    ratio=1.0e+12;
  } else {
    ratio= -volf_in / volf_out;
  }

  scale_outlet_velocity(ubo,ratio);

#if 0
  /* debug */
  volf_out = volf_bct( BndType::outlet(), &aox, &aoy, &aoz );
  ubo = volf_out / (aox + aoy + aoz);
  if(volf_out==0.0 && volf_in==0.0){
    ratio=1.0;
  } else if(volf_out==0.0 && volf_in!=0.0){
    ratio=1.0e+12;
  } else {
    ratio= -volf_in / volf_out;
  }
  boil::oout<<"Momentum::outlet: "<<time->current_time()<<" "<<volf_in<<" "<<volf_out<<" "<<ratio<<boil::endl;
  //exit(0);
#endif

  return;
}
