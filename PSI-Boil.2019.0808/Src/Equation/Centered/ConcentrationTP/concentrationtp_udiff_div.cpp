#include "concentrationtp.h"

/***************************************************************************//**
*  \brief Computes divergence of the diffusive velocity.
*******************************************************************************/
void ConcentrationTP::compute_udiff_div(Scalar & p, const Property * mu_fluid,
                                        const Vector & udiff) {
  real divu, um, up;
  real dxc,dyc,dzc;
  real gm, gp;

  real mu;
  Comp m;

  for_ijk(i,j,k) {
    divu = 0.0;
    
    /* i-direction */
    m = Comp::u();
    um = -udiff[m][i  ][j][k];
    up = -udiff[m][i+1][j][k];
    dxc= p.dxc(i);

    gm = surface_color(Sign::neg(),m,i,j,k);
    gp = surface_color(Sign::pos(),m,i,j,k);
    if(matter_sig==Sign::neg()) {
      gm = 1.-gm;
      gp = 1.-gp;
    }

    divu += (gp*up-gm*um)/dxc;

    /* j-direction */
    m = Comp::v();
    um = -udiff[m][i][j  ][k];
    up = -udiff[m][i][j+1][k];
    dyc= p.dyc(j);

    gm = surface_color(Sign::neg(),m,i,j,k);
    gp = surface_color(Sign::pos(),m,i,j,k);
    if(matter_sig==Sign::neg()) {
      gm = 1.-gm;
      gp = 1.-gp;
    }

    divu += (gp*up-gm*um)/dyc;

    /* k-direction */
    m = Comp::w();
    um = -udiff[m][i][j][k  ];
    up = -udiff[m][i][j][k+1];
    dzc= p.dzc(k);
                 
    gm = surface_color(Sign::neg(),m,i,j,k);
    gp = surface_color(Sign::pos(),m,i,j,k);
    if(matter_sig==Sign::neg()) {
      gm = 1.-gm;
      gp = 1.-gp;
    }

    divu += (gp*up-gm*um)/dzc;
   
    /* properly re-scale */
    mu = mu_fluid->value(i,j,k);
    divu = 2.0*mu*divu/3.0;
 
    /* add */
    p[i][j][k] += divu;
  }

  p.bnd_update();  /* questionable */
  p.exchange_all();

  return;
}
