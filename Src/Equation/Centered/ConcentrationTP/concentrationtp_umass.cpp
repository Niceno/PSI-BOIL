#include "concentrationtp.h"

/***************************************************************************//**
*  \brief Computes mass average velocity from the divergence-free uvol.
*******************************************************************************/
void ConcentrationTP::compute_umass(Vector & umass, const Vector & uvol,
                                    const Matter * gas,
                                    const Scalar * diff_eddy) {

  boil::timer.start("concentrationtp umass");

  Comp m = Comp::u();
  for_wvmijk(umass,m,i,j,k) {  /* don't use vmijk. wall will be skipped! */
    real phip, phim;
    phip = phi[i  ][j][k]; 
    phim = phi[i-1][j][k];

    /* contains density of diffusing species */
    real dc = dcoef->value(m,i,j,k);
    if(diff_eddy){ /* should be rescaled by rho_species/rho -> nonlinear! */
      dc += ((*diff_eddy)[i-1][j][k] +
             (*diff_eddy)[i  ][j][k] )/(2.0*turbS);
    }

    const real r = gas->rho(m,i,j,k);
    const real r_dif = gas->rho(1);
    const real r_env = gas->rho(0);
    real clrf = surface_color(m,i,j,k);
    if(matter_sig==Sign::neg()) clrf = 1.-clrf;
    
    umass[m][i][j][k] = uvol[m][i][j][k]
                      - clrf * dc/r_dif * (r_dif - r_env)/r
                      * (phip-phim) / (umass.dxc(m,i));
  }

  m = Comp::v();
  for_wvmijk(umass,m,i,j,k) {  /* don't use vmijk. wall will be skipped! */
    real phip, phim;
    phip = phi[i][j  ][k];
    phim = phi[i][j-1][k];

    /* contains density of diffusing species */
    real dc = dcoef->value(m,i,j,k);
    if(diff_eddy){ /* should be rescaled by rho_species/rho -> nonlinear! */
      dc += ((*diff_eddy)[i][j-1][k] +
             (*diff_eddy)[i][j  ][k] )/(2.0*turbS);
    }

    const real r = gas->rho(m,i,j,k);
    const real r_dif = gas->rho(1);
    const real r_env = gas->rho(0);
    real clrf = surface_color(m,i,j,k);
    if(matter_sig==Sign::neg()) clrf = 1.-clrf;

    umass[m][i][j][k] = uvol[m][i][j][k]
                      - clrf * dc/r_dif * (r_dif - r_env)/r
                      * (phip-phim) / (umass.dyc(m,j));
  }

  m = Comp::w();
  for_wvmijk(umass,m,i,j,k) {  /* don't use vmijk. wall will be skipped! */
    real phip, phim;
    phip = phi[i][j][k  ];
    phim = phi[i][j][k-1];

    /* contains density of diffusing species */
    real dc = dcoef->value(m,i,j,k);
    if(diff_eddy){ /* should be rescaled by rho_species/rho -> nonlinear! */
      dc += ((*diff_eddy)[i][j][k  ] +
             (*diff_eddy)[i][j][k-1] )/(2.0*turbS);
    }

    const real r = gas->rho(m,i,j,k);
    const real r_dif = gas->rho(1);
    const real r_env = gas->rho(0);
    real clrf = surface_color(m,i,j,k);
    if(matter_sig==Sign::neg()) clrf = 1.-clrf;

    umass[m][i][j][k] = uvol[m][i][j][k]
                      - clrf * dc/r_dif * (r_dif - r_env)/r
                      * (phip-phim) / (umass.dzc(m,k));
  }

  umass.exchange_all();

  boil::timer.stop("concentrationtp umass");
}
