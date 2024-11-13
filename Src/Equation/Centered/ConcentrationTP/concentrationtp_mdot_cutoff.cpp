#include "concentrationtp.h"

/***************************************************************************//**
*  cut off mass transfer rate if volume-fraction-of-vapor is smaller than crit
*******************************************************************************/
void ConcentrationTP::mdot_cutoff(Scalar & mdot) { 

  real sum_mdot_pos = 0.0;
  real sum_mdot_neg = 0.0;
  for_ijk(i,j,k) {
    real col = vfval(i,j,k);
    if(matter_sig==Sign::neg()) col = 1.0 - col;
    // col: volume fraction of vapor
    if (col<col_crit) {
      //if( mdot[i][j][k]!=0.0 ) {
      //  std::cout<<"mdot_cutoff: "<<i<<" "<<j<<" "<<k<<" "<<vfval(i,j,k)
      //           <<" "<<mdot[i][j][k]<<"\n";
      //}
      mdot[i][j][k]=0.0;
    }
    if (mdot[i][j][k]>0.0) {
      sum_mdot_pos += mdot[i][j][k]*dV(i,j,k);
    } else {
      sum_mdot_neg += mdot[i][j][k]*dV(i,j,k);
    }
  }
  mdot.exchange();

  boil::cart.sum_real(&sum_mdot_pos);
  boil::cart.sum_real(&sum_mdot_neg);

  boil::oout<<"mdot_cutoff: "<< time->current_time()
            <<" sum_mdot_pos[kg/s]= "<<sum_mdot_pos
            <<" sum_mdot_neg[kg/s]= "<<sum_mdot_neg<<"\n";
}
