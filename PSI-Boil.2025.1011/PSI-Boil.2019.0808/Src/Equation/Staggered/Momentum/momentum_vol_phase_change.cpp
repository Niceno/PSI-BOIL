#include "momentum.h"
#include <iomanip>

/******************************************************************************/
real Momentum::vol_phase_change(Scalar * psrc_old) {

  v_phase_change=0.0;
  for_vijk((*psrc_old),i,j,k)
    v_phase_change+=(*psrc_old)[i][j][k] * time->dt();

  boil::cart.sum_real(&v_phase_change);

  std::cout.setf(std::ios_base::scientific);
  std::cout<< std::setprecision(16);
  boil::oout<<"v_phase_change= "<<time->current_time()<<" "<<v_phase_change<<"\n";
  std::cout<< std::setprecision(8);
  std::cout.unsetf(std::ios_base::floatfield);
  return(v_phase_change);
}
