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
  boil::oout<<"v_phase_change= "<<v_phase_change<<"\n";
  std::cout<< std::setprecision(8);
  std::cout.unsetf(std::ios_base::floatfield);
  return(v_phase_change);
}

/*-----------------------------------------------------------------------------+
 '$Id: momentum_vol_phase_change.cpp,v 1.3 2014/02/03 14:37:21 sato Exp $'/
+-----------------------------------------------------------------------------*/
