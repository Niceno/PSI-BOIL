#include "centered.h"

/******************************************************************************/
void Centered::set_active_flag(ScalarInt & activeflag) { 
/***************************************************************************//**
*  Only the local value should be accessed! (no bnd_update or exchange)
*  \brief Set value of the active flag needed by the multigrid solver.
*******************************************************************************/

  for_ijk(i,j,k) {
    activeflag[i][j][k] = activeflag.domain()->ibody().on_p(i,j,k);
  }

  return;
}
