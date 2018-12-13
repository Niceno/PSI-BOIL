#include "enthalpytif.h"

/***************************************************************************//**
 *  Call to tifmodel
******************************************************************************/
void EnthalpyTIF::tint_field(const real factor, const bool iter) {
  tifmodel.tint_field(factor, iter); 

  return;
}
