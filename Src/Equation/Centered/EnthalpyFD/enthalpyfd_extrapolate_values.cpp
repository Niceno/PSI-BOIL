#include "enthalpyfd.h"

/***************************************************************************//**
*  \brief Extrapolate value stencil for gradient calculation.
*******************************************************************************/
void EnthalpyFD::extrapolate_values(std::vector<StencilPoint> & stencil,
                                    const StencilPoint & ctp,
                                    const StencilPoint & ctm) {

  /* extrapolate west */
  if(ctm.idx>0) {
  
    std::vector<StencilPoint> cut_stencil;
  }

  return;
}
