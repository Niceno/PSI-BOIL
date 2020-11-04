#include "commonheattransfer.h"

/******************************************************************************/
CommonHeatTransfer::CommonHeatTransfer(const Scalar & TPR,
                                       Topology * TOPO, 
                                       const TIF & tintmodel,
                                       Matter * F, Matter * S) :
  tpr(&TPR),
  topo(TOPO),
  tifmodel(tintmodel),
  flu(F),
  sol(S),
  bndtpr_sol( *TOPO->get_fs().domain() ),
  bndtpr_flu( *TOPO->get_fs().domain() )
{
  for_m(m) bndtpr_sol(m) = topo->get_fs()(m).shape(); /* a mistake? */
  for_m(m) bndtpr_flu(m) = topo->get_fs()(m).shape(); /* a mistake? */
  
  val_rhol = fluid()->rho(1);
  val_rhov = fluid()->rho(0);
  val_cpl  = fluid()->cp(1);
  val_cpv  = fluid()->cp(0);
  val_lambdal = fluid()->lambda(1);
  val_lambdav = fluid()->lambda(0);
  turbP=0.9;
  dirac_wall_source_val = 0.0;
  wall_resistance_val = 0.0;
  int_resistance_vap_val = 0.0;
  int_resistance_liq_val = 0.0;
}

/******************************************************************************/
CommonHeatTransfer::~CommonHeatTransfer() {
}

/******************************************************************************/
