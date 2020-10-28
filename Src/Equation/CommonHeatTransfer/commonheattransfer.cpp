#include "commonheattransfer.h"

/******************************************************************************/
CommonHeatTransfer::CommonHeatTransfer(const Scalar & TPR,
                                       Topology * TOPO, 
                                       const TIF & tintmodel,
                                       Matter * F, Matter * S,
                                       HTWallModel * HTM) :
  tpr(&TPR),
  topo(TOPO),
  tifmodel(tintmodel),
  flu(F),
  sol(S),
  htwallmodel(HTM),
  bndtpr( *TOPO->get_fs().domain() )
{
  for_m(m) bndtpr(m) = topo->get_fs()(m).shape(); /* a mistake? */
  
  val_rhol = fluid()->rho(1);
  val_rhov = fluid()->rho(0);
  val_cpl  = fluid()->cp(1);
  val_cpv  = fluid()->cp(0);
  val_lambdal = fluid()->lambda(1);
  val_lambdav = fluid()->lambda(0);
  turbP=0.9;

  /* heat transfer wall model should be equal to the one of enthalpy,
   * that's why it is a pointer. however, the default argument is a
   * nullptr, we need to fix that */
  if(!htwallmodel) {
    default_value_for_htwallmodel = true;
    htwallmodel = new HTWallModel();
  } else {
    default_value_for_htwallmodel = false;
  }
}

/******************************************************************************/
CommonHeatTransfer::~CommonHeatTransfer() {
  if(default_value_for_htwallmodel)
    delete htwallmodel;
}

/******************************************************************************/
